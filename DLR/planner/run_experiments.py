"""Experiment: does a wind-aware (energy) planner beat the shortest path on the DLR campus?

For every solved wind direction (multiples of 30 deg only) at its EDMO mean speed AND at 8 m/s, plan
N seeded random start/goal pairs over the campus strip with
    shortest path           (cost = distance)
    energy, fixed 10 m/s    (cost = J, ground speed 10 m/s)
    energy, fixed 15 m/s    (sensitivity: closer to the max-range airspeed)
    energy, best speed      (the drone picks its ground speed on every edge)
and score every path in the SAME true wind field with dense sampling (evaluate_path). The baseline for
each energy planner is the shortest path scored under that planner's flight regime.
Both with raw graph paths and after line-of-sight smoothing (each planner smoothed with its own cost).

    python run_experiments.py [--n 40] [--seed 0] [--strong 8] [--no-figs]
Outputs: results/plan_results.csv, results/plan_summary.json, fig/plan_*.png
"""
import argparse, csv, json, sys, time
from pathlib import Path
import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "atlas"))
from wind_atlas import WindAtlas                                           # noqa: E402
from planner import Environment, WindField, PlanningGraph, evaluate_path   # noqa: E402
from energy import max_range_speed                                         # noqa: E402

V_MR = max_range_speed()      # still-air max-range speed (~18.3 m/s)

# EDMO airfield 2020-2025, FROM-direction sectors: (% of hours, mean speed m/s). Calm (2.9 %) is excluded.
EDMO = {0: (3.6, 2.95), 30: (7.1, 3.15), 60: (8.3, 3.34), 90: (6.1, 2.94), 120: (3.1, 2.00), 150: (6.8, 2.07),
        180: (8.2, 2.11), 210: (15.8, 3.02), 240: (22.9, 5.15), 270: (10.3, 4.31), 300: (4.1, 3.39), 330: (3.7, 3.08)}
# planner tag -> (graph cost mode, v_g, evaluation mode)
PLANNERS = {"fix10": ("energy", 10.0, "energy"), "fix15": ("energy", 15.0, "energy"),
            "best": ("energy_bestspeed", None, "energy_bestspeed")}
SHORT = "distance"


def metrics(P, wind, env, tag):
    mode, vg, ev = PLANNERS[tag]
    return evaluate_path(P, wind, env, ev, v_g=vg or 10.0)


def pick_pairs(env, G, n, seed, min_sep=250.0, near=60.0, agl=(15.0, 30.0)):
    """Seeded random endpoint pairs: free nodes within `near` m of a building, 15-30 m above terrain."""
    bmask = (~env.valid) & ((env.z[None, None, :] - env.terrain[:, :, None]) > 3.0)   # solid well above bare ground
    bdist = distance_transform_edt(~bmask, sampling=env.h)
    f = G.node_fine
    h = G.node_xyz[:, 2] - env.terrain[f[:, 0], f[:, 1]]
    cand = np.nonzero((h >= agl[0]) & (h <= agl[1]) & (bdist[f[:, 0], f[:, 1], f[:, 2]] <= near))[0]
    rng = np.random.default_rng(seed)
    pairs, seen = [], set()
    while len(pairs) < n * 4:                                   # oversample; connectivity filter happens later
        a = int(rng.choice(cand))
        far = cand[np.linalg.norm(G.node_xyz[cand] - G.node_xyz[a], axis=1) >= min_sep]
        if far.size == 0:
            continue
        b = int(rng.choice(far))
        if (a, b) not in seen and (b, a) not in seen:
            seen.add((a, b)); pairs.append((a, b))
    return pairs, len(cand)


def path_cost(A, ids):
    return float(np.asarray(A[ids[:-1], ids[1:]]).sum())


COMPUTE_VERSION = "4"          # bump when the experiment logic changes; part of the cache key


def code_hash():
    import hashlib
    return hashlib.sha1((HERE / "planner.py").read_bytes() + COMPUTE_VERSION.encode() + repr(PLANNERS).encode()).hexdigest()[:12]


def compute_case(theta, case, atlas, env, G, pairs, short_ids, speed):
    """Plan every pair/leg for one (direction, speed case); returns rows, the most interesting pair, timing."""
    wind = WindField.from_atlas(atlas, theta, speed)
    t = time.time()
    mats = {tag: G.matrix(wind, m, vg if vg else 10.0) for tag, (m, vg, _) in PLANNERS.items()}
    timing = dict(matrices_s=time.time() - t); t_solve = []
    rows, best = [], None
    for pi, ((a, b), sid) in enumerate(zip(pairs, short_ids)):
        legs = {"AB": (a, b, sid), "BA": (b, a, sid[::-1])}
        plan_ids = {leg: {} for leg in legs}
        for tag, A in mats.items():
            ub = {leg: path_cost(A, s_ids) for leg, (s, g, s_ids) in legs.items()}   # shortest path re-priced = valid upper bound
            t = time.time()
            dist, pred = G.solve(A, [a, b], limit=max(ub.values()) * (1 + 1e-9))
            t_solve.append(time.time() - t)
            for k, leg in enumerate(legs):
                s, g, _ = legs[leg]
                plan_ids[leg][tag] = G.extract(pred[k], s, g)
        for leg, (s, g, s_ids) in legs.items():
            if any(plan_ids[leg][tag] is None for tag in PLANNERS):
                continue
            paths = {"short": G.node_xyz[s_ids], **{tag: G.node_xyz[plan_ids[leg][tag]] for tag in PLANNERS}}
            # smoothing: the SHORTEST path gets plain geometric string pulling (the naive "fly straight" drone);
            # every energy path gets ENERGY-AWARE string pulling (shortcut only if free AND not costlier in the true wind,
            # dense 2 m sampling, same speed policy as the planner)
            sm = {"short": G.smooth(paths["short"], mode=SHORT)}
            for tag, (m, vg, _) in PLANNERS.items():
                sm[tag] = G.smooth(paths[tag], wind, m, vg if vg else 10.0)
            P0, P1 = G.node_xyz[s], G.node_xyz[g]
            row = dict(theta=theta, case=case, wind_speed=round(speed, 2), pair=pi, leg=leg,
                       sep_m=round(float(np.linalg.norm(P0 - P1)), 1),
                       **{f"a{k}": round(float(P0[i]), 1) for i, k in enumerate("xyz")},
                       **{f"b{k}": round(float(P1[i]), 1) for i, k in enumerate("xyz")})
            for lab, pp in (("raw", paths), ("sm", sm)):
                for tag in PLANNERS:
                    ms, mp = metrics(pp["short"], wind, env, tag), metrics(pp[tag], wind, env, tag)
                    key = f"{lab}_{tag}_"
                    row[key + "E_short"], row[key + "E_plan"] = ms["energy_J"], mp["energy_J"]
                    row[key + "len_short"], row[key + "len_plan"] = ms["length_m"], mp["length_m"]
                    row[key + "t_short"], row[key + "t_plan"] = ms["time_s"], mp["time_s"]
                    row[key + "hw_short"], row[key + "hw_plan"] = ms["mean_headwind"], mp["mean_headwind"]
                    row[key + "clear_plan"] = mp["min_clearance_m"]
                    row[key + "save_pct"] = 100 * (1 - mp["energy_J"] / ms["energy_J"])
                    row[key + "extra_len_pct"] = 100 * (mp["length_m"] / ms["length_m"] - 1)
                    row[key + "dt_s"] = mp["time_s"] - ms["time_s"]
            # straight path at a CONSTANT still-air max-range ground speed (isolates wind-aware speed scheduling)
            mr = evaluate_path(sm["short"], wind, env, "energy", v_g=V_MR)
            row["sm_mr_E_short"], row["sm_mr_t_short"] = mr["energy_J"], mr["time_s"]
            rows.append(row)
            if case == "edmo" and leg == "AB" and (best is None or row["sm_best_save_pct"] > best["save"]):
                best = dict(save=row["sm_best_save_pct"], pair=pi, speed=speed,
                            paths={k: v.tolist() for k, v in paths.items()}, sm={k: v.tolist() for k, v in sm.items()})
    timing["dijkstra_call_s"] = float(np.mean(t_solve))
    return rows, best, timing


def case_stats(R):
    """Aggregate one (direction, case): raw / smoothed(fair) saving statistics and the separation of effects."""
    q = lambda v, p: float(np.percentile(v, p))
    g = lambda k: np.array([r[k] for r in R])
    ent = dict(wind_speed=R[0]["wind_speed"], n_legs=len(R))
    byp = {}
    for r in R:
        byp.setdefault(r["pair"], {})[r["leg"]] = r
    both = [d for d in byp.values() if len(d) == 2]
    for lab in ("raw", "sm"):
        for tag in PLANNERS:
            s = g(f"{lab}_{tag}_save_pct")
            e = dict(save_mean=float(s.mean()), save_median=float(np.median(s)), save_p90=q(s, 90), save_min=float(s.min()),
                     save_max=float(s.max()), frac_positive=float((s > 1e-9).mean()),
                     extra_len_mean=float(g(f"{lab}_{tag}_extra_len_pct").mean()), dt_mean_s=float(g(f"{lab}_{tag}_dt_s").mean()),
                     E_short_mean_J=float(g(f"{lab}_{tag}_E_short").mean()))
            rt = [100 * (1 - (d["AB"][f"{lab}_{tag}_E_plan"] + d["BA"][f"{lab}_{tag}_E_plan"]) /
                         (d["AB"][f"{lab}_{tag}_E_short"] + d["BA"][f"{lab}_{tag}_E_short"])) for d in both]
            asym = [100 * abs(d["AB"][f"{lab}_{tag}_E_short"] - d["BA"][f"{lab}_{tag}_E_short"]) /
                    (0.5 * (d["AB"][f"{lab}_{tag}_E_short"] + d["BA"][f"{lab}_{tag}_E_short"])) for d in both]
            e["roundtrip_save_mean"] = float(np.mean(rt))
            e["shortest_AB_vs_BA_energy_diff_pct_mean"] = float(np.mean(asym))
            ent[f"{lab}_{tag}"] = e
    for tag in PLANNERS:                     # headline: smoothed shortest (straight) vs energy-aware smoothed energy path
        ent[f"fair_{tag}"] = ent[f"sm_{tag}"]
    Es10, Ep10, Esb, Epb, Emr = (g(k) for k in ("sm_fix10_E_short", "sm_fix10_E_plan", "sm_best_E_short", "sm_best_E_plan", "sm_mr_E_short"))
    ent["effects"] = dict(          # all relative to the naive drone: straight path at fixed 10 m/s (mean % per leg)
        route_only_pct=float(np.mean(100 * (1 - Ep10 / Es10))),          # energy-aware route, still 10 m/s
        speed_only_pct=float(np.mean(100 * (1 - Esb / Es10))),           # straight path, best speed per segment (mostly NOT wind: 10 -> ~18 m/s)
        speed_choice_still_air_pct=float(np.mean(100 * (1 - Emr / Es10))),   # straight path at constant max-range speed vs 10 m/s: no wind logic
        wind_speed_adaptation_pct=float(np.mean(100 * (1 - Esb / Emr))),     # per-segment best speed in the wind vs constant max-range speed: what wind-aware speed scheduling buys
        route_and_speed_pct=float(np.mean(100 * (1 - Epb / Es10))),      # both
        route_given_best_speed_pct=float(np.mean(100 * (1 - Epb / Esb))),  # route benefit when both fly best speed
        direction_asym_fix10_pct=ent["sm_fix10"]["shortest_AB_vs_BA_energy_diff_pct_mean"],   # straight path A->B vs B->A
        direction_asym_best_pct=ent["sm_best"]["shortest_AB_vs_BA_energy_diff_pct_mean"])
    return ent


def build_summary(cases, dirs, G, args, rejected, n_cand, ratio, timing_now):
    flat = lambda e: {**{f"{lab}_{tag}": e[f"{lab}_{tag}"]["save_mean"] for lab in ("raw", "sm") for tag in PLANNERS},
                      **{f"fair_{tag}": e[f"fair_{tag}"]["save_mean"] for tag in PLANNERS},
                      **{f"effect_{k}": v for k, v in e["effects"].items()}}
    summary = dict(directions_available=dirs, n_pairs=args.n, pairs_rejected_unreachable=rejected, seed=args.seed,
                   strong_speed=args.strong, headline=("fair_* = sm_* on purpose: the energy paths are smoothed with ENERGY-AWARE string pulling (shortcut accepted only if collision-free and not costlier in the true wind, dense 2 m sampling, same speed policy), the shortest path with plain geometric pulling; both scored in the same field. "
                             "effects: route_only, speed_only (mostly the 10 -> 18 m/s efficiency gain, not wind), speed_choice_still_air, wind_speed_adaptation (per-segment best speed vs constant max-range speed %.1f m/s), direction_asym" % V_MR),
                   graph=dict(nodes=G.n, edges=G.n_edges, stride=G.stride, margin_m=G.margin, band_agl_m=[G.h_min, G.h_max]),
                   graph_shortest_over_euclid=ratio, timing=timing_now, per_direction={}, weighted={})
    for (theta, case), c in sorted(cases.items()):
        summary["per_direction"][f"{theta}_{case}"] = c["stats"]
    for case in ("edmo", "strong"):
        wts = {th: EDMO[th][0] for th in dirs if (th, case) in cases}
        tot = sum(wts.values())
        keys = flat(cases[(next(iter(wts)), case)]["stats"]).keys()
        summary["weighted"][case] = dict(
            directions=list(wts), weights={str(k): v / tot for k, v in wts.items()},
            expected_save_pct={k: sum(wts[th] * flat(cases[(th, case)]["stats"])[k] for th in wts) / tot for k in keys})
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--strong", type=float, default=8.0)
    ap.add_argument("--only-missing", action=argparse.BooleanOptionalAction, default=True,
                    help="reuse results/cache/<theta>_<case>.json when the atlas file, seed, n, strong and code are unchanged (default)")
    ap.add_argument("--force", action="store_true", help="recompute every direction")
    ap.add_argument("--no-figs", action="store_true")
    args = ap.parse_args()
    force = args.force or not args.only_missing
    (HERE / "results" / "cache").mkdir(parents=True, exist_ok=True); (HERE / "fig").mkdir(exist_ok=True)
    T0 = time.time(); timing = {}

    atlas = WindAtlas()
    dirs = sorted(int(round(d)) for d in atlas.directions if abs(d % 30) < 1e-9 and int(round(d)) in EDMO)
    print("directions used (multiples of 30 deg):", dirs, "| all loaded:", [float(d) for d in atlas.directions])
    t = time.time(); env = Environment.from_atlas(atlas); timing["environment_s"] = time.time() - t
    t = time.time(); G = PlanningGraph(env); timing["graph_build_s"] = time.time() - t
    print(f"graph: {G.n} nodes, {G.n_edges} directed edges, build {timing['graph_build_s']:.1f} s "
          f"(+ distance field {timing['environment_s']:.1f} s)")
    Ad = G.matrix(mode=SHORT)

    # ---- endpoints + distance paths (identical for every direction: same seed) ----------------
    cand_pairs, ncand = pick_pairs(env, G, args.n, args.seed)
    pairs, short_ids, rejected = [], [], 0
    for a, b in cand_pairs:
        if len(pairs) == args.n:
            break
        dist, pred = G.solve(Ad, [a])
        ids = G.extract(pred[0], a, b)
        if ids is None:
            rejected += 1; continue
        pairs.append((a, b)); short_ids.append(ids)
    seps = [np.linalg.norm(G.node_xyz[a] - G.node_xyz[b]) for a, b in pairs]
    rr = [np.linalg.norm(np.diff(G.node_xyz[i], axis=0), axis=1).sum() / s for i, s in zip(short_ids, seps)]
    ratio = dict(mean=float(np.mean(rr)), max=float(np.max(rr)))
    print(f"{len(pairs)} pairs (rejected unreachable: {rejected}), {ncand} candidate endpoints; separation "
          f"{min(seps):.0f}-{max(seps):.0f} m; graph-shortest/Euclid mean {ratio['mean']:.3f} max {ratio['max']:.3f}")

    # ---- per (direction, case): load from cache or compute -----------------------------------
    cases, computed = {}, 0
    ch = code_hash()
    for theta in dirs:
        f = atlas.runs[float(theta)]["file"]
        npz = HERE.parent / "atlas" / "grid.noindex" / f
        s = npz.stat()
        for case in ("edmo", "strong"):
            speed = EDMO[theta][1] if case == "edmo" else args.strong
            key = dict(npz=f, mtime_ns=s.st_mtime_ns, size=s.st_size, seed=args.seed, n=args.n, speed=speed,
                       stride=list(G.stride), margin=G.margin, band=[G.h_min, G.h_max], code=ch)
            cf = HERE / "results" / "cache" / f"{theta}_{case}.json"
            c = None
            if cf.exists() and not force:
                try:
                    c = json.loads(cf.read_text())
                    if c["key"] != key:
                        c = None
                except Exception:
                    c = None
            if c is None:
                t = time.time()
                rows, best, tm = compute_case(theta, case, atlas, env, G, pairs, short_ids, speed)
                c = dict(key=key, rows=rows, best=best, timing=tm, stats=case_stats(rows))
                cf.write_text(json.dumps(c)); computed += 1
                print(f"theta {theta:3d} {case:6s} {speed:5.2f} m/s COMPUTED in {time.time()-t:.0f} s "
                      f"(matrices {tm['matrices_s']:.0f} s, Dijkstra {tm['dijkstra_call_s']*1e3:.0f} ms/call)", flush=True)
            else:
                print(f"theta {theta:3d} {case:6s} {speed:5.2f} m/s cached", flush=True)
            cases[(theta, case)] = c
    allt = [c["timing"] for c in cases.values()]
    timing.update(energy_matrices_3modes_per_wind_s=float(np.mean([t["matrices_s"] for t in allt])),
                  dijkstra_per_call_2_sources_s=float(np.mean([t["dijkstra_call_s"] for t in allt])),
                  fields_computed_this_run=computed, run_s=time.time() - T0)
    print(f"timing: graph build {timing['graph_build_s']:.1f} s | cost matrices (3 modes) {timing['energy_matrices_3modes_per_wind_s']:.1f} s/wind | "
          f"Dijkstra {timing['dijkstra_per_call_2_sources_s']*1e3:.0f} ms/call (2 sources, pruned) | {computed} fields computed")

    # ---- csv + summary from the whole cache --------------------------------------------------
    rows = [r for k in sorted(cases) for r in cases[k]["rows"]]
    with open(HERE / "results" / "plan_results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    summary = build_summary(cases, dirs, G, args, rejected, ncand, ratio, timing)
    json.dump(summary, open(HERE / "results" / "plan_summary.json", "w"), indent=1)

    print("\nenergy saving % vs the shortest path (80 legs): FAIR = smoothed shortest vs energy-aware smoothed | raw graph paths")
    for tag in PLANNERS:
        print(f"-- planner {tag}")
        for k, e in summary["per_direction"].items():
            r, s = e[f"raw_{tag}"], e[f"fair_{tag}"]
            print(f"  {k:12s} {e['wind_speed']:5.2f} m/s n={e['n_legs']:3d}  FAIR mean {s['save_mean']:6.2f} med {s['save_median']:6.2f} p90 {s['save_p90']:6.2f} "
                  f"len {s['extra_len_mean']:+5.2f}% | raw {r['save_mean']:6.2f}/{r['save_median']:6.2f} len {r['extra_len_mean']:+5.2f}%")
    print("separation of effects (vs straight path at 10 m/s, smoothed; mean %):")
    for k, e in summary["per_direction"].items():
        print("  ", k, {a: round(b, 2) for a, b in e["effects"].items()})
    for case, w in summary["weighted"].items():
        print("weighted expected saving", case, {k: round(v, 2) for k, v in w["expected_save_pct"].items()})

    if not args.no_figs:
        best_pair = {th: (c["best"]["save"], c["best"]["pair"], c["best"]["speed"],
                          {k: np.array(v) for k, v in c["best"]["paths"].items()}, {k: np.array(v) for k, v in c["best"]["sm"].items()})
                     for (th, case), c in cases.items() if case == "edmo" and c["best"]}
        make_figures(atlas, env, rows, best_pair, dirs, summary)
    print("done in %.0f s" % (time.time() - T0))


# ============================================================================================
# figures
# ============================================================================================
def make_figures(atlas, env, rows, best_pair, dirs, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    COL = {"short": "#222222", "fix10": "#0072B2", "fix15": "#009E73", "best": "#D55E00"}
    NAME = {"short": "shortest", "fix10": "energy, 10 m/s", "fix15": "energy, 15 m/s", "best": "energy, best speed"}
    ground = atlas.ground()
    bld = np.nan_to_num(ground - env.terrain, nan=0.0) > 3.0                     # building footprint (top > 3 m above ground)
    X, Y = np.meshgrid(env.x, env.y, indexing="ij")

    for theta, (sv, pi, speed, paths, sm) in best_pair.items():
        wind = WindField.from_atlas(atlas, theta, speed)
        # --- map -------------------------------------------------------------------------
        k = np.clip(np.rint((env.terrain + 30.0 - env.z[0]) / env.h[2]).astype(int), 0, len(env.z) - 1)
        U30 = np.take_along_axis(wind.U, k[:, :, None, None].repeat(3, 3), axis=2)[:, :, 0, :]
        spd = np.linalg.norm(U30, axis=-1)
        fig, ax = plt.subplots(figsize=(9.5, 8))
        im = ax.pcolormesh(env.x, env.y, np.where(bld, np.nan, spd).T, cmap="viridis", shading="auto")
        ax.contourf(X, Y, bld.astype(float), levels=[0.5, 1.5], colors="none", hatches=["////"])
        ax.contour(X, Y, bld.astype(float), levels=[0.5], colors="k", linewidths=0.6)
        raw = paths
        for tag in ("fix10", "best"):                      # faint raw graph paths for reference
            ax.plot(raw[tag][:, 0], raw[tag][:, 1], color=COL[tag], lw=1.0, alpha=0.35, zorder=3)
        paths = sm                                         # headline = smoothed (shortest: geometric, energy: energy-aware)
        for tag in ("short", "fix10", "best"):
            P = paths[tag]
            if tag == "short":
                m = {t: metrics(P, wind, env, t) for t in ("fix10", "best")}
                lab = f"shortest: {m['fix10']['length_m']:.0f} m, {m['fix10']['energy_J']/1e3:.2f} kJ @10 m/s, {m['best']['energy_J']/1e3:.2f} kJ best-speed"
                ax.plot(P[:, 0], P[:, 1], "--", color="white", lw=3.2, zorder=4); ax.plot(P[:, 0], P[:, 1], "--", color=COL[tag], lw=1.8, label=lab, zorder=5)
            else:
                mp, ms = metrics(P, wind, env, tag), metrics(paths["short"], wind, env, tag)
                s = 100 * (1 - mp["energy_J"] / ms["energy_J"])
                lab = f"{NAME[tag]}: {mp['length_m']:.0f} m, {mp['energy_J']/1e3:.2f} kJ ({-s:+.1f}%)"
                ax.plot(P[:, 0], P[:, 1], color=COL[tag], lw=2.0, label=lab, zorder=6)
        s0, g0 = paths["short"][0], paths["short"][-1]
        ax.plot(*s0[:2], "o", ms=11, mfc="#FFD700", mec="k", zorder=8, label="start"); ax.plot(*g0[:2], "*", ms=17, mfc="#FFD700", mec="k", zorder=8, label="goal")
        ax.annotate("", xy=(-300 + 60 * -np.sin(np.radians(theta)), 380 + 60 * -np.cos(np.radians(theta))), xytext=(-300, 380),
                    arrowprops=dict(arrowstyle="-|>", color="k", lw=2))
        ax.text(-380, 400, f"wind from {theta} deg, {speed:.2f} m/s\n(speed plotted at 30 m above ground)", fontsize=9, va="bottom")
        ax.set_aspect("equal"); ax.set_xlabel("x east (m)"); ax.set_ylabel("y north (m)")
        ax.set_title(f"Pair {pi} (largest fair best-speed saving of {len(set(r['pair'] for r in rows))}): smoothed paths (faint = raw graph path)", fontsize=10)
        plt.colorbar(im, ax=ax, shrink=0.7, label="wind speed (m/s)")
        ax.legend(loc="lower right", fontsize=8, framealpha=0.92)
        fig.tight_layout(); fig.savefig(HERE / "fig" / f"plan_map_{theta}.png", dpi=130); plt.close(fig)

        # --- altitude profile -------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(10, 4.2))
        paths = sm
        P = paths["short"]
        sd = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
        dense = np.vstack([np.linspace(u, v, 20, endpoint=False) for u, v in zip(P[:-1], P[1:])] + [P[-1:]])
        dd = np.r_[0, np.cumsum(np.linalg.norm(np.diff(dense, axis=0), axis=1))]
        gi = np.stack([(dense[:, 0] - env.x[0]) / env.h[0], (dense[:, 1] - env.y[0]) / env.h[1]])
        gz = map_coordinates(np.nan_to_num(ground, nan=0.0), gi, order=0, mode="nearest")
        ax.fill_between(dd, 0, gz, color="0.75", step="mid", label="ground / buildings under the shortest path")
        for tag in ("short", "fix10", "best"):
            Q = paths[tag]; d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))]
            ax.plot(d, Q[:, 2], color=COL[tag], lw=2, ls="--" if tag == "short" else "-", label=NAME[tag])
        ax.set_ylim(0, max(60, paths["short"][:, 2].max() + 20)); ax.set_xlabel("distance along path (m)"); ax.set_ylabel("altitude z (m, atlas frame)")
        ax.set_title(f"Altitude profile, wind from {theta} deg at {speed:.2f} m/s (same pair as the map)"); ax.legend(fontsize=8, loc="upper right")
        fig.tight_layout(); fig.savefig(HERE / "fig" / f"plan_profile_{theta}.png", dpi=130); plt.close(fig)

    # --- savings distribution ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), sharey="row")
    rng = np.random.default_rng(0)
    for r_i, lab in enumerate(("raw", "sm")):
        for c_i, tag in enumerate(PLANNERS):
            ax = axes[r_i, c_i]
            for off, case, color in ((-0.2, "edmo", "#0072B2"), (0.2, "strong", "#E69F00")):
                data, pos = [], []
                for i, th in enumerate(dirs):
                    v = [r[f"{lab}_{tag}_save_pct"] for r in rows if r["theta"] == th and r["case"] == case]
                    if v:
                        data.append(v); pos.append(i + off)
                if not data:
                    continue
                bp = ax.boxplot(data, positions=pos, widths=0.32, patch_artist=True, showfliers=False)
                for b in bp["boxes"]:
                    b.set(facecolor=color, alpha=0.35, edgecolor=color)
                for m in bp["medians"]:
                    m.set(color="k")
                for p_, v in zip(pos, data):
                    ax.plot(p_ + rng.uniform(-0.07, 0.07, len(v)), v, ".", ms=2.5, color=color, alpha=0.6)
            ax.axhline(0, color="k", lw=0.8)
            ax.set_xticks(range(len(dirs))); ax.set_xticklabels([f"{d}" for d in dirs])
            ax.set_xlabel("wind from (deg)"); ax.grid(axis="y", alpha=0.3)
            ax.set_title(f"{'smoothed (fair)' if lab == 'sm' else 'raw graph path'} - planner: {NAME[tag]}", fontsize=10)
            if c_i == 0:
                ax.set_ylabel("energy saving vs shortest path (%)")
    axes[0, 0].legend(handles=[Line2D([], [], color="#0072B2", lw=6, alpha=0.5, label="EDMO mean speed"),
                               Line2D([], [], color="#E69F00", lw=6, alpha=0.5, label=f"{summary['strong_speed']:g} m/s")], fontsize=8, loc="best")
    fig.suptitle(f"Energy saving of the wind-aware planner over the shortest path ({summary['n_pairs']} pairs x 2 legs per box)")
    fig.tight_layout(); fig.savefig(HERE / "fig" / "plan_savings.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
