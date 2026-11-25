"""Planner tables and one-paragraph reading for slides.tex, from DLR/planner/results/plan_summary.json.
build(plan) -> (tab_plan: one row per direction, tab_levers: effect decomposition, sentence)."""


def _label(per):
    return "fair" if any(k.startswith("fair_") for v in per.values() for k in v) else "raw"


def build(plan):
    per = plan.get("per_direction", {})
    if not per:
        miss = "\\textit{planner results not available yet}\n"
        return miss, miss, "--"
    lab = _label(per)
    by_dir = {}
    for key, v in per.items():
        th, kind = key.split("_")
        by_dir.setdefault(int(th), {})[kind] = v
    t = ["\\footnotesize\\begin{tabular}{@{}rr rrr rrr@{}}", "\\toprule",
         " & & \\multicolumn{3}{c}{typical wind (EDMO mean)} & \\multicolumn{3}{c}{8\\,m/s} \\\\",
         "\\cmidrule(lr){3-5}\\cmidrule(l){6-8}",
         "from & reports & $U_{10}$ & 10\\,m/s & best speed & 10\\,m/s & best speed & p90 \\\\",
         " & (\\%) & (m/s) & \\multicolumn{2}{c}{energy saved (\\%)} & \\multicolumn{2}{c}{energy saved (\\%)} & best \\\\", "\\midrule"]
    edmo = {0: 3.6, 30: 7.1, 60: 8.3, 90: 6.1, 120: 3.1, 150: 6.8, 180: 8.2, 210: 15.8, 240: 22.9, 270: 10.3, 300: 4.1, 330: 3.7}
    for th in sorted(by_dir):
        e, s = by_dir[th].get("edmo"), by_dir[th].get("strong")
        if not e or not s:
            continue
        f = lambda v, k: f"{v[f'{lab}_{k}']['save_mean']:.1f}"
        t.append(f"{th}$^\\circ$ & {edmo.get(th, 0):.1f} & {e['wind_speed']:.1f} & {f(e, 'fix10')} & {f(e, 'best')} & "
                 f"{f(s, 'fix10')} & {f(s, 'best')} & {s[f'{lab}_best']['save_p90']:.0f} \\\\")
    t += ["\\bottomrule", "\\end{tabular}"]
    w = plan.get("weighted", {})
    ew, es = w.get("edmo", {}).get("expected_save_pct", {}), w.get("strong", {}).get("expected_save_pct", {})
    lv = ["\\small\\begin{tabular}{@{}lrr@{}}", "\\toprule",
          "Lever (climate-weighted, \\% of energy) & typical & 8\\,m/s \\\\", "\\midrule"]
    if "effect_wind_speed_adaptation_pct" in ew:
        rows = [("efficient speed, not wind: 18.3 vs.\\ 10\\,m/s, still air", "effect_speed_choice_still_air_pct", 0),
                ("wind-adaptive speed: best per segment vs.\\ 18.3\\,m/s", "effect_wind_speed_adaptation_pct", 1),
                ("wind-aware route, both at best speed", "effect_route_given_best_speed_pct", 1),
                ("wind-aware route, fixed 10\\,m/s", "effect_route_only_pct", 1),
                ("(not a saving: zero on a round trip) energy gap A$\\to$B vs.\\ B$\\to$A", "effect_direction_asym_best_pct", 0)]
        for n, k, nd in rows:
            lv.append(f"{n} & {ew[k]:.{nd}f} & {es[k]:.{nd}f} \\\\")
    lv += ["\\bottomrule", "\\end{tabular}"]
    ndir = len(plan.get("directions_available", []))
    s = (f"{plan.get('n_pairs', '--')} random start--goal pairs ($\\geq$\\,250\\,m apart, both legs), {ndir} wind directions. "
         "Shortest path (smoothed) vs.\\ energy path (smoothed only where it lowers energy), both re-evaluated in the same wind. ")
    if "effect_route_given_best_speed_pct" in ew:
        s += (f"\\textbf{{Reading:}} the wind-aware \\emph{{route}} saves about {ew['effect_route_given_best_speed_pct']:.0f}\\,\\% at typical "
              f"winds and {es['effect_route_given_best_speed_pct']:.0f}\\,\\% at 8\\,m/s (some trips 15--25\\,\\%). Flying at the right \\emph{{speed}} "
              f"is worth more ({ew['effect_wind_speed_adaptation_pct']:.0f}--{es['effect_wind_speed_adaptation_pct']:.0f}\\,\\% from adapting it to the wind, "
              f"{ew['effect_speed_choice_still_air_pct']:.0f}\\,\\% from leaving 10\\,m/s), and the \\emph{{direction}} of a leg changes its cost by "
              f"{ew['effect_direction_asym_best_pct']:.0f}--{es['effect_direction_asym_best_pct']:.0f}\\,\\% (A to B vs.\\ back). So the cost function pays off mostly in "
              "speed scheduling and leg ordering; the route matters on windy days and in the shelter of buildings.")
    return "\n".join(t) + "\n", "\n".join(lv) + "\n", s
