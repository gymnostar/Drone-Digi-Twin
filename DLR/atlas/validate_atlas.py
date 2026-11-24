#!/usr/bin/env python3
"""Validate the DLR wind atlas: is each shortcut in the lookup good enough, and what does it cost?

Run from the repo root:  .venv/bin/python DLR/atlas/validate_atlas.py
Writes DLR/atlas/results/validation.json, validation.md and figures in DLR/atlas/results/.

Checks (all errors are vector errors |U_pred - U_cfd|, divided by the reference speed Uref, so
"0.10" = 10 % of the 10 m wind; evaluated in the AIR over the campus, at fixed heights above ground):
1. Convergence of every run: probe drift over the last 50 iterations and final residuals.
2. Speed scaling (Reynolds independence): the 10 m/s run at 240 deg against 2 x the 5 m/s run.
3. Direction blending, leave-one-out: hide a solved direction, predict it from its neighbours
   60 deg apart, compare with the CFD. Upper bound on the error of the real atlas (30 deg apart).
   Compared with two cheaper models: nearest neighbour direction, and an empty-field log-law
   (the inflow profile with no buildings, i.e. what you get without any CFD).
4. Lookup cost: time per point query and per whole-field call.
"""
import json, re, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from wind_atlas import WindAtlas  # noqa: E402

CASES = HERE / "cases.noindex"
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
HEIGHTS = (10.0, 20.0, 40.0, 60.0)       # m above bare ground
KAPPA, Z0, ZREF = 0.41, 0.1, 10.0       # log-law constants used by the CFD inlet


# ---------- 1. convergence --------------------------------------------------------------------
PROBE_AGL = np.array([float(l.split(",")[4]) for l in open(HERE.parent / "meshStudy/probes_index.csv").read().split("\n")[1:] if l])


def convergence(case: Path):
    """Probe drift over the last 50 iterations, as |change of |U|| / Uref (local speeds can be ~0 in
    recirculation pockets, so dividing by the local speed would exaggerate tiny absolute changes)."""
    f = case / "postProcessing/probes/0/U"
    if not f.exists():
        return None
    uref = 5.0
    if (case / "RUN_INFO").exists():
        uref = float(re.search(r"Uref_10m\s+([0-9.]+)", (case / "RUN_INFO").read_text()).group(1))
    rows = [[float(v) for v in re.findall(r"[-+0-9.eE]+", l)] for l in open(f) if not l.startswith("#")]
    a = np.array(rows)
    it, m = a[:, 0], np.linalg.norm(a[:, 1:].reshape(len(a), -1, 3), axis=2)
    j = np.searchsorted(it, it[-1] - 50)
    d = np.abs(m[-1] - m[j]) / uref
    res = {}
    log = case / "log.simpleFoam"
    if log.exists():
        txt = log.read_text()
        for fld in ("Ux", "Uy", "p", "k", "epsilon"):
            r = re.findall(rf"Solving for {fld}, Initial residual = ([0-9.eE+-]+)", txt)
            if r:
                res[fld] = float(r[-1])
    return {"iterations": int(it[-1]),
            "drift_last50_over_Uref_all_probes": float(d.max()),
            "drift_last50_over_Uref_probes_ge10m": float(d[PROBE_AGL >= 10].max()),
            "final_initial_residuals": res}


# ---------- helpers ---------------------------------------------------------------------------
def campus_mask(atlas, radius_m=60.0):
    """Columns within radius_m of a building: where the CFD differs from the empty-field profile."""
    from scipy.ndimage import binary_dilation
    bld = (atlas.ground() - atlas.terrain()) > 3.0
    r = int(round(radius_m / (atlas.x[1] - atlas.x[0])))
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    return binary_dilation(bld, structure=(xx**2 + yy**2) <= r * r)


def at_heights(atlas, F):
    """Values of F (nx, ny, nz, ...) at HEIGHTS above bare terrain -> dict h -> (nx, ny, ...)."""
    t = atlas.terrain()
    dz, z0 = atlas.z[1] - atlas.z[0], atlas.z[0]
    I, J = np.indices(t.shape)
    out = {}
    for h in HEIGHTS:
        k = np.clip(np.round((t + h - z0) / dz).astype(int), 0, len(atlas.z) - 1)
        out[h] = F[I, J, k]
    return out


def loglaw_field(atlas, theta, speed):
    """Empty-field model: the CFD's inflow log-law profile everywhere, no buildings."""
    t = np.radians(theta)
    fdir = np.array([-np.sin(t), -np.cos(t), 0.0])
    ustar = KAPPA * speed / np.log((ZREF + Z0) / Z0)
    h = np.maximum(atlas.z[None, None, :] - atlas.terrain()[:, :, None], 0.0)
    mag = ustar / KAPPA * np.log((h + Z0) / Z0)
    U = mag[..., None] * fdir
    U[~atlas.valid] = np.nan
    return U


def errors(pred, true, mask, uref):
    """Normalised vector-error statistics per height, over mask columns where both are defined."""
    P, T = at_heights(ATLAS, pred), at_heights(ATLAS, true)
    out = {}
    for h in HEIGHTS:
        e = np.linalg.norm(P[h] - T[h], axis=-1) / uref
        ok = mask & np.isfinite(e)
        e = e[ok]
        out[f"{h:g}m"] = {"rmse": float(np.sqrt(np.mean(e**2))), "mean": float(e.mean()),
                          "p95": float(np.percentile(e, 95)), "n": int(e.size)}
    return out


if __name__ == "__main__":
    ATLAS = WindAtlas()
    dirs = [int(d) for d in ATLAS.directions]
    report = {"directions": dirs, "extra_runs": [f"{int(t)}@{u:g}" for t, u in ATLAS.extra]}
    mask = campus_mask(ATLAS)
    report["campus_columns"] = int(mask.sum())

    # 1. convergence
    conv = {}
    for d in sorted(CASES.glob("d*_u*")):
        c = convergence(d) if (d / "DONE").exists() else None
        if c:
            conv[d.name] = c
    c270 = convergence(HERE.parent / "meshStudy/coarse")
    if c270:
        conv["d270_u05 (mesh study coarse)"] = c270
    report["convergence"] = conv

    # 2. Reynolds independence
    if (240.0, 10.0) in ATLAS.extra and 240.0 in ATLAS.runs:
        U10 = ATLAS.extra[(240.0, 10.0)]["U"]           # already divided by 10
        U5 = ATLAS.runs[240.0]["U"]                       # divided by 5
        report["reynolds_240"] = errors(U5, U10, mask, 1.0)   # fields are already normalised
        k10, k5 = ATLAS.extra[(240.0, 10.0)]["k"], ATLAS.runs[240.0]["k"]
        kk = np.isfinite(k10) & np.isfinite(k5) & (k10 > 1e-6)
        report["reynolds_240_k_rel_median"] = float(np.median(np.abs(k5[kk] - k10[kk]) / k10[kk]))

    # 3. leave-one-out on direction blending (needs theta-30 and theta+30 solved). 255 is a
    #    validation-only run: always left out here so every test blends the 30-deg-grid neighbours.
    loo = {}
    for th in [d for d in dirs if d % 30 == 0]:
        a, b = (th - 30) % 360, (th + 30) % 360
        if a in dirs and b in dirs:
            sub = WindAtlas(exclude=(th, 255))
            true = ATLAS.field(th, 5.0)
            blend = sub.field(th, 5.0)
            loo[th] = {
                "blend_60deg": errors(blend, true, mask, 5.0),
                "nearest_60deg": {h: {"rmse": 0.5 * (errors(sub.field(a, 5.0), true, mask, 5.0)[h]["rmse"]
                                                     + errors(sub.field(b, 5.0), true, mask, 5.0)[h]["rmse"])}
                                  for h in [f"{x:g}m" for x in HEIGHTS]},
                "loglaw_no_buildings": errors(loglaw_field(ATLAS, th, 5.0), true, mask, 5.0),
            }
    report["leave_one_out"] = loo

    # 3b. the REAL atlas spacing: predict the extra 255-deg run from 240 and 270 (30 deg apart)
    if 255 in dirs and 240 in dirs and 270 in dirs:
        sub = WindAtlas(exclude=(255,))
        true = ATLAS.field(255, 5.0)
        report["blend_30deg_255"] = {
            "blend_30deg": errors(sub.field(255, 5.0), true, mask, 5.0),
            "nearest_30deg": {h: {"rmse": 0.5 * (errors(sub.field(240, 5.0), true, mask, 5.0)[h]["rmse"]
                                                 + errors(sub.field(270, 5.0), true, mask, 5.0)[h]["rmse"])}
                              for h in [f"{x:g}m" for x in HEIGHTS]},
            "loglaw_no_buildings": errors(loglaw_field(ATLAS, 255, 5.0), true, mask, 5.0),
        }

    # 4. lookup cost
    pts = np.c_[np.random.default_rng(0).uniform(-300, 350, 1000), np.random.default_rng(1).uniform(-300, 380, 1000),
                np.full(1000, 50.0)]
    t0 = time.perf_counter(); ATLAS.field(250, 6.0); t1 = time.perf_counter()
    ATLAS.query(pts, 250, 6.0); t2 = time.perf_counter()
    from scipy.interpolate import RegularGridInterpolator
    F = ATLAS.field(250, 6.0)
    it = RegularGridInterpolator((ATLAS.x, ATLAS.y, ATLAS.z), F, bounds_error=False, fill_value=np.nan)
    t3 = time.perf_counter(); it(pts); t4 = time.perf_counter(); it(pts[:1]); t5 = time.perf_counter()
    report["timing"] = {"field_blend_s": t1 - t0, "query_1000pts_incl_blend_s": t2 - t1,
                        "query_1000pts_prebuilt_s": t4 - t3, "query_1pt_prebuilt_ms": (t5 - t4) * 1e3}

    (OUT / "validation.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1)[:6000])
