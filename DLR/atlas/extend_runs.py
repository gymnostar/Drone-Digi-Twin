"""--list: names of atlas runs whose probe drift exceeds the criterion (for extend_runs.sh).
--compare <run>: |dU| / Uref between the shorter (prev_<run>.npz) and the extended field -> results/extension.json."""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
CRIT = 0.005          # 0.5 % of Uref over the last 50 iterations, probes >= 10 m above ground

if sys.argv[1] == "--list":
    val = json.loads((HERE / "results/validation.json").read_text())
    names = [k for k, c in val["convergence"].items()
             if c["drift_last50_over_Uref_probes_ge10m"] > CRIT and "(" not in k]
    print(" ".join(names))
else:
    import validate_atlas as V
    from wind_atlas import WindAtlas
    run = sys.argv[2]
    a, b = np.load(HERE / f"grid.noindex/prev_{run}.npz"), np.load(HERE / f"grid.noindex/{run}.npz")
    A = WindAtlas(); V.ATLAS = A
    mask, t = V.campus_mask(A), A.terrain()
    I, J = np.indices(t.shape)
    d = np.linalg.norm(b["U"] - a["U"], axis=-1) / float(a["Uref"])
    out = {"run": run}
    for h in V.HEIGHTS:
        k = np.clip(np.round((t + h - A.z[0]) / (A.z[1] - A.z[0])).astype(int), 0, len(A.z) - 1)
        e = d[I, J, k][mask]; e = e[np.isfinite(e)]
        out[f"{h:g}m"] = {"rms": float(np.sqrt((e**2).mean())), "p95": float(np.percentile(e, 95)), "max": float(e.max())}
    f = HERE / "results/extension.json"
    rec = json.loads(f.read_text()) if f.exists() else {}
    rec[run] = out
    f.write_text(json.dumps(rec, indent=1))
    print(json.dumps(out, indent=1))
