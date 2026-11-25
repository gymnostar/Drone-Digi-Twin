"""How much does a 350-iteration atlas field change if the solver keeps going?
usage: iteration_check.py <run> <new iterations> <old iterations>   (called by iteration_check.sh)"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from wind_atlas import WindAtlas  # noqa: E402
import validate_atlas as V  # noqa: E402

run, new, old = sys.argv[1:4]
a = np.load(HERE / f"grid.noindex/{run}.npz")
b = np.load(HERE / f"grid.noindex/iter_{run}.npz")
uref = float(a["Uref"])
A = WindAtlas()
V.ATLAS = A
dU = np.linalg.norm(b["U"] - a["U"], axis=-1) / uref
out = {"run": run, "iterations_old": int(old), "iterations_new": int(new), "uref": uref}
mask = V.campus_mask(A)
t = A.terrain()
for h in V.HEIGHTS:
    k = np.clip(np.round((t + h - A.z[0]) / (A.z[1] - A.z[0])).astype(int), 0, len(A.z) - 1)
    I, J = np.indices(t.shape)
    e = dU[I, J, k][mask]
    e = e[np.isfinite(e)]
    out[f"{h:g}m"] = {"rms": float(np.sqrt((e**2).mean())), "p95": float(np.percentile(e, 95)), "max": float(e.max())}
e = dU[np.isfinite(dU)]
out["whole_domain_rms"] = float(np.sqrt((e**2).mean()))
(HERE / "results/iteration_check.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
