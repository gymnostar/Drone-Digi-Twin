"""Map building heights + ground around the DLR weather station, in CASE coordinates (m).
case = scene_local - shift, shift = (-335, 279)  (see windCase/system/caseParams)."""
import numpy as np, trimesh, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, sys
C = "DLR/windCase/constant/triSurface/"
bld = trimesh.load(C + "buildings.stl"); ter = trimesh.load(C + "terrain.stl")
st = np.array([-213 + 335, 485 - 279])            # station in case coords
x = np.arange(-340, 380, 2.0); y = np.arange(-310, 400, 2.0)
X, Y = np.meshgrid(x, y); P = np.c_[X.ravel(), Y.ravel()]
def top(mesh):
    o = np.c_[P, np.full(len(P), 500.0)]; d = np.tile([0, 0, -1.0], (len(P), 1))
    loc, idx, _ = mesh.ray.intersects_location(o, d, multiple_hits=False)
    h = np.full(len(P), np.nan); h[idx] = loc[:, 2]; return h.reshape(X.shape)
hb, hg = top(bld), top(ter)
np.savez("DLR/meshStudy/site_map.npz", x=x, y=y, hb=hb, hg=hg, st=st)
fig, ax = plt.subplots(figsize=(9, 9))
ax.pcolormesh(x, y, hb - hg, cmap="Greys", shading="auto")
ax.add_patch(plt.Rectangle(st - 50, 100, 100, fill=False, ec="r", ls="--"))
ax.plot(*st, "r*", ms=15); ax.set_aspect("equal"); ax.grid(alpha=.3)
ax.set_title("buildings (height above ground), case coords; red = station ±50 m; wind from left")
fig.savefig("DLR/meshStudy/site_map.png", dpi=110, bbox_inches="tight")
print("station case xy", st, "ground there", hg[np.abs(y-st[1]).argmin(), np.abs(x-st[0]).argmin()])
