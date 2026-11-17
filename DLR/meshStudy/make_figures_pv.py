"""pvpython helper for make_figures.py: cut ONE mesh-study case and save the cut cells as .npz.

    pvpython make_figures_pv.py <case dir or .foam file> <out.npz> [--fields]

Two cuts through the internal mesh (case coordinates, m):
    h : horizontal plane z = 24  (about 7 m above ground), kept for x 20..265, y 115..290
    v : vertical plane   y = 206 (through the weather station), kept for x 20..265, z 10..80
Each cut cell becomes one polygon (the slice is NOT triangulated, so polygon = cell outline).
With --fields the wind speed |U| of each cut cell is saved too (latest time only).
Holes in the cut = solid (buildings, ground). Nothing is written inside the case.
make_figures.py draws the pictures from the .npz with matplotlib.
"""
import os, sys
import numpy as np
from paraview.simple import OpenFOAMReader, Slice, MergeBlocks, Delete, servermanager
from vtkmodules.util.numpy_support import vtk_to_numpy

CUTS = {   # name: (origin, normal, the 2 axes kept for drawing, (a0 min, a0 max, a1 min, a1 max))
    "h": ((0.0, 0.0, 24.0), (0, 0, 1), (0, 1), (20.0, 265.0, 115.0, 290.0)),
    "v": ((0.0, 206.0, 0.0), (0, 1, 0), (0, 2), (20.0, 265.0, 10.0, 80.0)),
}
MARGIN = 20.0   # m kept around the drawn window so cells cut by the frame are complete

src, out = sys.argv[1], sys.argv[2]
fields = "--fields" in sys.argv[3:]
foam = src if src.endswith(".foam") else os.path.join(src, "case.foam")

r = OpenFOAMReader(FileName=foam)
r.CaseType = "Reconstructed Case"
r.MeshRegions = ["internalMesh"]
r.CellArrays = ["U"] if fields else []
r.PointArrays = []
r.Createcelltopointfiltereddata = 0     # no point interpolation: less memory, cell values only
r.SkipZeroTime = 1
r.UpdatePipelineInformation()
times = list(r.TimestepValues) if hasattr(r.TimestepValues, "__len__") else [r.TimestepValues]
t = max(times) if times else 0.0
if fields and not times:
    sys.exit("no solved time step found in " + foam)
r.UpdatePipeline(t)
ncells = r.GetDataInformation().GetNumberOfCells()
print(f"{foam}: {ncells} cells, time {t if times else 'mesh only'}", flush=True)

res = {"ncells": np.int64(ncells), "time": np.float64(t), "has_fields": np.bool_(fields)}
for name, (origin, normal, axes, lim) in CUTS.items():
    s = Slice(Input=r)
    s.SliceType = "Plane"
    s.SliceType.Origin = list(origin)
    s.SliceType.Normal = list(normal)
    s.Triangulatetheslice = 0
    m = MergeBlocks(Input=s)
    m.UpdatePipeline(t)
    d = servermanager.Fetch(m)
    pts = vtk_to_numpy(d.GetPoints().GetData()).astype(np.float64)
    ca = d.GetCells()
    off = vtk_to_numpy(ca.GetOffsetsArray()).astype(np.int64)
    conn = vtk_to_numpy(ca.GetConnectivityArray()).astype(np.int64)
    n = np.diff(off)
    ok = n >= 3
    cen = np.add.reduceat(pts[conn], off[:-1], axis=0) / np.maximum(n, 1)[:, None]
    a0, a1 = axes
    keep = ok & (cen[:, a0] > lim[0] - MARGIN) & (cen[:, a0] < lim[1] + MARGIN) \
              & (cen[:, a1] > lim[2] - MARGIN) & (cen[:, a1] < lim[3] + MARGIN)
    kept = np.flatnonzero(keep)
    # flatten the kept polygons: vertex coordinates (2D) + offsets
    sel = np.concatenate([conn[off[i]:off[i + 1]] for i in kept]) if len(kept) else np.zeros(0, np.int64)
    res[name + "_xy"] = pts[sel][:, [a0, a1]].astype(np.float32)
    res[name + "_off"] = np.concatenate([[0], np.cumsum(n[kept])]).astype(np.int64)
    if fields:
        U = vtk_to_numpy(d.GetCellData().GetArray("U"))
        res[name + "_speed"] = np.linalg.norm(U, axis=1)[kept].astype(np.float32)
    print(f"  cut {name}: {len(kept)} cells kept of {len(n)}", flush=True)
    Delete(m); Delete(s)
    del d
Delete(r)
np.savez_compressed(out, **res)
print("wrote", out, flush=True)
