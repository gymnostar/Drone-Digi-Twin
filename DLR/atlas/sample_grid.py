"""pvpython: resample one solved atlas run onto the shared regular grid and save it as .npz.

    pvpython sample_grid.py <case dir> <out.npz>

The CFD mesh is unstructured (4 m cells on campus, 2 m at walls, 16 m far away). Planners, the RF
side and AirSim all want "value at (x, y, z)" in O(1), so every run is sampled ONCE onto the same
regular grid (case coordinates, m):
    x -400..440, y -360..440 every 5 m (horizontal), z 10..130 every 2.5 m (absolute; ground ~13-25 m)
Points inside buildings or below the terrain are outside the CFD mesh: valid = 0 there, U = NaN.
Arrays are stored [ix, iy, iz] (x fastest in VTK, transposed here so indexing reads naturally).
"""
import sys
import numpy as np
from paraview.simple import OpenFOAMReader, ResampleToImage, servermanager
from vtkmodules.util.numpy_support import vtk_to_numpy

X0, X1, DX = -400.0, 440.0, 5.0
Y0, Y1, DY = -360.0, 440.0, 5.0
Z0, Z1, DZ = 10.0, 130.0, 2.5
nx, ny, nz = int(round((X1 - X0) / DX)) + 1, int(round((Y1 - Y0) / DY)) + 1, int(round((Z1 - Z0) / DZ)) + 1

case, out = sys.argv[1], sys.argv[2]
r = OpenFOAMReader(FileName=case + "/case.foam")
r.CaseType = "Reconstructed Case"
r.MeshRegions = ["internalMesh"]
r.CellArrays = ["U", "k", "p", "nut"]
r.Createcelltopointfiltereddata = 1     # smooth point values -> trilinear-quality sampling
r.SkipZeroTime = 1
r.UpdatePipelineInformation()
t = max(r.TimestepValues) if hasattr(r.TimestepValues, "__len__") else r.TimestepValues
rs = ResampleToImage(Input=r)
rs.UseInputBounds = 0
rs.SamplingBounds = [X0, X1, Y0, Y1, Z0, Z1]
rs.SamplingDimensions = [nx, ny, nz]
rs.UpdatePipeline(t)
img = servermanager.Fetch(rs)
pd = img.GetPointData()

def arr(name, ncomp=1):
    a = vtk_to_numpy(pd.GetArray(name)).astype(np.float32)
    a = a.reshape(nz, ny, nx, ncomp) if ncomp > 1 else a.reshape(nz, ny, nx)
    return np.ascontiguousarray(np.moveaxis(a, (0, 1, 2), (2, 1, 0)))   # -> [ix, iy, iz, ...]

valid = arr("vtkValidPointMask") > 0.5
U, k, p, nut = arr("U", 3), arr("k"), arr("p"), arr("nut")
for a in (U, k, p, nut):
    a[~valid] = np.nan
info = open(case + "/RUN_INFO").read().split("\n")
meta = dict(l.split(" ", 1) for l in info if l.strip())
np.savez_compressed(
    out, x=np.linspace(X0, X1, nx), y=np.linspace(Y0, Y1, ny), z=np.linspace(Z0, Z1, nz),
    U=U, k=k, p=p, nut=nut, valid=valid, time=np.float64(t),
    theta_from_deg=np.float64(meta["theta_from_deg"]), Uref=np.float64(meta["Uref_10m"]))
print(f"{case}: t={t}, grid {nx}x{ny}x{nz}, valid {valid.mean()*100:.1f} %, "
      f"|U| max {np.nanmax(np.linalg.norm(U, axis=-1)):.2f} m/s -> {out}")
