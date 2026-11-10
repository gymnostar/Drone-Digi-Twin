#!/usr/bin/env python3
"""Extend the DLR terrain with flat ground so make_case.py can fit a large domain.

Why: make_case.py fits the largest circle (centred on the terrain bounding box) that has
ground everywhere under it. Arthur's terrain is a diagonal strip, so that circle was only
~123 m. Here we rebuild the terrain as ONE continuous square sheet (no gaps for air to leak
through), centred on the same point so local coordinates stay the same:

  * inside the original strip: real heights (linear interpolation of the original terrain)
  * outside: flat ground at the median original height
  * a smooth blend over BLEND_M metres so there is no step at the old edge

The flat part is invented, not measured. It is roughly right here (open airfield), but a
real DGM1 download would be the proper replacement.

Usage (from DLR/windCase):
  python3 ../scripts/extend_terrain.py terrain_mm.stl terrain_ext_mm.stl --side 1800 --dx 4
Output stays in mm, like the input, so make_case.py is called with --scale 0.001 as before.
"""
import argparse
import numpy as np
import trimesh
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("src")
ap.add_argument("dst")
ap.add_argument("--side", type=float, default=1800.0, help="square side [m]")
ap.add_argument("--dx", type=float, default=4.0, help="grid spacing [m]")
ap.add_argument("--blend", type=float, default=60.0, help="blend distance to flat [m]")
ap.add_argument("--scale", type=float, default=0.001, help="input unit -> m (mm: 0.001)")
A = ap.parse_args()

m = trimesh.load(A.src, force="mesh")
v = np.unique(m.vertices * A.scale, axis=0)              # metres
lo, hi = v.min(0), v.max(0)
cx, cy = (lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2        # same centre make_case.py will use
z_flat = float(np.median(v[:, 2]))

# heights inside the strip: interpolate the original triangulation's vertices
interp = LinearNDInterpolator(v[:, :2], v[:, 2])
tree = cKDTree(v[:, :2])

n = int(np.ceil(A.side / A.dx)) + 1
xs = cx + np.linspace(-A.side / 2, A.side / 2, n)
ys = cy + np.linspace(-A.side / 2, A.side / 2, n)
X, Y = np.meshgrid(xs, ys)
P = np.column_stack([X.ravel(), Y.ravel()])

dist, idx = tree.query(P)                                 # distance to nearest real terrain point
z_real = interp(P)
z_edge = np.where(np.isnan(z_real), v[idx, 2], z_real)    # outside the hull: nearest real height
w = np.clip((dist - A.dx) / A.blend, 0.0, 1.0)            # 0 on the real terrain, 1 far away
Z = (1 - w) * z_edge + w * z_flat

# triangulate the regular grid: two triangles per square, normals pointing up (+z)
verts = np.column_stack([P, Z]) / A.scale                 # back to input units (mm)
i = np.arange(n * n).reshape(n, n)
a, b, c, d = i[:-1, :-1].ravel(), i[:-1, 1:].ravel(), i[1:, 1:].ravel(), i[1:, :-1].ravel()
faces = np.vstack([np.column_stack([a, b, c]), np.column_stack([a, c, d])])
out = trimesh.Trimesh(verts, faces, process=False)
out.export(A.dst)

print(f"centre (m): {cx:.1f} {cy:.1f}   flat height: {z_flat:.2f} m")
print(f"grid: {n} x {n} points, {len(faces):,} triangles, side {A.side:.0f} m, dx {A.dx} m")
print(f"real terrain kept where within ~{A.dx:.0f} m of original points; blend {A.blend:.0f} m")
