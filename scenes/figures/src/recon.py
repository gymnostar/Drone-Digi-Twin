"""Reconstruct the post-snap sibling-octet merge (agglomerate.C) on the exported
422,894-cell common N mesh, so the merged mesh can be drawn. Verified against
the recorded per-level octet counts. Writes recon.npz next to this script."""
import sys, numpy as np, vtk, pyvista as pv
from vtk.util.numpy_support import vtk_to_numpy as v2n
from collections import defaultdict

from style import SCENES, FIGURES
R = str(SCENES / 'review-export') + '/'
OUT = sys.argv[1]
m = pv.read(R + 'joint-baseline.vtu')
pts = np.asarray(m.points, dtype=float)
types = np.asarray(m.celltypes)
conn = v2n(m.GetCells().GetConnectivityArray())
offs = v2n(m.GetCells().GetOffsetsArray())
nC = m.n_cells

HEX_FACES = [(0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)]
WEDGE_FACES = [(0,1,2),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5)]
cell_faces = [None]*nC
pf = m.GetPolyhedronFaces(); pfl = m.GetPolyhedronFaceLocations()
pf_conn = v2n(pf.GetConnectivityArray()); pf_off = v2n(pf.GetOffsetsArray())
fl_conn = v2n(pfl.GetConnectivityArray()); fl_off = v2n(pfl.GetOffsetsArray())
for c in range(nC):
    ids = conn[offs[c]:offs[c+1]]
    t = types[c]
    if t == 12: cell_faces[c] = [tuple(ids[list(f)]) for f in HEX_FACES]
    elif t == 13: cell_faces[c] = [tuple(ids[list(f)]) for f in WEDGE_FACES]
    else:
        fs = fl_conn[fl_off[c]:fl_off[c+1]]
        cell_faces[c] = [tuple(pf_conn[pf_off[f]:pf_off[f+1]]) for f in fs]
face_cells = defaultdict(list)
for c in range(nC):
    for f in cell_faces[c]: face_cells[tuple(sorted(f))].append(c)
bad = sum(1 for v in face_cells.values() if len(v) > 2)
patch_point = np.zeros(len(pts), bool)
nb = 0
for k, v in face_cells.items():
    if len(v) == 1:
        patch_point[list(k)] = True; nb += 1
print('faces', len(face_cells), 'boundary faces', nb, 'faces>2 cells', bad)

origin = np.array([-600., -1200., 0.]); tol = float(sys.argv[2]) if len(sys.argv) > 2 else 0.18
pmin = np.array([-26., -292., 1.5]); pmax = np.array([26., -268., 66.])
level = np.full(nC, -1)
for c in np.where(types == 12)[0]:
    v = pts[conn[offs[c]:offs[c+1]]]
    w = (v.max(0) - v.min(0)).max()
    L = int(round(np.log2(20.0 / w)))
    level[c] = L
# true cell levels for non-hex cells are approximated from bbox (used only for 2:1 check)
for c in np.where(types != 12)[0]:
    v = pts[conn[offs[c]:offs[c+1]]]
    w = (v.max(0) - v.min(0)).max()
    level[c] = int(round(np.log2(20.0 / max(w, 1e-9))))
groups = {}
for c in np.where(types == 12)[0]:
    L = level[c]
    if L < 0 or L > 3: continue
    ids = conn[offs[c]:offs[c+1]]; v = pts[ids]
    width = 20.0 * 2.0**-L
    cen = v.mean(0)  # approximates OpenFOAM centroid for near-cubic cells
    g = np.floor((cen - origin) / width).astype(int)
    if (g < 0).any(): continue
    low = origin + width * g
    if patch_point[ids].any(): continue
    t = (v - low) / width
    corners = 0; ok = True
    for row in t:
        corner = 0
        for a in range(3):
            if abs(row[a] - 1) < tol: corner |= 1 << a
            elif abs(row[a]) >= tol: ok = False; break
        if not ok: break
        corners |= 1 << corner
    if not ok or corners != 255: continue
    high = low + width
    if (low < pmax).all() and (high > pmin).all(): continue
    key = (L, g[0]//2, g[1]//2, g[2]//2)
    slot = g[0] % 2 + 2*(g[1] % 2) + 4*(g[2] % 2)
    o = groups.setdefault(key, [[-1]*8, True])
    if o[0][slot] != -1: o[1] = False
    o[0][slot] = c
master = np.full(nC, -1); group_of = np.full(nC, -1)
count = [0]*4; parents = []
for key in sorted(groups, reverse=True):
    children, valid = groups[key]
    if not valid or -1 in children: continue
    L = key[0]; good = True
    for slot, ch in enumerate(children):
        conns = 0
        for f in cell_faces[ch]:
            cs = face_cells[tuple(sorted(f))]
            if len(cs) != 2 or len(f) != 4: good = False; break
            other = cs[0] if cs[1] == ch else cs[1]
            if other in children:
                bit = slot ^ children.index(other)
                if bit not in (1, 2, 4) or (conns & bit): good = False; break
                conns |= bit
            else:
                ol = level[other] - (master[other] >= 0)
                if abs(ol - (L - 1)) > 1: good = False; break
        good = good and conns == 7
        if not good: break
    if not good: continue
    count[L] += 1
    mst = min(children)
    for ch in children: master[ch] = mst; group_of[ch] = len(parents)
    parents.append(key)
print('octets by level', count, 'total', sum(count), 'cells after', nC - 7*sum(count))
print('recorded: [7772, 2395, 1304, 15957] total 27428 cells 230898')
np.savez_compressed(OUT, level=level, group_of=group_of, types=types, count=np.array(count))
