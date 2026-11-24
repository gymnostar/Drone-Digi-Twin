"""Mesh before/after from the exported common N mesh plus the reconstructed
sibling-octet merge (recon.npz, verified against the recorded per-level counts)."""
import sys, numpy as np, pyvista as pv
from collections import defaultdict
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle, Patch
from vtk.util.numpy_support import vtk_to_numpy as v2n
from style import *

S = sys.argv[1]
R = str(SCENES / 'review-export') + '/'
m = pv.read(R + 'joint-baseline.vtu')
d = np.load(S + '/recon.npz')
level, group = d['level'], d['group_of']
m.cell_data['cid'] = np.arange(m.n_cells)
sizes = [1.25, 2.5, 5, 10, 20, 40]
cols = ['#0d366b', '#1c5cab', '#3987e5', '#86b6ef', '#cde2fb', '#f7f9fc']
cmap = ListedColormap(cols); norm = BoundaryNorm(np.arange(7) - 0.5, 6)

def section(normal, origin, axes, lim):
    sl = m.slice(normal=normal, origin=origin)
    cid = np.asarray(sl.cell_data['cid']); P = np.asarray(sl.points)[:, axes]
    conn = v2n(sl.GetPolys().GetConnectivityArray()); off = v2n(sl.GetPolys().GetOffsetsArray())
    polys = [P[conn[off[i]:off[i+1]]] for i in range(len(off) - 1)]
    (x0, x1), (y0, y1) = lim
    keep = [i for i, p in enumerate(polys) if x0 - 45 < p[:, 0].mean() < x1 + 45 and y0 - 45 < p[:, 1].mean() < y1 + 45]
    polys = [polys[i] for i in keep]; cid = cid[keep]
    def si(p):
        w = (p.max(0) - p.min(0)).max(); return int(np.argmin([abs(np.log2(w / s)) for s in sizes]))
    before = [(p, si(p)) for p in polys]
    box = defaultdict(lambda: [np.inf, np.inf, -np.inf, -np.inf]); gsize = {}
    after, flag = [], []
    for p, c in zip(polys, cid):
        g = group[c]
        if g < 0:
            after.append((p, si(p))); flag.append((p, 0))
        else:
            b = box[g]; b[0] = min(b[0], p[:, 0].min()); b[1] = min(b[1], p[:, 1].min())
            b[2] = max(b[2], p[:, 0].max()); b[3] = max(b[3], p[:, 1].max()); gsize[g] = si(p) + 1
    for g, b in box.items():
        q = np.array([[b[0], b[1]], [b[2], b[1]], [b[2], b[3]], [b[0], b[3]]])
        after.append((q, gsize[g])); flag.append((q, 1))
    return before, after, flag

def draw(ax, items, lim, cm=cmap, nm=norm, lw=0.3):
    pc = PolyCollection([p for p, _ in items], array=np.array([c for _, c in items]), cmap=cm, norm=nm,
                        edgecolors='#3a4654', linewidths=lw)
    ax.add_collection(pc); ax.set_xlim(*lim[0]); ax.set_ylim(*lim[1]); ax.set_aspect('equal')
    return pc

# ---- Figure A: vertical north-south section through the rotor (x = 0)
lim = ((-470, 220), (0, 130))
before, after, flag = section('x', (0, 0, 0), [1, 2], lim)
fig, axs = plt.subplots(2, 1, figsize=(12.8, 7.6), constrained_layout=True)
pc = draw(axs[0], before, lim); draw(axs[1], after, lim)
axs[0].set_title('Strict reference mesh: 422,894 cells', fontsize=19)
axs[1].set_title('Production mesh after merging groups of 8 cells: 230,898 cells', fontsize=19)
for a in axs:
    a.set_ylabel('height (m)')
    a.plot([-280, -280], [32.1 - 13.5, 32.1 + 13.5], color=ORANGE, lw=4, solid_capstyle='butt', zorder=6)
    a.add_patch(Rectangle((-292, 1.5), 24, 64.5, fill=False, ec=ORANGE, lw=2, ls='--', zorder=6))
    a.annotate('rotor + protected box', (-268, 60), (-230, 105), fontsize=15, color=INK,
               bbox=dict(fc='white', ec='none', alpha=0.85, pad=2),
               arrowprops=dict(arrowstyle='-', color=INK, lw=1.2), zorder=7)
    a.annotate('', (130, 118), (205, 118), arrowprops=dict(arrowstyle='-|>', color=INK, lw=2.5), zorder=7)
    a.text(167, 122, 'wind (from N)', ha='center', va='bottom', fontsize=15, zorder=7,
           bbox=dict(fc='white', ec='none', alpha=0.85, pad=1))
axs[1].set_xlabel('north (m)   (section along x = 0, through the rotor; the wind blows right to left)')
cb = fig.colorbar(pc, ax=axs, orientation='vertical', fraction=0.03, pad=0.01, ticks=range(6))
cb.ax.set_yticklabels([f'{s:g} m' for s in sizes]); cb.set_label('cell edge length')
save(fig, 'fig03-mesh-before-after.png')

# ---- Figure B: near-ground plan view (z = 6 m): merged vs unchanged + whole-mesh counts
limp = ((-215, 215), (-470, 230))
_, _, flag = section('z', (0, 0, 6.0), [0, 1], limp)
fig = plt.figure(figsize=(12.8, 7.6), constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[0.9, 1])
ax = fig.add_subplot(gs[0])
draw(ax, flag, limp, ListedColormap(['#c9c8c3', AQUA]), BoundaryNorm([-.5, .5, 1.5], 2), lw=0.2)
ax.set_xlabel('east (m)'); ax.set_ylabel('north (m)')
ax.set_title('Merged cells at 6 m height', fontsize=19)
ax.plot([-13.5, 13.5], [-280, -280], color=ORANGE, lw=4, solid_capstyle='butt', zorder=6)
ax.add_patch(Rectangle((-26, -292), 52, 24, fill=False, ec=ORANGE, lw=2, ls='--', zorder=6))
ax.legend(handles=[Patch(fc=AQUA, ec='#55606b', label='merged (8 → 1)'),
                   Patch(fc='#c9c8c3', ec='#55606b', label='unchanged'),
                   Patch(fc='white', ec=ORANGE, ls='--', lw=2, label='rotor, protected')],
          loc='upper left', bbox_to_anchor=(0, -0.09), ncol=2, fontsize=14)
cnt = d['count']
bef = np.zeros(6, int)
for L in range(0, 5): bef[4 - L] = int((level == L).sum())
aft = bef.copy()
for L in range(4):
    aft[4 - L] -= 8 * int(cnt[L]); aft[5 - L] += int(cnt[L])
ax2 = fig.add_subplot(gs[1]); y = np.arange(6); h = 0.38
ax2.barh(y + h/2 + 0.03, bef, h, color=BLUE, label='strict reference (422,894)')
ax2.barh(y - h/2 - 0.03, aft, h, color=AQUA, label='production (230,898)')
ax2.set_yticks(y); ax2.set_yticklabels([f'{s:g} m' for s in sizes]); ax2.set_ylabel('cell edge length')
ax2.set_xlabel('cells in the whole domain (thousands)')
ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v/1000:.0f}'))
for yi, (b, a) in enumerate(zip(bef, aft)):
    ax2.text(b + 2500, yi + h/2 + 0.03, f'{b:,}', va='center', fontsize=14, color=INK2)
    ax2.text(a + 2500, yi - h/2 - 0.03, f'{a:,}', va='center', fontsize=14, color=INK2)
ax2.set_xlim(0, 290000); ax2.set_title('Cells by size, whole domain', fontsize=19)
ax2.legend(loc='lower right', fontsize=14); ax2.grid(axis='x'); ax2.set_axisbelow(True)
save(fig, 'fig04-mesh-what-changed.png')
print('before', bef.tolist(), bef.sum(), 'after', aft.tolist(), aft.sum())
