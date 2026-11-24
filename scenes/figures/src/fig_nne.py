import json, numpy as np, pyvista as pv
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from style import *
W = str(SCENES) + '/'
d = json.load(open(W + 'joint/results/nne-residual-diagnostic.json'))
sc = d['cases'][0]; assert sc['coupling'] == 'SIMPLEC'
it = sc['iterations']['902']
stl = pv.read(W + 'review-export/buildings.stl').triangulate()
P = np.asarray(stl.points); F = np.asarray(stl.faces).reshape(-1, 4)[:, 1:]
tris = P[F][:, :, :2]
u = tris[:, 1] - tris[:, 0]; v = tris[:, 2] - tris[:, 0]; area = np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
tris = tris[area > 1e-3]
fig = plt.figure(figsize=(12.8, 7.6), constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1])
ax = fig.add_subplot(gs[0])
ax.add_collection(PolyCollection(tris, facecolors='#c9c8c3', edgecolors='#c9c8c3', linewidths=0.3))
# bins (25 m cubes): sum fractions per plan column
def columns(field):
    col = {}
    for b in it[field]['world_25m_bins_by_L1']:
        x, y, z = b['lower_corner_m']; col.setdefault((x, y), []).append((z, b['L1_fraction']))
    return col
for (x, y), zs in columns('Ux').items():
    f = sum(v for _, v in zs)
    ax.add_patch(Rectangle((x, y), 25, 25, fc=ORANGE, alpha=0.35, ec=ORANGE, lw=2.5, zorder=3))
    zr = f'{min(z for z, _ in zs):.0f}-{max(z for z, _ in zs) + 25:.0f} m'
    if f < 0.03: continue
    tx, ty = (-120, 75) if x == -100 else (5, -60)
    ax.annotate(f'{100*f:.1f}% of crosswind\nresidual (z {zr})', (x + 12.5, y + 12.5), (tx, ty), ha='center', fontsize=13, zorder=6,
                bbox=dict(fc='white', ec='none', alpha=0.9, pad=1), arrowprops=dict(arrowstyle='-', color=INK, lw=1.2))
pc = it['p']['world_25m_bins_by_L1'][1]; assert pc['lower_corner_m'] == [-175, 0, 0]
ax.add_patch(Rectangle((-175, 0), 25, 25, fc=BLUE, alpha=0.3, ec=BLUE, lw=2.5, zorder=3))
ax.text(-162.5, -12, f'{100*pc["L1_fraction"]:.1f}% of pressure\nresidual (z 0-25 m)', ha='center', va='top', fontsize=13, zorder=6,
        bbox=dict(fc='white', ec='none', alpha=0.8, pad=1))
for c in it['Ux']['top_3_cells']:
    ax.scatter(*c['world_m'][:2], marker='D', s=110, color=ORANGE, edgecolor=INK, lw=1, zorder=5)
for c in it['p']['top_3_cells'][:2]:
    ax.scatter(*c['world_m'][:2], marker='o', s=110, color=BLUE, edgecolor=INK, lw=1, zorder=5)
a = np.radians(22.5); flow = -np.array([np.sin(a), np.cos(a)]); perp = np.array([np.cos(a), -np.sin(a)])
hub = np.array([0, -280]); e = np.array([hub - 13.5 * perp, hub + 13.5 * perp])
ax.plot(e[:, 0], e[:, 1], color=INK, lw=5, solid_capstyle='butt', zorder=5)
ax.text(-8, -262, 'rotor: 0.02% of the\ncrosswind residual', ha='right', fontsize=14)
s0 = np.array([40, 60]); ax.annotate('', s0 + 45 * flow, s0, arrowprops=dict(arrowstyle='-|>', lw=2.5, color=INK))
ax.text(s0[0] + 4, s0[1] + 4, 'wind from 22.5°', fontsize=14)
ax.set_xlim(-215, 75); ax.set_ylim(-305, 110); ax.set_aspect('equal')
ax.set_xlabel('east (m)'); ax.set_ylabel('north (m)')
ax.set_title('Where the NNE solve is stuck (plan view)', fontsize=19)
ax.legend(handles=[Line2D([], [], marker='s', ls='', ms=14, mfc=ORANGE, alpha=0.5, mec=ORANGE, label='25 m blocks with the most crosswind-momentum residual'),
                   Line2D([], [], marker='D', ls='', ms=10, mfc=ORANGE, mec=INK, label='3 largest crosswind cells'),
                   Line2D([], [], marker='o', ls='', ms=10, mfc=BLUE, mec=INK, label='2 largest pressure cells (building-ground corner)')],
          loc='upper left', bbox_to_anchor=(0, -0.1), fontsize=13)
ax2 = fig.add_subplot(gs[1])
reg = it['Ux']['regions']
cols = columns('Ux'); top_col = max(cols.values(), key=lambda zs: sum(v for _, v in zs))
vals = [('2 stacked 25 m blocks\nat one building edge', 100 * sum(v for _, v in top_col)),
        (f'cells touching buildings\n({reg["buildings"]["mesh_cells"]:,} cells)', 100 * reg['buildings']['L1_fraction']),
        (f'cells touching ground\n({reg["ground"]["mesh_cells"]:,} cells)', 100 * reg['ground']['L1_fraction']),
        (f'rotor force cells\n({reg["rotor_source_cells"]["mesh_cells"]:,} cells)', 100 * reg['rotor_source_cells']['L1_fraction'])]
y = np.arange(len(vals))
ax2.barh(y, [v for _, v in vals], 0.6, color=[ORANGE, MUTED, MUTED, INK])
for yi, (_, v) in zip(y, vals):
    ax2.text(v + 0.2, yi, f'{v:.2f}%' if v < 0.1 else f'{v:.1f}%', va='center', fontsize=16)
ax2.set_yticks(y); ax2.set_yticklabels([l for l, _ in vals], fontsize=15); ax2.invert_yaxis()
ax2.set_xlim(0, 13); ax2.set_xlabel(f'% of crosswind residual\n(whole mesh: {d["cells"]:,} cells)')
ax2.set_title('Share of remaining imbalance', fontsize=19)
ax2.grid(axis='x'); ax2.set_axisbelow(True)
save(fig, 'fig09-nne-residual-location.png')
