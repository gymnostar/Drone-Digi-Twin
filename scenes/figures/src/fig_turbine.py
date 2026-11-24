import numpy as np, pyvista as pv
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from vtk.util.numpy_support import vtk_to_numpy as v2n
from style import *
R = str(SCENES / 'review-export') + '/'
cmap = LinearSegmentedColormap.from_list('b', ['#f4f8fd', '#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'])
nrm = Normalize(4.5, 10)
fig, axs = plt.subplots(2, 1, figsize=(12.8, 7.4), constrained_layout=True, sharex=True)
for ax, (name, title) in zip(axs, [('turbine-reference', 'Reference mesh: 359,334 cells, 131.1 s to generate the map'),
                                   ('turbine-cheap', 'Cheap mesh: 112,206 cells, 22.1 s (wake cells 4x wider)')]):
    t = pv.read(R + name + '.vtu'); t.cell_data['speed'] = np.linalg.norm(t.cell_data['U'], axis=1)
    sl = t.clip_box((-60, 340, -5, 5, 0, 85), invert=False).slice(normal='y', origin=(0, 0, 0))
    P = np.asarray(sl.points)[:, [0, 2]]
    conn = v2n(sl.GetPolys().GetConnectivityArray()); off = v2n(sl.GetPolys().GetOffsetsArray())
    polys = [P[conn[off[i]:off[i+1]]] for i in range(len(off) - 1)]
    pc = PolyCollection(polys, array=np.asarray(sl.cell_data['speed']), cmap=cmap, norm=nrm, edgecolors='#3a4654', linewidths=0.25)
    ax.add_collection(pc); ax.set_xlim(-40, 330); ax.set_ylim(0, 80); ax.set_aspect('equal')
    ax.plot([0, 0], [18.6, 45.6], color=ORANGE, lw=4, solid_capstyle='butt')
    ax.set_title(title, fontsize=19); ax.set_ylabel('height (m)')
axs[0].annotate('rotor, 27 m', (0, 45.6), (25, 70), fontsize=15, arrowprops=dict(arrowstyle='-', color=INK, lw=1.2),
                bbox=dict(fc='white', ec='none', alpha=0.85, pad=2))
axs[1].set_xlabel('distance downstream of the rotor (m)   (vertical section through the hub; wind blows left to right)')
cb = fig.colorbar(pc, ax=axs, fraction=0.03, pad=0.01); cb.set_label('wind speed (m/s)')
save(fig, 'fig02-turbine-alone-mesh.png')
