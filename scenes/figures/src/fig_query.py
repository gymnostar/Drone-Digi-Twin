"""What the lookup returns: real calls to scenes/joint/lookup.py on the released
N and NE tables, plotted as a wake transect and a vertical profile."""
import json, sys
import numpy as np
from matplotlib.lines import Line2D
from style import *

sys.path.insert(0, str(SCENES / 'joint'))
import lookup  # noqa: E402

F = {d: lookup.WindField.for_direction(d) for d in ('N', 'NE')}
timing = json.loads((SCENES / 'joint/lookup/query-timing.json').read_text())
COL = {'N': BLUE, 'NE': ORANGE}; MK = {'N': 'o', 'NE': 'D'}
spd = lambda v: np.nan if v is None else float(np.linalg.norm(v))
HUB, YT = 32.1, -350.0                     # transect: 70 m (2.6 rotor diameters) south of the hub
P = {'N': (0.0, -350.0), 'NE': (-50.0, -330.0)}  # profile points ~70 m downwind of the hub for each heading

fig = plt.figure(figsize=(12.8, 7.4))
ax = fig.add_axes([0.07, 0.42, 0.56, 0.51]); axp = fig.add_axes([0.72, 0.1, 0.26, 0.83])

# ---- transect across the wakes at hub height
xf = np.arange(-200, 200.01, 0.5); xn = np.arange(-200, 201, 5.0)
for d, f in F.items():
    ax.plot(xf, [spd(f.interpolate(x, YT, HUB)) for x in xf], color=COL[d], lw=2, label=f'wind from {d}')
    ax.scatter(xn, [spd(f.at(x, YT, HUB)) for x in xn], s=22, marker=MK[d], color=COL[d], edgecolor='white', lw=0.6, zorder=3)
for d, x, tx, ty, what in [('N', 0.0, 30, 3.4, 'N rotor wake'), ('N', -55.0, -195, 3.4, 'N: 60 m tower wake'),
                           ('NE', -70.0, -30, 10.1, 'NE rotor wake')]:
    v = F[d].at(x, YT, HUB)
    ax.annotate(f'{what}\nat({x:g}, {YT:g}, {HUB}) → {spd(v):.2f} m/s', (x, spd(v)), (tx, ty), fontsize=13, ha='left', va='center', arrowprops=dict(arrowstyle='-', color=INK, lw=1.2),
                bbox=dict(fc='white', ec='none', alpha=0.9, pad=1))
ax.set_xlim(-200, 200); ax.set_ylim(0, 11)
ax.set_xlabel(f'east (m), at north {YT:g} m and height {HUB} m')
ax.set_ylabel('wind speed |U| (m/s)')
ax.set_title('Transect across the wakes, 70 m south of the hub, at hub height', fontsize=18)
ax.grid(True); ax.set_axisbelow(True)
ax.legend(handles=[Line2D([], [], color=COL[d], marker=MK[d], ms=6, lw=2, label=f'wind from {d}') for d in F] +
          [Line2D([], [], color=MUTED, marker='o', ms=5, ls='', label='markers: grid nodes, at()'),
           Line2D([], [], color=MUTED, lw=2, label='line: interpolate() every 0.5 m')],
          loc='lower right', ncol=2, fontsize=13)

# ---- vertical profile 70 m downwind of the rotor, per heading
zf = np.arange(5, 80.01, 0.5); zn = lookup.HEIGHTS
for d, f in F.items():
    x, y = P[d]
    axp.plot([spd(f.interpolate(x, y, z)) for z in zf], zf, color=COL[d], lw=2)
    axp.scatter([spd(f.at(x, y, z)) for z in zn], zn, s=60, marker=MK[d], color=COL[d], edgecolor='white', lw=1, zorder=3,
                label=f'from {d} at ({x:g}, {y:g})')
    v = f.at(x, y, HUB)
    axp.annotate(f'{spd(v):.2f}', (spd(v), HUB), (spd(v) + (0.6 if d == 'NE' else -0.6), HUB + (6 if d == 'NE' else -6)),
                 fontsize=13, ha='left' if d == 'NE' else 'right', arrowprops=dict(arrowstyle='-', color=INK, lw=1))
axp.axhspan(HUB - 13.5, HUB + 13.5, color=GRID, zorder=0)
axp.text(0.3, HUB + 11.5, 'rotor height range', fontsize=13, color=INK2)
axp.set_xlim(0, 11); axp.set_ylim(0, 82); axp.set_yticks(zn); axp.set_yticklabels([f'{z:g}' for z in zn])
axp.set_xlabel('wind speed |U| (m/s)'); axp.set_ylabel('height above ground (m)')
axp.set_title('Profile 70 m downwind', fontsize=18)
axp.legend(loc='upper left', fontsize=13); axp.grid(True); axp.set_axisbelow(True)

# ---- the calls themselves, with what came back
n = F['N']; fmt = lambda v: 'None  (masked: inside a building)' if v is None else '(' + ', '.join(f'{c:+.2f}' for c in v) + ') m/s'
_, basis = lookup.WindField.for_heading(22.5)
picked = lookup.WindField.for_heading(22.5)[0].metadata['direction']
lines = [('N.at(0, -350, 32.1)', fmt(n.at(0, -350, 32.1))),
         ('N.interpolate(0, -352.5, 33)', fmt(n.interpolate(0, -352.5, 33))),
         ('N.at(-60, 0, 32.1)', fmt(n.at(-60, 0, 32.1))),
         ('WindField.for_heading(22.5)', f'{picked} field, basis "{basis}"')]
fig.text(0.015, 0.305, 'Real calls to scenes/joint/lookup.py (velocity = east, north, up)', fontsize=15, weight='bold', va='top')
for i, (call, out) in enumerate(lines):
    yy = 0.245 - i * 0.042
    fig.text(0.015, yy, call, fontsize=13, family='monospace', va='top')
    fig.text(0.27, yy, '→ ' + out, fontsize=13, va='top', color=INK2)
fig.text(0.015, 0.008, f'Measured query time: {timing["direct_query"]["median_batch_mean_us_per_call"]:.2f} µs per grid-node call, '
         f'{timing["interpolated_query"]["median_batch_mean_us_per_call"]:.2f} µs interpolated\n'
         '(median of 5 warm batches of 10,000 Python calls, one Mac)', fontsize=13, va='bottom', color=INK2)
save(fig, 'figquery-lookup-example.png')
for call, out in lines: print(call, '->', out)
for d in F:
    print(d, 'transect min speed', np.nanmin([spd(F[d].at(x, YT, HUB)) for x in xn]), 'profile', [round(spd(F[d].at(*P[d], z)), 2) for z in zn])
