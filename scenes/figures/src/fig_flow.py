"""Pipeline from geometry to a wind query. Numbers are read from the recorded
results and the lookup files, so the diagram cannot drift from them."""
import json
from matplotlib.patches import FancyBboxPatch
from style import *

J = SCENES / 'joint'
def load(p): return json.loads((J / p).read_text())
m = {k: load(f'results/{c}/metrics.json') for k, c in [
    ('Nref', 'poc-N-reference-tight-transport'), ('N', 'poc-N-agglomerated-siblings'),
    ('NEref', 'table-NE-reference-retimed'), ('NE', 'table-NE-simplec-retimed')]}
idx = load('lookup/index.json')['fields']
timing = load('lookup/query-timing.json')
with __import__('gzip').open(J / 'lookup' / idx['N']['file'], 'rt') as f:
    points = json.loads(f.readline()[2:])['unique_positions']
def worst(k):
    r = m[k]['comparison']['regions']
    return max(100 * r[z]['vector_RMSE_over_Uref'] for z in ('city', 'building_wake', 'rotor_wake', 'flight'))
def cover(k):
    r = m[k]['comparison']['regions']
    return min(r[z]['coverage'] for z in ('city', 'building_wake', 'rotor_wake', 'flight'))
t = {k: m[k]['complete_map_generation_elapsed_s'] for k in m}
it = {k: m[k]['solver_final_iteration'] for k in m}

fig = plt.figure(figsize=(12.8, 7.2))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.set_xlim(0, 128); ax.set_ylim(0, 72)
W, H = 23.5, 16


def box(x, y, title, body, fc='#f2f6fc', ec=BLUE):
    ax.add_patch(FancyBboxPatch((x, y), W, H, boxstyle='round,pad=0,rounding_size=1.2', fc=fc, ec=ec, lw=2))
    ax.text(x + 1.2, y + H - 1.4, title, fontsize=15, weight='bold', va='top', color=INK)
    ax.text(x + 1.2, y + H - 5.0, body, fontsize=12.5, va='top', color=INK2, linespacing=1.35)


def arrow(a, b, color=INK, **kw):
    ax.annotate('', b, a, arrowprops=dict(arrowstyle='-|>', color=color, lw=2, shrinkA=0, shrinkB=0, **kw))


xs = [1.5 + 25.5 * i for i in range(5)]; y1, y2 = 41, 12
ax.text(1.5, 68, 'From geometry to a wind value: one pipeline per wind heading', fontsize=20, weight='bold', va='center')
box(xs[0], y1, '1  Geometry', 'Niigata buildings\n(STL, full scale)\nSWiFT rotor: 27 m disk,\nhub 32.1 m (force region)')
box(xs[1], y1, '2  Mesh', 'blockMesh +\nsnappyHexMesh, strict:\n' f'{m["Nref"]["cells"]:,} cells (N)\n{m["NEref"]["cells"]:,} cells (NE)')
box(xs[2], y1, '3  Merge (production)', 'groups of 8 interior cells\nbecome 1; walls, ground,\nrotor-force cells untouched\n' f'{m["N"]["cells"]:,} (N), {m["NE"]["cells"]:,} (NE)',
    fc='#e7f5ee', ec=AQUA)
box(xs[3], y1, '4  Solve', 'steady RANS, OpenFOAM\nv2306, realizable k-ε\n' f'strict: {it["Nref"]} / {it["NEref"]} it. (N / NE)\nproduction: {it["N"]} / {it["NE"]} it.')
box(xs[4], y1, '5  Sample', 'velocity at fixed points\n(OpenFOAM cellPoint\ninterpolation); same\npoints in both runs')
for a, b in zip(xs[:-1], xs[1:]):
    arrow((a + W, y1 + H / 2), (b, y1 + H / 2))
# the strict reference skips the merge
arrow((xs[1] + W - 4, y1 + H), (xs[3] + 4, y1 + H), color=MUTED, ls='--', connectionstyle='arc3,rad=-0.25')
ax.text(xs[2] + W / 2, y1 + H + 4.3, 'strict reference skips step 3', ha='center', va='bottom', fontsize=13, color=INK2,
        bbox=dict(fc='white', ec='none', pad=1))
# 5 -> 6 wraps down between the rows
ym = (y1 + y2 + H) / 2
ax.plot([xs[4] + W / 2, xs[4] + W / 2, xs[0] + W / 2], [y1, ym, ym], color=INK, lw=2, solid_joinstyle='miter')
arrow((xs[0] + W / 2, ym), (xs[0] + W / 2, y2 + H))

r2 = [xs[0], xs[0] + 32, xs[0] + 64]
box(r2[0], y2, '6  Compare', 'production vs strict,\n4 zones: RMSE ≤ 5%,\nP95 ≤ 15%, coverage ≥ 95%\n' f'worst zone: N {worst("N"):.2f}%,\nNE {worst("NE"):.2f}%; coverage {100 * min(cover("N"), cover("NE")):.0f}%',
    fc='#fdf1ea', ec=ORANGE)
box(r2[1], y2, '7  Lookup field', f'production field, 5 m grid,\n7 heights (5 to 80 m):\n{points:,} points per heading\n(east, north, up velocity)')
box(r2[2], y2, '8  Query', 'query(x, y, z, heading)\n→ velocity, or "no value"\nin a building\n'
    f'{timing["direct_query"]["median_batch_mean_us_per_call"]:.2f} µs per grid-point call\n'
    f'{timing["interpolated_query"]["median_batch_mean_us_per_call"]:.2f} µs interpolated')
arrow((r2[0] + W, y2 + H / 2), (r2[1], y2 + H / 2))
ax.text(r2[0] + W + 4.25, y2 + H / 2 + 1, 'passes', ha='center', fontsize=12.5, color=INK2)
arrow((r2[1] + W, y2 + H / 2), (r2[2], y2 + H / 2))

x = r2[2] + W + 8
ax.text(x, y2 + H - 1.4, 'Cold time, steps 2–5', fontsize=15, weight='bold', va='top')
ax.text(x, y2 + H - 5.0, f'N:  {t["Nref"]:.0f} → {t["N"]:.0f} s  ({t["Nref"] / t["N"]:.2f}×)\n'
        f'NE: {t["NEref"]:.0f} → {t["NE"]:.0f} s  ({t["NEref"] / t["NE"]:.2f}×)\n'
        'strict → production,\none timing each (Studio,\n3 CPUs, 5 GiB)', fontsize=12.5, va='top', color=INK2, linespacing=1.35)
ax.text(1.5, 4.5, 'Each wind heading needs its own pass through steps 1–7. Released so far: N (0°) and NE (45°).',
        fontsize=13.5, color=INK2)
save(fig, 'figflow-pipeline.png')
print({k: (m[k]['cells'], round(t[k], 2), it[k]) for k in m}, 'worst', worst('N'), worst('NE'), 'coverage', cover('N'), cover('NE'), 'points', points)
