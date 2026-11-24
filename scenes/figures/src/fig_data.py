import json, numpy as np
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from style import *

W = str(SCENES) + '/'
J = W + 'joint/results/'
def load(p): return json.load(open(p))
def metrics(c): return load(J + c + '/metrics.json')

# ------------------------------------------------------------ 06 timing
runs = [('N strict reference', 'poc-N-reference-tight-transport'),
        ('N same mesh, looser stopping', 'poc-N-baseline-retimed'),
        ('N production (merged mesh)', 'poc-N-agglomerated-siblings'),
        None,
        ('NE strict reference (retimed)', 'table-NE-reference-retimed'),
        ('NE production, SIMPLE factors', 'table-NE-production'),
        ('NE production, SIMPLEC 0.9 (retimed)', 'table-NE-simplec-retimed')]
stage_cols = [('mesh', 'mesh', BLUE), ('solve', 'solve', ORANGE),
              ('sample', 'extract samples', AQUA), ('other', 'setup + checks', YELLOW)]
fig, ax = plt.subplots(figsize=(12.8, 7.2), constrained_layout=True)
ylab, ypos, y = [], [], 0
rows = {}
for r in runs:
    if r is None: y += 0.6; continue
    name, case = r; m = metrics(case); c = m['cost']
    st = {'mesh': c['mesh']['elapsed_s'], 'solve': c['solve']['elapsed_s'], 'sample': c['sample']['elapsed_s']}
    st['other'] = m['complete_map_generation_elapsed_s'] - sum(st.values())
    rows[case] = (st, m)
    left = 0
    for k, _, col in stage_cols:
        ax.barh(y, st[k], left=left, color=col, height=0.68, edgecolor='white', linewidth=2)
        if st[k] > 35: ax.text(left + st[k] / 2, y, f'{st[k]:.0f}', ha='center', va='center', fontsize=15, color='white' if col in (BLUE, ORANGE) else INK)
        left += st[k]
    tot = m['complete_map_generation_elapsed_s']; it = m['solver_final_iteration']
    tx = tot + 6 if not (case == 'poc-N-agglomerated-siblings') else 272
    ax.text(tx, y, f'{tot:.1f} s  ({m["cells"]:,} cells, {it} it.)', va='center', fontsize=15, color=INK)
    ylab.append(name); ypos.append(y); y += 1
# 2x lines: half of each reference
nref = rows['poc-N-reference-tight-transport'][1]['complete_map_generation_elapsed_s']
eref = rows['table-NE-reference-retimed'][1]['complete_map_generation_elapsed_s']
for half, (y0, y1) in [(nref / 2, (ypos[0] - .5, ypos[2] + .5)), (eref / 2, (ypos[3] - .5, ypos[5] + .5))]:
    ax.plot([half, half], [y0, y1], color=INK, lw=2, ls=(0, (4, 3)))
    ax.text(half + 4, y0 + 0.02, f'2x target: {half:.1f} s', ha='left', va='bottom', fontsize=14, color=INK,
            bbox=dict(fc='white', ec='none', pad=1))
ax.set_yticks(ypos); ax.set_yticklabels(ylab)
ax.set_xlabel('cold complete map generation time (s), one timing each, Studio, 3 CPUs, 5 GiB')
ax.set_xlim(0, 780); ax.set_ylim(ypos[-1] + 0.7, ypos[0] - 0.95); ax.grid(axis='x'); ax.set_axisbelow(True)
ax.legend(handles=[Patch(color=c, label=l) for _, l, c in stage_cols], ncol=4, loc='lower right',
          bbox_to_anchor=(1.0, 1.0), fontsize=15)
save(fig, 'fig05-timing-breakdown.png')

# ------------------------------------------------------------ 07 verdict scatter
fig, ax = plt.subplots(figsize=(12.8, 7.2), constrained_layout=True)
ax.add_patch(Rectangle((2, 0), 8, 5, color='#e7f5ee', zorder=0))
ax.text(9.8, 0.25, 'target: at least 2x faster and worst zone at most 5%', ha='right', fontsize=14, color=INK2)
ax.axvline(2, color=INK, lw=1.6, ls=(0, (4, 3))); ax.axhline(5, color=INK, lw=1.6, ls=(0, (4, 3)))
pts = []
for i in (1, 2, 3):
    d = load(W + f'turbine/results/poc-studio-holdout{i}-verified-comparison.json')
    pts.append(('turbine', d['cost']['speedup'], 100 * d['numerical_regions']['wakeROI']['vector_RMSE_over_Uref']))
for x, yv in [(p[1], p[2]) for p in pts]:
    ax.scatter(x, yv, s=150, color=BLUE, edgecolor='white', linewidth=2, zorder=3)
ax.text(6.0, 2.25, 'turbine alone\n3 held-out conditions\n5.97-6.09x, 1.38-1.63%', ha='center', fontsize=15)
city = load(W + 'buildings/results/poc-study-summary.json')['provisional_comparison']
cx, cy = city['complete_elapsed_ratio_descriptive'], 100 * city['regions']['city']['vector_RMSE_over_Uref']
ax.scatter(cx, cy, s=150, facecolor='white', edgecolor=ORANGE, linewidth=3, zorder=3)
ax.text(cx - 0.25, cy, 'city alone: 6.42%\n(8.6x only vs an unconverged\nreference: no certified speed-up)', ha='right', va='center', fontsize=14)
def worst(case):
    r = metrics(case)['comparison']['regions']
    return max(100 * r[z]['vector_RMSE_over_Uref'] for z in ('city', 'building_wake', 'rotor_wake', 'flight'))
cpl = [('poc-N-baseline-retimed', 'N, same mesh, looser stopping', 's', 'left', (-0.55, 0.35)),
       ('poc-N-agglomerated-siblings', 'N production (released)', 'o', 'left', (0.12, 0.2)),
       ('table-NE-production', 'NE first release (SIMPLE)', 'D', 'left', (-0.64, -1.0)),
       ('table-NE-simplec-retimed', 'NE released (SIMPLEC)', '^', 'left', (0.12, 0.0))]
for case, lab, mk, ha, (dx, dy) in cpl:
    m = metrics(case); x = m['comparison']['speedup']; yv = worst(case)
    ax.scatter(x, yv, s=170, marker=mk, color=AQUA, edgecolor=INK, linewidth=1.2, zorder=4)
    sep = ': ' if 'SIMPLEC' in lab else '\n'
    ax.text(x + dx, yv + dy, f'{lab}{sep}{x:.2f}x, {yv:.2f}%', ha=ha, va='center' if sep == ': ' else 'baseline', fontsize=14, bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.5))
ax.set_xlim(1, 10); ax.set_ylim(0, 7.5)
ax.set_xlabel('speed-up vs strict reference (reference time / cheap time)')
ax.set_ylabel('vector RMSE, worst zone (% of inflow)')
ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=12, color=BLUE, label='turbine alone'),
                   Line2D([], [], marker='o', ls='', ms=12, mfc='white', mec=ORANGE, mew=3, label='city alone'),
                   Line2D([], [], marker='o', ls='', ms=12, color=AQUA, mec=INK, label='buildings + turbine')],
          loc='upper left', bbox_to_anchor=(0.32, 1.0), fontsize=15)
ax.grid(True); ax.set_axisbelow(True)
save(fig, 'fig06-verdict-speedup-vs-error.png')

# ------------------------------------------------------------ 08 per-zone errors
zones = [('city', 'city'), ('building_wake', 'building wake'), ('rotor_wake', 'rotor wake'), ('flight', 'flight region')]
cases = [('poc-N-agglomerated-siblings', 'N production', BLUE),
         ('table-NE-production', 'NE first release (SIMPLE)', ORANGE),
         ('table-NE-simplec-retimed', 'NE released (SIMPLEC)', AQUA)]
fig, axs = plt.subplots(1, 2, figsize=(12.8, 7.0), constrained_layout=True)
for ax, key, lim, title in [(axs[0], 'vector_RMSE_over_Uref', 5, 'Vector RMSE (limit 5%)'),
                            (axs[1], 'P95_over_Uref', 15, 'P95 vector error (limit 15%)')]:
    x = np.arange(4); w = 0.26
    for i, (case, lab, col) in enumerate(cases):
        r = metrics(case)['comparison']['regions']
        v = [100 * r[z][key] for z, _ in zones]
        b = ax.bar(x + (i - 1) * w, v, w - 0.03, color=col, label=lab)
        for xi, vi in zip(x + (i - 1) * w, v):
            ax.text(xi, vi + lim * 0.012, f'{vi:.1f}', ha='center', va='bottom', fontsize=12.5, color=INK2)
    ax.axhline(lim, color=INK, lw=2, ls=(0, (4, 3)))
    ax.text(3.45, lim * 1.01, f'limit {lim}%', ha='right', va='bottom', fontsize=15)
    ax.set_xticks(x); ax.set_xticklabels([n.replace(' ', '\n') for _, n in zones])
    ax.set_ylim(0, lim * 1.12); ax.set_title(title, fontsize=19)
    ax.set_ylabel('% of 8.7 m/s inflow'); ax.grid(axis='y'); ax.set_axisbelow(True)
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc='outside lower center', ncol=3, fontsize=15)
save(fig, 'fig07-errors-by-zone.png')

# ------------------------------------------------------------ 09 residual vs limit for strict solves
def rv(d, key='last_initial_residuals'):
    r = d[key]; u = max(r.get(k, 0) for k in ('Ux', 'Uy', 'Uz', 'max_U'))
    return r['p'], u
rows = []
m = metrics('poc-N-reference-fixed'); rows.append(('N strict, original settings', *rv(m), int(m['solver_final_iteration']), False))
att = {a['case']: a for a in load(J + 'optimization-summary.json')['attempts']}
d = att['poc-N-reference-damped']; rows.append(('N strict, damped momentum (0.45)', *rv(d, 'last_complete_initial_residuals'), d['last_complete_iteration'], False))
m = metrics('poc-N-reference-tight-transport'); rows.append(('N strict, tighter inner tolerance', *rv(m), int(m['solver_final_iteration']), True))
m = metrics('table-NE-reference-retimed'); rows.append(('NE strict', *rv(m), int(m['solver_final_iteration']), True))
for case, lab, f in [('table-NNE-corners', 'NNE 22.5° strict, SIMPLEC', 'solve-verdict.json'),
                     ('table-NNE-simple', 'NNE, standard SIMPLE', 'solve-verdict.json'),
                     ('nne-one-corrector', 'NNE, one fewer pressure correction', 'trial.json'),
                     ('nne-simplec-tutorial', 'NNE, SIMPLEC relaxation 0.9', 'trial.json'),
                     ('table-NbE-corner', 'NbE 11.25° strict', 'solve-verdict.json'),
                     ('table-NW-oriented', 'NW 315° strict', 'solve-verdict.json')]:
    d = load(J + case + '/' + f); rows.append((lab, *rv(d, 'last_complete_max_initial_residuals'), d['last_complete_iteration'], False))
fig, ax = plt.subplots(figsize=(12.8, 7.6), constrained_layout=True)
ys = np.arange(len(rows))
ax.axvspan(0.05, 1, color='#e7f5ee', zorder=0)
ax.axvline(1, color=INK, lw=2, ls=(0, (4, 3)))
for y, (lab, p, u, it, ok) in zip(ys, rows):
    pr, ur = p / 1e-4, u / 1e-5
    ax.plot([min(pr, ur), max(pr, ur)], [y, y], color=GRID, lw=3, zorder=1)
    ax.scatter(pr, y, s=170, marker='o', color=BLUE, edgecolor='white', linewidth=1.5, zorder=3)
    ax.scatter(ur, y, s=190, marker='D', color=ORANGE, edgecolor='white', linewidth=1.5, zorder=3)
    tag = 'converged' if ok else 'stopped at 600 s'
    ax.text(40, y, f'{it} it., {tag}', va='center', fontsize=14, color=INK if ok else INK2, weight='bold' if ok else 'normal')
ax.set_xscale('log'); ax.set_xlim(0.05, 400)
ax.set_xticks([0.1, 1, 10, 30]); ax.set_xticklabels(['0.1x', '1x = limit', '10x', '30x'])
ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows]); ax.invert_yaxis()
ax.set_xlabel('final residual divided by its strict limit (log scale; left of the dashed line passes)')
ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=12, color=BLUE, label='pressure (limit 1e-4)'),
                   Line2D([], [], marker='D', ls='', ms=12, color=ORANGE, label='velocity, largest component (limit 1e-5)')],
          loc='lower center', bbox_to_anchor=(0.4, 1.0), ncol=2, fontsize=15)
ax.grid(axis='x'); ax.set_axisbelow(True)
save(fig, 'fig08-strict-residuals-vs-limit.png')
for r in rows: print('resid', r[0], f'p/lim={r[1]/1e-4:.2f} U/lim={r[2]/1e-5:.2f}', r[3])
