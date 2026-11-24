from style import *
groups = [
 ('Mesh changes', [
  ('N: coarsen the city floor one level', 'rotor-force cells kept, but their corners moved', 'rejected before solve'),
  ('N: merge only perfectly box-shaped cells', 'only 2.7% of cells qualified; boundary check stopped it', 'rejected before solve'),
  ('N: merge groups of 8 after snapping', '0.55 mm audit flag was a false alarm; exact check passed', 'accepted: 230,898 cells'),
  ('NE: first strict mesh', '531,607 cells, above the 500,000 cap', 'fixed: 425,006 cells'),
  ('NNE, NbE: first strict meshes', 'a few distorted cells at building corners', 'fixed: local refinement'),
  ('NNE: refine edge before snapping', 'rotor-force corners moved up to 0.65 µm', 'rejected before solve'),
  ('NNE: split edge cells after snapping', 'one cell had 6 eligible corners; the splitter needs 8', 'stopped by OpenFOAM'),
 ]),
 ('Solver changes', [
  ('N strict: original settings', '600 s limit hit at 995 it.; pressure 2.4x its limit', 'timed out'),
  ('N strict: damp momentum to 0.45', '600 s hit; worse (pressure 8.9x, velocity 13x)', 'timed out'),
  ('N strict: inner tolerance 1e-7 to 1e-10', 'converged in 639 it., 525.27 s', 'accepted as reference'),
  ('N: SIMPLEC on the fine mesh', 'wake monitors moved 0.247% > 0.20% limit', 'rejected'),
  ('NNE (4 settings), NbE, NW strict', 'all hit 600 s above limits (worst: velocity 27x)', 'no checked field'),
  ('NE production: SIMPLEC 0.9', '136 vs 277 it.; 200.40 s, needed 199.43 s', '1.99x: 2x missed'),
 ])]
fig = plt.figure(figsize=(12.8, 7.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off')
X = [0.02, 0.33, 0.77]; y = 0.955; dy = 0.054
for i, h in enumerate(['What was tried', 'Measured result', 'Outcome']):
    ax.text(X[i], y, h, fontsize=16, color=INK2, weight='bold', va='center')
y -= dy * 0.9
for g, rows in groups:
    ax.text(X[0], y - 0.012, g, fontsize=17, weight='bold', va='center', color=INK)
    y -= dy
    for a, b, c in rows:
        ax.plot([0.02, 0.98], [y + dy / 2, y + dy / 2], color=GRID, lw=1)
        ax.text(X[0], y, a, fontsize=14.5, va='center')
        ax.text(X[1], y, b, fontsize=14.5, va='center', color=INK2)
        ok = c.startswith('accepted') or c.startswith('fixed')
        ax.text(X[2], y, ('✓ ' if ok else '✗ ') + c, fontsize=14.5, va='center', weight='bold', color=INK)
        y -= dy
    y -= dy * 0.25
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
save(fig, 'fig10-what-failed-and-why.png')
