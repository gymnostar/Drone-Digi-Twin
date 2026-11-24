"""Hub-height wind from N and from NE, read from the released lookup tables,
plus how much the two fields differ at every grid point."""
import sys
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Patch
from style import *

sys.path.insert(0, str(SCENES / 'joint'))
import lookup  # noqa: E402  (the released lookup API, with its checksum checks)

Z = 32.1
xs = np.arange(-200, 201, 5.0); ys = np.arange(-450, 201, 5.0)
fields = {d: lookup.WindField.for_direction(d) for d in ('N', 'NE')}
U = {d: np.array([[f.at(x, y, Z) or (np.nan,) * 3 for x in xs] for y in ys]) for d, f in fields.items()}
speed = {d: np.linalg.norm(u, axis=-1) for d, u in U.items()}
vdiff = np.linalg.norm(U['NE'] - U['N'], axis=-1)
sdiff = speed['NE'] - speed['N']
ok = ~np.isnan(vdiff)
rot = 2 * fields['N'].metadata['Uref_m_s'] * np.sin(np.radians(22.5))  # two equal winds 45° apart

BUILDING = '#c9c8c3'
blues = LinearSegmentedColormap.from_list('b', ['#f4f8fd'] + SEQ); blues.set_bad(BUILDING)
oranges = LinearSegmentedColormap.from_list('o', ['#fdf1ea', '#f9cbb4', '#f3a17c', ORANGE, '#c44c1c', '#8f3510', '#5c2109'])
oranges.set_bad(BUILDING)
div = LinearSegmentedColormap.from_list('d', ['#b0272a', RED, '#f2a7a6', '#f0efec', '#9ec5f4', BLUE, '#184f95'])
div.set_bad(BUILDING)

ex = (xs[0] - 2.5, xs[-1] + 2.5, ys[0] - 2.5, ys[-1] + 2.5)
fig, axs = plt.subplots(1, 4, figsize=(12.8, 6.9), constrained_layout=True, sharey=True)


def frame(ax, title, heading):
    ax.set_title(title, fontsize=17)
    ax.set_xlim(ex[:2]); ax.set_ylim(ex[2:]); ax.set_aspect('equal')
    ax.set_xticks([-150, 0, 150]); ax.set_xlabel('east (m)', fontsize=15)
    for h in heading:  # rotor bar across the wind, at the hub
        a = np.radians(h); perp = np.array([np.cos(a), -np.sin(a)])
        e = np.array([[0, -280] - 13.5 * perp, [0, -280] + 13.5 * perp])
        ax.plot(e[:, 0], e[:, 1], color=INK, lw=4 if len(heading) == 1 else 2, solid_capstyle='butt', zorder=5)


def arrow(ax, h, label, xy=(120, 115)):
    a = np.radians(h); d = -np.array([np.sin(a), np.cos(a)]) * 55
    ax.annotate('', np.array(xy) + d / 2, np.array(xy) - d / 2, arrowprops=dict(arrowstyle='-|>', color=INK, lw=2.5), zorder=7)
    ax.text(xy[0], xy[1] + 38, label, ha='center', fontsize=14, zorder=7, bbox=dict(fc='white', ec='none', alpha=0.85, pad=1))


vmax = 10
for ax, d, h in [(axs[0], 'N', 0), (axs[1], 'NE', 45)]:
    im = ax.imshow(np.ma.masked_invalid(speed[d]), origin='lower', extent=ex, cmap=blues, norm=Normalize(0, vmax), interpolation='nearest')
    frame(ax, f'From {d}: speed at {Z} m', [h]); arrow(ax, h, f'wind from {d}')
cb = fig.colorbar(im, ax=axs[:2], orientation='horizontal', fraction=0.04, pad=0.01, aspect=30)
cb.set_label('wind speed (m/s), one scale for N and NE')

im = axs[2].imshow(np.ma.masked_invalid(vdiff), origin='lower', extent=ex, cmap=oranges, norm=Normalize(0, 11), interpolation='nearest')
frame(axs[2], 'Vector difference', [0, 45])
cb = fig.colorbar(im, ax=axs[2], orientation='horizontal', fraction=0.04, pad=0.01, aspect=15)
cb.set_label('vector diff. (m/s)'); cb.ax.axvline(rot, color=INK, lw=2)
axs[2].text(0, 170, f'median {np.median(vdiff[ok]):.1f} m/s;\nthe 45° turn of an\n8.7 m/s wind alone\ngives {rot:.1f} m/s (tick)', ha='center', va='top', fontsize=13,
            zorder=7, bbox=dict(fc='white', ec='none', alpha=0.85, pad=2))

lim = np.nanmax(np.abs(sdiff))
im = axs[3].imshow(np.ma.masked_invalid(sdiff), origin='lower', extent=ex, cmap=div, norm=TwoSlopeNorm(0, -8, 8), interpolation='nearest')
frame(axs[3], 'Speed difference', [0, 45])
cb = fig.colorbar(im, ax=axs[3], orientation='horizontal', fraction=0.04, pad=0.01, aspect=15)
cb.set_label('NE − N speed (m/s)')
axs[3].text(0, 170, 'blue: faster from NE\nred: faster from N', ha='center', va='top', fontsize=13, zorder=7,
            bbox=dict(fc='white', ec='none', alpha=0.85, pad=2))

axs[0].set_ylabel('north (m)')
axs[0].annotate('rotor (hub 32.1 m)', (13, -280), (40, -400), fontsize=13, arrowprops=dict(arrowstyle='-', color=INK, lw=1.2),
                bbox=dict(fc='white', ec='none', alpha=0.85, pad=1), zorder=8)
fig.legend(handles=[Patch(fc=BUILDING, label='no value: inside a building (masked cell)')], loc='outside lower right', fontsize=14)
save(fig, 'figNE-heading-difference.png')

print(f'masked cells at {Z} m: N {np.isnan(speed["N"]).sum()}, NE {np.isnan(speed["NE"]).sum()}, both-valid {ok.sum()}')
print(f'speed max N {np.nanmax(speed["N"]):.2f} NE {np.nanmax(speed["NE"]):.2f}')
print(f'|dU| median {np.median(vdiff[ok]):.2f} P95 {np.percentile(vdiff[ok], 95):.2f} max {vdiff[ok].max():.2f}; rotation-only {rot:.2f}')
print(f'speed diff min {np.nanmin(sdiff):.2f} max {np.nanmax(sdiff):.2f}; |dspeed|>1 m/s share {(np.abs(sdiff[ok]) > 1).mean():.3f}')
