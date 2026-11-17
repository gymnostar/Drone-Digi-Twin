#!/usr/bin/env python3
"""Pictures for DLR/meshStudy/explainer.html. Light: numpy + matplotlib only, no OpenFOAM.
Run from the repo root:  .venv/bin/python DLR/results/explainer_fig/make_figs.py

fig_site_campus.png, fig_site_zoom.png   top view: buildings, cluster box, station +-50 m, probe lines, wind
fig_section.png   vertical cut along the wind through the station: ground, buildings, box, probe heights
fig_cells.png     ILLUSTRATION (made-up flow): what 4 m / 2 m / 1 m cells "see" behind a building
fig_expected.png  SKETCH (invented curves): what the morning plots should look like
"""
import logging
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
HERE = Path(__file__).resolve().parent
MS = HERE.parent.parent / "meshStudy"

# ---- facts copied from DLR/meshStudy/setup.py (keep in sync) ----
BOX_MIN, BOX_MAX = (50.0, 145.0, 12.0), (235.0, 262.0, 60.0)
LOC = {  # name: (x, y, heights above ground)
    "upstream": (-150.0, 206.0, [2, 5, 10, 20, 40]),
    "gap":      (172.0, 182.0, [2, 5, 10, 20, 40]),
    "wake":     (140.0, 210.0, [2, 5, 10, 20, 40]),
    "station":  (122.0, 206.0, [15, 20, 25, 40]),
}
COL = {"upstream": "#059669", "gap": "#d97706", "wake": "#7c3aed", "station": "#dc2626"}
# the 5 highly skewed faces of the coarse mesh (centres, from coarse/postProcessing/constant/skewFaces/skewFaces.vtp)
SKEW = [(-145, -70), (-217, -153), (150, 290), (-193, -255)]   # last location holds 2 faces

site = np.load(MS / "site_map.npz")
x, y, hb, hg, st = site["x"], site["y"], site["hb"], site["hg"], site["st"]
H = np.where(np.isfinite(hb), hb - hg, 0.0)
H = np.where(np.isfinite(H), H, 0.0)
H = ndi.grey_closing(H, size=(3, 3))          # the ray-cast map is speckled: close the small holes
Hm = np.ma.masked_less(H, 1.0)


def ground(px, py):
    return float(hg[np.abs(y - py).argmin(), np.abs(x - px).argmin()])


def wind_arrow(ax, x0, y0, L, **kw):
    ax.add_patch(FancyArrowPatch((x0, y0), (x0 + L, y0), arrowstyle="simple,head_width=10,head_length=12",
                                 mutation_scale=1.6, color="#2563eb", **kw))


# ------------------------------------------------------------------ fig_site_campus / fig_site_zoom
WB = dict(fc="white", alpha=.85, ec="none", pad=1.5)


def base_map(ax, ms):
    im = ax.pcolormesh(x, y, Hm, cmap="Greys", vmin=0, vmax=24, shading="auto")
    ax.add_patch(Rectangle(BOX_MIN[:2], BOX_MAX[0] - BOX_MIN[0], BOX_MAX[1] - BOX_MIN[1],
                           fc="#3b82f6", alpha=.13, ec="#1d4ed8", lw=2, zorder=3))
    ax.add_patch(Rectangle((st[0] - 50, st[1] - 50), 100, 100, fill=False, ec="#dc2626", ls="--", lw=1.6, zorder=4))
    ax.set_aspect("equal"); ax.grid(alpha=.25)
    ax.set_xlabel("case x (m)  →  east"); ax.set_ylabel("case y (m)  →  north")
    for name, (px, py, hs) in LOC.items():
        ax.plot(px, py, "o", ms=ms, mfc=COL[name], mec="white", mew=1.5, zorder=6)
    return im


fig, ax = plt.subplots(figsize=(8.6, 9.2))
im = base_map(ax, 9)
ax.set_xlim(-225, 390); ax.set_ylim(-320, 410)
ax.text(BOX_MAX[0] + 6, BOX_MIN[1] - 6, "cluster box\n(refined 4 / 2 / 1 m)", color="#1d4ed8", fontsize=10,
        va="top", weight="bold", bbox=WB, zorder=9)
ax.text(st[0] - 50, st[1] + 54, "station ±50 m", color="#dc2626", fontsize=9, bbox=WB, zorder=9)
ax.annotate("upstream line\n(control, outside the box)", (-150, 206), (-215, 290), fontsize=10, color=COL["upstream"],
            arrowprops=dict(arrowstyle="-", color=COL["upstream"]), bbox=WB)
for sx, sy in SKEW:
    ax.plot(sx, sy, "x", color="#ea580c", ms=10, mew=2.5, zorder=7)
ax.text(-215, -305, "×  the 5 badly skewed faces of the coarse mesh\n    (all outside the box, far from every probe)",
        color="#ea580c", fontsize=9.5, bbox=WB)
wind_arrow(ax, -215, 60, 80)
ax.text(-215, 74, "wind from the west, 5 m/s at 10 m", color="#2563eb", fontsize=10, bbox=WB)
ax.set_title("Whole campus, top view (case coordinates)", fontsize=12)
cb = fig.colorbar(im, ax=ax, shrink=.6, pad=.02); cb.set_label("building height above ground (m)")
fig.savefig(HERE / "fig_site_campus.png", dpi=110, bbox_inches="tight"); plt.close(fig)

fig, ax = plt.subplots(figsize=(11, 8.4))
im = base_map(ax, 12)
ax.set_xlim(35, 255); ax.set_ylim(122, 288)
ax.text(BOX_MIN[0] + 2, BOX_MIN[1] - 3, "cluster box: x 50…235, y 145…262, z 12…60 m  (185 × 117 × 48 m)",
        color="#1d4ed8", fontsize=10.5, va="top", weight="bold", bbox=WB, zorder=9)
ax.plot(*st, "*", ms=24, mfc="#dc2626", mec="white", zorder=7)
ax.text(st[0] - 49, st[1] + 51.5, "station position: known only ±50 m (dashed square)", color="#dc2626",
        fontsize=9.5, bbox=WB, zorder=9)
lab = {"station": ("station line\n15, 20, 25, 40 m", (-64, 22)),
       "wake": ("wake line\n2, 5, 10, 20, 40 m", (12, 22)),
       "gap": ("gap line\n2, 5, 10, 20, 40 m", (14, -30))}
for name, (txt, (dx, dy)) in lab.items():
    px, py, _ = LOC[name]
    ax.annotate(txt, (px, py), (px + dx, py + dy), fontsize=10.5, color=COL[name], weight="bold",
                arrowprops=dict(arrowstyle="-", color=COL[name], lw=1.4), zorder=8, bbox=WB)
ax.text(88, 176, "station building (about 10–13 m tall)", fontsize=9.5, color="#111827", rotation=-29,
        bbox=WB, zorder=9)
ax.text(200, 272, "taller block\n(about 16–21 m)", fontsize=9.5, color="#111827", bbox=WB, zorder=9)
wind_arrow(ax, 40, 279, 30)
ax.text(73, 276, "wind", color="#2563eb", fontsize=10.5)
ax.set_title("Zoom on the cluster box: where the 4 probe lines are", fontsize=12)
cb = fig.colorbar(im, ax=ax, shrink=.7, pad=.02); cb.set_label("building height above ground (m)")
fig.savefig(HERE / "fig_site_zoom.png", dpi=110, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ fig_section
i0 = np.abs(y - 206).argmin()
xs = x.copy()
hrow = H[i0 - 1:i0 + 2].max(0)                     # building height along the cut (±2 m)
grow = hg[i0].copy(); ok = np.isfinite(grow)
grow = np.interp(xs, xs[ok], grow[ok])             # fill the ray-cast holes in the ground
fig, ax = plt.subplots(figsize=(15, 5.6))
ax.fill_between(xs, 0, grow, color="#a3a380", lw=0, zorder=2)
top = grow + hrow
ax.fill_between(xs, grow, top, where=hrow > 1, color="#6b7280", step="mid", lw=0, zorder=2)
ax.add_patch(Rectangle((BOX_MIN[0], BOX_MIN[2]), BOX_MAX[0] - BOX_MIN[0], BOX_MAX[2] - BOX_MIN[2],
                       fc="#3b82f6", alpha=.12, ec="#1d4ed8", lw=2, zorder=1))
ax.text(BOX_MIN[0] + 2, BOX_MAX[2] + 0.8, "cluster box (z 12…60 m): 4 m cells (coarse), 2 m (medium), 1 m (fine)",
        color="#1d4ed8", fontsize=10, va="bottom", weight="bold")
ax.text(-110, 50, "outside the box: the SAME mesh in all 3 cases\n(4 m cells near the ground, bigger higher up)",
        fontsize=9.5, color="#374151")
for name, (px, py, hs) in LOC.items():
    g = ground(px, py)
    proj = name == "gap"
    ax.plot([px, px], [g, g + max(hs)], color=COL[name], lw=1.2, ls=":" if proj else "-", zorder=4)
    ax.plot([px] * len(hs), [g + h for h in hs], "o", ms=8, mfc="white" if proj else COL[name], mec=COL[name],
            mew=2, zorder=5)
    for h in hs:
        ax.text(px + 3, g + h, f"{h} m", fontsize=8, color=COL[name], va="center", zorder=6,
                bbox=dict(fc="white", alpha=.7, ec="none", pad=.5))
lt = {"upstream": ("upstream", -150, 63.5, "center"), "station": ("station", 120, 65.5, "right"),
      "wake": ("wake", 142, 65.5, "left"), "gap": ("gap line\n(24 m south of\nthis cut, drawn here)", 182, 47, "left")}
for name, (txt, tx, tz, ha) in lt.items():
    ax.text(tx, tz, txt, ha=ha, fontsize=10, color=COL[name], weight="bold", va="bottom")
ax.text(-212, 13, "ground (terrain)", fontsize=9, color="#3f3f2a")
ax.text(80, 34, "station\nbuilding", fontsize=9, color="#374151", ha="left")
ax.text(186, 40, "taller block", fontsize=9, color="#374151", ha="left")
# log-law inflow profile sketch at the left
zz = np.linspace(0.3, 44, 60); g0 = ground(-200, 206)
uu = 5 * np.log(zz / 0.1) / np.log(10 / 0.1)
ax.plot(-205 + 4.5 * uu, g0 + zz, color="#2563eb", lw=2)
for zq in (2, 5, 10, 20, 40):
    uq = 5 * np.log(zq / 0.1) / np.log(100)
    ax.annotate("", (-205 + 4.5 * uq, g0 + zq), (-205, g0 + zq), arrowprops=dict(arrowstyle="->", color="#2563eb"))
ax.text(-212, g0 + 47, "inflow: 5 m/s at 10 m,\nmore wind higher up", color="#2563eb", fontsize=9.5)
ax.set_xlim(-215, 262); ax.set_ylim(8, 72)
ax.set_xlabel("case x (m), along the wind (west → east)"); ax.set_ylabel("height z (m)")
ax.set_title("Side view along the wind, through the weather station (cut at y = 206 m). "
             "Heights next to the dots = metres above the local ground.  Vertical scale stretched ×3.", fontsize=11)
ax.set_aspect(3.0); ax.grid(alpha=.2)
fig.savefig(HERE / "fig_section.png", dpi=110, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ fig_cells (illustration)
# Made-up but plausible flow behind a 12 m building (side view). Not simulation output.
d = 0.125
X0, X1, Z1, HB = -8.0, 48.0, 28.0, 12.0
xf = np.arange(X0, X1, d) + d / 2; zf = np.arange(0, Z1, d) + d / 2
XF, ZF = np.meshgrid(xf, zf)
U0 = 5 * np.log(np.maximum(ZF, 0.05) / 0.1) / np.log(100); U0 = np.maximum(U0, 0)
hbub = HB * (1 + 0.012 * np.clip(XF, 0, None))
A = 1.35 * np.exp(-(np.clip(ZF, 0, None) / hbub) ** 4) * np.exp(-np.clip(XF, 0, None) / 20.0)
A += 0.25 * np.exp(-((ZF - HB * 1.05) / 2.2) ** 2) * np.exp(-np.clip(XF, 0, None) / 10.0)   # shear layer
UF = U0 * (1 - A) + 0.6 * np.sin(XF / 2.3) * np.exp(-((ZF - 6) / 3) ** 2) * np.exp(-XF / 14.0) * (XF > 0)
inside = (XF < 0) & (ZF < HB)
UF[inside] = np.nan


def coarsen(h):
    n = int(round(h / d))
    nz, nx = UF.shape[0] // n, UF.shape[1] // n
    blk = UF[:nz * n, :nx * n].reshape(nz, n, nx, n)
    with np.errstate(all="ignore"):
        return np.nanmean(blk, axis=(1, 3))


fig = plt.figure(figsize=(15, 7.6))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.5], wspace=.18, hspace=.32)
maps = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
prof = fig.add_subplot(gs[:, 2])
panels = [("Smooth 'real' flow (made up)", None), ("4 m cells: coarse box", 4), ("2 m cells: medium box", 2),
          ("1 m cells: fine box (drone size)", 1)]
kw = dict(cmap="viridis", vmin=-1.5, vmax=6.2)
for ax, (title, h) in zip(maps, panels):
    if h is None:
        ax.pcolormesh(xf - d / 2, zf - d / 2, UF, shading="auto", **kw)
    else:
        C = coarsen(h)
        xe = X0 + h * np.arange(C.shape[1] + 1); ze = h * np.arange(C.shape[0] + 1)
        ax.pcolormesh(xe, ze, C, shading="flat", edgecolors=(1, 1, 1, .35 if h > 1 else .18),
                      linewidth=.5 if h > 1 else .3, **kw)
    ax.add_patch(Rectangle((X0, 0), -X0, HB, fc="#6b7280", ec="k", zorder=3))
    ax.add_patch(Rectangle((14, 7), 1, 1, fc="none", ec="#ff2d55", lw=2.2, zorder=5))
    ax.annotate("drone\n1 × 1 m", (15, 8), (21, 17), color="#ff2d55", fontsize=10, weight="bold",
                arrowprops=dict(arrowstyle="-", color="#ff2d55", lw=1.5), zorder=6)
    ax.axvline(12, color="w", ls="--", lw=1)
    ax.set_xlim(X0, X1 - 8); ax.set_ylim(0, Z1); ax.set_aspect("equal")
    ax.set_title(title, fontsize=11)
for ax in maps[2:]:
    ax.set_xlabel("distance behind the building (m)")
for ax in maps[::2]:
    ax.set_ylabel("height (m)")
maps[0].text(-7.5, 13, "building", fontsize=9, color="#374151")
maps[0].text(12.5, 25.5, "profile line →", color="w", fontsize=8.5)
# vertical profile at x = 12 m
ax = prof
j = np.abs(xf - 12.0).argmin()
ax.plot(UF[:, j], zf, color="k", lw=2, label="smooth flow")
for h, c in ((4, "#9ca3af"), (2, "#2563eb"), (1, "#dc2626")):
    C = coarsen(h); jj = int((12.0 - X0) // h)
    ze = h * np.arange(C.shape[0] + 1)
    ax.stairs(C[:, jj], ze, orientation="horizontal", color=c, lw=2, label=f"{h} m cells", baseline=None)
ax.axvline(0, color="#999", lw=.8)
ax.set_ylim(0, Z1); ax.set_xlabel("wind speed (m/s)"); ax.set_title("Same line, seen by each mesh", fontsize=11)
ax.legend(fontsize=9, loc="upper left"); ax.grid(alpha=.3); ax.set_ylabel("height (m)")
fig.suptitle("ILLUSTRATION with a made-up flow (not our simulation): each cell keeps ONE average value.\n"
             "Big cells blur the slow, reversed air behind the building; the drone (red square) is smaller than a 4 m cell.",
             fontsize=12, y=.99)
fig.colorbar(maps[0].collections[0], ax=maps, shrink=.8, pad=.02,
             label="wind speed along x (m/s); below 0 = air flowing back")
fig.savefig(HERE / "fig_cells.png", dpi=110, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ fig_expected (sketch)
with plt.xkcd(scale=.6, length=120, randomness=2):
    plt.rcParams["font.family"] = ["DejaVu Sans"]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    c3 = {"coarse": "#9ca3af", "medium": "#2563eb", "fine": "#dc2626"}
    z = np.linspace(1.5, 40, 50)
    up = 5 * np.log(z / 0.1) / np.log(100)
    ax = axs[0]
    for k, (c, col) in enumerate(c3.items()):
        ax.plot(up + 0.03 * (k - 1), z, color=col, lw=2.5, label=c)
    ax.set_title("UPSTREAM: 3 curves on top of each other"); ax.set_xlabel("wind speed"); ax.set_ylabel("height above ground")
    ax.legend(loc="upper left"); ax.set_xticks([]); ax.set_yticks([0, 10, 20, 30, 40])
    ax = axs[1]
    base = 5.6 * (1 - np.exp(-z / 13)) - 0.4
    for c, col, s in (("coarse", c3["coarse"], 1.0), ("medium", c3["medium"], 0.35), ("fine", c3["fine"], 0.0)):
        ax.plot(base + s * 1.4 * np.exp(-z / 6), z, color=col, lw=2.5, label=c)
    ax.annotate("medium and fine close:\nthe answer has settled", (2.2, 7), (3.2, 16), fontsize=10,
                arrowprops=dict(arrowstyle="->"))
    ax.annotate("coarse differs most\nnear the ground", (1.7, 2.6), (2.8, 3.0), fontsize=10,
                arrowprops=dict(arrowstyle="->"))
    ax.set_title("WAKE: curves split near the ground"); ax.set_xlabel("wind speed"); ax.set_xticks([])
    ax.set_yticks([0, 10, 20, 30, 40])
    ax = axs[2]
    it = np.arange(0, 1000)
    for c, col, v in (("coarse", c3["coarse"], 3.9), ("medium", c3["medium"], 3.45), ("fine", c3["fine"], 3.35)):
        ax.plot(it, v + 1.8 * np.exp(-it / 90) * np.cos(it / 40), color=col, lw=2.2, label=c)
    ax.annotate("flat at the end\n= converged", (900, 3.4), (560, 4.9), fontsize=10, arrowprops=dict(arrowstyle="->"))
    ax.set_title("PROBE vs ITERATION: settles, then flat"); ax.set_xlabel("iteration"); ax.set_ylabel("wind speed")
    ax.set_yticks([])
    fig.suptitle("SKETCH of what I expect (invented curves, NOT results)", fontsize=13)
    fig.tight_layout()
    fig.savefig(HERE / "fig_expected.png", dpi=100, bbox_inches="tight"); plt.close(fig)
print("wrote", *sorted(p.name for p in HERE.glob("fig_*.png")))
