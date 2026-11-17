#!/usr/bin/env python3
"""Pictures of the mesh-study meshes, their refinement box and (when solved) the wind.
Run from anywhere:   .venv/bin/python DLR/meshStudy/make_figures.py  [--force]
Works with partial results: it skips cases that are not ready and says so.

Writes DLR/results/mesh_fig/ (other documents link to these exact names):
    box_topview.png                          buildings, refinement box, probes, station +-50 m, wind, scale
    mesh_slice_{coarse,medium,fine}.png      horizontal cut z = 24 m (~7 m above ground), same window
    mesh_vertical_fine.png                   vertical cut y = 206 m (also _coarse/_medium if available)
    wind_slice_{coarse,medium,fine}.png      wind speed on the z = 24 m cut, same window and colour range

How: make_figures_pv.py (run with ParaView's pvpython, one case at a time, never in parallel) cuts
the mesh and saves the cut cells to DLR/meshStudy/.figcache/<case>_cuts.npz; this script draws them.
A cut is redone only when its case changed (or with --force). ParaView reads <case>/case.foam;
the empty file is created if missing (the only thing ever written inside a case folder).

Readiness:  mesh pictures  -> constant/polyMesh/points and log.checkMesh exist
            wind pictures  -> run finished: 1000/U exists and processor* folders are gone
"""
import csv, os, subprocess, sys, time
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle, FancyArrow, Patch
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results" / "mesh_fig"
CACHE = HERE / ".figcache"
CASES = ["coarse", "medium", "fine"]
BOX_CELL = {"coarse": 4, "medium": 2, "fine": 1}            # m, inside the cluster box
BOX_MIN, BOX_MAX = (50.0, 145.0, 12.0), (235.0, 262.0, 60.0)  # = setup.py
STATION = (122.0, 206.0)
END_TIME = "1000"
H_WIN = (20.0, 265.0, 115.0, 290.0)    # x0 x1 y0 y1 of the horizontal pictures (= make_figures_pv.py)
V_WIN = (20.0, 265.0, 10.0, 80.0)      # x0 x1 z0 z1 of the vertical picture
Z_CUT, Y_CUT = 24.0, 206.0
VMAX = float(os.environ.get("WIND_VMAX", 6.0))   # fixed colour range 0..VMAX m/s for all wind pictures
PVPYTHON = os.environ.get("PVPYTHON", "/Applications/ParaView-6.2.0.app/Contents/bin/pvpython")
FORCE = "--force" in sys.argv

SOLID = "#9a9a9a"        # holes in a cut = buildings / ground
FACE = "#f3f5f9"         # mesh cell fill
EDGE = "#1d2533"         # mesh cell edges
BOXC = "#d62828"         # refinement box
STC = "#e07a00"          # station

IDX = list(csv.DictReader(open(HERE / "probes_index.csv")))
PROBE_XY = {}
for r in IDX:
    PROBE_XY.setdefault(r["location"], (float(r["x"]), float(r["y"]), []))[2].append(float(r["height_above_ground"]))


# ---------------------------------------------------------------- case status
def mesh_file(c):
    p = HERE / c / "constant" / "polyMesh"
    return next((f for f in (p / "points", p / "points.gz") if f.exists()), None)

def finished(c):
    d = HERE / c / END_TIME
    return ((d / "U").exists() or (d / "U.gz").exists()) and not list((HERE / c).glob("processor*"))

def mesh_ready(c):
    return mesh_file(c) is not None and ((HERE / c / "log.checkMesh").exists() or finished(c))

def field_file(c):
    d = HERE / c / END_TIME
    return d / "U" if (d / "U").exists() else d / "U.gz"


# ---------------------------------------------------------------- ParaView part
def extract(c):
    """Run the pvpython cutter for case c if its cache is missing or stale. Returns the npz path or None."""
    npz = CACHE / f"{c}_cuts.npz"
    want_fields = finished(c)
    src_time = max(mesh_file(c).stat().st_mtime, (HERE / c / "log.checkMesh").stat().st_mtime
                   if (HERE / c / "log.checkMesh").exists() else 0)
    if want_fields:
        src_time = max(src_time, field_file(c).stat().st_mtime)
    if npz.exists() and not FORCE and npz.stat().st_mtime > src_time:
        with np.load(npz) as z:
            if bool(z["has_fields"]) or not want_fields:
                print(f"  {c}: cut cache up to date ({npz.name})")
                return npz
    if not Path(PVPYTHON).exists():
        print(f"  {c}: pvpython not found at {PVPYTHON} (set PVPYTHON=...); using old cache if any")
        return npz if npz.exists() else None
    foam = HERE / c / "case.foam"
    if not foam.exists():
        foam.touch()                       # empty file the OpenFOAM reader needs (allowed)
        print(f"  {c}: created empty {foam.relative_to(HERE.parent.parent)}")
    CACHE.mkdir(exist_ok=True)
    tmp = CACHE / f"{c}_cuts.tmp.npz"
    cmd = [PVPYTHON, str(HERE / "make_figures_pv.py"), str(HERE / c), str(tmp)] + (["--fields"] if want_fields else [])
    print(f"  {c}: cutting with pvpython ({'mesh + wind' if want_fields else 'mesh only'}) ...", flush=True)
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except subprocess.TimeoutExpired:
        print(f"  {c}: pvpython took more than 30 min, stopped; using old cache if any")
        return npz if npz.exists() else None
    if p.returncode != 0 or not tmp.exists():
        print(f"  {c}: pvpython FAILED (exit {p.returncode}):\n" + (p.stdout + p.stderr)[-2000:])
        return npz if npz.exists() else None
    tmp.replace(npz)
    print(f"  {c}: done in {time.time() - t0:.0f} s")
    return npz


# ---------------------------------------------------------------- drawing helpers
def wrote(f):
    print("  wrote", os.path.relpath(f, HERE.parent))

def polys(z, key):
    xy, off = z[key + "_xy"], z[key + "_off"]
    return [xy[off[i]:off[i + 1]] for i in range(len(off) - 1)]

def box_outline(ax, vertical=False):
    if vertical:
        r = Rectangle((BOX_MIN[0], BOX_MIN[2]), BOX_MAX[0] - BOX_MIN[0], BOX_MAX[2] - BOX_MIN[2])
    else:
        r = Rectangle(BOX_MIN[:2], BOX_MAX[0] - BOX_MIN[0], BOX_MAX[1] - BOX_MIN[1])
    r.set(fill=False, ec=BOXC, lw=1.8, ls=(0, (5, 3)), zorder=5)
    ax.add_patch(r)

def probe_marks(ax, fs):
    """Probe locations inside the horizontal window, with short labels."""
    for loc, (x, y, _) in PROBE_XY.items():
        if not (H_WIN[0] < x < H_WIN[1] and H_WIN[2] < y < H_WIN[3]):
            continue
        if loc == "station":
            ax.plot(x, y, marker="*", ms=15, mfc=STC, mec="k", mew=0.8, zorder=7, ls="none")
        else:
            ax.plot(x, y, marker="o", ms=7, mfc="white", mec="k", mew=1.2, zorder=7, ls="none")
        dx, dy = {"station": (-8, 9), "wake": (7, 7), "gap": (7, -11)}.get(loc, (6, 6))
        ax.annotate(loc, (x, y), (x + dx, y + dy), fontsize=fs, zorder=8, ha="right" if dx < 0 else "left",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

def style_axes(ax, fs, vertical=False):
    ax.set_aspect("equal")
    ax.set_xlabel("x (m), wind blows towards +x", fontsize=fs)
    ax.set_ylabel("z (m)" if vertical else "y (m)", fontsize=fs)
    ax.tick_params(labelsize=fs - 1)
    ax.set_facecolor(SOLID)


def mesh_slice(c, z):
    fs = 13
    fig, ax = plt.subplots(figsize=(6.4, 5.0), dpi=250, layout="constrained")
    ax.add_collection(PolyCollection(polys(z, "h"), facecolors=FACE, edgecolors=EDGE, linewidths=0.22))
    box_outline(ax); probe_marks(ax, fs - 2)
    ax.set_xlim(H_WIN[:2]); ax.set_ylim(H_WIN[2:]); style_axes(ax, fs)
    ax.set_title(f"{c}: {BOX_CELL[c]} m cells in the box   ({int(z['ncells']):,} cells)\n"
                 f"horizontal cut z = {Z_CUT:g} m (about 7 m above ground)", fontsize=fs)
    ax.legend(handles=[Line2D([], [], color=BOXC, ls=(0, (5, 3)), lw=1.8, label="refinement box"),
                       Patch(fc=SOLID, label="solid (building)")],
              loc="lower right", fontsize=fs - 3, framealpha=0.92)
    f = OUT / f"mesh_slice_{c}.png"; fig.savefig(f); plt.close(fig); wrote(f)


def mesh_vertical(c, z):
    fs = 12
    fig, ax = plt.subplots(figsize=(10.0, 3.9), dpi=160, layout="constrained")
    ax.add_collection(PolyCollection(polys(z, "v"), facecolors=FACE, edgecolors=EDGE, linewidths=0.3))
    box_outline(ax, vertical=True)
    ax.axhline(Z_CUT, color="#2563eb", lw=1.2, ls=":", zorder=6)
    ax.text(V_WIN[1] - 2, Z_CUT + 1, f"horizontal cut z = {Z_CUT:g} m", color="#2563eb", fontsize=fs - 2,
            ha="right", va="bottom", zorder=8, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    st = [float(r["z"]) for r in IDX if r["location"] == "station"]
    ax.plot([STATION[0]] * len(st), st, marker="*", ms=12, mfc=STC, mec="k", mew=0.7, ls="none", zorder=7)
    ax.annotate("station probes", (STATION[0], st[-1]), (STATION[0] + 6, st[-1] + 6), fontsize=fs - 2, zorder=8,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax.set_xlim(V_WIN[:2]); ax.set_ylim(V_WIN[2:]); style_axes(ax, fs, vertical=True)
    ax.set_title(f"{c} mesh: vertical cut y = {Y_CUT:g} m through the weather station "
                 f"({BOX_CELL[c]} m cells in the box, {int(z['ncells']):,} cells)", fontsize=fs)
    ax.legend(handles=[Line2D([], [], color=BOXC, ls=(0, (5, 3)), lw=1.8, label="refinement box (z 12–60 m)"),
                       Patch(fc=SOLID, label="solid (ground, building)")],
              loc="upper left", fontsize=fs - 3, framealpha=0.92)
    f = OUT / f"mesh_vertical_{c}.png"; fig.savefig(f); plt.close(fig); wrote(f)


def wind_slice(c, z):
    fs = 13
    cmap, norm = plt.get_cmap("viridis"), plt.Normalize(0.0, VMAX)
    fig, ax = plt.subplots(figsize=(6.4, 4.9), dpi=250, layout="constrained")
    pc = PolyCollection(polys(z, "h"), array=z["h_speed"], cmap=cmap, norm=norm, edgecolors="face", linewidths=0.15)
    ax.add_collection(pc)
    box_outline(ax); probe_marks(ax, fs - 2)
    ax.set_xlim(H_WIN[:2]); ax.set_ylim(H_WIN[2:]); style_axes(ax, fs)
    cb = fig.colorbar(pc, ax=ax, shrink=0.85, pad=0.02, extend="max")
    cb.set_label("wind speed |U| (m/s)", fontsize=fs - 1); cb.ax.tick_params(labelsize=fs - 2)
    ax.set_title(f"{c} ({BOX_CELL[c]} m cells in the box): wind speed\nat z = {Z_CUT:g} m (about 7 m above ground), "
                 f"iteration {z['time']:g}", fontsize=fs)
    f = OUT / f"wind_slice_{c}.png"; fig.savefig(f); plt.close(fig); wrote(f)


def box_topview():
    s = np.load(HERE / "site_map.npz")
    x, y, h = s["x"], s["y"], s["hb"] - s["hg"]
    # the ray-cast map has pinholes (rays slipping between STL triangles): fill a hole when most neighbours are building
    hp = np.pad(h, 1, constant_values=np.nan)
    nb = np.stack([hp[1 + i:hp.shape[0] - 1 + i, 1 + j:hp.shape[1] - 1 + j]
                   for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0)])
    fill = np.isnan(h) & (np.sum(nb > 0.5, axis=0) >= 5)
    h = h.copy()
    h[fill] = np.nanmedian(np.where(nb > 0.5, nb, np.nan)[:, fill], axis=0)
    win = (-205.0, 300.0, 60.0, 350.0)
    fs = 11
    fig, ax = plt.subplots(figsize=(10.0, 6.3), dpi=160, layout="constrained")
    ix = (x >= win[0] - 4) & (x <= win[1] + 4); iy = (y >= win[2] - 4) & (y <= win[3] + 4)
    hh = np.ma.masked_invalid(np.where(h[np.ix_(iy, ix)] > 0.5, h[np.ix_(iy, ix)], np.nan))
    pm = ax.pcolormesh(x[ix], y[iy], hh, cmap="Greys", vmin=0, vmax=30, shading="auto", zorder=1)
    cb = fig.colorbar(pm, ax=ax, shrink=0.8, pad=0.01)
    cb.set_label("building height above ground (m)", fontsize=fs)
    # windows of the other pictures
    ax.add_patch(Rectangle((H_WIN[0], H_WIN[2]), H_WIN[1] - H_WIN[0], H_WIN[3] - H_WIN[2], fill=False,
                           ec="#6b7280", lw=1.0, ls=":", zorder=3))
    ax.text(H_WIN[0] + 2, H_WIN[3] + 3, "window of the mesh and wind pictures", fontsize=fs - 2, color="#4b5563")
    ax.plot([V_WIN[0], V_WIN[1]], [Y_CUT, Y_CUT], color="#2563eb", lw=1.0, ls="-.", zorder=3)
    ax.text(V_WIN[1] + 3, Y_CUT, f"vertical cut\ny = {Y_CUT:g} m", fontsize=fs - 2, color="#2563eb", va="center")
    # refinement box
    box_outline(ax)
    ax.text(BOX_MIN[0] + 2, BOX_MIN[1] + 3, "refinement box, z 12–60 m\n4 / 2 / 1 m cells (coarse / medium / fine)",
            fontsize=fs - 1, color=BOXC, va="bottom", zorder=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    # station and its +-50 m uncertainty
    ax.add_patch(Rectangle((STATION[0] - 50, STATION[1] - 50), 100, 100, fill=False, ec=STC, lw=1.4, ls="--", zorder=4))
    ax.text(STATION[0] - 48, STATION[1] + 48, "station position\n±50 m", fontsize=fs - 2, color=STC, zorder=8, va="top",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    # probes
    lab = {"upstream": ((-150, 236), "upstream (control)"), "gap": ((182, 168), "gap"),
           "wake": ((150, 222), "wake"), "station": ((76, 192), "station (above roof)")}
    for loc, (px, py, hs) in PROBE_XY.items():
        if loc == "station":
            ax.plot(px, py, marker="*", ms=17, mfc=STC, mec="k", mew=0.8, ls="none", zorder=7)
        else:
            ax.plot(px, py, marker="o", ms=8, mfc="white", mec="k", mew=1.3, ls="none", zorder=7)
        (tx, ty), name = lab[loc]
        ax.annotate(f"{name}\n{min(hs):g}–{max(hs):g} m above ground", (px, py), (tx, ty), fontsize=fs - 1,
                    ha="right" if tx < px - 5 else ("center" if abs(tx - px) <= 5 else "left"), va="center", zorder=8,
                    arrowprops=dict(arrowstyle="-", color="k", lw=0.8, shrinkB=5),
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6", lw=0.5))
    # wind arrow and scale bar
    ax.add_patch(FancyArrow(-195, 325, 60, 0, width=4, head_width=14, head_length=14, color="#1f4e9c", zorder=8))
    ax.text(-195, 308, "wind from the west (270°)\n5 m/s at 10 m above ground", fontsize=fs - 1, color="#1f4e9c",
            va="top", zorder=8)
    ax.plot([190, 290], [72, 72], color="k", lw=3, solid_capstyle="butt", zorder=8)
    for xx in (190, 240, 290):
        ax.plot([xx, xx], [69, 75], color="k", lw=1.2, zorder=8)
    ax.text(240, 77, "100 m", ha="center", va="bottom", fontsize=fs, zorder=8,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))
    ax.set_xlim(win[:2]); ax.set_ylim(win[2:]); ax.set_aspect("equal")
    ax.set_xlabel("x (m, case coordinates)", fontsize=fs); ax.set_ylabel("y (m)", fontsize=fs)
    ax.grid(alpha=0.25, zorder=0)
    ax.set_title("DLR Oberpfaffenhofen campus: refinement box, probe locations and weather station (top view)",
                 fontsize=fs + 1)
    f = OUT / "box_topview.png"; fig.savefig(f); plt.close(fig); wrote(f)


# ---------------------------------------------------------------- main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("box top view:")
    box_topview()
    done, skipped = [], []
    for c in CASES:
        print(f"{c}:")
        if not mesh_ready(c):
            why = "no mesh yet" if mesh_file(c) is None else "mesh not checked yet (no log.checkMesh)"
            print(f"  SKIPPED {c}: {why}"); skipped.append(f"{c} (all pictures: {why})"); continue
        npz = extract(c)                     # one pvpython process at a time
        if npz is None:
            skipped.append(f"{c} (cut failed)"); continue
        with np.load(npz) as z:
            mesh_slice(c, z)
            mesh_vertical(c, z)
            if bool(z["has_fields"]):
                wind_slice(c, z); done.append(f"{c} (mesh + wind)")
            else:
                print(f"  SKIPPED wind_slice_{c}.png: run not finished (no {END_TIME}/U yet, or not reconstructed)")
                skipped.append(f"{c} (wind picture: run not finished)"); done.append(f"{c} (mesh only)")
    print("\nSummary: rendered " + (", ".join(done) or "nothing"))
    if skipped:
        print("Skipped: " + "; ".join(skipped))
    print(f"Pictures in {OUT}")


if __name__ == "__main__":
    main()
