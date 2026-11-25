#!/usr/bin/env python3
"""Figures for the wind atlas (repo root):  .venv/bin/python DLR/atlas/make_figures.py

results/atlas_maps_<h>m.png   one small map per solved direction: wind speed at h m above ground
                              divided by the free-stream speed at that height (1 = undisturbed)
results/atlas_reynolds.png    10 m/s run vs 2 x the 5 m/s run (240 deg): scatter + error map
results/atlas_loo.png         leave-one-out error of direction blending vs nearest vs log-law
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from wind_atlas import WindAtlas  # noqa: E402
import validate_atlas as V  # noqa: E402

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "figure.dpi": 130})


def arrow(ax, theta, x=0.1, y=0.88, L=0.12):
    """Arrow showing where the air goes (wind FROM theta)."""
    f = np.array([-np.sin(np.radians(theta)), -np.cos(np.radians(theta))])
    ax.annotate("", xy=(x + L / 2 * f[0], y + L / 2 * f[1]), xytext=(x - L / 2 * f[0], y - L / 2 * f[1]),
                xycoords="axes fraction", arrowprops=dict(arrowstyle="-|>", color="k", lw=1.4))


def maps(atlas, h=30.0):
    V.ATLAS = atlas
    dirs = [d for d in atlas.directions if round(d) % 30 == 0]   # 255 is validation-only
    n = len(dirs)
    cols = 4 if n > 4 else n
    rows = int(np.ceil(n / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(3.0 * cols + 0.8, 2.75 * rows + 0.5), squeeze=False, constrained_layout=True)
    bld = (atlas.ground() - atlas.terrain()) > 3.0
    ustar = V.KAPPA * 5.0 / np.log((V.ZREF + V.Z0) / V.Z0)
    free = ustar / V.KAPPA * np.log((h + V.Z0) / V.Z0)
    for ax, th in zip(axs.flat, dirs):
        F = V.at_heights(atlas, atlas.field(th, 5.0))[h] if h in V.HEIGHTS else None
        if F is None:
            V.HEIGHTS = tuple(sorted(set(V.HEIGHTS) | {h}))
            F = V.at_heights(atlas, atlas.field(th, 5.0))[h]
        s = np.linalg.norm(F, axis=-1) / free
        im = ax.pcolormesh(atlas.x, atlas.y, s.T, cmap="RdBu_r", shading="auto",
                           norm=TwoSlopeNorm(vmin=0.3, vcenter=1.0, vmax=1.3))
        ax.contour(atlas.x, atlas.y, bld.T.astype(float), levels=[0.5], colors="k", linewidths=0.4)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"from {int(th)}°")
        arrow(ax, th)
    for ax in list(axs.flat)[n:]:
        ax.axis("off")
    cb = fig.colorbar(im, ax=axs, shrink=0.7, pad=0.01)
    cb.set_label(f"|U| / free-stream |U| at {h:g} m above ground")
    fig.suptitle(f"DLR campus wind atlas: speed-up (red) and shelter (blue) {h:g} m above ground, "
                 f"steady RANS k-ε, 5 m/s at 10 m")
    f = OUT / f"atlas_maps_{h:g}m.png"
    fig.savefig(f, bbox_inches="tight"); plt.close(fig)
    return f


def reynolds(atlas):
    if (240.0, 10.0) not in atlas.extra or 240.0 not in atlas.runs:
        return None
    V.ATLAS = atlas
    U10 = V.at_heights(atlas, atlas.extra[(240.0, 10.0)]["U"])
    U5 = V.at_heights(atlas, atlas.runs[240.0]["U"])
    mask = V.campus_mask(atlas)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.3))
    for h, c in zip(V.HEIGHTS, plt.cm.viridis(np.linspace(0, 0.9, len(V.HEIGHTS)))):
        a = np.linalg.norm(U5[h], axis=-1)[mask]; b = np.linalg.norm(U10[h], axis=-1)[mask]
        ok = np.isfinite(a) & np.isfinite(b)
        axs[0].scatter(a[ok][::7], b[ok][::7], s=2, color=c, label=f"{h:g} m", alpha=0.5)
    axs[0].plot([0, 1.6], [0, 1.6], "k--", lw=0.8)
    axs[0].set_xlabel("|U| / Uref, 5 m/s run"); axs[0].set_ylabel("|U| / Uref, 10 m/s run")
    axs[0].legend(markerscale=4, title="above ground"); axs[0].set_title("Same pattern at both speeds?")
    e = np.linalg.norm(U5[20.0] - U10[20.0], axis=-1)
    im = axs[1].pcolormesh(atlas.x, atlas.y, e.T, cmap="magma_r", vmin=0, vmax=0.05, shading="auto")
    axs[1].set_aspect("equal"); axs[1].set_title("|ΔU| / Uref at 20 m above ground")
    fig.colorbar(im, ax=axs[1])
    f = OUT / "atlas_reynolds.png"
    fig.tight_layout(); fig.savefig(f, bbox_inches="tight"); plt.close(fig)
    return f


def loo_plot():
    rep = json.loads((OUT / "validation.json").read_text())
    loo = rep.get("leave_one_out", {})
    if not loo:
        return None
    hs = ["10m", "20m", "40m", "60m"]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    models = [("loglaw_no_buildings", "log-law, no buildings (no CFD)", "#b8b8b8"),
              ("nearest_unrotated_60deg", "nearest solved run, as is", "#e8a35c"),
              ("nearest_rotated_60deg", "nearest solved run, rotated", "#c96f1f"),
              ("blend_linear_60deg", "linear blend of raw vectors", "#8fb3dd"),
              ("blend_rotate_60deg", "rotate, then blend (atlas default)", "#2c6fbb")]
    w = 0.8 / len(models)
    for i, (key, lab, col) in enumerate(models):
        vals = [np.mean([loo[t][key][h]["rmse"] for t in loo]) for h in hs]
        ax.bar(np.arange(len(hs)) + (i - (len(models) - 1) / 2) * w, vals, w, label=lab, color=col)
    ax.set_xticks(range(len(hs))); ax.set_xticklabels([h.replace("m", " m") for h in hs])
    ax.set_xlabel("height above ground"); ax.set_ylabel("RMS vector error / Uref")
    ax.set_title(f"Predicting a hidden direction (mean over {len(loo)} directions: "
                 f"{', '.join(str(t) for t in loo)}°) — campus columns")
    b30 = rep.get("blend_30deg_255", {})
    if b30:
        ax.plot(np.arange(len(hs)) + 2 * w, [b30["blend_rotate_30deg"][h]["rmse"] for h in hs], "k*", ms=10,
                label="rotate+blend at the real 30° gap (255° held out)")
    ax.legend(fontsize=7.5, ncol=2); ax.grid(axis="y", alpha=0.3)
    f = OUT / "atlas_loo.png"
    fig.tight_layout(); fig.savefig(f, bbox_inches="tight"); plt.close(fig)
    return f


if __name__ == "__main__":
    A = WindAtlas()
    for h in (10.0, 30.0):
        print(maps(A, h))
    print(reynolds(A))
    print(loo_plot())
