"""
sionna_cost_map.py
==================
Converts a Sionna RF simulation (radio maps) into a 3D cost map for the drone
path planner, so the planner prefers voxels where RF localization is accurate.

Pipeline
--------
Sionna radio maps (path gain per transmitter, one horizontal slice per altitude)
  -> SNR per transmitter per voxel (link budget)
  -> range-measurement Fisher information with real geometry (unit vectors to
     each transmitter) -> position CRLB (RMSE bound, metres) per voxel
  -> normalised cost multiplier >= 1 per voxel  (inf where blocked)

Because every cost is >= 1, a straight-line Euclidean heuristic in A* stays
admissible.

Assumption: the Sionna transmitters are the RF anchors/sources, and the metric is
the position CRLB from range (TOA-style) measurements to the usable ones.

Typical use
-----------
    # 1) Inside your Sionna environment (once), export the radio maps:
    from sionna_cost_map import RadioMapStack
    stack = RadioMapStack.from_sionna(scene, altitudes=[10, 20, 30, 40, 50],
                                      center_xy=(0, 0), size_xy=(400, 400),
                                      cell_size=(2.0, 2.0))
    stack.to_npz("radio_maps.npz")

    # 2) In the planner:
    from sionna_cost_map import RadioMapStack, GridSpec, RFConfig, build_cost_map
    stack = RadioMapStack.from_npz("radio_maps.npz")
    grid  = GridSpec(origin=(-200, -200, 10), shape=(200, 200, 9), resolution=(2, 2, 5))
    cmap  = build_cost_map(stack, grid, RFConfig(tx_power_dbm=23, bandwidth_hz=20e6),
                           obstacles=my_obstacle_mask)      # bool (nx, ny, nz), optional
    step_cost = cmap.make_cost_fn(weight=1.0)               # f(idx_a, idx_b) -> float

Grid index order is (x, y, z). If your planner indexes [z][y][x], transpose with
cmap.cost.transpose(2, 1, 0).

Requires: numpy, scipy.  Optional: matplotlib (plots), sionna (live export).
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Callable, Optional, Tuple

import numpy as np

C_LIGHT = 299_792_458.0
GAIN_FLOOR_DB = -200.0


# --------------------------------------------------------------------------- #
# Planner grid description
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GridSpec:
    """Planner voxel grid. `origin` is the world position (m) of voxel (0,0,0)."""
    origin: Tuple[float, float, float]
    shape: Tuple[int, int, int]              # (nx, ny, nz)
    resolution: Tuple[float, float, float]   # (dx, dy, dz) in metres

    def axes(self):
        return tuple(self.origin[a] + np.arange(self.shape[a]) * self.resolution[a]
                     for a in range(3))

    def centers(self) -> np.ndarray:
        X, Y, Z = np.meshgrid(*self.axes(), indexing="ij")
        return np.stack([X, Y, Z], axis=-1)  # (nx, ny, nz, 3)

    @staticmethod
    def from_bounds(min_corner, max_corner, resolution) -> "GridSpec":
        res = np.broadcast_to(np.asarray(resolution, float), (3,))
        lo, hi = np.asarray(min_corner, float), np.asarray(max_corner, float)
        shape = tuple(int(np.floor((hi - lo)[a] / res[a])) + 1 for a in range(3))
        return GridSpec(tuple(lo), shape, tuple(res))


# --------------------------------------------------------------------------- #
# Sionna radio-map container
# --------------------------------------------------------------------------- #
def _to_np(x) -> np.ndarray:
    """Convert a Sionna / DrJit / TensorFlow tensor to numpy."""
    if hasattr(x, "numpy"):
        x = x.numpy()
    return np.array(x)


@dataclass
class RadioMapStack:
    """
    path_gain : (n_tx, n_z, n_y, n_x) LINEAR path gain (not dB)
    x, y, z   : 1D cell-centre coordinates (m) of the radio-map layers
    tx_pos    : (n_tx, 3) transmitter positions (m)
    """
    path_gain: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    tx_pos: np.ndarray

    # ---- persistence ----
    def to_npz(self, path: str) -> None:
        np.savez_compressed(path, path_gain=self.path_gain, x=self.x, y=self.y,
                            z=self.z, tx_pos=self.tx_pos)

    @classmethod
    def from_npz(cls, path: str) -> "RadioMapStack":
        d = np.load(path)
        pg = d["path_gain"]
        if pg.ndim == 3:                       # (n_tx, n_y, n_x): single altitude
            pg = pg[:, None]
        return cls(pg, d["x"], d["y"], np.atleast_1d(d["z"]), d["tx_pos"])

    # ---- live export from a Sionna scene ----
    @classmethod
    def from_sionna(cls, scene, altitudes, center_xy, size_xy, cell_size=(2.0, 2.0),
                    samples_per_tx: int = 1_000_000, max_depth: int = 5,
                    **solver_kwargs) -> "RadioMapStack":
        """
        Run one horizontal radio map per altitude. Supports Sionna RT >= 1.0
        (RadioMapSolver) and the legacy 0.x API (scene.coverage_map).
        Check argument names against your installed version.
        """
        cx, cy = center_xy
        layers = []
        try:
            from sionna.rt import RadioMapSolver
            solver = RadioMapSolver()
            for z in altitudes:
                rm = solver(scene, center=[cx, cy, float(z)], orientation=[0, 0, 0],
                            size=list(size_xy), cell_size=list(cell_size),
                            samples_per_tx=samples_per_tx, max_depth=max_depth,
                            **solver_kwargs)
                layers.append(_to_np(rm.path_gain))
        except ImportError:                    # legacy Sionna 0.x
            for z in altitudes:
                cm = scene.coverage_map(cm_center=[cx, cy, float(z)],
                                        cm_orientation=[0, 0, 0],
                                        cm_size=list(size_xy),
                                        cm_cell_size=list(cell_size),
                                        num_samples=samples_per_tx,
                                        max_depth=max_depth, **solver_kwargs)
                layers.append(_to_np(cm.path_gain))

        pg = np.stack(layers, axis=1)          # (n_tx, n_z, n_y, n_x)
        n_y, n_x = pg.shape[-2:]
        x = cx - size_xy[0] / 2 + (np.arange(n_x) + 0.5) * cell_size[0]
        y = cy - size_xy[1] / 2 + (np.arange(n_y) + 0.5) * cell_size[1]
        tx_pos = np.array([_to_np(tx.position).reshape(-1)[:3]
                           for tx in scene.transmitters.values()], dtype=float)
        return cls(pg.astype(float), x, y, np.asarray(altitudes, float), tx_pos)


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
@dataclass
class RFConfig:
    tx_power_dbm: float | tuple = 20.0   # scalar or one value per transmitter
    bandwidth_hz: float = 20e6
    noise_figure_db: float = 6.0
    rms_bandwidth_hz: Optional[float] = None   # default: bandwidth / sqrt(12)
    min_snr_db: float = 0.0              # a link below this is not usable
    min_anchors: int = 3                 # fewer usable links -> unlocalizable
    dims: str = "3d"                     # "3d" or "2d" (horizontal only)
    sigma_floor_m: float = 0.05          # lower bound on single-range std dev
    unlocalizable_rmse_m: float = 1e3
    crlb_good_m: float = 1.0             # RMSE at/below this -> no penalty
    crlb_bad_m: float = 50.0             # RMSE at/above this -> max penalty
    loc_weight: float = 4.0              # cost = 1 + loc_weight * normalised badness
    reg: float = 1e-6                    # relative FIM regularisation


# --------------------------------------------------------------------------- #
# Resampling radio maps onto the planner grid
# --------------------------------------------------------------------------- #
def _resample(vol_db: np.ndarray, src_axes, grid: GridSpec) -> np.ndarray:
    """Trilinear resample a (nz, ny, nx) volume to the planner grid -> (nx, ny, nz)."""
    from scipy.ndimage import map_coordinates

    def frac_idx(target, src):
        if len(src) == 1:
            return np.zeros_like(target, dtype=float)
        step = src[1] - src[0]
        return np.clip((target - src[0]) / step, 0, len(src) - 1)

    xs, ys, zs = src_axes
    gx, gy, gz = grid.axes()
    FX, FY, FZ = np.meshgrid(frac_idx(gx, xs), frac_idx(gy, ys), frac_idx(gz, zs),
                             indexing="ij")
    return map_coordinates(vol_db, np.stack([FZ, FY, FX]), order=1, mode="nearest")


def _xy_footprint(grid: GridSpec, xs, ys) -> np.ndarray:
    gx, gy, _ = grid.axes()
    hx = (xs[1] - xs[0]) / 2 if len(xs) > 1 else 0.0
    hy = (ys[1] - ys[0]) / 2 if len(ys) > 1 else 0.0
    mx = (gx >= xs[0] - hx) & (gx <= xs[-1] + hx)
    my = (gy >= ys[0] - hy) & (gy <= ys[-1] + hy)
    return (mx[:, None] & my[None, :])[:, :, None]   # broadcast over z


# --------------------------------------------------------------------------- #
# Cost map
# --------------------------------------------------------------------------- #
@dataclass
class RFCostMap:
    grid: GridSpec
    cost: np.ndarray            # (nx,ny,nz) multiplier >= 1, inf where blocked
    crlb_rmse_m: np.ndarray     # position CRLB RMSE bound (m)
    snr_db: np.ndarray          # (n_tx, nx, ny, nz)
    n_usable: np.ndarray        # usable anchors per voxel
    blocked: np.ndarray         # bool

    def world_to_index(self, p) -> Tuple[int, int, int]:
        idx = np.rint((np.asarray(p, float) - np.asarray(self.grid.origin))
                      / np.asarray(self.grid.resolution)).astype(int)
        return tuple(int(i) for i in np.clip(idx, 0, np.asarray(self.grid.shape) - 1))

    def make_cost_fn(self, weight: float = 1.0) -> Callable:
        """
        Returns f(idx_a, idx_b) -> total step cost (metres * RF multiplier).
        weight=0 ignores RF; weight=1 uses the map as-is; >1 emphasises RF.
        Replace the distance term with your planner's own base step cost
        (e.g. the wind/energy cost) if needed.
        """
        res = np.asarray(self.grid.resolution, float)

        def f(a, b) -> float:
            a, b = tuple(a), tuple(b)
            d = float(np.linalg.norm((np.asarray(b) - np.asarray(a)) * res))
            c = 0.5 * (self.cost[a] + self.cost[b])
            return d * (1.0 + weight * (c - 1.0))
        return f

    def save(self, path: str) -> None:
        np.savez_compressed(path, origin=self.grid.origin, shape=self.grid.shape,
                            resolution=self.grid.resolution, cost=self.cost,
                            crlb_rmse_m=self.crlb_rmse_m, snr_db=self.snr_db,
                            n_usable=self.n_usable, blocked=self.blocked)

    @classmethod
    def load(cls, path: str) -> "RFCostMap":
        d = np.load(path)
        g = GridSpec(tuple(d["origin"]), tuple(int(s) for s in d["shape"]),
                     tuple(d["resolution"]))
        return cls(g, d["cost"], d["crlb_rmse_m"], d["snr_db"], d["n_usable"], d["blocked"])

    def plot_slice(self, k: int = 0, what: str = "cost", show: bool = True):
        import matplotlib.pyplot as plt
        data = {"cost": self.cost, "crlb": self.crlb_rmse_m}[what][:, :, k].astype(float)
        data = np.where(np.isfinite(data), data, np.nan)
        gx, gy, gz = self.grid.axes()
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(data.T, origin="lower", extent=[gx[0], gx[-1], gy[0], gy[-1]],
                       cmap="viridis_r" if what == "crlb" else "magma")
        ax.set_title(f"{what} @ z={gz[k]:.1f} m")
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
        fig.colorbar(im, ax=ax)
        if show:
            plt.show()
        return fig


def build_cost_map(stack: RadioMapStack, grid: GridSpec, cfg: RFConfig = RFConfig(),
                   obstacles: Optional[np.ndarray] = None) -> RFCostMap:
    n_tx = stack.path_gain.shape[0]
    k = 3 if cfg.dims.lower() == "3d" else 2

    # --- link budget: path gain -> SNR on the planner grid ---
    noise_dbm = -174.0 + 10 * np.log10(cfg.bandwidth_hz) + cfg.noise_figure_db
    p_tx = np.broadcast_to(np.atleast_1d(np.asarray(cfg.tx_power_dbm, float)), (n_tx,))
    footprint = _xy_footprint(grid, stack.x, stack.y)
    src_axes = (stack.x, stack.y, stack.z)

    snr_db = np.empty((n_tx,) + grid.shape)
    for t in range(n_tx):
        gain_db = 10 * np.log10(np.maximum(stack.path_gain[t], 10 ** (GAIN_FLOOR_DB / 10)))
        g = _resample(gain_db, src_axes, grid)
        snr = p_tx[t] + g - noise_dbm
        snr_db[t] = np.where(footprint, snr, GAIN_FLOOR_DB)

    # --- per-link ranging variance (Cramer-Rao for TOA) ---
    beta = cfg.rms_bandwidth_hz or cfg.bandwidth_hz / np.sqrt(12.0)
    snr_lin = 10 ** (snr_db / 10)
    sigma2 = C_LIGHT ** 2 / (8 * np.pi ** 2 * beta ** 2 * snr_lin)
    sigma2 = np.maximum(sigma2, cfg.sigma_floor_m ** 2)
    usable = snr_db >= cfg.min_snr_db
    w = np.where(usable, 1.0 / sigma2, 0.0)                       # (n_tx, nx, ny, nz)

    # --- geometry: unit vectors voxel -> transmitter ---
    centers = grid.centers()                                       # (nx,ny,nz,3)
    d = stack.tx_pos[:, None, None, None, :] - centers[None]
    r = np.maximum(np.linalg.norm(d, axis=-1, keepdims=True), 1e-3)
    u = (d / r)[..., :k]                                           # (n_tx,nx,ny,nz,k)

    # --- Fisher information and CRLB ---
    fim = np.einsum("t...,t...a,t...b->...ab", w, u, u)            # (nx,ny,nz,k,k)
    tr = np.trace(fim, axis1=-2, axis2=-1)
    fim = fim + (cfg.reg * np.maximum(tr, 1e-12))[..., None, None] * np.eye(k)
    rmse = np.sqrt(np.trace(np.linalg.inv(fim), axis1=-2, axis2=-1))

    n_usable = usable.sum(axis=0)
    rmse = np.where(n_usable >= cfg.min_anchors, rmse, cfg.unlocalizable_rmse_m)

    # --- normalised cost multiplier ---
    lo, hi = np.log10(cfg.crlb_good_m), np.log10(cfg.crlb_bad_m)
    bad = np.clip((np.log10(np.maximum(rmse, 1e-6)) - lo) / (hi - lo), 0.0, 1.0)
    cost = 1.0 + cfg.loc_weight * bad

    blocked = np.zeros(grid.shape, bool) if obstacles is None else obstacles.astype(bool)
    cost = np.where(blocked, np.inf, cost)
    return RFCostMap(grid, cost, rmse, snr_db, n_usable, blocked)


# --------------------------------------------------------------------------- #
# Synthetic demo (no Sionna needed) and CLI
# --------------------------------------------------------------------------- #
def synthetic_stack(freq_hz=3.5e9, extent=400.0, cell=4.0,
                    altitudes=(10, 30, 50)) -> RadioMapStack:
    """Free-space radio maps from 4 ground anchors, for testing the pipeline."""
    tx = np.array([[-150, -150, 5], [150, -150, 5], [0, 160, 5], [100, 100, 20]], float)
    ax = np.arange(-extent / 2 + cell / 2, extent / 2, cell)
    lam = C_LIGHT / freq_hz
    Z, Y, X = np.meshgrid(np.asarray(altitudes, float), ax, ax, indexing="ij")
    pg = []
    for t in tx:
        dist = np.maximum(np.sqrt((X - t[0]) ** 2 + (Y - t[1]) ** 2 + (Z - t[2]) ** 2), 1.0)
        pg.append((lam / (4 * np.pi * dist)) ** 2)
    return RadioMapStack(np.stack(pg), ax, ax, np.asarray(altitudes, float), tx)


def main():
    p = argparse.ArgumentParser(description="Sionna radio maps -> planner cost map")
    p.add_argument("--input", help="radio_maps.npz produced by RadioMapStack.to_npz")
    p.add_argument("--out", default="rf_cost_map.npz")
    p.add_argument("--grid-min", type=float, nargs=3, default=[-200, -200, 10])
    p.add_argument("--grid-max", type=float, nargs=3, default=[200, 200, 50])
    p.add_argument("--res", type=float, nargs=3, default=[4, 4, 10])
    p.add_argument("--demo", action="store_true", help="use a synthetic free-space map")
    p.add_argument("--plot", action="store_true")
    a = p.parse_args()

    stack = synthetic_stack() if a.demo or not a.input else RadioMapStack.from_npz(a.input)
    grid = GridSpec.from_bounds(a.grid_min, a.grid_max, a.res)
    cmap = build_cost_map(stack, grid)
    cmap.save(a.out)
    fin = cmap.cost[np.isfinite(cmap.cost)]
    print(f"grid {grid.shape} | cost min/mean/max = {fin.min():.2f}/{fin.mean():.2f}/{fin.max():.2f}"
          f" | best CRLB {cmap.crlb_rmse_m.min():.2f} m | saved {a.out}")
    if a.plot:
        cmap.plot_slice(k=0, what="crlb")


if __name__ == "__main__":
    main()
