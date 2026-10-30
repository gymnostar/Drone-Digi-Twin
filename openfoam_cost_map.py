"""
openfoam_cost_map.py
====================
Convert an OpenFOAM CFD result into a voxel cost map for a 3D A* drone planner.

Pipeline
--------
1. Load the OpenFOAM case (.foam / case directory) or any VTK file via PyVista.
2. Resample the (usually unstructured) CFD mesh onto a regular voxel grid that
   matches the planner's grid (origin + voxel size).
3. Derive per-voxel quantities: wind vector, wind speed, turbulence (from `k`
   if present, otherwise from velocity shear), and a blocked mask for solid /
   outside-domain voxels (buildings, terrain).
4. Combine them into
     * a static per-voxel cost volume  (CostMap.static_cost), and
     * a directional edge cost         (CostMap.edge_cost) that models the
       extra energy of flying against / across the wind.

Dependencies
------------
    pip install numpy pyvista            # required
    pip install scipy                    # optional (safety margin around obstacles)

Typical use with the planner
----------------------------
    from openfoam_cost_map import CostMap

    cm = CostMap.from_openfoam(r"C:\\cfd\\myCase", voxel_size=2.0, time="latest")
    cm.save("costmap.npz")                  # reuse later without re-reading OpenFOAM

    cm = CostMap.load("costmap.npz")
    blocked = cm.blocked                    # bool (nx, ny, nz) -> planner obstacles
    edge_cost = cm.make_edge_cost_fn()      # f(node_a, node_b) -> float

    # inside your A* neighbour loop, replace the plain distance with:
    #     step_cost = edge_cost(current, neighbour)
    # and skip neighbours where cm.is_blocked(neighbour).

Run `python openfoam_cost_map.py --demo` to test without any CFD data.
Run `python openfoam_cost_map.py CASE_DIR --voxel 2 --out costmap.npz` to convert.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Sequence, Tuple

import numpy as np

try:
    import pyvista as pv
except ImportError:  # only needed for from_openfoam()
    pv = None

Index = Tuple[int, int, int]


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
@dataclass
class CostWeights:
    """Tuning knobs. All penalty terms are normalised to ~0..1 before weighting,
    and the result is in 'calm-air-equivalent metres' per metre flown."""
    energy: float = 1.0        # directional wind energy cost (headwind/crosswind)
    turbulence: float = 1.5    # penalty for gusty / turbulent air
    shear: float = 1.0         # penalty for strong velocity gradients
    speed: float = 0.25        # penalty for high absolute wind speed
    cruise_speed: float = 10.0     # m/s ground speed of the drone
    max_airspeed: float = 20.0     # m/s; edges needing more are infeasible
    safety_margin_m: float = 2.0   # inflate obstacles by this distance (needs scipy)
    percentile: float = 95.0       # robust normalisation percentile


# --------------------------------------------------------------------------- #
# Cost map
# --------------------------------------------------------------------------- #
class CostMap:
    """Voxel cost map. Array axes are (x, y, z) = (i, j, k)."""

    def __init__(
        self,
        origin: Sequence[float],
        voxel_size: float,
        wind: np.ndarray,          # (nx, ny, nz, 3) m/s, zeros where blocked
        blocked: np.ndarray,       # (nx, ny, nz) bool
        turbulence: np.ndarray,    # (nx, ny, nz) m/s (u' ~ sqrt(2k/3)) or shear proxy
        shear: np.ndarray,         # (nx, ny, nz) 1/s
        weights: Optional[CostWeights] = None,
    ):
        self.origin = np.asarray(origin, dtype=float)
        self.voxel = float(voxel_size)
        self.wind = wind.astype(np.float32)
        self.blocked = blocked.astype(bool)
        self.turbulence = turbulence.astype(np.float32)
        self.shear = shear.astype(np.float32)
        self.weights = weights or CostWeights()
        self.shape = blocked.shape
        self.speed = np.linalg.norm(self.wind, axis=-1)
        self._build_static_cost()

    # ---- static (direction independent) cost ------------------------------ #
    def _norm(self, a: np.ndarray) -> np.ndarray:
        free = ~self.blocked
        if not free.any():
            return np.zeros_like(a)
        ref = np.percentile(a[free], self.weights.percentile)
        return np.clip(a / ref, 0.0, 2.0) if ref > 1e-9 else np.zeros_like(a)

    def _build_static_cost(self) -> None:
        w = self.weights
        pen = (
            w.turbulence * self._norm(self.turbulence)
            + w.shear * self._norm(self.shear)
            + w.speed * self._norm(self.speed)
        )
        self.penalty = pen.astype(np.float32)               # per metre flown
        self.static_cost = np.where(self.blocked, np.inf, 1.0 + pen).astype(np.float32)

    # ---- coordinate helpers ---------------------------------------------- #
    def world_to_index(self, xyz: Sequence[float]) -> Index:
        idx = np.rint((np.asarray(xyz, float) - self.origin) / self.voxel).astype(int)
        return tuple(int(v) for v in idx)  # type: ignore[return-value]

    def index_to_world(self, idx: Index) -> np.ndarray:
        return self.origin + np.asarray(idx, float) * self.voxel

    def in_bounds(self, idx: Index) -> bool:
        return all(0 <= idx[d] < self.shape[d] for d in range(3))

    def is_blocked(self, idx: Index) -> bool:
        return (not self.in_bounds(idx)) or bool(self.blocked[idx])

    # ---- directional edge cost ------------------------------------------- #
    def edge_cost(self, a: Index, b: Index) -> float:
        """Cost of flying from voxel a to voxel b (any pair, normally neighbours).

        Energy model: the drone holds ground velocity v_g = cruise * direction.
        Required airspeed = v_g - wind. Power ~ airspeed^3, time = L / cruise, so
        energy / (calm-air energy) = L * (airspeed / cruise)^3.
        Turbulence / shear / speed penalties are added per metre flown.
        """
        if self.is_blocked(a) or self.is_blocked(b):
            return float("inf")
        w = self.weights
        d = (np.asarray(b, float) - np.asarray(a, float)) * self.voxel
        L = float(np.linalg.norm(d))
        if L == 0.0:
            return 0.0
        v_ground = w.cruise_speed * d / L
        wind = 0.5 * (self.wind[a] + self.wind[b])
        airspeed = float(np.linalg.norm(v_ground - wind))
        if airspeed > w.max_airspeed:
            return float("inf")
        energy = L * (airspeed / w.cruise_speed) ** 3
        penalty = L * 0.5 * (float(self.penalty[a]) + float(self.penalty[b]))
        return w.energy * energy + penalty

    def make_edge_cost_fn(self) -> Callable[[Index, Index], float]:
        return self.edge_cost

    # ---- persistence ------------------------------------------------------ #
    def save(self, path: str) -> None:
        np.savez_compressed(
            path, origin=self.origin, voxel=self.voxel, wind=self.wind,
            blocked=self.blocked, turbulence=self.turbulence, shear=self.shear,
            weights=np.array(list(vars(self.weights).values()), dtype=float),
        )

    @classmethod
    def load(cls, path: str, weights: Optional[CostWeights] = None) -> "CostMap":
        z = np.load(path)
        if weights is None:
            weights = CostWeights(*[float(v) for v in z["weights"]])
        return cls(z["origin"], float(z["voxel"]), z["wind"], z["blocked"],
                   z["turbulence"], z["shear"], weights)

    # ---- construction from CFD ------------------------------------------- #
    @classmethod
    def from_openfoam(
        cls,
        case: str,
        voxel_size: float,
        time: str | float = "latest",
        velocity_field: str = "U",
        tke_field: str = "k",
        bounds: Optional[Sequence[float]] = None,   # (xmin,xmax,ymin,ymax,zmin,zmax)
        origin: Optional[Sequence[float]] = None,   # align to planner grid origin
        weights: Optional[CostWeights] = None,
    ) -> "CostMap":
        """Read an OpenFOAM case and resample onto a regular voxel grid.

        `case` may be a case directory, a *.foam file, or a VTK/VTU file
        (e.g. produced by `foamToVTK`).
        """
        if pv is None:
            raise ImportError("pyvista is required: pip install pyvista")
        weights = weights or CostWeights()
        mesh = _read_cfd(case, time)

        # Point data gives smooth interpolation on resampling.
        if velocity_field in mesh.cell_data and velocity_field not in mesh.point_data:
            mesh = mesh.cell_data_to_point_data()
        if velocity_field not in mesh.point_data:
            raise KeyError(
                f"Field '{velocity_field}' not found. Available: "
                f"{list(mesh.point_data.keys())}"
            )

        b = np.array(bounds if bounds is not None else mesh.bounds, dtype=float)
        org = np.array(origin if origin is not None else b[[0, 2, 4]], dtype=float)
        dims = np.floor((b[[1, 3, 5]] - org) / voxel_size).astype(int) + 1
        dims = np.maximum(dims, 1)

        grid_cls = getattr(pv, "ImageData", None) or pv.UniformGrid
        grid = grid_cls(dimensions=tuple(int(d) for d in dims),
                        spacing=(voxel_size,) * 3, origin=tuple(org))
        sampled = grid.sample(mesh)

        def to_vol(arr: np.ndarray) -> np.ndarray:  # VTK x-fastest -> (nx,ny,nz[,c])
            shape = tuple(int(d) for d in dims)
            if arr.ndim == 1:
                return arr.reshape(shape, order="F")
            return arr.reshape(shape + (arr.shape[1],), order="F")

        valid = to_vol(np.asarray(sampled["vtkValidPointMask"])).astype(bool)
        wind = to_vol(np.asarray(sampled[velocity_field], dtype=np.float64))
        wind = np.where(valid[..., None], wind, 0.0)

        # Turbulence: u' = sqrt(2k/3) when k exists, else fall back to shear.
        shear = _shear_magnitude(wind, voxel_size)
        if tke_field in sampled.point_data:
            k = to_vol(np.asarray(sampled[tke_field], dtype=np.float64))
            turb = np.sqrt(np.maximum(2.0 * k / 3.0, 0.0))
            turb = np.where(valid, turb, 0.0)
        else:
            turb = shear * voxel_size  # velocity change across one voxel (m/s)

        blocked = ~valid
        margin = int(np.ceil(weights.safety_margin_m / voxel_size))
        if margin > 0:
            blocked = _dilate(blocked, margin)

        return cls(org, voxel_size, wind, blocked, turb, shear, weights)

    # ---- synthetic data for testing -------------------------------------- #
    @classmethod
    def demo(cls, voxel_size: float = 2.0, weights: Optional[CostWeights] = None) -> "CostMap":
        """60x60x30 m domain, 6 m/s wind in +x, one building with a turbulent wake."""
        n = (60, 60, 30)
        x, y, z = np.meshgrid(*[np.arange(s) * voxel_size for s in n], indexing="ij")
        wind = np.zeros(n + (3,))
        wind[..., 0] = 6.0 * np.clip(z / 20.0, 0.3, 1.0)           # boundary-layer profile
        blocked = (abs(x - 40) < 8) & (abs(y - 60) < 16) & (z < 30)   # building
        wake = np.exp(-(((x - 62) / 20) ** 2 + ((y - 60) / 14) ** 2 + (z / 35) ** 2))
        wind[..., 0] *= 1.0 - 0.8 * wake
        turb = 2.5 * wake
        wind[blocked] = 0.0
        shear = _shear_magnitude(wind, voxel_size)
        return cls((0, 0, 0), voxel_size, wind, blocked, turb, shear, weights)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _read_cfd(case: str, time: str | float):
    p = Path(case)
    if p.is_file() and p.suffix.lower() not in (".foam", ".openfoam"):
        return pv.read(str(p))                     # .vtk / .vtu / .vtm etc.
    if p.is_dir():                                 # make the .foam stub ParaView-style
        foam = p / f"{p.name}.foam"
        if not foam.exists():
            foam.touch()
    else:
        foam = p
    reader = pv.OpenFOAMReader(str(foam))
    times = list(reader.time_values)
    if times:
        t = times[-1] if time == "latest" else float(time)
        reader.set_active_time_value(min(times, key=lambda v: abs(v - t)))
    data = reader.read()
    mesh = data["internalMesh"] if "internalMesh" in data.keys() else data.combine()
    return mesh


def _shear_magnitude(wind: np.ndarray, h: float) -> np.ndarray:
    """Frobenius norm of grad(U) via central differences (1/s)."""
    total = np.zeros(wind.shape[:3])
    for comp in range(3):
        for axis in range(3):
            if wind.shape[axis] > 1:
                total += np.gradient(wind[..., comp], h, axis=axis) ** 2
    return np.sqrt(total)


def _dilate(mask: np.ndarray, iterations: int) -> np.ndarray:
    try:
        from scipy.ndimage import binary_dilation
        return binary_dilation(mask, iterations=iterations)
    except ImportError:
        print("scipy not installed: skipping obstacle safety margin")
        return mask


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _main() -> None:
    ap = argparse.ArgumentParser(description="OpenFOAM -> drone planner cost map")
    ap.add_argument("case", nargs="?", help="OpenFOAM case dir, .foam file or VTK file")
    ap.add_argument("--voxel", type=float, default=2.0, help="voxel size in metres")
    ap.add_argument("--time", default="latest", help="time step or 'latest'")
    ap.add_argument("--out", default="costmap.npz")
    ap.add_argument("--demo", action="store_true", help="use synthetic wind field")
    ap.add_argument("--plot", action="store_true", help="show a z-slice of the cost")
    args = ap.parse_args()

    if args.demo:
        cm = CostMap.demo(args.voxel)
    elif args.case:
        cm = CostMap.from_openfoam(args.case, args.voxel, args.time)
    else:
        ap.error("provide a case path or --demo")

    cm.save(args.out)
    free = ~cm.blocked
    print(f"Saved {args.out}: grid {cm.shape}, voxel {cm.voxel} m, "
          f"{free.mean() * 100:.1f}% free, wind max {cm.speed.max():.2f} m/s")

    a, b = (2, 2, 10), (cm.shape[0] - 3, 2, 10)
    print(f"edge cost {a}->{b} (downwind-ish): {cm.edge_cost(a, (a[0] + 1, a[1], a[2])):.3f}")
    print(f"edge cost reverse (upwind):        {cm.edge_cost((a[0] + 1, a[1], a[2]), a):.3f}")

    if args.plot:
        import matplotlib.pyplot as plt
        k = min(10, cm.shape[2] - 1)
        sl = np.where(cm.blocked[:, :, k], np.nan, cm.static_cost[:, :, k])
        plt.imshow(sl.T, origin="lower", cmap="viridis")
        plt.colorbar(label="static cost per metre")
        plt.title(f"Cost map, z-slice k={k}")
        plt.show()


if __name__ == "__main__":
    _main()
