"""Wind lookup for the DLR campus: any wind direction and speed, from a handful of CFD runs.

    from wind_atlas import WindAtlas
    atlas = WindAtlas()                       # loads DLR/atlas/grid.noindex/*.npz
    U = atlas.field(theta=250, speed=6.0)     # (nx, ny, nz, 3) m/s on the shared grid, NaN in solids
    u = atlas.query(pts, theta=250, speed=6)  # (N, 3) at arbitrary points (trilinear)

How a query is answered (three approximations, each one tested in DLR/atlas/validate_atlas.py):
1. SPACE: each CFD run was resampled once onto a regular 5 x 5 x 2.5 m grid (sample_grid.py);
   between grid points, trilinear interpolation.
2. SPEED: Reynolds independence. U scales linearly with the 10 m reference speed and k with its
   square, so a 5 m/s run serves every speed. Checked against a 10 m/s run at 240 deg.
3. DIRECTION: linear blend of the two nearest solved directions (30 deg apart). This smears
   wakes that sit at different angles, so it is the weakest of the three; leave-one-out tests
   (predict a solved direction from its two neighbours 60 deg apart) give an upper bound on its error.

theta = meteorological FROM direction (deg clockwise from north); speed = wind at 10 m above
ground upstream (what EDMO reports). Coordinates: case frame, x east, y north, z up (m).
"""
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator

GRID_DIR = Path(__file__).resolve().parent / "grid.noindex"


class WindAtlas:
    def __init__(self, grid_dir=GRID_DIR, uref=5.0, exclude=()):
        """Load every <grid_dir>/dTTT_uSS.npz with reference speed `uref` (other speeds are kept
        in self.extra for validation). `exclude` = directions to leave out (leave-one-out tests)."""
        self.runs, self.extra = {}, {}
        for f in sorted(Path(grid_dir).glob("d*_u*.npz")):
            z = np.load(f)
            theta, u = float(z["theta_from_deg"]) % 360, float(z["Uref"])
            if not hasattr(self, "x"):
                self.x, self.y, self.z = z["x"], z["y"], z["z"]
                self.valid = z["valid"]
            entry = {"U": z["U"] / u, "k": z["k"] / u**2, "file": f.name}  # normalised by Uref
            if abs(u - uref) < 1e-6 and round(theta) not in exclude:
                self.runs[theta] = entry
            else:
                self.extra[(theta, u)] = entry
        if not self.runs:
            raise FileNotFoundError(f"no {uref:g} m/s runs in {grid_dir}")
        self.directions = np.array(sorted(self.runs))

    # ---- direction blending -------------------------------------------------------------
    def _weights(self, theta):
        """Two nearest solved directions around theta (wrapping at 360) and their linear weights."""
        th = np.asarray(self.directions)
        theta = theta % 360
        if len(th) == 1:
            return [(th[0], 1.0)]
        d = (theta - th) % 360                     # how far theta is past each solved direction
        a = th[np.argmin(d)]                       # nearest solved direction at or below theta
        b_gap = (th - theta) % 360
        b = th[np.argmin(b_gap)]                   # nearest solved direction at or above theta
        if a == b:
            return [(a, 1.0)]
        span = (b - a) % 360
        wb = ((theta - a) % 360) / span
        return [(a, 1 - wb), (b, wb)]

    def field(self, theta, speed=5.0, var="U"):
        """Whole-grid field for wind FROM theta at reference speed `speed` (m/s at 10 m).
        var = "U" -> (nx, ny, nz, 3) m/s; var = "k" -> (nx, ny, nz) m^2/s^2. NaN inside solids."""
        out = sum(w * self.runs[t][var] for t, w in self._weights(theta))
        return out * (speed if var == "U" else speed**2)

    def query(self, pts, theta, speed=5.0, var="U"):
        """Trilinear lookup at arbitrary points (N, 3). NaN if a neighbouring grid node is solid."""
        f = self.field(theta, speed, var)
        it = RegularGridInterpolator((self.x, self.y, self.z), f, bounds_error=False, fill_value=np.nan)
        return it(np.asarray(pts, dtype=float))

    # ---- geometry helpers ---------------------------------------------------------------
    def ground(self):
        """Height (absolute z, m) of the first air node in each column: top of terrain/building."""
        first = np.argmax(self.valid, axis=2)
        h = self.z[first].astype(float)
        h[~self.valid.any(axis=2)] = np.nan
        return h

    def terrain(self):
        """Approximate bare-ground height per column (smooth lower envelope, buildings removed):
        the minimum of ground() over a 60 m window. Used to express heights above ground."""
        from scipy.ndimage import minimum_filter
        return minimum_filter(np.nan_to_num(self.ground(), nan=1e9), size=13, mode="nearest")


if __name__ == "__main__":
    a = WindAtlas()
    print("solved directions:", a.directions, "| extra runs:", list(a.extra))
    print("grid", len(a.x), len(a.y), len(a.z), "air fraction", a.valid.mean().round(3))
    for th in (240, 255, 270):
        U = a.field(th, 5.0)
        print(f"theta {th}: mean |U| at z=60 m = {np.nanmean(np.linalg.norm(U[:, :, np.argmin(abs(a.z-60))], axis=-1)):.2f} m/s")
