"""Wind-aware 3D path planner on the DLR wind atlas: shortest path vs energy-optimal path.

    env   = Environment.from_atlas(atlas)              # geometry: solid mask, distance to solids, terrain
    wind  = WindField.from_atlas(atlas, theta, speed)  # U(x, y, z), the TRUE field used for all evaluation
    G     = PlanningGraph(env, stride=(2, 2, 2))       # coarse 26-connected directed graph (geometry only)
    A     = G.matrix(wind, mode="energy", v_g=10.0)    # sparse directed cost matrix (wind dependent)
    path, cost = G.plan(A, start_xyz, goal_xyz)        # exact Dijkstra on the discretised problem
    m     = evaluate_path(path, wind, env, mode="energy", v_g=10.0)   # dense, graph-independent metrics

Cost modes (edge i -> j, length L, unit direction d, wind w = mean of the two node winds):
  "distance"          L                                    (the shortest-path baseline)
  "energy"            L * energy_per_metre(d, v_g, w)      fixed cruise GROUND speed v_g
  "energy_bestspeed"  L * min_v energy_per_metre(d, v, w)  drone picks its speed on every edge
  "time"              L / v_ground, fixed AIRSPEED V_air   (edge dropped if it cannot be flown)

All coordinates are metres, x east, y north, z up (absolute, as in the atlas).
"""
from itertools import product
import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

from energy import energy_per_metre, best_energy_per_metre, segment_power, SPEEDS

MODES = ("distance", "energy", "energy_bestspeed", "time")


def _axis_step(a):
    return float(a[1] - a[0])


# --------------------------------------------------------------------------------------------
# geometry and wind on the fine grid
# --------------------------------------------------------------------------------------------
class Environment:
    """Fine-grid geometry: air mask, distance to the nearest solid node, bare-ground height."""

    def __init__(self, x, y, z, valid, terrain=None):
        self.x, self.y, self.z = (np.asarray(a, float) for a in (x, y, z))
        self.valid = np.asarray(valid, bool)
        self.h = (_axis_step(self.x), _axis_step(self.y), _axis_step(self.z))
        self.terrain = np.zeros(self.valid.shape[:2]) if terrain is None else np.asarray(terrain, float)
        if self.valid.all():
            self.dist = np.full(self.valid.shape, 1e6)
        else:  # distance (m) from every air node to the nearest solid node; 0 inside solids
            self.dist = distance_transform_edt(self.valid, sampling=self.h)

    @classmethod
    def from_atlas(cls, atlas):
        return cls(atlas.x, atlas.y, atlas.z, atlas.valid, atlas.terrain())

    def index(self, pts):
        """Fractional fine-grid index coordinates (3, N) of points (N, 3)."""
        pts = np.atleast_2d(pts)
        return np.stack([(pts[:, 0] - self.x[0]) / self.h[0],
                         (pts[:, 1] - self.y[0]) / self.h[1],
                         (pts[:, 2] - self.z[0]) / self.h[2]])

    def clearance(self, pts):
        """Trilinear distance (m) to the nearest solid grid node at points (N, 3)."""
        return map_coordinates(self.dist, self.index(pts), order=1, mode="nearest")

    def height_above_terrain(self, pts):
        pts = np.atleast_2d(pts)
        ij = self.index(pts)[:2]
        t = map_coordinates(self.terrain, ij, order=1, mode="nearest")
        return pts[:, 2] - t


class WindField:
    """Wind vectors U (nx, ny, nz, 3) on the fine grid, NaN in solids; trilinear sampling that
    stays finite right next to a solid (valid-corner-weighted average)."""

    def __init__(self, x, y, z, U):
        self.x, self.y, self.z = (np.asarray(a, float) for a in (x, y, z))
        self.U = U
        self.h = (_axis_step(self.x), _axis_step(self.y), _axis_step(self.z))
        ok = np.isfinite(U[..., 0])
        self._mask = ok.astype(np.float32)
        self._U0 = np.where(ok[..., None], U, 0.0).astype(np.float32)

    @classmethod
    def from_atlas(cls, atlas, theta, speed):
        return cls(atlas.x, atlas.y, atlas.z, atlas.field(theta, speed))

    @classmethod
    def uniform(cls, env, vec):
        U = np.broadcast_to(np.asarray(vec, np.float32), env.valid.shape + (3,)).copy()
        U[~env.valid] = np.nan
        return cls(env.x, env.y, env.z, U)

    def sample(self, pts):
        """Wind (N, 3) at arbitrary points. Points deep inside solids return 0."""
        pts = np.atleast_2d(pts)
        c = np.stack([(pts[:, 0] - self.x[0]) / self.h[0], (pts[:, 1] - self.y[0]) / self.h[1],
                      (pts[:, 2] - self.z[0]) / self.h[2]])
        m = map_coordinates(self._mask, c, order=1, mode="nearest")
        out = np.stack([map_coordinates(self._U0[..., k], c, order=1, mode="nearest") for k in range(3)], -1)
        return np.where(m[:, None] > 1e-6, out / np.maximum(m, 1e-6)[:, None], 0.0)


# --------------------------------------------------------------------------------------------
# dense polyline evaluation (shared by evaluate_path, smoothing)
# --------------------------------------------------------------------------------------------
def _subdivide(P, ds):
    """Cut every segment of polyline P (m, 3) into pieces <= ds. Returns piece midpoints, unit
    directions, piece lengths and the index of the parent segment."""
    P = np.asarray(P, float)
    d = np.diff(P, axis=0)
    L = np.linalg.norm(d, axis=1)
    keep = L > 1e-9
    seg = np.nonzero(keep)[0]
    n = np.maximum(1, np.ceil(L[seg] / ds).astype(int))
    parent = np.repeat(seg, n)
    first = np.cumsum(n) - n
    k = np.arange(parent.size) - np.repeat(first, n)             # piece number inside its segment
    frac = (k + 0.5) / np.repeat(n, n)
    mid = P[parent] + frac[:, None] * d[parent]
    return mid, d[parent] / L[parent][:, None], L[parent] / np.repeat(n, n), parent


def _piece_cost(dirs, dl, W, mode, v_g, V_air, speeds=SPEEDS):
    """Per-piece (cost, time, speed) of flying unit direction dirs for length dl in wind W."""
    if mode in ("distance", "energy"):
        e = energy_per_metre(dirs, v_g, W)
        return e * dl, dl / v_g, np.full(dl.shape, float(v_g))
    if mode == "energy_bestspeed":
        e, v = best_energy_per_metre(dirs, W, speeds=speeds)
        return e * dl, dl / v, v
    if mode == "time":
        vg = _airspeed_ground_speed(dirs, W, V_air)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = dl / vg
            e = segment_power(dirs, np.where(np.isfinite(vg), vg, 1.0), W) * t
        return e, t, vg
    raise ValueError(mode)


def _airspeed_ground_speed(dirs, W, V_air):
    """Ground speed along unit dirs when the drone holds airspeed MAGNITUDE V_air (3D): solves
    |v d - w| = V_air for v > 0. NaN where the wind is too strong to make headway."""
    dw = np.einsum("...k,...k->...", dirs, W)
    disc = dw**2 - np.einsum("...k,...k->...", W, W) + V_air**2
    with np.errstate(invalid="ignore"):
        v = dw + np.sqrt(disc)
    return np.where((disc >= 0) & (v > 1e-3), v, np.nan)


def evaluate_path(path_xyz, wind, env=None, mode="energy", v_g=10.0, V_air=10.0, ds=2.0):
    """Metrics of a polyline in the TRUE wind field, sampled every `ds` m (trilinear), independent of
    the planning grid. mode: "energy" (fixed ground speed v_g), "energy_bestspeed", "time" (fixed
    airspeed V_air); "distance" is reported like "energy". The same function scores every path.
    Returns length_m, energy_J, time_s, mean_headwind (m/s, + = against the drone, length-weighted),
    max_wind (m/s), min_clearance_m (to the nearest solid grid node), min_height_agl_m."""
    P = np.asarray(path_xyz, float)
    mid, dirs, dl, _ = _subdivide(P, ds)
    W = wind.sample(mid)
    cost, t, v = _piece_cost(dirs, dl, W, mode, v_g, V_air)
    e = cost                                      # for the three energy-type modes cost IS energy (J)
    L = dl.sum()
    pts = np.vstack([P, mid])
    out = dict(length_m=float(L), energy_J=float(np.nansum(e)), time_s=float(np.nansum(t)),
               mean_headwind=float(-(np.einsum("ij,ij->i", W, dirs) * dl).sum() / L) if L > 0 else 0.0,
               max_wind=float(np.linalg.norm(W, axis=1).max()),
               min_clearance_m=float(env.clearance(pts).min()) if env is not None else np.nan,
               min_height_agl_m=float(env.height_above_terrain(pts).min()) if env is not None else np.nan)
    return out


# --------------------------------------------------------------------------------------------
# planning graph
# --------------------------------------------------------------------------------------------
class PlanningGraph:
    """Coarse, directed, 26-connected graph over the free nodes of the planning band.

    stride    (sx, sy, sz): coarse nodes are every sx-th fine node (default 2,2,2 = 10 m x 10 m x 5 m).
    margin    obstacle inflation: a node is free only if no solid fine node is within `margin` m.
    h_min/max planning band above the local bare terrain (m); also capped by the atlas ceiling.
    The graph structure (edges, lengths, directions) is wind independent; costs are built per wind.
    """

    OFFSETS = np.array([o for o in product((-1, 0, 1), repeat=3) if o != (0, 0, 0)], dtype=int)

    def __init__(self, env, stride=(2, 2, 2), margin=5.0, h_min=10.0, h_max=100.0, speeds=SPEEDS[::2]):
        self.speeds = speeds     # candidate ground speeds for "energy_bestspeed" edges (1 m/s steps; evaluation uses 0.5)
        self.env, self.stride, self.margin = env, tuple(stride), margin
        self.h_min, self.h_max = h_min, h_max
        sx, sy, sz = self.stride
        agl = env.z[None, None, :] - env.terrain[:, :, None]
        band = (agl >= h_min - 1e-9) & (agl <= h_max + 1e-9)
        self.ok_fine = (env.dist > margin) & band               # inflated-obstacle + altitude-band mask
        free = self.ok_fine[::sx, ::sy, ::sz]
        self.shape = free.shape
        self.ids = np.full(free.shape, -1, np.int32)
        self.n = int(free.sum())
        self.ids[free] = np.arange(self.n, dtype=np.int32)
        ijk = np.argwhere(free)                                   # coarse indices
        self.node_fine = ijk * np.array(self.stride)              # fine indices of the nodes
        self.node_xyz = np.stack([env.x[self.node_fine[:, 0]], env.y[self.node_fine[:, 1]],
                                  env.z[self.node_fine[:, 2]]], axis=1)
        self.step = np.array([env.h[0] * sx, env.h[1] * sy, env.h[2] * sz])
        self._build_edges(free)
        self._tree = None

    def _build_edges(self, free):
        sx, sy, sz = self.stride
        shp = self.shape
        src, dst, off = [], [], []
        vec = self.OFFSETS * self.step
        self.off_len = np.linalg.norm(vec, axis=1)
        self.off_dir = vec / self.off_len[:, None]
        for k, o in enumerate(self.OFFSETS):
            s_sl, d_sl, s_ix, d_ix = [], [], [], []
            for a in range(3):
                n = shp[a]
                if o[a] == 0:
                    sl_s = sl_d = slice(0, n)
                elif o[a] == 1:
                    sl_s, sl_d = slice(0, n - 1), slice(1, n)
                else:
                    sl_s, sl_d = slice(1, n), slice(0, n - 1)
                s_sl.append(sl_s); d_sl.append(sl_d)
                s_ix.append(np.arange(n)[sl_s] * self.stride[a]); d_ix.append(np.arange(n)[sl_d] * self.stride[a])
            ok = free[tuple(s_sl)] & free[tuple(d_sl)]
            mid = self.ok_fine[np.ix_(*[(s + d) // 2 for s, d in zip(s_ix, d_ix)])]   # midpoint must be free too
            ok &= mid
            src.append(self.ids[tuple(s_sl)][ok]); dst.append(self.ids[tuple(d_sl)][ok])
            off.append(np.full(int(ok.sum()), k, np.int8))
        self.src = np.concatenate(src); self.dst = np.concatenate(dst); self.off = np.concatenate(off)
        self.length = self.off_len[self.off]
        self.n_edges = self.src.size

    # ---- costs ------------------------------------------------------------------------
    def costs(self, wind, mode="energy", v_g=10.0, V_air=10.0, chunk=40000):
        """Directed edge costs (n_edges,) for this wind field; np.inf = edge cannot be flown."""
        if mode == "distance":
            return self.length.astype(float).copy()
        f = self.node_fine
        Un = wind.U[f[:, 0], f[:, 1], f[:, 2]].astype(float)
        out = np.empty(self.n_edges)
        for a in range(0, self.n_edges, chunk):
            s, d, o = self.src[a:a + chunk], self.dst[a:a + chunk], self.off[a:a + chunk]
            W = 0.5 * (Un[s] + Un[d])
            c, _, _ = _piece_cost(self.off_dir[o], self.off_len[o], W, mode, v_g, V_air, speeds=self.speeds)
            out[a:a + chunk] = np.where(np.isfinite(c), c, np.inf)
        return out

    def matrix(self, wind=None, mode="distance", v_g=10.0, V_air=10.0):
        """Sparse directed cost matrix (csr). A[i, j] = cost of flying i -> j."""
        c = self.costs(wind, mode, v_g, V_air)
        good = np.isfinite(c) & (c > 0)
        return csr_matrix((c[good], (self.src[good], self.dst[good])), shape=(self.n, self.n))

    # ---- search -----------------------------------------------------------------------
    def snap(self, xyz):
        """Index of the nearest free node to each point (N, 3)."""
        if self._tree is None:
            self._tree = cKDTree(self.node_xyz)
        return self._tree.query(np.atleast_2d(xyz))[1]

    @staticmethod
    def solve(A, sources, limit=np.inf):
        """Dijkstra from each source node: (dist, predecessors), each (len(sources), n)."""
        return dijkstra(A, directed=True, indices=np.atleast_1d(sources), return_predecessors=True, limit=limit)

    def extract(self, pred_row, source, goal):
        """Node-id path source -> goal from a predecessor row, or None if unreachable."""
        if source == goal:
            return np.array([source])
        if pred_row[goal] < 0:
            return None
        p = [goal]
        while p[-1] != source:
            p.append(pred_row[p[-1]])
        return np.array(p[::-1])

    def plan(self, A, start_xyz, goal_xyz):
        """Optimal path (xyz (m, 3), total cost) between two points (snapped to free nodes), or (None, inf)."""
        s, g = (int(i) for i in self.snap(np.vstack([start_xyz, goal_xyz])))
        dist, pred = self.solve(A, s)
        ids = self.extract(pred[0], s, g)
        return (None, np.inf) if ids is None else (self.node_xyz[ids], float(dist[0, g]))

    # ---- smoothing (line-of-sight string pulling) -------------------------------------
    def segments_free(self, p, Q, spacing=None):
        """For one point p (3,) and candidate endpoints Q (J, 3): is the straight segment p -> Q[j] inside
        the free (inflated, in-band) region? Checked every `spacing` m against the fine mask."""
        spacing = spacing or 0.5 * min(self.env.h)
        Lmax = np.linalg.norm(Q - p, axis=1).max()
        t = np.linspace(0, 1, int(np.ceil(Lmax / spacing)) + 1)
        pts = p + t[None, :, None] * (Q - p)[:, None, :]            # (J, K, 3)
        e = self.env
        ix = np.rint((pts[..., 0] - e.x[0]) / e.h[0]).astype(int)
        iy = np.rint((pts[..., 1] - e.y[0]) / e.h[1]).astype(int)
        iz = np.rint((pts[..., 2] - e.z[0]) / e.h[2]).astype(int)
        inb = (ix >= 0) & (ix < len(e.x)) & (iy >= 0) & (iy < len(e.y)) & (iz >= 0) & (iz < len(e.z))
        okp = np.zeros(inb.shape, bool)
        okp[inb] = self.ok_fine[ix[inb], iy[inb], iz[inb]]
        return okp.all(axis=1)

    def smooth(self, path_xyz, wind=None, mode="distance", v_g=10.0, V_air=10.0, ds=2.0):
        """Greedy line-of-sight shortcutting: from each vertex jump to the farthest later vertex that is
        visible (straight segment stays free) AND, for energy modes, whose straight segment costs no more
        than the path piece it replaces (in `wind`). The same algorithm serves every cost mode."""
        P = np.asarray(path_xyz, float)
        n = len(P)
        if n < 3:
            return P.copy()
        cum = np.zeros(n)
        if mode != "distance":
            mid, dirs, dl, par = _subdivide(P, ds)
            c = _piece_cost(dirs, dl, wind.sample(mid), mode, v_g, V_air)[0]
            cum[1:] = np.cumsum(np.bincount(par, weights=c, minlength=n - 1))
        out, i = [0], 0
        while i < n - 1:
            js = np.arange(n - 1, i + 1, -1)                        # candidates beyond the next vertex
            vis = self.segments_free(P[i], P[js]) if js.size else np.zeros(0, bool)
            nxt = i + 1
            for j, v in zip(js, vis):
                if not v:
                    continue
                if mode != "distance":
                    m2, d2, l2, _ = _subdivide(P[[i, j]], ds)
                    cs = _piece_cost(d2, l2, wind.sample(m2), mode, v_g, V_air)[0].sum()
                    if cs > (cum[j] - cum[i]) * (1 + 1e-9):
                        continue
                nxt = int(j)
                break
            out.append(nxt); i = nxt
        return P[out]
