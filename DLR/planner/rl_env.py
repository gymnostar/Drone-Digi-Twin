"""Reinforcement-learning environment for the wind-aware drone planner, with an exact Dijkstra oracle.

    from rl_env import WindDroneEnv
    env = WindDroneEnv(atlas, stride=(4, 4, 2))                       # coarse 26-connected graph, 20 x 20 x 5 m
    obs, info = env.reset(seed=0, options={"theta": 240, "speed": 8.0})
    obs, r, terminated, truncated, info = env.step(env.oracle_action(info["node"]))

Gymnasium API (reset / step with a 5-tuple), implemented here WITHOUT importing gymnasium so that the
prototype runs in the plain project venv. If gymnasium is installed its Discrete/Box spaces are used.

THE POINT. For one fixed, known wind field the planner's Dijkstra is EXACT on the discretised problem,
so it is the oracle any learning agent is measured against:

    regret = (agent energy - optimal energy) / optimal energy        >= 0, and 0 only for a perfect agent

Reward = -(edge energy in J) / E_REF, no goal bonus  ->  undiscounted return == -(total energy) / E_REF.
Maximising return is therefore EXACTLY minimising energy, and `-return * E_REF` can be compared with
`optimal_cost()` with no unit or shaping conversion. See RL_ENV.md.

Observation (float32, length 36):  [0:3] position in the graph bounding box, in [0, 1]
                                   [3:6] (goal - position) / box extent, in [-1, 1]
                                   [6:9] local wind (u, v, w) / 10 m/s       (x east, y north, z up)
                                   [9]   straight-line distance to goal / box diagonal
                                   [10:36] 1.0 where move k (OFFSETS[k]) is a valid edge, else 0.0
Action: Discrete(26), move k = PlanningGraph.OFFSETS[k] (dx, dy, dz in coarse steps);
        with speed_choice=True: Discrete(26 * n_speeds), action = k * n_speeds + speed_index.
"""
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "atlas"))
from planner import Environment, PlanningGraph, WindField, _piece_cost    # noqa: E402
from energy import energy_per_metre                                        # noqa: E402

try:                                                    # optional: use real Gymnasium spaces when present
    from gymnasium import spaces as _gym_spaces
except Exception:                                       # pragma: no cover
    _gym_spaces = None

E_REF = 1000.0          # J. Reward = -energy / E_REF, i.e. one reward unit is one kJ.
P_INVALID = 1.0         # reward units (1 kJ): penalty for a move with no edge. A single edge costs ~0.1-0.5 kJ,
                        # so bumping into a wall is never cheaper than flying a real edge, yet it is small enough
                        # not to swamp the energy signal.
P_TRUNC = 30.0          # reward units (30 kJ): charged once when the episode is cut by max_steps. About 2-3x the
                        # energy of a typical campus crossing, so "wandering until the clock ends" always loses.

# EDMO airfield 2020-2025, FROM-direction sectors: (% of hours, mean speed m/s). Calm (2.9 %) is excluded.
EDMO = {0: (3.6, 2.95), 30: (7.1, 3.15), 60: (8.3, 3.34), 90: (6.1, 2.94), 120: (3.1, 2.00), 150: (6.8, 2.07),
        180: (8.2, 2.11), 210: (15.8, 3.02), 240: (22.9, 5.15), 270: (10.3, 4.31), 300: (4.1, 3.39), 330: (3.7, 3.08)}


# --------------------------------------------------------------------------------------------
# minimal spaces (used only when gymnasium is not importable)
# --------------------------------------------------------------------------------------------
class _Discrete:
    def __init__(self, n):
        self.n, self.shape, self.dtype = int(n), (), np.int64
        self._rng = np.random.default_rng()

    def seed(self, seed=None):
        self._rng = np.random.default_rng(seed)

    def sample(self, mask=None):
        if mask is None:
            return int(self._rng.integers(self.n))
        return int(self._rng.choice(np.nonzero(np.asarray(mask))[0]))

    def contains(self, x):
        return bool(np.ndim(x) == 0 and float(x) == int(x) and 0 <= int(x) < self.n)

    def __repr__(self):
        return f"Discrete({self.n})"


class _Box:
    def __init__(self, low, high, shape, dtype=np.float32):
        self.low, self.high, self.shape, self.dtype = low, high, tuple(shape), np.dtype(dtype)
        self._rng = np.random.default_rng()

    def seed(self, seed=None):
        self._rng = np.random.default_rng(seed)

    def sample(self):
        return self._rng.normal(size=self.shape).astype(self.dtype)

    def contains(self, x):
        x = np.asarray(x)
        return bool(x.shape == self.shape and np.all(x >= self.low) and np.all(x <= self.high))

    def __repr__(self):
        return f"Box({self.low}, {self.high}, {self.shape}, {self.dtype})"


def make_discrete(n):
    return _gym_spaces.Discrete(n) if _gym_spaces else _Discrete(n)


def make_box(shape):
    return (_gym_spaces.Box(-np.inf, np.inf, shape, np.float32) if _gym_spaces
            else _Box(-np.inf, np.inf, shape, np.float32))


def uniform_wind_factory(env):
    """(theta, speed) -> uniform WindField blowing FROM theta (met. convention): handy for tests and toys."""
    def make(theta, speed):
        a = np.radians(theta)
        return WindField.uniform(env, [-speed * np.sin(a), -speed * np.cos(a), 0.0])
    return make


class _Case:
    """Everything that depends on the wind (theta, speed): edge costs, local wind, oracle cache."""
    pass


class WindDroneEnv:
    """Drone navigation on the planner's coarse 26-connected graph in one frozen wind field per episode.

    Parameters
    atlas          WindAtlas (default: load DLR/atlas/grid.noindex). Ignored if `environment` and `wind_factory` given.
    environment    planner.Environment (default: Environment.from_atlas(atlas)).
    wind_factory   callable (theta, speed) -> WindField (default: WindField.from_atlas(atlas, ...)).
    stride         PlanningGraph coarse stride, (4, 4, 2) = 20 x 20 x 5 m nodes. (2, 2, 2) is the planner default.
    cost_mode      "energy_bestspeed" (default; the drone flies each edge at its cheapest ground speed) or
                   "energy" (fixed ground speed v_g). Ignored when speed_choice=True (the agent picks the speed).
    speed_choice   action = (move, speed): 26 * len(speeds) discrete actions. Oracle = best speed on every edge.
    speeds         candidate ground speeds (m/s); default planner's.
    layer_z        if set, fly only in the horizontal layer nearest this absolute altitude (m): the 8 horizontal
                   moves are valid, the 18 vertical ones are not. A deliberately small sub-problem for tabular RL.
    max_steps      truncation horizon (counts invalid moves too).
    shaping        potential-based reward shaping (Ng, Harada & Russell 1999); default off; info["energy_J"] always holds the true energy.
    gamma          discount used ONLY inside the shaping term (1.0 = undiscounted, the benchmark setting).
    directions     wind directions to sample at reset (default: climate sectors present in the atlas).
    speed_jitter   log-normal sigma of the sampled speed around the sector mean (0 = always the mean).
    theta, speed   fix the wind for every episode (options at reset override).
    cache_size     number of wind cases (cost tables + oracle) kept; building one costs seconds on the full graph.
    """

    def __init__(self, atlas=None, environment=None, wind_factory=None, stride=(4, 4, 2), cost_mode="energy_bestspeed",
                 v_g=10.0, speed_choice=False, speeds=None, layer_z=None, max_steps=300, shaping=False, gamma=1.0,
                 directions=None, speed_jitter=0.35, theta=None, speed=None, min_sep=200.0, cache_size=6,
                 margin=5.0, h_min=10.0, h_max=100.0, seed=None):
        if cost_mode not in ("energy_bestspeed", "energy"):
            raise ValueError("cost_mode must be 'energy_bestspeed' or 'energy' (rewards are energies)")
        if environment is None or wind_factory is None:
            if atlas is None:
                from wind_atlas import WindAtlas
                atlas = WindAtlas()
            if environment is None:
                environment = Environment.from_atlas(atlas)
            if wind_factory is None:
                wind_factory = lambda th, sp: WindField.from_atlas(atlas, th, sp)   # noqa: E731
        self.atlas, self.environment, self.wind_factory = atlas, environment, wind_factory
        kw = {} if speeds is None else {"speeds": np.asarray(speeds, float)}
        self.G = G = PlanningGraph(environment, stride=stride, margin=margin, h_min=h_min, h_max=h_max, **kw)
        self.cost_mode, self.v_g, self.speed_choice = cost_mode, float(v_g), bool(speed_choice)
        self.speeds = np.asarray(G.speeds, float)
        self.n_speeds = len(self.speeds) if self.speed_choice else 1
        self.max_steps, self.shaping, self.gamma, self.min_sep = int(max_steps), bool(shaping), float(gamma), min_sep
        self.speed_jitter, self.fixed_wind = float(speed_jitter), (theta, speed)
        self.cache_size = int(cache_size)

        # static structure: which of the 26 moves exists from each node (wind independent)
        n = G.n
        off = G.OFFSETS
        keep = np.ones(G.n_edges, bool)
        self.layer_z = None
        if layer_z is not None:
            zs = np.unique(G.node_xyz[:, 2])
            self.layer_z = float(zs[np.argmin(abs(zs - layer_z))])
            keep = (off[G.off, 2] == 0) & (G.node_xyz[G.src, 2] == self.layer_z)
        self.e_src, self.e_dst, self.e_off = G.src[keep], G.dst[keep], G.off[keep].astype(np.int64)
        self.next_node = np.full((n, 26), -1, np.int32)
        self.next_node[self.e_src, self.e_off] = self.e_dst
        self.edge_index = np.full((n, 26), -1, np.int64)
        self.edge_index[self.e_src, self.e_off] = np.arange(self.e_src.size)
        touched = np.zeros(n, bool)
        touched[self.e_src] = True
        touched[self.e_dst] = True
        self.active = np.nonzero(touched)[0]                       # nodes with at least one edge in the sub-problem
        self._tree = cKDTree(G.node_xyz[self.active])
        lo, hi = G.node_xyz.min(0), G.node_xyz.max(0)
        self._lo, self._ext = lo, np.maximum(hi - lo, 1.0)
        self._diag = float(np.linalg.norm(hi - lo))
        self._pos = ((G.node_xyz - lo) / self._ext).astype(np.float32)       # normalised positions (n, 3)

        # sampling of the climate (only sectors that exist in the atlas)
        self.directions, self.dir_p, self.dir_speed = self._climate(directions)

        self.n_actions = 26 * self.n_speeds
        self.action_space = make_discrete(self.n_actions)
        self.observation_space = make_box((36,))
        self._cases = OrderedDict()
        self.case = None
        self._rng = np.random.default_rng(seed)
        if seed is not None:
            self.action_space.seed(seed)
        self.node = self.goal = None
        self.steps = 0
        self.energy_J = 0.0

    # ----------------------------------------------------------------------------------
    # climate and wind cases
    # ----------------------------------------------------------------------------------
    def _climate(self, directions):
        have = None
        if directions is not None:
            have = [float(d) % 360 for d in directions]
        elif self.atlas is not None:
            have = [float(d) for d in self.atlas.directions]
        if have is None:                                           # custom factory: all 12 sectors
            sect = dict(EDMO)
        else:
            sect = {s: v for s, v in EDMO.items() if any(abs((s - d + 180) % 360 - 180) < 0.5 for d in have)}
            if not sect:                                           # no overlap with the climate: uniform, mean 4 m/s
                sect = {d: (1.0, 4.0) for d in have}
        ths = np.array(sorted(sect), float)
        p = np.array([sect[int(t)][0] for t in ths], float)
        return ths, p / p.sum(), np.array([sect[int(t)][1] for t in ths], float)

    def sample_wind(self):
        i = int(self._rng.choice(len(self.directions), p=self.dir_p))
        sp = self.dir_speed[i]
        if self.speed_jitter > 0:
            sp = sp * np.exp(self.speed_jitter * self._rng.standard_normal() - 0.5 * self.speed_jitter**2)
            sp = float(np.clip(np.round(sp * 2) / 2, 1.0, 20.0))   # 0.5 m/s bins so cases can be re-used
        return float(self.directions[i]), float(sp)

    def set_wind(self, theta, speed):
        """Select (and build, or fetch from the cache) the wind case for FROM-direction theta, speed m/s."""
        key = (round(float(theta), 3), round(float(speed), 3))
        c = self._cases.get(key)
        if c is None:
            c = self._build_case(*key)
            self._cases[key] = c
            while len(self._cases) > self.cache_size:
                self._cases.popitem(last=False)
        else:
            self._cases.move_to_end(key)
        self.case = c
        return c

    def _edge_costs(self, Un):
        """Energy (J) of every active edge; shape (E,) or (E, n_speeds) when speed_choice. Same formulas as
        PlanningGraph.costs (mean of the two node winds), restricted to the active edges."""
        G, s, d, o = self.G, self.e_src, self.e_dst, self.e_off
        W = 0.5 * (Un[s] + Un[d])
        dirs, ln = G.off_dir[o], G.off_len[o]
        if self.speed_choice:
            e = energy_per_metre(dirs[:, None, :], self.speeds[None, :], W[:, None, :]) * ln[:, None]
            return np.where(np.isfinite(e), e, np.inf)
        c, _, _ = _piece_cost(dirs, ln, W, self.cost_mode, self.v_g, 10.0, speeds=G.speeds)
        return np.where(np.isfinite(c), c, np.inf)

    def _build_case(self, theta, speed):
        G = self.G
        wind = self.wind_factory(theta, speed)
        f = G.node_fine
        Un = wind.U[f[:, 0], f[:, 1], f[:, 2]].astype(float)
        c = _Case()
        c.theta, c.speed = theta, speed
        c.node_wind = np.nan_to_num(Un, nan=0.0)
        c.wind10 = (c.node_wind / 10.0).astype(np.float32)
        edge = self._edge_costs(Un)
        if self.speed_choice:
            c.edge_cost_s = edge
            edge_best = edge.min(axis=1)
        else:
            edge_best = edge
        c.edge_cost = edge_best
        good = np.isfinite(edge_best) & (edge_best > 0)
        c.cost = np.full((G.n, 26), np.inf)
        c.cost[self.e_src[good], self.e_off[good]] = edge_best[good]
        c.valid = np.isfinite(c.cost)
        if self.speed_choice:
            c.cost_s = np.full((G.n, 26, self.n_speeds), np.inf)
            c.cost_s[self.e_src[good], self.e_off[good]] = c.edge_cost_s[good]
        c.csr_T = csr_matrix((edge_best[good], (self.e_dst[good], self.e_src[good])), shape=(G.n, G.n))
        c.csr_len = csr_matrix((G.off_len[self.e_off[good]], (self.e_src[good], self.e_dst[good])), shape=(G.n, G.n))
        c.min_epm = float((edge_best[good] / G.off_len[self.e_off[good]]).min())    # J/m, lower bound of any edge
        c.oracle = {}
        return c

    # ----------------------------------------------------------------------------------
    # oracle
    # ----------------------------------------------------------------------------------
    def node_of(self, x):
        """Node index of an int or of the nearest ACTIVE node to xyz (3,)."""
        if np.ndim(x) == 0:
            return int(x)
        return int(self.active[self._tree.query(np.asarray(x, float))[1]])

    def oracle_policy(self, goal=None):
        """(action (n,), cost_to_go (n,) in J) for the current wind: one reverse-graph Dijkstra from `goal`
        on the transposed cost matrix. action[i] is the optimal move (with speed if speed_choice) from node i;
        -1 at the goal and at nodes that cannot reach it (cost_to_go = inf). Cached per (wind, goal)."""
        c = self._need_case()
        goal = self.goal if goal is None else self.node_of(goal)
        hit = c.oracle.get(goal)
        if hit is not None:
            return hit
        dist, pred = dijkstra(c.csr_T, directed=True, indices=goal, return_predecessors=True)
        # reverse graph: pred[i] is the NEXT node after i on the optimal forward path i -> goal
        match = (self.next_node == pred[:, None]) & c.valid
        move = np.argmax(match, axis=1)
        has = match.any(axis=1)
        if self.speed_choice:
            sidx = np.argmin(c.cost_s[np.arange(self.G.n), move], axis=1)
            act = move * self.n_speeds + sidx
        else:
            act = move
        act = np.where(has, act, -1).astype(np.int64)
        if len(c.oracle) > 32:
            c.oracle.clear()
        c.oracle[goal] = (act, dist)
        return act, dist

    def oracle_action(self, node=None):
        """Optimal action from `node` (default: the current node) towards the current goal; -1 if at the goal."""
        node = self.node if node is None else self.node_of(node)
        return int(self.oracle_policy()[0][node])

    def optimal_cost(self, start, goal):
        """Minimum total energy (J) from start to goal (node ids or xyz) in the current wind; inf if unreachable."""
        return float(self.oracle_policy(self.node_of(goal))[1][self.node_of(start)])

    def _move_between(self, i, j):
        k = np.nonzero((self.next_node[i] == j) & self._need_case().valid[i])[0]
        if k.size == 0:
            raise ValueError(f"no edge {i} -> {j}")
        return int(k[0])

    def path_cost(self, nodes):
        """Energy (J) of a node path in the current wind (best speed on every edge if speed_choice)."""
        c = self._need_case()
        return float(sum(c.cost[a, self._move_between(a, b)] for a, b in zip(nodes[:-1], nodes[1:])))

    def shortest_path(self, start, goal):
        """Geometric (distance-only) shortest path start -> goal on the same graph: (nodes, energy in J in the
        current wind). This is the 'ignore the wind' baseline the energy planner is meant to beat."""
        c = self._need_case()
        s, g = self.node_of(start), self.node_of(goal)
        d, pred = dijkstra(c.csr_len, directed=True, indices=s, return_predecessors=True)
        if not np.isfinite(d[g]):
            return None, np.inf
        p = [g]
        while p[-1] != s:
            p.append(int(pred[p[-1]]))
        p = p[::-1]
        return np.array(p), self.path_cost(p)

    def _need_case(self):
        if self.case is None:
            raise RuntimeError("no wind selected yet: call reset() or set_wind(theta, speed)")
        return self.case

    # ----------------------------------------------------------------------------------
    # Gymnasium API
    # ----------------------------------------------------------------------------------
    def reset(self, seed=None, options=None):
        """Start an episode. options (all optional): theta, speed (wind FROM direction deg, m/s at 10 m),
        start, goal (node id or xyz in m). Anything not given is sampled (wind from the climate, goal from
        the active nodes, start from nodes that can reach the goal and are >= min_sep m away)."""
        if seed is not None:
            self._rng = np.random.default_rng(seed)
            self.action_space.seed(seed)
        o = dict(options or {})
        th, sp = o.get("theta", self.fixed_wind[0]), o.get("speed", self.fixed_wind[1])
        if th is None or sp is None:
            sth, ssp = self.sample_wind()
            th = sth if th is None else th
            sp = ssp if sp is None else sp
        self.set_wind(th, sp)
        xyz = self.G.node_xyz
        for _ in range(50):
            goal = self.node_of(o["goal"]) if "goal" in o else int(self._rng.choice(self.active))
            _, dist = self.oracle_policy(goal)
            if "start" in o:
                start = self.node_of(o["start"])
                break
            far = np.linalg.norm(xyz[self.active] - xyz[goal], axis=1) >= self.min_sep
            cand = self.active[np.isfinite(dist[self.active]) & far]
            if cand.size:
                start = int(self._rng.choice(cand))
                break
        else:
            raise RuntimeError("could not sample a reachable start; relax min_sep")
        if not np.isfinite(dist[start]):
            raise ValueError("goal is unreachable from the requested start in this wind")
        self.goal, self.node, self.steps, self.energy_J = goal, start, 0, 0.0
        self._gdist = (np.linalg.norm(xyz - xyz[goal], axis=1) / self._diag).astype(np.float32)
        self.start = start
        info = self._info(extra={"optimal_cost_J": float(dist[start]), "theta": th, "speed": sp, "start": start})
        return self._obs(), info

    def step(self, action):
        if self.node is None:
            raise RuntimeError("call reset() before step()")
        try:
            a = int(action)
            ok = 0 <= a < self.n_actions and float(action) == a
        except (TypeError, ValueError):
            ok = False
        if not ok:
            raise ValueError(f"action {action!r} outside Discrete({self.n_actions})")
        move, sidx = divmod(a, self.n_speeds)
        c, here = self.case, self.node
        invalid = not c.valid[here, move]
        if self.speed_choice and not invalid:
            invalid = not np.isfinite(c.cost_s[here, move, sidx])
        if invalid:
            cost_J, reward_raw = 0.0, -P_INVALID
            new = here
        else:
            cost_J = float(c.cost_s[here, move, sidx] if self.speed_choice else c.cost[here, move])
            reward_raw = -cost_J / E_REF
            new = int(self.next_node[here, move])
        self.node = new
        self.steps += 1
        self.energy_J += cost_J
        terminated = new == self.goal
        truncated = (not terminated) and self.steps >= self.max_steps
        if truncated:
            reward_raw -= P_TRUNC
        reward = reward_raw
        if self.shaping:
            reward += self.gamma * self._phi(new, terminated) - self._phi(here)
        info = self._info(extra={"invalid": invalid, "cost_J": cost_J, "reward_raw": reward_raw,
                                 "is_success": bool(terminated)})
        if terminated or truncated:
            opt = float(self.oracle_policy()[1][self.start])
            info["optimal_cost_J"] = opt
            info["regret"] = (self.energy_J - opt) / opt if terminated and opt > 0 else float("nan")
        return self._obs(), float(reward), bool(terminated), bool(truncated), info

    # ----------------------------------------------------------------------------------
    def _phi(self, node, terminal=False):
        """Potential = -(straight-line distance to goal x cheapest J/m of any edge) / E_REF. The factor is a LOWER
        bound on energy per metre, so Phi is an admissible (optimistic) estimate of -cost-to-go; Phi(goal) = 0."""
        if terminal:
            return 0.0
        return -float(self._gdist[node]) * self._diag * self.case.min_epm / E_REF

    def _obs(self):
        n, c = self.node, self.case
        o = np.empty(36, np.float32)
        o[0:3] = self._pos[n]
        o[3:6] = self._pos[self.goal] - self._pos[n]          # (g - p) / ext
        o[6:9] = c.wind10[n]
        o[9] = self._gdist[n]
        o[10:] = c.valid[n]
        return o

    def action_mask(self, node=None):
        """Boolean mask over the n_actions discrete actions that are valid edges (and speeds) from `node`."""
        node = self.node if node is None else node
        v = self.case.valid[node]
        if not self.speed_choice:
            return v.copy()
        return np.repeat(v, self.n_speeds) & np.isfinite(self.case.cost_s[node]).reshape(-1)

    def _info(self, extra=None):
        d = {"node": int(self.node), "goal": int(self.goal), "steps": self.steps, "energy_J": self.energy_J,
             "action_mask": self.action_mask()}
        if self.shaping:
            d["shaping"] = True
        d.update(extra or {})
        return d
