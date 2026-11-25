"""Fast tests of the wind-aware planner on small synthetic grids (no atlas needed)."""
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.sparse.csgraph import dijkstra

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from planner import Environment, WindField, PlanningGraph, evaluate_path  # noqa: E402
from energy import energy_per_metre  # noqa: E402


def make_env(n=31, nz=13, block=None, dx=5.0):
    x = np.arange(n) * dx
    z = np.arange(nz) * dx
    valid = np.ones((n, n, nz), bool)
    if block is not None:
        (x0, x1), (y0, y1), (z0, z1) = block
        valid[x0:x1, y0:y1, z0:z1] = False
    return Environment(x, x.copy(), z, valid, terrain=np.zeros((n, n)))


def graph(env, **kw):
    kw = dict(stride=(1, 1, 1), margin=5.0, h_min=5.0, h_max=60.0) | kw
    return PlanningGraph(env, **kw)


def plan(G, A, a, b):
    return G.plan(A, np.array(a, float), np.array(b, float))


def test_zero_wind_energy_path_is_a_shortest_path():
    env = make_env(); G = graph(env); wind = WindField.uniform(env, [0, 0, 0])
    a, b = [10, 10, 30], [140, 90, 30]
    pd, _ = plan(G, G.matrix(wind, "distance"), a, b)
    pe, ce = plan(G, G.matrix(wind, "energy", v_g=10.0), a, b)
    L = lambda p: np.linalg.norm(np.diff(p, axis=0), axis=1).sum()
    assert L(pe) == pytest.approx(L(pd), rel=1e-9)
    # level flight in still air: cost = (J/m at 10 m/s) x length
    epm = energy_per_metre(np.array([1.0, 0, 0]), 10.0, np.zeros(3))
    assert ce == pytest.approx(float(epm) * L(pe), rel=1e-9)
    assert evaluate_path(pe, wind, env, "energy")["energy_J"] == pytest.approx(ce, rel=1e-6)


def test_uniform_wind_energy_path_stays_straightish_and_beats_shortest():
    env = make_env(); G = graph(env); wind = WindField.uniform(env, [3, 2, 0])
    a, b = [10, 10, 30], [140, 120, 30]
    pd, _ = plan(G, G.matrix(wind, "distance"), a, b)
    pe, _ = plan(G, G.matrix(wind, "energy"), a, b)
    straight = np.linalg.norm(np.array(b) - np.array(a))
    L = lambda p: np.linalg.norm(np.diff(p, axis=0), axis=1).sum()
    assert L(pe) <= 1.09 * straight                      # 26-connectivity: at most ~8% detour
    Ee = evaluate_path(pe, wind, env)["energy_J"]; Ed = evaluate_path(pd, wind, env)["energy_J"]
    assert Ee <= Ed * (1 + 1e-9)                         # the energy path is optimal over all graph paths


def test_paths_never_enter_inflated_obstacles():
    env = make_env(block=((10, 20), (5, 25), (0, 9)))    # a 50 x 100 x 45 m wall in the middle
    G = graph(env)
    wind = WindField.uniform(env, [-4, 1, 0])
    a, b = [10, 70, 30], [145, 70, 30]
    for mode in ("distance", "energy", "energy_bestspeed"):
        p, _ = plan(G, G.matrix(wind, mode), a, b)
        assert p is not None
        idx = np.rint(p / 5.0).astype(int)
        assert (env.dist[idx[:, 0], idx[:, 1], idx[:, 2]] > G.margin).all()
        ps = G.smooth(p, wind, mode)
        m = evaluate_path(ps, wind, env, "energy" if mode != "energy_bestspeed" else mode)
        assert m["min_clearance_m"] > 0                 # smoothed path never touches a solid node
        dense = np.vstack([np.linspace(u, v, 60) for u, v in zip(ps[:-1], ps[1:])])
        i = np.rint(dense / 5.0).astype(int)
        assert env.valid[i[:, 0], i[:, 1], i[:, 2]].all()


def test_directed_costs_downwind_cheaper_than_upwind():
    env = make_env(); G = graph(env); wind = WindField.uniform(env, [4, 0, 0])   # blowing toward +x
    a, b = [10, 70, 30], [140, 70, 30]
    # NB: at v_g = 10 m/s (~ the max-endurance airspeed, P(V) is flat there) the sign is not guaranteed;
    # at a max-range-like cruise of 18 m/s the headwind side of the power curve is steep.
    A = G.matrix(wind, "energy", v_g=18.0)
    _, c_down = plan(G, A, a, b)
    _, c_up = plan(G, A, b, a)
    assert c_down < c_up
    s, t = G.snap(np.array([a, b]))
    assert not np.isclose(A[s, t], A[t, s]) or A[s, t] == 0   # edge costs themselves are asymmetric
    Ab = G.matrix(wind, "energy_bestspeed")
    assert plan(G, Ab, a, b)[1] < plan(G, Ab, b, a)[1]


def test_bestspeed_never_costs_more_than_fixed_speed():
    env = make_env(); G = graph(env); wind = WindField.uniform(env, [-3, 2, 0.5])
    assert (G.costs(wind, "energy_bestspeed") <= G.costs(wind, "energy", v_g=10.0) + 1e-9).all()


def test_dijkstra_matches_bruteforce_on_tiny_graph():
    env = make_env(n=5, nz=4, block=((2, 3), (1, 4), (0, 2)))
    G = graph(env, margin=2.0, h_min=0.0, h_max=100.0)
    rng = np.random.default_rng(1)
    U = rng.normal(0, 3, env.valid.shape + (3,)).astype(np.float32); U[~env.valid] = np.nan
    wind = WindField(env.x, env.y, env.z, U)
    A = G.matrix(wind, "energy", v_g=8.0)
    n = G.n
    D = np.full((n, n), np.inf); np.fill_diagonal(D, 0)
    s, t = A.nonzero(); D[s, t] = np.asarray(A[s, t]).ravel()
    for k in range(n):                                   # Floyd-Warshall
        D = np.minimum(D, D[:, [k]] + D[[k], :])
    d = dijkstra(A, directed=True, indices=np.arange(n))
    assert np.allclose(d, D)
    assert not np.allclose(D, D.T)                       # directed: costs are not symmetric


def test_evaluate_path_uniform_wind_matches_analytic_and_sampling_is_independent_of_ds():
    env = make_env(); wind = WindField.uniform(env, [3, -1, 0])
    P = np.array([[10, 10, 30], [100, 50, 30], [140, 140, 40.0]])
    m1 = evaluate_path(P, wind, env, "energy", ds=1.0); m2 = evaluate_path(P, wind, env, "energy", ds=7.0)
    assert m1["energy_J"] == pytest.approx(m2["energy_J"], rel=1e-9)
    d = np.diff(P, axis=0); L = np.linalg.norm(d, axis=1)
    e = sum(energy_per_metre(di / li, 10.0, np.array([3, -1, 0.0])) * li for di, li in zip(d, L))
    assert m1["energy_J"] == pytest.approx(float(e), rel=1e-9)
    assert m1["length_m"] == pytest.approx(L.sum())


def test_smoothing_does_not_lengthen_and_stays_free():
    env = make_env(block=((10, 20), (5, 25), (0, 9))); G = graph(env)
    wind = WindField.uniform(env, [0, 0, 0])
    p, _ = plan(G, G.matrix(wind, "distance"), [10, 70, 30], [145, 70, 30])
    ps = G.smooth(p, wind, "distance")
    L = lambda q: np.linalg.norm(np.diff(q, axis=0), axis=1).sum()
    assert L(ps) <= L(p) + 1e-9 and len(ps) <= len(p)
    assert np.allclose(ps[0], p[0]) and np.allclose(ps[-1], p[-1])
