"""Fast tests of the RL environment on a small synthetic grid with a wall (no atlas, no CFD needed)."""
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from planner import Environment                                                          # noqa: E402
from rl_env import (WindDroneEnv, uniform_wind_factory, E_REF, P_INVALID, P_TRUNC, EDMO)  # noqa: E402

SPEEDS = np.array([5.0, 8.0, 10.0, 12.0, 15.0])


def make_env(n=25, nz=9, dx=10.0, wall=True):
    x = np.arange(n) * dx
    z = np.arange(nz) * dx
    valid = np.ones((n, n, nz), bool)
    if wall:
        valid[8:12, 4:20, 0:6] = False                                                    # a wall across the middle
    return Environment(x, x.copy(), z, valid, terrain=np.zeros((n, n)))


def make(**kw):
    env = make_env()
    base = dict(environment=env, wind_factory=uniform_wind_factory(env), stride=(2, 2, 1), margin=5.0, h_min=10.0,
                h_max=60.0, speeds=SPEEDS, min_sep=60.0, max_steps=200, speed_jitter=0.0)
    base.update(kw)
    return WindDroneEnv(**base)


WIND = {"theta": 250.0, "speed": 7.0}


def play_nodes(env, nodes):
    """Fly a node path with the env API; return (return, last info)."""
    ret, info = 0.0, None
    for a, b in zip(nodes[:-1], nodes[1:]):
        k = env._move_between(a, b)
        _, r, te, tr, info = env.step(k)
        ret += r
    return ret, info


def test_api_shapes_and_spaces():
    env = make()
    assert env.action_space.n == 26
    obs, info = env.reset(seed=0, options=WIND)
    assert obs.shape == (36,) and obs.dtype == np.float32 and np.isfinite(obs).all()
    assert env.observation_space.contains(obs)
    assert {"node", "goal", "optimal_cost_J", "action_mask"} <= set(info)
    assert info["action_mask"].shape == (26,) and (obs[10:] == info["action_mask"]).all()
    out = env.step(env.action_space.sample())
    assert len(out) == 5 and out[0].shape == (36,) and isinstance(out[1], float)
    assert isinstance(out[2], bool) and isinstance(out[3], bool)
    assert (obs[3:6] == env._pos[env.goal] - env._pos[info["node"]]).all()
    assert make(speed_choice=True).action_space.n == 26 * len(SPEEDS)
    with pytest.raises(ValueError):
        env.step(26)
    with pytest.raises(ValueError):
        env.step(-1)


def test_seed_reproducibility():
    runs = []
    for _ in range(2):
        env = make()
        seq = []
        for ep in range(3):
            obs, info = env.reset(seed=11 if ep == 0 else None)
            seq.append((info["theta"], info["speed"], info["start"], info["goal"], obs.copy()))
            for _ in range(5):
                obs, r, *_ = env.step(env.action_space.sample())
                seq.append((r, obs.copy()))
        runs.append(seq)
    for a, b in zip(*runs):
        for u, v in zip(a, b):
            assert np.array_equal(u, v)
    e1, e2 = make(), make()
    a, b = e1.reset(seed=1)[1], e2.reset(seed=2)[1]
    assert (a["start"], a["goal"], a["theta"], a["speed"]) != (b["start"], b["goal"], b["theta"], b["speed"])


@pytest.mark.parametrize("kw", [dict(), dict(cost_mode="energy", v_g=10.0), dict(speed_choice=True),
                                dict(layer_z=30.0)])
def test_oracle_policy_reaches_goal_with_return_minus_optimal(kw):
    env = make(**kw)
    for seed in range(6):
        obs, info = env.reset(seed=seed, options=WIND)
        opt = info["optimal_cost_J"]
        assert np.isfinite(opt) and opt > 0
        assert opt == pytest.approx(env.optimal_cost(info["start"], info["goal"]))
        ret, done = 0.0, False
        while not done:
            a = env.oracle_action()
            assert a >= 0 and env.action_mask()[a]
            obs, r, te, tr, info = env.step(a)
            ret += r
            done = te or tr
        assert te and not tr and info["is_success"]
        assert ret == pytest.approx(-opt / E_REF, abs=1e-6)
        assert env.energy_J == pytest.approx(opt, rel=1e-9)
        assert info["regret"] == pytest.approx(0.0, abs=1e-9)


def test_optimum_matches_planner_dijkstra():
    env = make()
    env.reset(seed=3, options=WIND)
    G = env.G
    A = G.matrix(env.wind_factory(**{"theta": 250.0, "speed": 7.0}), "energy_bestspeed")
    for seed in range(4):
        _, info = env.reset(seed=seed, options=WIND)
        d = G.solve(A, info["start"])[0][0, info["goal"]]
        assert info["optimal_cost_J"] == pytest.approx(d, rel=1e-9)


def test_oracle_is_a_lower_bound_for_other_paths_and_beats_shortest_path_in_wind():
    env = make()
    gaps = []
    for seed in range(6):
        _, info = env.reset(seed=seed, options={"theta": 250.0, "speed": 12.0})
        path, e_sp = env.shortest_path(info["start"], info["goal"])
        assert e_sp >= info["optimal_cost_J"] * (1 - 1e-12)
        gaps.append(e_sp / info["optimal_cost_J"] - 1)
    assert max(gaps) > 0.0                                   # the wind-blind path loses somewhere


def test_random_policy_is_no_better_than_oracle():
    env = make(max_steps=60)
    for seed in range(5):
        _, info = env.reset(seed=seed, options=WIND)
        oracle_ret = -info["optimal_cost_J"] / E_REF
        ret, done = 0.0, False
        while not done:
            _, r, te, tr, info = env.step(env.action_space.sample())
            ret += r
            done = te or tr
        assert ret <= oracle_ret + 1e-9
        assert info["energy_J"] >= (info["optimal_cost_J"] if te else 0.0) - 1e-9


def test_invalid_move_is_penalised_and_state_unchanged():
    env = make()
    obs, info = env.reset(seed=0, options=WIND)
    # walk to a node that has at least one invalid move (ceiling or floor of the band always do)
    bad = np.nonzero(~env.action_mask())[0]
    assert bad.size > 0
    node, steps, energy = env.node, env.steps, env.energy_J
    obs2, r, te, tr, info = env.step(int(bad[0]))
    assert r == pytest.approx(-P_INVALID) and info["invalid"] and not te
    assert env.node == node and env.energy_J == energy and env.steps == steps + 1
    assert np.array_equal(obs2, obs)


def test_truncation_and_its_cost():
    env = make(max_steps=4)
    _, info = env.reset(seed=0, options=WIND)
    bad = int(np.nonzero(~env.action_mask())[0][0])
    out = None
    for i in range(4):
        out = env.step(bad)
        assert out[3] == (i == 3) and not out[2]
    assert out[1] == pytest.approx(-P_INVALID - P_TRUNC) and out[4]["steps"] == 4
    assert np.isnan(out[4]["regret"])


def test_reaching_goal_on_the_last_allowed_step_is_termination_not_truncation():
    env = make()
    env.reset(seed=2, options=WIND)
    n = 0
    while True:                                              # length of the optimal path
        a = env.oracle_action()
        if a < 0:
            break
        _, _, te, tr, _ = env.step(a)
        n += 1
        if te:
            break
    env2 = make(max_steps=n)
    env2.reset(seed=2, options=WIND)
    for _ in range(n):
        out = env2.step(env2.oracle_action())
    assert out[2] and not out[3]


def test_goal_gives_no_bonus_reward_is_pure_energy():
    env = make()
    env.reset(seed=0, options=WIND)
    done, rs, costs = False, [], []
    while not done:
        _, r, te, tr, info = env.step(env.oracle_action())
        rs.append(r), costs.append(info["cost_J"])
        done = te or tr
    assert np.allclose(rs, -np.array(costs) / E_REF, atol=1e-12)


def test_potential_shaping_changes_return_by_a_policy_independent_constant():
    raw, shp = make(shaping=False), make(shaping=True)
    diffs = []
    for seed in range(4):
        _, i0 = raw.reset(seed=seed, options=WIND)
        shp.reset(seed=seed, options=WIND)
        assert i0["start"] == shp.start
        # two different complete paths: the oracle's and the wind-blind shortest path
        nodes_o = [i0["start"]]
        while nodes_o[-1] != i0["goal"]:
            nodes_o.append(int(raw.next_node[nodes_o[-1], raw.oracle_action(nodes_o[-1])]))
        nodes_s, _ = raw.shortest_path(i0["start"], i0["goal"])
        vals = []
        for nodes in (nodes_o, list(nodes_s)):
            raw.reset(options=dict(WIND, start=i0["start"], goal=i0["goal"]))
            shp.reset(options=dict(WIND, start=i0["start"], goal=i0["goal"]))
            vals.append(play_nodes(shp, nodes)[0] - play_nodes(raw, nodes)[0])
        assert vals[0] == pytest.approx(vals[1], abs=1e-9)           # same offset for different policies
        phi0 = -float(shp._gdist[i0["start"]]) * shp._diag * shp.case.min_epm / E_REF
        assert vals[0] == pytest.approx(-phi0, abs=1e-9)             # telescoping: sum F = Phi(goal) - Phi(start)
        diffs.append(vals[0])
    assert min(diffs) > 0


def test_layer_mode_allows_only_horizontal_moves():
    env = make(layer_z=30.0)
    _, info = env.reset(seed=0, options=WIND)
    from planner import PlanningGraph
    vertical = PlanningGraph.OFFSETS[:, 2] != 0
    assert not (info["action_mask"] & vertical).any() and info["action_mask"].sum() > 0
    assert env.G.node_xyz[env.node, 2] == 30.0


def test_sampling_uses_only_atlas_directions_and_renormalises():
    env = make(directions=[210, 240, 270])
    assert list(env.directions) == [210.0, 240.0, 270.0] and env.dir_p.sum() == pytest.approx(1.0)
    w = np.array([EDMO[d][0] for d in (210, 240, 270)])
    assert env.dir_p == pytest.approx(w / w.sum())
    seen = {env.reset(seed=s)[1]["theta"] for s in range(30)}
    assert seen <= {210.0, 240.0, 270.0} and len(seen) >= 2
    env_j = make(directions=[240], speed_jitter=0.35)
    sp = {env_j.reset(seed=s)[1]["speed"] for s in range(20)}
    assert all(1.0 <= s <= 20.0 and (2 * s) == int(2 * s) for s in sp)


def test_speed_choice_oracle_picks_the_best_speed_and_agrees_with_bestspeed_mode():
    a, b = make(speed_choice=True), make()
    _, ia = a.reset(seed=1, options=WIND)
    _, ib = b.reset(seed=1, options=WIND)
    assert ia["optimal_cost_J"] == pytest.approx(ib["optimal_cost_J"], rel=1e-9)
    assert a.oracle_action() % len(SPEEDS) in range(len(SPEEDS))
    assert a.action_mask().shape == (26 * len(SPEEDS),)
