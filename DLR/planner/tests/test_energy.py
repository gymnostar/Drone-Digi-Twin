"""Fast tests for the propulsion (energy.py) and battery (battery.py) models."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import battery as B  # noqa: E402
import energy as E  # noqa: E402

X = np.array([1.0, 0.0, 0.0])


def test_hover_power_zeng():
    assert E.power_level(0.0) == pytest.approx(168.5, abs=0.5)
    assert E.ZENG2019.P0 + E.ZENG2019.Pi == pytest.approx(168.5, abs=0.5)


def test_power_is_u_shaped_with_expected_speeds():
    V = np.linspace(0.0, 30.0, 3001)
    P = E.power_level(V)
    i = np.argmin(P)
    assert 0 < i < len(V) - 1                       # interior minimum
    assert np.all(np.diff(P[: i + 1]) <= 1e-9) and np.all(np.diff(P[i:]) >= -1e-9)
    assert 9.0 < V[i] < 11.5                        # max endurance ~10 m/s
    assert 17.0 < V[1 + np.argmin(P[1:] / V[1:])] < 19.5   # max range ~18 m/s
    assert E.max_endurance_speed() == pytest.approx(V[i], abs=0.1)


def test_wind_effect_on_best_energy_and_speed():
    e_still, v_still = E.best_energy_per_metre(X, np.zeros(3))
    e_head, v_head = E.best_energy_per_metre(X, np.array([-5.0, 0, 0]))   # wind against motion
    e_tail, v_tail = E.best_energy_per_metre(X, np.array([5.0, 0, 0]))
    assert e_head > e_still > e_tail
    # Physics: optimum AIRSPEED is ~constant (max-range speed). Ground speed = airspeed -/+ wind, so the best
    # GROUND speed is LOWER into a headwind and HIGHER with a tailwind; airspeed is higher into a headwind.
    assert v_head < v_still < v_tail
    assert v_head + 5 > v_still > v_tail - 5 - 1e-9


def test_drag_force_antiparallel_and_magnitude():
    d = E.ZENG2019
    v = np.array([3.0, -4.0, 12.0])
    f = E.drag_force(v)
    assert np.dot(f, v) < 0
    cross = np.cross(f, v)
    assert np.allclose(cross, 0, atol=1e-9)
    assert np.linalg.norm(f) == pytest.approx(0.5 * d.rho * d.d0 * d.s * d.A * np.linalg.norm(v) ** 2)
    assert np.allclose(E.drag_force(np.zeros(3)), 0)


def test_drag_power_equals_parasite_term():
    V = 15.0
    _, _, par = E.power_components(V)
    assert -np.dot(E.drag_force(np.array([V, 0, 0])), [V, 0, 0]) == pytest.approx(float(par))


def test_climb_costs_weight_times_rate():
    d = E.ZENG2019
    assert E.power(12.0, 2.0) - E.power(12.0, 0.0) == pytest.approx(d.W * 2.0)
    assert E.power(12.0, -3.0) == pytest.approx(E.power(12.0, 0.0))   # no credit for descent


def test_vectorised_shapes_broadcast():
    V = np.linspace(0, 20, 7)
    assert E.power_level(V).shape == (7,)
    assert E.power(V[:, None], np.array([0.0, 1.0, 2.0])).shape == (7, 3)
    dirs = np.tile(X, (4, 1))
    winds = np.zeros((4, 3))
    assert E.segment_power(dirs, 10.0, winds).shape == (4,)
    e, v = E.best_energy_per_metre(dirs, winds)
    assert e.shape == v.shape == (4,)
    assert E.power_momentum(np.zeros((5, 2))).shape == (5,)


def test_airsim_set_and_momentum_check():
    a = E.AIRSIM_GENERIC_QUAD
    assert a.W == pytest.approx(9.81)
    assert 40 < E.power_level(0.0, a) < 80
    # first-principles momentum theory and Zeng closed form agree within 10% at moderate speeds
    V = np.linspace(0, 25, 26)
    for d in (E.ZENG2019, a):
        m = E.power_momentum(np.stack([V, 0 * V], -1), d)
        assert np.all(np.abs(m / E.power_level(V, d) - 1) < 0.10)


# ----------------------------------------------------------------------------------- battery
def test_cell_energy_not_less_than_electrical_and_loss_positive():
    for P in (100, 200, 400):
        r = B.energy_constant_power(P, 300)
        assert r.e_cell >= r.e_elec
        assert r.e_loss > 0
        assert r.e_elec == pytest.approx(P / B.PACK_6S_5000.eta * 300, rel=1e-9)


def test_loss_grows_superlinearly_with_power():
    # same duration; loss(2P) > 2 loss(P)
    l1 = B.energy_constant_power(150, 120).e_loss
    l2 = B.energy_constant_power(300, 120).e_loss
    assert l2 > 2 * l1


def test_soc_monotonic_decreasing_and_bounded():
    rng = np.random.default_rng(0)
    p = 100 + 200 * rng.random(300)
    r = B.simulate(p, 1.0)
    assert np.all(np.diff(r.soc) < 0)
    assert 0 < r.soc[-1] < 1 and r.soc_used == pytest.approx(1 - r.soc[-1])
    assert not r.cutoff


def test_cutoff_flag_for_absurd_power():
    r = B.energy_constant_power(5000, 30)
    assert r.infeasible and r.cutoff
    assert r.t_cutoff > 0
    assert not B.energy_constant_power(150, 60).cutoff


def test_cutoff_triggers_when_pack_is_nearly_empty():
    r = B.energy_constant_power(250, 60, soc0=0.03)
    assert r.cutoff


def test_cold_pack_loses_more_and_rc_model_runs():
    warm = B.energy_constant_power(250, 200)
    cold = B.energy_constant_power(250, 200, B.PACK_COLD)
    assert cold.e_loss > 2 * warm.e_loss
    rc = B.energy_constant_power(250, 200, B.PACK_6S_5000_RC)
    assert rc.e_cell >= rc.e_elec


def test_peukert_option_increases_soc_use():
    p0 = B.Pack(peukert=1.0)
    p1 = B.Pack(peukert=1.05)
    s0 = B.energy_constant_power(400, 120, p0).soc_used
    s1 = B.energy_constant_power(400, 120, p1).soc_used
    assert s1 > s0
