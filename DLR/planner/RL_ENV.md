# RL environment for the wind-aware drone planner

`rl_env.py` wraps the planner's coarse 26-connected graph as a Gymnasium-style environment (no `gymnasium`
import needed). `rl_demo.py` trains a tabular agent on it. The classical planner is the **oracle**.

## Why reward = -cost, and why no goal bonus
Reward per step = -(edge energy in J) / `E_REF` (1000 J). The undiscounted return of an episode is then
`-(total energy)/E_REF`, so maximising return is exactly minimising energy, and the planner's cost function is
the reward function (Arthur's "use the cost functions as rewards").

A goal bonus (+100 on arrival, say) would add a constant to every successful episode and break the benchmark:
return would no longer equal -cost, so you could not compare it with `optimal_cost()` without knowing the bonus,
and a large bonus drowns small energy differences (the 5-10 % we care about). Termination at the goal already
tells the agent the episode is over; the penalty for not arriving is the energy it keeps burning.

Two documented extra charges, both outside normal flight: invalid move (no edge) = position unchanged,
reward `-P_INVALID` (1 kJ, above any single edge); truncation after `max_steps` = one-off `-P_TRUNC` (30 kJ,
about 3x a campus crossing). Treat truncation as a time limit, not a terminal state, when you bootstrap
(the demo removes the penalty from the TD target for this reason).

## Why Dijkstra is the oracle, and what regret means
For one fixed, known wind the environment is a deterministic shortest-path problem with non-negative edge
costs. Dijkstra solves it exactly: `env.oracle_policy(goal)` runs one reverse Dijkstra from the goal on the
transposed cost matrix and returns the best move and the optimal cost-to-go of every node.
No agent can beat it on this graph, so

    regret = (agent energy - optimal energy) / optimal energy  >= 0

is the right metric: 0 = perfect, 0.13 = 13 % more energy than necessary. Always report regret, not raw return.
(Exact on the *discretised* problem; the graph itself approximates continuous flight.)

## What RL could add beyond Dijkstra, and what it cannot
Cannot beat: the static, known-wind case above. Do not claim "RL finds better routes" there.

Could add where Dijkstra's assumptions break:
- unknown wind: learn from experience or from onboard measurements instead of a CFD atlas;
- time-varying wind (gusts, direction shifts): a policy conditions on the current observation; re-planning is costly;
- partial observability: only the local wind is measured (this is already the observation);
- multi-drone coordination and RF-link constraints (SWARM-LINK): cost is no longer a sum of edge costs;
- continuous actions in AirSim: a controller instead of 26 grid moves.
Whenever you build such a variant, keep this env as the sanity floor: the agent must match the oracle in the easy case.

## Tables
Observation, `float32[36]`:

| index | content |
|---|---|
| 0:3 | position in the graph bounding box, [0, 1] |
| 3:6 | (goal - position) / box extent |
| 6:9 | local wind (u, v, w) / 10 m/s (x east, y north, z up) |
| 9 | straight-line distance to goal / box diagonal |
| 10:36 | 1 if move k is a valid edge from here, else 0 |

Action: `Discrete(26)`, move k = `PlanningGraph.OFFSETS[k]`. With `speed_choice=True`: `Discrete(26*n_speeds)`,
action = k * n_speeds + speed index. `info["action_mask"]` marks valid actions.

| step outcome | reward | terminated | truncated |
|---|---|---|---|
| valid move | -edge J / 1000 | node == goal | steps >= max_steps |
| invalid move | -1 (state unchanged) | no | as above |
| truncation step | extra -30 | no | yes |

Wind: sampled on `reset` from the EDMO climate (12 sectors, only those present in the atlas, renormalised;
speed = sector mean with log-normal jitter, 0.5 m/s bins), or fixed by `options={"theta", "speed", "start", "goal"}`.
Cost mode `energy_bestspeed` (default) or `energy` (fixed ground speed). Optional `layer_z` flies one altitude only.

## How to run
    cd DLR/planner
    ../../.venv/bin/python rl_demo.py                       # about 1 min, writes fig/rl_learning_curve.png
    ../../.venv/bin/python -m pytest tests/test_rl_env.py -q

    env = WindDroneEnv(atlas, stride=(4, 4, 2))             # first reset of a new wind costs ~10 s (cost table)
    obs, info = env.reset(seed=0, options={"theta": 240, "speed": 8.0})
    a = env.oracle_action(); env.optimal_cost(info["start"], info["goal"]); env.shortest_path(start, goal)

Demo result (240 deg, 8 m/s, one layer at z = 40 m, 460 nodes): the wind-blind shortest path costs +11.5 %; plain
Q-learning reaches the oracle exactly after about 700 episodes; with potential shaping after about 400.

## Potential-based shaping (flag `shaping=True`, default off)
Adds F(s, s') = gamma*Phi(s') - Phi(s) to the reward, with Phi = -(distance to goal x cheapest J/m)/E_REF.
Ng, Harada & Russell (1999): any shaping of this form leaves the optimal policy unchanged, because along any
full path the F terms telescope to Phi(goal) - Phi(start), a constant independent of the route (the tests check
this). It only gives the learner a sense of direction, so it learns faster. `info["energy_J"]` stays the true
energy, so score regret from that, not from the shaped return.

## Swapping in PPO (stable-baselines3) later
SB3 needs a `gymnasium.Env`; the class is duck-typed, so wrap it (about 10 lines, after `pip install gymnasium
stable-baselines3`):

    class GymWindDrone(gymnasium.Env):
        def __init__(self, **kw):
            self.e = WindDroneEnv(**kw)
            self.action_space, self.observation_space = self.e.action_space, self.e.observation_space
        def reset(self, *, seed=None, options=None): return self.e.reset(seed=seed, options=options)
        def step(self, a): return self.e.step(a)

Use `MaskablePPO` (sb3-contrib) with `info["action_mask"]` or rely on `-P_INVALID`. Train on the climate wind,
then evaluate the regret per (start, goal, wind) against `optimal_cost`. Start with the demo's single layer.

## Limits
- Tabular Q-learning memorises one start/goal/wind; it shows the benchmark works, it does not generalise.
- Reward resolution is the graph (20 m x 20 m x 5 m at stride (4,4,2)); the oracle is optimal only on that graph.
- Wind is frozen within an episode and perfectly known to the oracle; the agent sees only the local wind.
- Wind cases are cached (6 kept); sampled speeds are binned so cases repeat. `energy_bestspeed` on the full graph is slow
  to build (seconds); use `layer_z` or a coarser stride for quick experiments.
- Energy model and its assumptions are in ENERGY_MODEL.md; tuned constants `P_INVALID`, `P_TRUNC` are design choices.
