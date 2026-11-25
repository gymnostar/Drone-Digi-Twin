"""Tabular Q-learning vs the Dijkstra oracle on a deliberately small sub-problem.

Problem: one horizontal layer (z = 40 m) of the DLR campus graph at 20 m spacing, one frozen wind (240 deg from,
8 m/s), fixed start and goal ~470 m apart. State = graph node (about 1,700), actions = the 8 horizontal moves,
reward = -energy/1000 J per edge, no goal bonus. The agent never sees the wind or the map: it learns
Q(node, move) from experience only. The oracle (one reverse Dijkstra, rl_env.WindDroneEnv.oracle_policy) gives the
exact optimum of the same discretised problem, so the final gap to the oracle is a true regret.

    python rl_demo.py                      # ~1-2 min, one process
Outputs: fig/rl_learning_curve.png, results/rl_demo.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from rl_env import WindDroneEnv, E_REF, P_TRUNC     # noqa: E402

# ---- problem ----------------------------------------------------------------------------------------------------
THETA, SPEED = 240.0, 8.0                    # wind FROM 240 deg at 8 m/s (10 m reference speed)
LAYER_Z = 40.0                               # fly in the horizontal layer at this altitude (m, absolute)
START_XYZ = (240.0, 260.0, LAYER_Z)          # chosen from a scan: the wind-blind shortest path costs +13 % here
GOAL_XYZ = (-180.0, 40.0, LAYER_Z)

# ---- hyperparameters (what each one does) -----------------------------------------------------------------------
STRIDE = (8, 8, 2)     # graph coarseness: 40 m x 40 m nodes in the layer (about 460 states). (4, 4, 2) works too but
                       # a tabular agent then needs ~5x more episodes: the table grows with the number of nodes.
EPISODES = 1500        # training episodes (one episode = one flight from start until goal or MAX_STEPS)
MAX_STEPS = 120        # time limit per episode. The optimum needs ~12 moves; 120 lets the agent wander a lot first.
ALPHA = 0.5            # learning rate: how far Q moves toward the new target. Environment is deterministic, so a
                       # large alpha is safe and fast; with stochastic winds you would decay it.
GAMMA = 1.0            # discount. 1.0 = undiscounted, so Q is exactly "minus remaining kJ" and matches the benchmark
                       # (any gamma < 1 would make the agent trade energy against path length).
EPS_START, EPS_END = 0.20, 0.02   # epsilon-greedy exploration: probability of a random valid move, decayed linearly
EPS_DECAY_EPISODES = 1000         # from EPS_START to EPS_END over this many episodes, then held
Q_INIT = 0.0           # OPTIMISTIC initialisation. Every reward is negative, so Q=0 over-promises on every untried
                       # move; the greedy agent is drawn to try them all (cheap directed exploration).
EVAL_EVERY = 25        # greedy-policy evaluation period (episodes)
WALL_BUDGET_S = 150.0  # safety stop; the demo must stay well under 3 minutes
SEED = 0


def greedy_rollout(env, Q, start, goal, max_steps=MAX_STEPS):
    """Energy (J) of the greedy policy argmax_a Q(s, a) over valid moves, from start; nan if it fails to reach goal."""
    env.reset(options={"theta": THETA, "speed": SPEED, "start": start, "goal": goal})
    seen = set()
    for _ in range(max_steps):
        node = env.node
        if node in seen:                     # greedy policy cycles: it has not learned a route yet
            return np.nan, False
        seen.add(node)
        q = np.where(env.action_mask(), Q[node], -np.inf)
        _, _, term, trunc, info = env.step(int(np.argmax(q)))
        if term:
            return env.energy_J, True
        if trunc:
            break
    return np.nan, False


def trace(env, Q, start, goal):
    """Node path of the greedy policy (for plotting)."""
    env.reset(options={"theta": THETA, "speed": SPEED, "start": start, "goal": goal})
    path = [env.node]
    for _ in range(MAX_STEPS):
        q = np.where(env.action_mask(), Q[env.node], -np.inf)
        _, _, term, trunc, _ = env.step(int(np.argmax(q)))
        path.append(env.node)
        if term or trunc:
            break
    return np.array(path)


def train(env, start, goal, optimal, shaped, t0):
    """Tabular Q-learning. shaped=True adds the potential-based shaping term (same optimal policy, see RL_ENV.md);
    the greedy evaluation always scores TRUE energy (env.energy_J), so the curves are comparable."""
    env.shaping = shaped
    rng = np.random.default_rng(SEED)
    _va = {}

    def valid_moves(n):                          # cached action mask of node n as an index array
        v = _va.get(n)
        if v is None:
            v = _va[n] = np.nonzero(env.case.valid[n])[0]
        return v

    Q = np.full((env.G.n, 26), Q_INIT)
    curve = []
    t_train = time.time()
    ep = 0
    for ep in range(EPISODES + 1):
        if ep % EVAL_EVERY == 0:
            e, ok = greedy_rollout(env, Q, start, goal)
            curve.append((ep, e))
            if ep % 500 == 0:
                print(f"  [{'shaped' if shaped else 'plain '}] ep {ep:5d}  greedy energy {e:8.0f} J  "
                      f"gap {100 * (e / optimal - 1):6.2f} %  [{time.time() - t_train:5.1f} s]")
            if time.time() - t0 > WALL_BUDGET_S:
                print("  wall budget reached, stopping early")
                break
        if ep == EPISODES:
            break
        eps = max(EPS_END, EPS_START + (EPS_END - EPS_START) * ep / EPS_DECAY_EPISODES)
        env.reset(options={"theta": THETA, "speed": SPEED, "start": start, "goal": goal})
        node, done = env.node, False
        while not done:
            va = valid_moves(node)                                          # the moves that exist from this node
            if rng.random() < eps:
                a = int(va[rng.integers(len(va))])                          # explore
            else:
                q = Q[node, va]
                a = int(rng.choice(va[q == q.max()]))                       # exploit (random tie-break)
            _, r, term, trunc, info = env.step(a)
            r_td = r + (P_TRUNC if trunc else 0.0)       # the time-limit penalty is a benchmark charge, not part of the
            nxt = info["node"]                           # MDP: do not let it leak into Q (we bootstrap on truncation)
            target = r_td if term else r_td + GAMMA * Q[nxt, valid_moves(nxt)].max()
            Q[node, a] += ALPHA * (target - Q[node, a])                     # temporal-difference update
            node, done = nxt, term or trunc
    seconds = time.time() - t_train
    env.shaping = False
    final, ok = greedy_rollout(env, Q, start, goal)
    return dict(Q=Q, curve=np.array(curve, float), final=final, ok=ok, episodes=ep, seconds=seconds)


def main():
    t_all = time.time()
    env = WindDroneEnv(stride=STRIDE, layer_z=LAYER_Z, theta=THETA, speed=SPEED, max_steps=MAX_STEPS, seed=SEED)
    env.reset(seed=SEED)
    start, goal = env.node_of(START_XYZ), env.node_of(GOAL_XYZ)
    optimal = env.optimal_cost(start, goal)                                  # J, exact (Dijkstra)
    sp_nodes, shortest = env.shortest_path(start, goal)                      # J, geometric shortest path
    opt_path = [start]
    while opt_path[-1] != goal:
        a = env.oracle_policy(goal)[0][opt_path[-1]]
        opt_path.append(int(env.next_node[opt_path[-1], a]))
    print(f"start {start} goal {goal}  optimal {optimal:.0f} J  shortest-path {shortest:.0f} J "
          f"(+{100 * (shortest / optimal - 1):.1f} %)  nodes in layer {len(env.active)}")

    runs = {}
    for shaped in (False, True):
        runs[shaped] = train(env, start, goal, optimal, shaped, t_all)
    plain = runs[False]
    final, ok = plain["final"], plain["ok"]
    gap = max(0.0, 100 * (final / optimal - 1)) if ok else float("nan")
    ep = plain["episodes"]
    learned_path = trace(env, plain["Q"], start, goal)
    wall = time.time() - t_all

    def summary(r):
        within = lambda tol: next((int(c[0]) for c in r["curve"] if np.isfinite(c[1]) and c[1] <= optimal * (1 + tol)), None)  # noqa: E731
        g = max(0.0, 100 * (r["final"] / optimal - 1)) if r["ok"] else None
        return dict(final_greedy_energy_J=float(r["final"]) if r["ok"] else None, final_gap_to_oracle_pct=g,
                    reached_goal=bool(r["ok"]), episodes=int(r["episodes"]), episodes_to_within_5pct=within(0.05),
                    episodes_to_within_1pct=within(0.01), episodes_to_exact=within(1e-9),
                    train_seconds=round(r["seconds"], 1))
    res = dict(wind={"theta_deg": THETA, "speed_mps": SPEED}, layer_z_m=LAYER_Z, stride=list(STRIDE),
               start_xyz=list(START_XYZ), goal_xyz=list(GOAL_XYZ), n_layer_nodes=int(len(env.active)),
               optimal_energy_J=optimal, shortest_path_energy_J=shortest,
               shortest_path_excess_pct=100 * (shortest / optimal - 1),
               final_gap_to_oracle_pct=float(gap), episodes=int(ep), wall_seconds=round(wall, 1),
               tabular_q_learning=summary(plain), tabular_q_learning_with_potential_shaping=summary(runs[True]),
               hyperparameters=dict(alpha=ALPHA, gamma=GAMMA, eps_start=EPS_START, eps_end=EPS_END,
                                    eps_decay_episodes=EPS_DECAY_EPISODES, q_init=Q_INIT, max_steps=MAX_STEPS,
                                    seed=SEED),
               note="gap = (greedy energy - optimal energy) / optimal energy; optimal is exact Dijkstra on the same graph. "
                    "Headline numbers are the plain agent (no shaping).")
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / "rl_demo.json").write_text(json.dumps(res, indent=2))
    print(f"final greedy {final:.0f} J, gap to oracle {gap:.2f} %, {ep} episodes, total {wall:.0f} s")
    print("with shaping:", res["tabular_q_learning_with_potential_shaping"])
    plot(env, runs, optimal, shortest, sp_nodes, np.array(opt_path), learned_path, gap)
    return res


def plot(env, runs, optimal, shortest, sp_nodes, opt_path, learned, gap):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    C_OPT, C_SP, C_RL, C_GRAY = "#1b7837", "#c0392b", "#2c6fbb", "#b5b5b5"
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw={"width_ratios": [1.25, 1]})
    top = shortest / 1000 * 1.35
    for shaped, color, lab in ((False, C_RL, "Q-learning (plain)"), (True, "#8e5ea2", "Q-learning + potential shaping")):
        curve = runs[shaped]["curve"]
        ok = np.isfinite(curve[:, 1])
        ax.plot(curve[ok, 0], curve[ok, 1] / 1000, color=color, lw=1.8, label=lab + ", greedy route")
        if (~ok).any():
            ax.plot(curve[~ok, 0], np.full((~ok).sum(), top * (0.985 if not shaped else 0.965)), "|", color=color,
                    alpha=.6, label=None if shaped else "tick row: greedy route does not reach the goal yet")
    ax.axhline(optimal / 1000, color=C_OPT, lw=2, label=f"oracle (Dijkstra) {optimal / 1000:.2f} kJ")
    ax.axhline(shortest / 1000, color=C_SP, lw=2, ls="--",
               label=f"shortest path, same wind {shortest / 1000:.2f} kJ (+{100 * (shortest / optimal - 1):.0f} %)")
    ax.set_xlabel("training episode"); ax.set_ylabel("energy of the greedy route (kJ)")
    ax.set_ylim(optimal / 1000 * 0.9, top)
    ax.set_title(f"Tabular Q-learning reaches the exact optimum (plain agent: final gap {gap:.2f} %)")
    ax.grid(alpha=.25); ax.legend(loc="upper right", fontsize=8.5, frameon=False)
    xy = env.G.node_xyz
    bx.scatter(xy[env.active, 0], xy[env.active, 1], s=3, color=C_GRAY, label="free nodes at z = 40 m")
    for p, c, ls, lw, lab in ((sp_nodes, C_SP, "--", 2, "shortest path"), (opt_path, C_OPT, "-", 3, "oracle"),
                              (learned, C_RL, ":", 2.5, "learned")):
        bx.plot(xy[p, 0], xy[p, 1], color=c, ls=ls, lw=lw, label=lab)
    bx.plot(*xy[opt_path[0], :2], "ko", ms=6); bx.plot(*xy[opt_path[-1], :2], "k*", ms=11)
    a = np.radians(env.case.theta)                                           # arrow: where the wind blows TO
    bx.annotate("", xy=(-380 + 90 * -np.sin(a), 400 + 90 * -np.cos(a)), xytext=(-380, 400),
                arrowprops=dict(arrowstyle="->", lw=2))
    bx.text(-330, 430, f"wind {env.case.theta:.0f} deg, {env.case.speed:g} m/s", fontsize=8.5)
    bx.set_aspect("equal"); bx.set_xlabel("x east (m)"); bx.set_ylabel("y north (m)")
    bx.set_title("Routes (circle = start, star = goal)"); bx.legend(loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout()
    (HERE / "fig").mkdir(exist_ok=True)
    fig.savefig(HERE / "fig" / "rl_learning_curve.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
