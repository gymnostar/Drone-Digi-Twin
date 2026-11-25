# Wind-aware path planner (DLR campus)

Files: `planner.py` (library), `run_experiments.py` (comparison), `tests/test_planner.py`,
`results/plan_results.csv`, `results/plan_summary.json`, `fig/plan_*.png`. Energy model: `energy.py`.

## 1. Problem

Arthur's question: a drone that only minimises path length vs a drone that uses the wind lookup table
(the CFD atlas). Both fly between the same two waypoints; both are scored in the SAME true wind. If the
energy-aware path really saves energy, the cost function is good, and later it becomes an RL reward.
## 2. Graph formulation

- Atlas grid (5 m, 5 m, 2.5 m) is **coarsened** to 10 x 10 x 5 m nodes (every 2nd fine node, so wind at a node
  is an exact atlas value). About 128k free nodes.
- **Free node** = no solid fine node within 5 m (Euclidean distance transform of the air mask = obstacle
  inflation) AND between 10 and 100 m above bare terrain (and below the atlas ceiling, 130 m).
- **26-connectivity**: every node links to its 3x3x3 neighbours (3.1 M directed edges). An edge exists only if
  both nodes and the midpoint are free. Edges are built by looping over the 26 offsets, not over nodes.
- Wind on an edge = mean of its two node winds. Costs are computed per wind field, structure is shared.
## 3. Why directed edges

Airspeed = ground velocity - wind and power depends on airspeed, so cost(i,j) != cost(j,i) (headwind vs tailwind):
the graph is a digraph and A->B, B->A are different problems. With zero wind it is symmetric = "shortest path".
## 4. Cost modes

| mode | edge cost | meaning |
|---|---|---|
| `distance` | L | shortest-path baseline |
| `energy` | L * E/m(dir, v_g, wind) | fixed ground speed (10 m/s default); Zeng 2019 power / v_g |
| `energy_bestspeed` | L * min_v E/m(dir, v, wind) | drone picks its speed per edge (1 m/s steps) |
| `time` | L / v_ground | fixed airspeed 10 m/s, ground speed from |v d - w| = V_air |

Climbing adds W*c, descent gets no credit (`energy.py`). `evaluate_path` re-scores any polyline densely (wind
every 2 m, trilinear): length, energy, time, headwind, clearance. It agrees with the graph cost within ~0.1 %.

## 5. Why Dijkstra is exact

Edge costs are >= 0. Dijkstra settles nodes in order of increasing cost, so when the goal is settled no cheaper
path can exist: it returns the global optimum **of the discretised problem** (scipy `csgraph.dijkstra`, with
`limit` = cost of the shortest path re-priced in the energy cost, a valid upper bound, to prune the search).
A* gives the same answer faster if the heuristic never overestimates: for distance use the Euclidean distance;
for energy use Euclidean distance x (minimum J/m over all directions and the largest wind in the field), which
is admissible but weak. It is worth it only for single queries; for 40 pairs x 3 planners Dijkstra is fine.
Energy is exact only on the grid: the continuous optimum can be a little lower.

## 6. Discretisation artifacts (read before trusting small savings)

- 26-connectivity only allows 26 headings, so graph paths are up to ~8 % longer than the Euclidean line
  (measured: mean 1.062, max 1.088 for the 40 pairs). They are staircases (see the maps).
- Many equal-length staircases exist; the energy planner picks the one that happens to suit the wind, so the **raw**
  saving partly measures a better staircase and overstates the wind effect (raw paths carry ~6 % zig-zag excess).
- **Smoothing** = greedy line-of-sight string pulling on the free mask. The baseline (shortest) path gets the plain
  geometric version ("fly straight", what a naive drone does). Energy paths get the **energy-aware** version: a
  shortcut is accepted only if collision-free AND not costlier in the true wind (dense 2 m sampling, same speed
  policy as the path). Plain geometric smoothing of an energy path would erase its wind-aware detours and
  understate it. **`fair_*` (headline) is therefore identical in value to `sm_*` on purpose**, not a copy:
  fair = smoothed shortest vs energy-aware smoothed energy path. `raw_*` is kept for transparency.

## 7. Results (N = 40 pairs, both legs A->B and B->A: 80 paths per row)

Endpoints: random free nodes 15-30 m above terrain within 60 m of a building, separation 251-752 m (seed 0, the same
pairs for every direction). Saving = 1 - E(energy path)/E(shortest path), both scored in the same true field in the
same flight regime. **Directions available when run: 210, 240, 270** (EDMO mean speed 3.02 / 5.15 / 4.31 m/s;
"8 m/s" = strong wind). Run `python run_experiments.py` again as more appear (results are cached per direction).

FAIR saving %, mean / median (p90), extra length %:

| wind | fixed 10 m/s | best speed per edge | extra len (best) |
|---|---|---|---|
| 210 @ 3.02 | 0.19 / 0.0 (1.1) | 0.48 / 0.0 (1.4) | -0.2 |
| 240 @ 5.15 | 0.45 / 0.0 (1.8) | 1.64 / 0.2 (4.9) | +0.5 |
| 270 @ 4.31 | 0.38 / 0.0 (1.5) | 0.65 / 0.0 (2.5) | +0.1 |
| 210 @ 8 | 2.32 / 1.6 (6.6) | 5.06 / 3.2 (12.8) | +2.2 |
| 240 @ 8 | 2.65 / 2.4 (5.9) | 5.33 / 3.4 (12.8) | +2.3 |
| 270 @ 8 | 1.67 / 0.0 (5.8) | 4.33 / 1.6 (12.7) | +2.0 |

Raw (unsmoothed) best-speed means: 1.7 / 3.8 / 2.4 % (EDMO), 7.5 / 7.6 / 6.8 % (8 m/s).
**EDMO-frequency-weighted expected saving** (weights renormalised over 210/240/270: 0.32/0.47/0.21):
fixed 10 m/s 0.35 % fair (1.5 % raw); best speed 1.06 % fair (2.8 % raw). At 8 m/s everywhere: 2.3 % / 5.0 % fair.

### Separating route, speed and travel direction (mean % per leg, straight path = smoothed shortest)

| effect (key in plan_summary.json) | EDMO weighted | 8 m/s weighted |
|---|---|---|
| `route_only`: energy-aware route vs straight, both at 10 m/s | 0.35 | 2.3 |
| `route_given_best_speed`: same, both at best speed | 1.06 | 5.0 |
| `speed_only`: straight path, best speed vs fixed 10 m/s (**mostly not wind**, see below) | 31.3 | 32.9 |
| `speed_choice_still_air`: straight path at constant max-range speed 18.3 m/s vs 10 m/s (no wind logic) | 30.1 | 29.6 |
| `wind_speed_adaptation`: best speed per segment in the wind vs constant 18.3 m/s | 1.9 | 5.8 |
| `route_and_speed`: energy route + best speed vs straight at 10 m/s | 32.1 | 37.1 |
| `direction_asym_fix10` / `_best`: straight path A->B vs B->A energy difference | 3.0 / 39.6 | 7.2 / 70.8 |

How to read it:
1. Route choice: wind-aware routing gives only ~0.3-1 % at typical EDMO speeds and 2-5 % at 8 m/s, with <= 2 % more
   length and no time penalty. It grows with wind speed. Small compared with the next point.
2. Speed: 10 m/s is near the max-endurance speed, where energy per metre is 12.6 J/m against 8.8 J/m at 18.3 m/s
   even in still air. The ~30 % "speed_only" gain is therefore the efficient-speed gain, nothing to do with the wind.
   The wind-related part of speed is `wind_speed_adaptation` (1.9 % EDMO, 5.8 % at 8 m/s).
3. Travel direction: at a fixed speed the straight path costs 3 % (EDMO) to 7 % (8 m/s) more one way than the other;
   with speed adaptation the A->B vs B->A gap reaches 40-70 %, because energy per metre in a tailwind is much
   lower than in a headwind. Round trips (A->B->A) save about as much as one-way legs (see `roundtrip_save_mean`).

## 8. Limits

Steady RANS wind (no gusts, no turbulence k in the cost yet); quasi-steady energy (no acceleration or turning costs, no battery model); coarse 10 m nodes with 5 m inflation (narrow gaps closed); only 3 of 12
directions at the time of writing (blending between solved directions is a weaker approximation); direction-
independent EDMO speeds; the strip pairs are mostly SW-NE, which is a specific geometry relative to the wind.

## 9. Next steps

1. Re-run `python run_experiments.py` when more directions exist (globs the atlas; multiples of 30 deg only; only
   new or changed directions are computed, `--force` recomputes all; cache in `results/cache/`, keyed on npz mtime,
   seed, n, strong speed and code; ~2-3 min per wind field under load).
2. Energy/battery term (`battery.py`), then an RF-source localisation term and a turbulence (k) penalty in the
   edge cost; weigh terms and re-check with the same comparison.
3. RL: reward per step = -(edge cost) with the same `energy`/`evaluate_path` code; the Dijkstra optimum is the
   benchmark an RL policy must approach (a clean upper bound on achievable return).
4. Non-greedy smoothing and A* search for real-time use.

Tests: `.venv/bin/python -m pytest DLR/planner/tests/test_planner.py`.
