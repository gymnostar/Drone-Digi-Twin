# Drone energy model: propulsion, drag, battery

Code: `energy.py` (propulsion + drag + AirSim set + momentum check), `battery.py` (LiPo pack),
`tests/test_energy.py`. Figures in `fig/energy_*.png`.

## 1. Propulsion power (Zeng, Xu & Zhang 2019, IEEE TWC 18(4), eq. 12)

Level flight at **airspeed** V (speed relative to the air, not the ground):

    P(V) = P0 (1 + 3V^2/U_tip^2)  +  Pi (sqrt(1 + V^4/(4 v0^4)) - V^2/(2 v0^2))^(1/2)  +  1/2 d0 rho s A V^3
           blade profile             induced (makes the lift)                              parasite (fights drag)

with P0 = delta/8 rho s A Omega^3 R^3, Pi = (1+k) W^1.5 / sqrt(2 rho A), v0 = sqrt(W/(2 rho A)).

- **Induced power falls with speed.** Hovering, the rotor must push a small mass of air down fast. Moving forward it
  meets fresh air, so it can push a larger mass down more slowly: same thrust, less kinetic energy (momentum theory).
  Hover is the worst case (Pi = 88.6 W); at 10 m/s it is 35 W.
- **Parasite power rises as V^3.** Drag is 1/2 rho CdA V^2 and power is drag x speed.
- **Profile power** is almost flat (the blades spin at constant rate).
- Sum: a U. Minimum power (max **endurance**) at 10.2 m/s; minimum power **per metre** P/V (max **range**) at 18.3 m/s,
  where the tangent from the origin touches the curve. Above that, each extra metre costs more than it saves in time.
- **Headwind/tailwind:** wind enters only through airspeed, v_air = v_ground - wind. The optimal AIRSPEED is
  still about 18 m/s, so the best GROUND speed is 15.8 m/s into a 5 m/s headwind (airspeed 20.8: fly faster through
  the air) and 21.6 m/s with a tailwind. Drifting with the wind at v_g = w costs hover power, not zero.
- Vertical extension (not in Zeng): climbing at air-relative rate c adds W*c; descent gets no credit (conservative).

| ZENG2019 (Table I) | value | | value |
|---|---|---|---|
| W (weight) | 20 N (~2 kg) | solidity s | 0.05 |
| rotor radius R / disc area A | 0.4 m / 0.503 m^2 | fuselage drag ratio d0 | 0.6 |
| blade rate Omega (U_tip) | 300 rad/s (120 m/s) | induced factor k | 0.1 |
| air density rho | 1.225 kg/m^3 | profile drag coeff. delta | 0.012 |

Results: hover 168.5 W (P0 79.9 + Pi 88.6); P(10.2) = 126.0 W; P(18.3) = 161.6 W; best 8.83 J/m (still air),
11.86 J/m (5 m/s headwind, at 15.8 m/s ground), 6.86 J/m (tailwind, at 21.6 m/s). See `energy_per_metre_wind.png`,
`energy_power_breakdown.png`.

## 2. Drag function

`drag_force(v_air)` returns D = -1/2 rho (d0 s A) |v| v (N, antiparallel to the airspeed vector); CdA = d0 s A
= 0.0151 m^2 for Zeng. Its power D.v is exactly the parasite term of P(V), so the drag function and the energy
cost are one model, not two. Drag is also what tilts the thrust vector in forward flight (T sin(theta) = D), which
is how a flight controller / AirSim would see it.

## 3. Battery model (`battery.py`)

Equivalent circuit: terminal V = V_oc(SOC) - V_rc - I R0, optional RC pair (R1, tau) for polarisation. Given
electrical power P (= P_prop / eta + p_aux) it solves P = (V_oc - I R0) I for I each step, integrates SOC, and
reports current, terminal voltage, SOC, energy delivered, energy drawn from the cells (V_oc I dt, includes I^2 R
losses), and flags: `cutoff` (V_term < 3.3 V/cell) and `infeasible` (no real root: demand above V_oc^2/4R0).

| parameter | default | source |
|---|---|---|
| pack | 6S, 5000 mAh, ~111 Wh, 22.2 V nom. | chosen: ~0.65 kg, 2 C at hover, ~31 min hover endurance |
| OCV(SOC) | 4.20 V (100%) ... 3.27 V (0%) per cell | widely circulated hobby-LiPo table, approximation |
| R0 | 10 mOhm/cell (60 mOhm pack) | typical 5 Ah 25C LiPo datasheets, 5-15 mOhm/cell, assumed |
| motor+ESC efficiency eta | 0.80 | assumed (~0.85 x 0.95) |
| cutoff | 3.3 V/cell | common low-voltage cutoff |
| cold / aged pack | R x3 (`PACK_COLD`) | rule of thumb 2-4x at ~0 C, assumed |
| Peukert exponent | 1.0 (off); LiPo 1.00-1.05 | Doerffel & Sharkh 2006; effect is small for LiPo |

Model structure follows Chen & Rincon-Mora (2006) and Tremblay & Dessaint (2009). Peukert (capacity shrinks at
high current) is minor for LiPo; the dominant rate effect is the resistive loss I^2 R0, which is in the model.

## 4. Does the battery change the answer? (Stefano's question)

1 km straight trip at constant ground speed, 6S/5 Ah, eta 0.8 (still air; kJ):

| ground speed | 4 | 8 | 12 | 16 | 20 | 24 m/s |
|---|---|---|---|---|---|---|
| propulsion | 37.7 | 16.1 | 10.6 | 9.0 | 8.9 | 9.7 |
| electrical | 47.1 | 20.2 | 13.3 | 11.3 | 11.1 | 12.1 |
| from cells | 48.0 | 20.5 | 13.5 | 11.5 | 11.4 | 12.4 |
| SOC used (%) | 10.7 | 4.5 | 3.0 | 2.5 | 2.5 | 2.8 |

Headwind 5 m/s: cells 40.3 / 20.7 / 16.0 / 15.2 / 16.0 / 18.0 kJ at the same speeds (best 15.2 kJ, 3.4% SOC);
tailwind 5 m/s: 53.3 / 25.1 / 14.0 / 10.0 / 8.8 / 8.9 kJ (best 8.7 kJ, 1.9% SOC).

- **Optimal speed shift** (0.1 m/s grid): the cell-energy optimum is 0.2 m/s lower than the propulsion optimum in all
  three wind cases (18.3 to 18.1, 15.8 to 15.6, 21.6 to 21.4 m/s). Using the propulsion-optimal speed costs only
  0.008-0.017% extra cell energy. I^2R is 1.8-2.3% of the cell energy at the optimum.
- **Same propulsion energy, different power profile** (still air, 1 km): flying 14 m/s (134 W, 71 s) vs 23.6 m/s
  (226 W, 42 s) both need 9.56 kJ of propulsion energy; cells deliver 12.15 vs 12.29 kJ (+1.2% for the fast one).
  12 vs 27 m/s (128 W vs 288 W): +2.1%. Pulsing the same mean power as +-50% (30 s blocks): +0.6% cell energy,
  +-100%: +2.4%.
- **When it matters** (`energy_battery_losses.png`, right): cold/aged pack (R x3) at 20% SOC: loss share rises to 8%
  at the optimum; the 14 vs 23.6 m/s pair then differs by 5.7%, 12 vs 27 m/s by 11%; the +-100% pulsing case
  costs 19% more than constant power. Also: higher C-rate (4S pack, small battery, heavy payload), SOC below ~20%
  (voltage sag toward cutoff, power-limited manoeuvres), aggressive acceleration.

**Conclusion:** for this vehicle (about 2 C), a warm pack and a speed-optimal planner, the battery is a nearly
linear 2% overhead; minimising propulsion energy already minimises the cell energy to within 0.02%, so Arthur's
statement holds. A battery term is worth adding only for low-SOC/cold/high-C-rate regimes, or as a hard
constraint: SOC reserve, cutoff and `max_power` check.

## 5. AirSim set (`AIRSIM_GENERIC_QUAD`)

AirSim SimpleFlight `setupFrameGenericQuad` (read from GitHub microsoft/AirSim main): mass 1.0 kg, arm 0.2275 m
(F450), body box 0.18 x 0.11 x 0.04 m, propeller 0.2286 m, C_T 0.109919, C_P 0.040164, max 6396.667 rpm,
rho 1.225, `linear_drag_coefficient` 0.325 (AirSim drag F = Cd A rho v^2 per face; face list from memory).
AirSim has no battery or power model, so the Zeng parameters are derived or assumed.

| Zeng param | value | how |
|---|---|---|
| W | 9.81 N | 1 kg x g |
| R / A (4 rotors) | 0.1143 m / 0.1642 m^2 | D/2; 4 pi R^2 |
| Omega | 513 rad/s (4900 rpm, 77% of max) | from T = C_T rho n^2 D^4 = W/4 |
| s | 0.12 | assumed (2 blades, ~22 mm chord) |
| d0 | 0.447 | so that d0 s A = CdA = 0.0088 m^2 (AirSim drag on frontal area) |
| k, delta | 0.1, 0.012 | assumed (Zeng) |

Hover 60.6 W (P0 7.3 + Pi 53.3); AirSim's own C_P gives 66.9 W shaft power (Zeng-model value is 10% lower, plausible).
Max endurance 10.8 m/s, max range 15.4 m/s at 2.96 J/m. This set is 4x lighter, so it is slower and 3x cheaper.

## 6. Momentum cross-check (`power_momentum`)

Thrust balances weight + drag (T = sqrt(W^2 + D^2)), Glauert induced velocity solved numerically,
P = (1+k) T v_i + D V + profile power (4.65 mu^2 term, helicopter value). Agreement with Zeng: within 0.6% below
10 m/s, 2.3% at most up to 30 m/s (`energy_momentum_check.png`); with Zeng's coefficient 3: within 0.7%. This is a
consistency check (Zeng derives from the same theory), not validation against flight data.

## 7. Assumptions and limits

- Quasi-steady level flight: no acceleration, turn, or tilt-dynamics cost. Climb/descent only via W*c.
- Zeng is a single-rotor equivalent (all four rotors lumped); fine for quads, ignores rotor-rotor interference.
- Air density 1.225 (sea level). The DLR site (~590 m) is ~1.16: hover -1% (166.7 W), best J/m -4% (8.47),
  best-range speed 18.3 to 18.5 m/s. Set `Drone(rho=1.16)` to use it.
- Battery: constant OCV curve, no temperature or ageing model beyond `r_scale`; efficiency eta constant; avionics
  `p_aux` default 0 (10 W moves the optimum +0.4 m/s, 30 W: +1.1 m/s).
- Time-varying wind during a leg is not modelled; the planner uses a mean wind per edge.

## 8. How the planner uses it

Edge cost = P(|v_g d - w|) / v_g x length, with d the unit direction, w the local wind from the CFD cell, v_g
the ground speed (fixed, or minimised over `SPEEDS` by `best_energy_per_metre`). Through `segment_power` this
adds climb and the horizontal/vertical airspeed split automatically. Optionally convert to SOC with
`battery.energy_constant_power(P, T)` per edge, or just check the reserve at the end of the route.
