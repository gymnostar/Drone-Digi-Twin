"""Rotary-wing drone power and drag as functions of AIRSPEED (the wind enters only through airspeed).

Core model: Zeng, Xu & Zhang (2019), "Energy minimization for wireless communication with
rotary-wing UAV", IEEE Trans. Wireless Commun. 18(4), eq. (12). It is the standard propulsion
model in the UAV-communications literature (the SWARM-LINK side), derived from helicopter
blade-element / momentum theory for level forward flight at airspeed V:

    P(V) = P0 (1 + 3V^2/U_tip^2)                              blade profile power
         + Pi (sqrt(1 + V^4/(4 v0^4)) - V^2/(2 v0^2))^(1/2)    induced power (lift)
         + 1/2 d0 rho s A V^3                                  parasite power = drag x V

    P0 = delta/8 rho s A Omega^3 R^3,  Pi = (1+k) W^1.5 / sqrt(2 rho A),  v0 = sqrt(W / (2 rho A))

The drag force behind the parasite term is D = 1/2 rho (d0 s A) V^2, opposite the airspeed.
Hover (V = 0) costs P0 + Pi. Induced power FALLS with speed (translational lift), parasite power
RISES as V^3, so P(V) is U-shaped: there is a speed that minimises energy per metre.

Wind: a drone moving at ground velocity v_g in wind w sees airspeed v_a = v_g - w. Headwind ->
larger |v_a| -> more power for the same ground speed, and more time per metre.

Vertical motion (extension, not in Zeng 2019): climbing at air-relative rate c = v_gz - w_z adds
W*c (potential-energy rate). Descent gets no credit (no regeneration; conservative).

All functions broadcast over numpy arrays.
"""
from dataclasses import dataclass, replace
import numpy as np

RHO = 1.225  # air density, kg/m^3 (ISA sea level; Oberpfaffenhofen ~590 m asl would be ~1.16)


@dataclass(frozen=True)
class Drone:
    """Parameters of Zeng et al. 2019, Table I (W = 20 N, i.e. a ~2 kg rotary-wing UAV)."""
    W: float = 20.0          # weight, N
    rho: float = RHO         # air density, kg/m^3
    R: float = 0.4           # rotor radius, m
    A: float = 0.503         # rotor disc area, m^2
    Omega: float = 300.0     # blade angular velocity, rad/s
    s: float = 0.05          # rotor solidity
    d0: float = 0.6          # fuselage drag ratio
    k: float = 0.1           # induced-power correction factor
    delta: float = 0.012     # profile drag coefficient

    @property
    def U_tip(self):          # blade tip speed, m/s
        return self.Omega * self.R

    @property
    def P0(self):             # blade profile power in hover, W
        return self.delta / 8 * self.rho * self.s * self.A * self.Omega**3 * self.R**3

    @property
    def Pi(self):             # induced power in hover, W
        return (1 + self.k) * self.W**1.5 / np.sqrt(2 * self.rho * self.A)

    @property
    def v0(self):             # mean rotor induced velocity in hover, m/s
        return np.sqrt(self.W / (2 * self.rho * self.A))

    @property
    def drag_area(self):      # equivalent flat-plate area d0*s*A = C_D * A_ref, m^2
        return self.d0 * self.s * self.A


ZENG2019 = Drone()


def power_level(V, d: Drone = ZENG2019):
    """Propulsion power (W) in level flight at horizontal AIRSPEED V (m/s)."""
    V = np.asarray(V, dtype=float)
    blade = d.P0 * (1 + 3 * V**2 / d.U_tip**2)
    induced = d.Pi * np.sqrt(np.sqrt(1 + V**4 / (4 * d.v0**4)) - V**2 / (2 * d.v0**2))
    parasite = 0.5 * d.d0 * d.rho * d.s * d.A * V**3
    return blade + induced + parasite


def power(v_air_h, v_air_z=0.0, d: Drone = ZENG2019):
    """Power (W) for horizontal airspeed magnitude v_air_h and vertical air-relative rate v_air_z
    (m/s, + = climbing through the air). Climb adds W*v_air_z; descent gets no credit."""
    return power_level(v_air_h, d) + d.W * np.maximum(np.asarray(v_air_z, dtype=float), 0.0)


def drag_force(v_air, d: Drone = ZENG2019):
    """Parasite drag force vector (N) for airspeed vector(s) v_air (..., 3): D = -1/2 rho CdA |v| v."""
    v_air = np.asarray(v_air, dtype=float)
    return -0.5 * d.rho * d.drag_area * np.linalg.norm(v_air, axis=-1, keepdims=True) * v_air


def segment_power(direction, v_ground, wind, d: Drone = ZENG2019):
    """Power (W) to fly with ground speed v_ground (m/s) along unit vector(s) `direction` (..., 3)
    through wind vector(s) `wind` (..., 3). v_ground broadcasts against direction[..., 0]."""
    v = np.asarray(v_ground, dtype=float)
    va = v[..., None] * np.asarray(direction, dtype=float) - np.asarray(wind, dtype=float)
    return power(np.hypot(va[..., 0], va[..., 1]), va[..., 2], d)


def energy_per_metre(direction, v_ground, wind, d: Drone = ZENG2019):
    """Energy per metre of ground track (J/m) = power / ground speed."""
    return segment_power(direction, v_ground, wind, d) / np.asarray(v_ground, dtype=float)


def power_components(V, d: Drone = ZENG2019):
    """(blade profile, induced, parasite) power in W at horizontal airspeed V; their sum is power_level."""
    V = np.asarray(V, dtype=float)
    blade = d.P0 * (1 + 3 * V**2 / d.U_tip**2)
    induced = d.Pi * np.sqrt(np.sqrt(1 + V**4 / (4 * d.v0**4)) - V**2 / (2 * d.v0**2))
    parasite = 0.5 * d.d0 * d.rho * d.s * d.A * V**3
    return blade, induced, parasite


_V_FINE = np.linspace(0.05, 40.0, 8000)


def max_endurance_speed(d: Drone = ZENG2019):
    """Airspeed (m/s) that minimises power P(V) (longest flight time)."""
    return float(_V_FINE[np.argmin(power_level(_V_FINE, d))])


def max_range_speed(d: Drone = ZENG2019):
    """Airspeed (m/s) that minimises energy per metre P(V)/V in still air (longest range)."""
    return float(_V_FINE[np.argmin(power_level(_V_FINE, d) / _V_FINE)])


SPEEDS = np.arange(3.0, 25.01, 0.5)   # candidate ground speeds, m/s


def best_energy_per_metre(direction, wind, d: Drone = ZENG2019, speeds=SPEEDS):
    """Minimum J/m over candidate ground speeds, and the speed that achieves it.
    A drone that may pick its speed flies faster into a headwind and slower with a tailwind."""
    direction = np.asarray(direction, dtype=float)[..., None, :]
    wind = np.asarray(wind, dtype=float)[..., None, :]
    e = energy_per_metre(direction, speeds, wind, d)                 # (..., n_speeds)
    i = np.argmin(e, axis=-1)
    return np.take_along_axis(e, i[..., None], axis=-1)[..., 0], speeds[i]


# ---------------------------------------------------------------------------------------------
# AirSim SimpleFlight default multirotor ("generic quad", F450-like frame)
# Source: microsoft/AirSim, AirLib/include/vehicles/multirotor/MultiRotorParams.hpp
# (setupFrameGenericQuad, Params) and RotorParams.hpp (read from GitHub main in Nov 2026).
# AirSim has NO battery/power model in SimpleFlight: it only needs thrust, so Zeng's blade/induced
# parameters below are mapped (or assumed) by us; the thrust side comes from AirSim's numbers.
#   mass 1.0 kg (setupFrameGenericQuad); 4 rotors; arm length 0.2275 m (F450); body box 0.18 x 0.11 x 0.04 m
#   propeller diameter 0.2286 m (DJI Phantom 2 prop); C_T = 0.109919, C_P = 0.040164 (UIUC, T = C_T rho n^2 D^4)
#   max_rpm 6396.667; air density 1.225; linear_drag_coefficient = 1.3/4 = 0.325 (Params default).
# Drag in AirSim's FastPhysicsEngine: per face F = -(Cd_lin * face_area) * rho * v_n^2 (no 1/2); faces are
# the six body-box sides plus the rotor discs and rotor cylinder cross-sections (1 cm high). The face list
# is from memory of MultiRotorPhysicsBody::createDragVertices, not re-verified line by line.
# ---------------------------------------------------------------------------------------------
_AS_MASS = 1.0                                   # kg, setupFrameGenericQuad
_AS_G = 9.81                                     # m/s^2, AirSim default gravity (assumed 9.81)
_AS_D = 0.2286                                   # m, RotorParams.propeller_diameter
_AS_CT = 0.109919                                # RotorParams.C_T
_AS_RPM_MAX = 6396.667                           # RotorParams.max_rpm
_AS_R = _AS_D / 2                                # m, rotor radius
_AS_A = 4 * np.pi * _AS_R**2                     # m^2, disc area of 4 rotors (single-rotor-equivalent)
_AS_N_HOVER = np.sqrt((_AS_MASS * _AS_G / 4) / (_AS_CT * RHO * _AS_D**4))   # rev/s per rotor at hover
_AS_OMEGA = 2 * np.pi * _AS_N_HOVER              # rad/s at hover (about 77% of max rpm)
_AS_S = 0.12                                     # solidity: ASSUMED (2 blades, ~22 mm mean chord, 9.4-inch prop)
# Frontal area for AirSim drag: body box front face 0.11*0.04 plus 4 rotor cylinders 0.2286*0.01 each.
_AS_FRONT = 0.11 * 0.04 + 4 * _AS_D * 0.01       # m^2 (from-memory face list, see above)
_AS_CDA = 2 * 0.325 * _AS_FRONT                  # m^2: AirSim F = 0.325*A*rho*v^2 == 1/2 rho (CdA) v^2
AIRSIM_GENERIC_QUAD = Drone(
    W=_AS_MASS * _AS_G,                          # 9.81 N, AirSim mass x gravity
    rho=RHO,                                     # AirSim air_density 1.225
    R=_AS_R,                                     # 0.1143 m, half of propeller_diameter
    A=_AS_A,                                     # 0.1645 m^2, 4 x pi R^2
    Omega=_AS_OMEGA,                             # ~513 rad/s, from C_T rho n^2 D^4 = W/4 (derived)
    s=_AS_S,                                     # assumed
    d0=_AS_CDA / (_AS_S * _AS_A),                # chosen so that d0*s*A = CdA of AirSim drag model (derived)
    k=0.1,                                       # induced-power factor: assumed, same as Zeng
    delta=0.012,                                 # profile drag coefficient: assumed, same as Zeng
)


def airsim_shaft_power_hover():
    """Aerodynamic shaft power (W) in hover implied by AirSim's own C_P (P = C_P rho n^3 D^5 per rotor).
    Independent sanity check of AIRSIM_GENERIC_QUAD's hover power (no profile-power assumptions)."""
    return 4 * 0.040164 * RHO * _AS_N_HOVER**3 * _AS_D**5


def power_momentum(v_air_vec, d: Drone = ZENG2019, mu_coeff=4.65):
    """First-principles level-flight power (W) for a quad: momentum theory with Glauert's forward-flight
    induced velocity. Horizontal air-relative velocity vector (..., 2 or 3; only the horizontal norm is used).

    Thrust balances weight and parasite drag: T = sqrt(W^2 + D^2), tilt theta = atan(D/W), D = 1/2 rho CdA V^2.
    Induced velocity: v_i = v_h^2 / sqrt((V cos th)^2 + (V sin th + v_i)^2), v_h = sqrt(T/(2 rho A)) (solved by
    bisection). Power: P = (1+k) T v_i + T V sin th (= D V, parasite) + profile power
    P0 (1 + mu_coeff mu^2), mu = V cos th / U_tip (4.65 is the helicopter value, Leishman 2006; Zeng uses 3)."""
    v = np.asarray(v_air_vec, dtype=float)
    V = np.hypot(v[..., 0], v[..., 1]) if v.shape[-1] >= 2 else np.abs(v[..., 0])
    D = 0.5 * d.rho * d.drag_area * V**2
    T = np.hypot(d.W, D)
    th = np.arctan2(D, d.W)
    vh = np.sqrt(T / (2 * d.rho * d.A))
    lo, hi = np.zeros_like(V), 2 * vh * np.ones_like(V)         # v_i in (0, v_h] always
    for _ in range(60):
        vi = 0.5 * (lo + hi)
        f = vi * np.sqrt((V * np.cos(th))**2 + (V * np.sin(th) + vi)**2) - vh**2   # increasing in vi
        lo, hi = np.where(f < 0, vi, lo), np.where(f < 0, hi, vi)
    vi = 0.5 * (lo + hi)
    induced = (1 + d.k) * T * vi
    parasite = T * V * np.sin(th)
    mu = V * np.cos(th) / d.U_tip
    profile = d.P0 * (1 + mu_coeff * mu**2)
    return induced + parasite + profile


if __name__ == "__main__":
    d = ZENG2019
    print(f"P0 = {d.P0:.2f} W, Pi = {d.Pi:.2f} W, v0 = {d.v0:.2f} m/s, hover = {power_level(0.0):.1f} W")
    for V in (0, 5, 10, 15, 20, 25):
        print(f"V = {V:4.1f} m/s  P = {power_level(V):6.1f} W  E/m = {power_level(V)/max(V, 1e-9):8.2f} J/m")
    e, v = best_energy_per_metre(np.array([1.0, 0, 0]), np.zeros(3))
    print(f"still air: best {e.item():.2f} J/m at {v.item()} m/s")
    for w in (-5, 5):
        e, v = best_energy_per_metre(np.array([1.0, 0, 0]), np.array([w, 0, 0.0]))
        print(f"{'tail' if w > 0 else 'head'}wind 5 m/s: best {e.item():.2f} J/m at ground speed {v.item()} m/s")
    print(f"max-endurance speed {max_endurance_speed():.1f} m/s, max-range speed {max_range_speed():.1f} m/s")
    a = AIRSIM_GENERIC_QUAD
    print(f"AirSim generic quad: W={a.W:.2f} N R={a.R:.4f} m A={a.A:.4f} m2 Omega={a.Omega:.0f} rad/s "
          f"(U_tip {a.U_tip:.1f} m/s) s={a.s} d0={a.d0:.3f} CdA={a.drag_area:.4f} m2")
    print(f"  P0={a.P0:.1f} W Pi={a.Pi:.1f} W hover={power_level(0.0, a):.1f} W "
          f"(AirSim C_P implied shaft power {airsim_shaft_power_hover():.1f} W)")
    print(f"  max-endurance {max_endurance_speed(a):.1f} m/s, max-range {max_range_speed(a):.1f} m/s, "
          f"min J/m {np.min(power_level(_V_FINE, a) / _V_FINE):.2f}")
