"""Equivalent-circuit LiPo pack model (Rint, optionally + one RC pair) for the drone energy model.

Question it answers (prof. Rini): "if two flights use the same propulsion energy but pull it from the
battery differently (fast and high-power vs slow and low-power), does the battery nonlinearity change
which one is better?"  Short answer for this drone: only by a few percent (see ENERGY_MODEL.md).

Model (Thevenin / "Rint" family; Chen & Rincon-Mora, IEEE Trans. Energy Conversion 21(2), 2006; Tremblay &
Dessaint, World Electric Vehicle J. 3, 2009 for the generic discharge-curve idea):

    V_t = V_oc(SOC) - V_rc - I R0,        dV_rc/dt = (I R1 - V_rc)/tau        (R1 = 0 -> pure Rint)
    dSOC/dt = -I_eff / (3600 Q),          I_eff = I (I / I_1C)^(k_p - 1)       (Peukert, k_p = 1 -> ideal)
    P_elec = V_t I   (ESC input power) ->  I = (V_oc' - sqrt(V_oc'^2 - 4 R0 P)) / (2 R0),  V_oc' = V_oc - V_rc

No real root (V_oc'^2 < 4 R0 P) means the pack cannot deliver P at all: flagged, current clamped to the
maximum-power point I = V_oc'/(2 R0).  Energy drawn from the cells = int V_oc I dt (the chemical energy,
including I^2 R0 and RC dissipation; energy temporarily stored in C1 is counted as "loss" until it relaxes).

Propulsion power -> electrical power: P_elec = P_prop / eta + p_aux, with eta the combined motor+ESC
efficiency (default 0.80: ~0.85 brushless motor x ~0.95 ESC, typical hobby-grade figures; ASSUMED) and p_aux a
constant avionics/payload draw (default 0).

OCV curve: per-cell resting voltage vs SOC, linearly interpolated from the widely circulated hobby-LiPo
table (4.20 V full ... 3.27 V empty). It is an APPROXIMATION of typical LiPo data, not a single peer-reviewed
curve; real cells differ by +-30 mV with chemistry and temperature.

Internal resistance: 10 mOhm per cell default (typical 5 Ah 25C-class LiPo is ~5-15 mOhm/cell when fresh and
warm; pack datasheets, e.g. Tattu/Turnigy 5000 mAh 6S ~ 40-80 mOhm pack). Resistance rises roughly 2-4x at
0 C and with age, and rises as SOC falls below ~20% (modelled with an optional multiplier `r_scale`).

Pack choice: 6S (22.2 V nominal, 25.2 V full), 5000 mAh = 111 Wh, ~0.65 kg, i.e. ~1/3 of the ~2 kg
ZENG2019 take-off weight; hover is then ~170/0.8 = 210 W electric = 9.6 A (about 1.9 C) and endurance
~31 min to the 3.3 V/cell cutoff (~43 min at the best-endurance speed), which is typical for a 2 kg survey/inspection quad. A 4S pack would carry 14 A at hover (same
energy, 1.5x the I^2R loss for the same cell resistance), so 6S is the more favourable/plausible choice.

Peukert's law (t = H (C/(I H))^k): lead-acid k ~ 1.1-1.3; LiPo/Li-ion k ~ 1.0-1.05 (Peukert 1897;
Doerffel & Sharkh, J. Power Sources 155, 2006 for the limits of the law). So rate-dependent capacity loss is
small for LiPo at <= 3C and the dominant rate-dependent effect is the resistive loss I^2 R0 (which
is already inside the Rint model: loss/energy = R0 I / V_oc, i.e. it grows LINEARLY with current and
quadratically with power). Peukert is provided as an option, default off.
"""
from dataclasses import dataclass
import numpy as np

# resting voltage per cell vs state of charge (0..1), hobby-LiPo table (approximation, see docstring)
_OCV_SOC = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                     0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
_OCV_V = np.array([3.27, 3.61, 3.69, 3.71, 3.73, 3.75, 3.77, 3.79, 3.80, 3.82, 3.84,
                   3.85, 3.87, 3.91, 3.95, 3.98, 4.02, 4.08, 4.11, 4.15, 4.20])


def ocv_cell(soc):
    """Open-circuit voltage of one LiPo cell (V) at state of charge soc in [0, 1]."""
    return np.interp(np.clip(soc, 0.0, 1.0), _OCV_SOC, _OCV_V)


@dataclass(frozen=True)
class Pack:
    n_series: int = 6              # cells in series (6S: 22.2 V nominal)
    capacity_Ah: float = 5.0       # Q, Ah (5000 mAh)
    r_cell: float = 0.010          # R0 per cell, ohm (assumed typical, see docstring)
    r1_cell: float = 0.0           # RC-pair resistance per cell, ohm (0 -> pure Rint)
    tau: float = 30.0              # RC time constant, s (typical 10-60 s for LiPo diffusion/polarisation)
    v_cutoff_cell: float = 3.3     # terminal-voltage cutoff per cell, V
    eta: float = 0.80              # motor + ESC efficiency (assumed)
    p_aux: float = 0.0             # constant avionics/payload power, W
    peukert: float = 1.0           # Peukert exponent (1.0 = off; LiPo 1.0-1.05)
    r_scale: float = 1.0           # resistance multiplier for cold / aged packs (e.g. 2.5)

    @property
    def R0(self):
        return self.n_series * self.r_cell * self.r_scale

    @property
    def R1(self):
        return self.n_series * self.r1_cell * self.r_scale

    @property
    def energy_Wh(self):           # nominal energy, using the mean OCV over the curve
        return self.capacity_Ah * self.n_series * float(np.trapezoid(_OCV_V, _OCV_SOC))

    def voc(self, soc):
        return self.n_series * ocv_cell(soc)


PACK_6S_5000 = Pack()
PACK_6S_5000_RC = Pack(r_cell=0.008, r1_cell=0.004, tau=30.0)   # same total resistance, split R0 + R1
PACK_COLD = Pack(r_scale=3.0)                                     # ~0 C or aged pack


@dataclass
class BatteryResult:
    t: np.ndarray              # time at end of each step, s
    current: np.ndarray        # A
    v_term: np.ndarray         # V
    soc: np.ndarray            # state of charge at end of each step (initial value not included)
    p_elec: np.ndarray         # demanded electrical power, W
    e_elec: float              # energy delivered to the ESC, J
    e_cell: float              # energy drawn from the cells (chemical), J, includes losses
    e_loss: float              # e_cell - e_elec, J
    soc_used: float
    cutoff: bool               # terminal voltage fell below the cutoff at some step
    infeasible: bool           # demand exceeded max deliverable power at some step (no real root)
    t_cutoff: float = float("nan")   # time of first cutoff/infeasible event, s

    @property
    def loss_fraction(self):
        return self.e_loss / self.e_cell if self.e_cell > 0 else 0.0


def simulate(p_prop, dt, pack: Pack = PACK_6S_5000, soc0=1.0, from_electrical=False) -> BatteryResult:
    """Discharge the pack with a propulsion power series p_prop (W, one value per step of length dt s).
    If from_electrical is True the series is taken as electrical power (eta and p_aux not applied)."""
    p_prop = np.atleast_1d(np.asarray(p_prop, dtype=float))
    p_el = p_prop if from_electrical else p_prop / pack.eta + pack.p_aux
    n = p_el.size
    I = np.zeros(n); Vt = np.zeros(n); S = np.zeros(n)
    soc, vrc = float(soc0), 0.0
    R0, R1 = pack.R0, pack.R1
    a = 1.0 - np.exp(-dt / pack.tau) if R1 > 0 else 0.0
    I_1C = pack.capacity_Ah
    e_cell = e_el = 0.0
    cutoff = infeasible = False
    t_event = float("nan")
    v_cut = pack.n_series * pack.v_cutoff_cell
    for i in range(n):
        voc = float(pack.voc(soc)) - vrc
        disc = voc * voc - 4.0 * R0 * p_el[i]
        if disc < 0.0:                                   # pack cannot deliver this power
            infeasible = True
            if np.isnan(t_event):
                t_event = (i + 1) * dt
            cur = voc / (2.0 * R0)
        else:
            cur = (voc - np.sqrt(disc)) / (2.0 * R0)
        vt = voc - cur * R0
        if vt < v_cut and not cutoff:
            cutoff = True
            if np.isnan(t_event):
                t_event = (i + 1) * dt
        i_eff = cur * (cur / I_1C) ** (pack.peukert - 1.0) if pack.peukert != 1.0 and cur > 0 else cur
        e_cell += float(pack.voc(soc)) * cur * dt        # chemical power = V_oc * I
        e_el += (vt * cur if disc < 0.0 else p_el[i]) * dt   # delivered (not demanded) power
        vrc += (cur * R1 - vrc) * a
        soc = max(soc - i_eff * dt / (3600.0 * pack.capacity_Ah), 0.0)
        I[i], Vt[i], S[i] = cur, vt, soc
    return BatteryResult(np.arange(1, n + 1) * dt, I, Vt, S, p_el, e_el, e_cell, e_cell - e_el,
                         soc0 - soc, cutoff, infeasible, t_event)


def energy_constant_power(P_prop, T, pack: Pack = PACK_6S_5000, soc0=1.0, dt=2.0) -> BatteryResult:
    """Battery result for constant propulsion power P_prop (W) held for T seconds."""
    n = max(int(np.ceil(T / dt)), 1)
    return simulate(np.full(n, float(P_prop)), T / n, pack, soc0)


def max_power(pack: Pack = PACK_6S_5000, soc=1.0):
    """Maximum electrical power (W) the pack can deliver at this SOC (V_oc^2 / 4 R0), no cutoff applied."""
    return float(pack.voc(soc)) ** 2 / (4 * pack.R0)


if __name__ == "__main__":
    p = PACK_6S_5000
    print(f"6S 5000 mAh: {p.n_series * 3.7 * p.capacity_Ah:.0f} Wh nominal, R0 = {p.R0 * 1e3:.0f} mOhm, "
          f"max power {max_power(p):.0f} W")
    for P in (126, 170, 250, 400):
        r = energy_constant_power(P, 600, p)
        print(f"P_prop {P:4d} W for 10 min: electrical {r.e_elec/1e3:6.1f} kJ, cells {r.e_cell/1e3:6.1f} kJ, "
              f"loss {100*r.loss_fraction:4.1f} %, SOC used {100*r.soc_used:4.1f} %, I = {r.current[0]:.1f} A")
