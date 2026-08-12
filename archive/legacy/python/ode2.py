"""
AAV PBPK + liver cellular fate model, smoother plotting version.

Main changes relative to the first draft:
1. Optional zero-order infusion instead of instantaneous bolus.
2. Non-uniform dense time grid near t = 0, so the fast early phase is not missed.
3. Log plotting uses NaN for zero/negative values instead of replacing zeros by 1.
4. Effective flow scaling is explicit. In a one-central-blood toy model, raw mouse CO
   can make vascular compartments equilibrate in seconds.
5. Figures are saved to ./aav_pbpk_outputs by default.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


ORGANS = ["liver", "spleen", "kidney", "heart", "muscle", "lung", "rest"]

IDX = {
    "A_blood": 0,
    "A_liver_v": 1,
    "A_liver_isf": 2,
    "A_spleen_v": 3,
    "A_spleen_isf": 4,
    "A_kidney_v": 5,
    "A_kidney_isf": 6,
    "A_heart_v": 7,
    "A_heart_isf": 8,
    "A_muscle_v": 9,
    "A_muscle_isf": 10,
    "A_lung_v": 11,
    "A_lung_isf": 12,
    "A_rest_v": 13,
    "A_rest_isf": 14,
    "B": 15,
    "EE": 16,
    "LE": 17,
    "CY": 18,
    "Ncap": 19,
    "Nss": 20,
    "Nds": 21,
    "Epi": 22,
    "M": 23,
    "P": 24,
    "Ab": 25,
}

COLORS = {
    "liver": "red",
    "spleen": "green",
    "kidney": "blue",
    "heart": "orange",
    "muscle": "purple",
    "lung": "cyan",
    "rest": "brown",
}


# ---------------------------------------------------------------------
# User-facing switches
# ---------------------------------------------------------------------
DOSE_VG = 1e12

# "infusion" gives a visually smoother early curve.
# Set ADMINISTRATION = "bolus" to recover the original instantaneous dose style.
ADMINISTRATION = "infusion"  # allowed: "infusion" or "bolus"
INFUSION_DURATION_MIN = 5.0

# In this one-central-blood toy model, using raw CO makes all vascular
# concentrations equilibrate within seconds. Q_SCALE converts CO into an
# effective central-organ exchange rate. For the original raw-flow behavior,
# set Q_SCALE = 1.0. For smoother iGEM-demo curves, try 0.01 to 0.10.
Q_SCALE = 0.05

SAVE_FIGURES = True
SHOW_FIGURES = False
OUTPUT_DIR = Path("aav_pbpk_outputs")


# ---------------------------------------------------------------------
# Model parameters
# ---------------------------------------------------------------------
def make_params() -> Dict[str, float | str]:
    co = 500.0  # mL/h, approximate cardiac output for a 25 g mouse
    q_scale = Q_SCALE

    return {
        # Administration
        "dose_vg": DOSE_VG,
        "administration": ADMINISTRATION,
        "T_inf_h": INFUSION_DURATION_MIN / 60.0,

        # Central blood
        "V_blood": 1.5,      # mL
        "CL_blood": 0.01,    # mL/h-like effective clearance

        # Organ vascular and interstitial effective volumes, mL
        "V_liver_v": 0.14,
        "V_liver_isf": 0.22,
        "V_spleen_v": 0.01,
        "V_spleen_isf": 0.02,
        "V_kidney_v": 0.05,
        "V_kidney_isf": 0.08,
        "V_heart_v": 0.012,
        "V_heart_isf": 0.025,
        "V_muscle_v": 0.4,
        "V_muscle_isf": 1.2,
        "V_lung_v": 0.08,
        "V_lung_isf": 0.04,
        "V_rest_v": 0.7,
        "V_rest_isf": 2.16,

        # Effective organ blood-flow-like exchange rates, mL/h
        "CO": co,
        "Q_scale": q_scale,
        "Q_lung": q_scale * co,
        "Q_liver": q_scale * 0.25 * co,
        "Q_spleen": q_scale * 0.06 * co,
        "Q_kidney": q_scale * 0.20 * co,
        "Q_heart": q_scale * 0.05 * co,
        "Q_muscle": q_scale * 0.15 * co,
        "Q_rest": q_scale * 0.14 * co,

        # Vascular-to-interstitial permeability terms, mL/h-like
        "PS_liver": 0.20,
        "PS_spleen": 0.15,
        "PS_kidney": 0.08,
        "PS_heart": 0.05,
        "PS_muscle": 0.03,
        "PS_lung": 0.10,
        "PS_rest": 0.02,

        # Tissue partition coefficients
        "Kp_liver": 1.5,
        "Kp_spleen": 1.2,
        "Kp_kidney": 0.8,
        "Kp_heart": 0.6,
        "Kp_muscle": 0.5,
        "Kp_lung": 0.7,
        "Kp_rest": 0.4,

        # Organ vascular RES / nonspecific clearance, 1/h
        "k_res_liver": 0.02,
        "k_res_spleen": 0.03,
        "k_res_kidney": 0.005,
        "k_res_heart": 0.002,
        "k_res_muscle": 0.001,
        "k_res_lung": 0.004,
        "k_res_rest": 0.001,

        # Interstitial degradation/loss, 1/h
        "k_deg_isf_liver": 0.001,
        "k_deg_isf_spleen": 0.001,
        "k_deg_isf_kidney": 0.001,
        "k_deg_isf_heart": 0.001,
        "k_deg_isf_muscle": 0.001,
        "k_deg_isf_lung": 0.001,
        "k_deg_isf_rest": 0.001,

        # Liver cell-surface binding
        "R_tot": 1e5,
        "k_on": 1e-6,
        "k_off": 0.05,
        "k_int": 0.2,

        # Liver intracellular trafficking, 1/h
        "k_ee_le": 0.3,
        "k_rec": 0.05,
        "k_deg_ee": 0.02,
        "k_escape": 0.005,
        "k_lys": 0.1,
        "k_nuc": 0.02,
        "k_uncoat_cyto": 0.005,
        "k_uncoat_nuc": 0.02,
        "k_deg_cyto": 0.01,
        "k_deg_ncap": 0.005,
        "k_ds": 0.01,
        "k_deg_ss": 0.02,
        "k_epi": 0.01,
        "k_deg_ds": 0.005,
        "k_loss_epi": 0.0005,
        "k_dil": 0.0,

        # Expression module
        "k_tx": 2.0,
        "h": 1.2,
        "EC50_tx": 100.0,
        "k_deg_m": 0.2,
        "k_tl": 5.0,
        "k_deg_p": 0.02,

        # Antibody module
        "k_neut": 1e-14,
        "k_ab_max": 0.05,
        "EC50_ab": 1e10,
        "k_deg_ab": 0.005,
    }


def make_initial_condition(p: Dict[str, float | str]) -> np.ndarray:
    y0 = np.zeros(26, dtype=float)

    if p["administration"] == "bolus":
        y0[IDX["A_blood"]] = float(p["dose_vg"])
    elif p["administration"] == "infusion":
        y0[IDX["A_blood"]] = 0.0
    else:
        raise ValueError("ADMINISTRATION must be 'bolus' or 'infusion'.")

    return y0


# ---------------------------------------------------------------------
# ODE right-hand side
# ---------------------------------------------------------------------
def dose_input_rate(t: float, p: Dict[str, float | str]) -> float:
    """Zero-order infusion into central blood, vg/h."""
    if p["administration"] != "infusion":
        return 0.0

    T_inf = float(p["T_inf_h"])
    if 0.0 <= t <= T_inf:
        return float(p["dose_vg"]) / T_inf
    return 0.0


def organ_fluxes(A_blood: float, A_v: float, A_isf: float, organ: str, p: Dict[str, float | str]) -> Tuple[float, float, float, float]:
    Cb = A_blood / float(p["V_blood"])
    Cv = A_v / float(p[f"V_{organ}_v"])
    Cisf = A_isf / float(p[f"V_{organ}_isf"])

    # Blood <-> organ vascular exchange
    J_blood_to_v = float(p[f"Q_{organ}"]) * (Cb - Cv)

    # Organ vascular <-> organ interstitial exchange
    J_v_to_isf = float(p[f"PS_{organ}"]) * (Cv - Cisf / float(p[f"Kp_{organ}"]))

    # Organ vascular loss and interstitial loss
    J_res = float(p[f"k_res_{organ}"]) * A_v
    J_deg_isf = float(p[f"k_deg_isf_{organ}"]) * A_isf

    return J_blood_to_v, J_v_to_isf, J_res, J_deg_isf


def rhs(t: float, y: np.ndarray, p: Dict[str, float | str]) -> list[float]:
    (
        A_blood,
        A_liver_v, A_liver_isf,
        A_spleen_v, A_spleen_isf,
        A_kidney_v, A_kidney_isf,
        A_heart_v, A_heart_isf,
        A_muscle_v, A_muscle_isf,
        A_lung_v, A_lung_isf,
        A_rest_v, A_rest_isf,
        B, EE, LE, CY, Ncap, Nss, Nds, Epi, M, P, Ab,
    ) = y

    A_v = {
        "liver": A_liver_v,
        "spleen": A_spleen_v,
        "kidney": A_kidney_v,
        "heart": A_heart_v,
        "muscle": A_muscle_v,
        "lung": A_lung_v,
        "rest": A_rest_v,
    }
    A_isf = {
        "liver": A_liver_isf,
        "spleen": A_spleen_isf,
        "kidney": A_kidney_isf,
        "heart": A_heart_isf,
        "muscle": A_muscle_isf,
        "lung": A_lung_isf,
        "rest": A_rest_isf,
    }

    J_blood_to_v = {}
    J_v_to_isf = {}
    J_res = {}
    J_deg_isf = {}

    for organ in ORGANS:
        J_blood_to_v[organ], J_v_to_isf[organ], J_res[organ], J_deg_isf[organ] = organ_fluxes(
            A_blood, A_v[organ], A_isf[organ], organ, p
        )

    # Liver cell-surface binding. Clamp only for rate evaluation to avoid
    # non-integer powers or negative receptor artifacts from tiny numerical noise.
    C_liver_isf = max(A_liver_isf, 0.0) / float(p["V_liver_isf"])
    B_eff = max(B, 0.0)
    R_free = max(float(p["R_tot"]) - B_eff, 0.0)
    J_bind = float(p["k_on"]) * C_liver_isf * R_free - float(p["k_off"]) * B_eff
    if B <= 0.0 and J_bind < 0.0:
        J_bind = 0.0

    Ab_eff = max(Ab, 0.0)
    J_neut_blood = float(p["k_neut"]) * Ab_eff * max(A_blood, 0.0)

    dA_blood = (
        dose_input_rate(t, p)
        - sum(J_blood_to_v.values())
        - float(p["CL_blood"]) * (A_blood / float(p["V_blood"]))
        - J_neut_blood
    )

    dA_v = {}
    dA_isf = {}
    for organ in ORGANS:
        dA_v[organ] = J_blood_to_v[organ] - J_v_to_isf[organ] - J_res[organ]
        dA_isf[organ] = J_v_to_isf[organ] - J_deg_isf[organ]

    dA_isf["liver"] -= J_bind

    dB = J_bind - float(p["k_int"]) * B
    dEE = float(p["k_int"]) * B - (float(p["k_ee_le"]) + float(p["k_rec"]) + float(p["k_deg_ee"])) * EE
    dLE = float(p["k_ee_le"]) * EE - (float(p["k_escape"]) + float(p["k_lys"])) * LE
    dCY = float(p["k_escape"]) * LE - (float(p["k_nuc"]) + float(p["k_uncoat_cyto"]) + float(p["k_deg_cyto"])) * CY
    dNcap = float(p["k_nuc"]) * CY - (float(p["k_uncoat_nuc"]) + float(p["k_deg_ncap"])) * Ncap
    dNss = float(p["k_uncoat_cyto"]) * CY + float(p["k_uncoat_nuc"]) * Ncap - (float(p["k_ds"]) + float(p["k_deg_ss"])) * Nss
    dNds = float(p["k_ds"]) * Nss - (float(p["k_epi"]) + float(p["k_deg_ds"])) * Nds
    dEpi = float(p["k_epi"]) * Nds - (float(p["k_loss_epi"]) + float(p["k_dil"])) * Epi

    Epi_eff = max(Epi, 0.0)
    h = float(p["h"])
    tx = float(p["k_tx"]) * (Epi_eff ** h) / (float(p["EC50_tx"]) ** h + Epi_eff ** h + 1e-30)
    dM = tx - float(p["k_deg_m"]) * M
    dP = float(p["k_tl"]) * M - float(p["k_deg_p"]) * P

    Ag = max(A_blood + sum(A_v.values()) + 0.5 * sum(A_isf.values()) + LE, 0.0)
    dAb = float(p["k_ab_max"]) * Ag / (float(p["EC50_ab"]) + Ag + 1e-30) - float(p["k_deg_ab"]) * Ab

    return [
        dA_blood,
        dA_v["liver"], dA_isf["liver"],
        dA_v["spleen"], dA_isf["spleen"],
        dA_v["kidney"], dA_isf["kidney"],
        dA_v["heart"], dA_isf["heart"],
        dA_v["muscle"], dA_isf["muscle"],
        dA_v["lung"], dA_isf["lung"],
        dA_v["rest"], dA_isf["rest"],
        dB, dEE, dLE, dCY, dNcap, dNss, dNds, dEpi, dM, dP, dAb,
    ]


# ---------------------------------------------------------------------
# Solving utilities
# ---------------------------------------------------------------------
def make_short_grid() -> np.ndarray:
    """0 to 2 h, dense early sampling to show smooth early kinetics."""
    return np.unique(np.r_[
        np.linspace(0.0, 2.0 / 60.0, 450, endpoint=False),  # first 2 min
        np.linspace(2.0 / 60.0, 30.0 / 60.0, 450, endpoint=False),  # 2 to 30 min
        np.linspace(30.0 / 60.0, 2.0, 500),  # 30 min to 2 h
    ])


def make_long_grid() -> np.ndarray:
    """0 to 56 days; also dense early enough to show administration phase."""
    return np.unique(np.r_[
        np.linspace(0.0, 2.0, 300, endpoint=False),       # first 2 h
        np.linspace(2.0, 24.0 * 7.0, 300, endpoint=False), # day 0 to 7
        np.linspace(24.0 * 7.0, 24.0 * 56.0, 500),         # day 7 to 56
    ])


def solve_single_interval(t_eval: np.ndarray, y0: np.ndarray, p: Dict[str, float | str], max_step: float) -> solve_ivp:
    sol = solve_ivp(
        lambda t, y: rhs(t, y, p),
        t_span=(float(t_eval[0]), float(t_eval[-1])),
        y0=y0,
        t_eval=t_eval,
        method="Radau",
        rtol=1e-6,
        atol=1e-9,
        max_step=max_step,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


class SimpleSolution:
    def __init__(self, t: np.ndarray, y: np.ndarray):
        self.t = t
        self.y = y
        self.success = True


def solve_model(t_eval: np.ndarray, y0: np.ndarray, p: Dict[str, float | str], post_infusion_max_step: float) -> SimpleSolution:
    """Solve with a separate small-step infusion segment if needed."""
    t_eval = np.asarray(t_eval, dtype=float)
    t_eval = np.unique(t_eval)

    if p["administration"] != "infusion":
        sol = solve_single_interval(t_eval, y0, p, max_step=post_infusion_max_step)
        return SimpleSolution(sol.t, sol.y)

    T_inf = float(p["T_inf_h"])
    if t_eval[-1] <= T_inf:
        sol = solve_single_interval(t_eval, y0, p, max_step=T_inf / 50.0)
        return SimpleSolution(sol.t, sol.y)

    # Segment 1: 0 to infusion end, use small steps.
    t1 = np.unique(np.r_[t_eval[t_eval <= T_inf], T_inf])
    sol1 = solve_single_interval(t1, y0, p, max_step=T_inf / 50.0)
    y_T = sol1.y[:, -1]

    # Segment 2: after infusion. Input is zero, larger steps are fine.
    t2 = t_eval[t_eval > T_inf]
    sol2 = solve_single_interval(np.r_[T_inf, t2], y_T, p, max_step=post_infusion_max_step)

    # Avoid duplicating T_inf.
    t = np.r_[sol1.t, sol2.t[1:]]
    y = np.c_[sol1.y, sol2.y[:, 1:]]
    return SimpleSolution(t, y)


# ---------------------------------------------------------------------
# Plotting utilities
# ---------------------------------------------------------------------
def log_safe(x: Iterable[float]) -> np.ndarray:
    """For log-axis plotting: do not draw zeros/negative values."""
    x = np.asarray(x, dtype=float)
    return np.where(x > 0.0, x, np.nan)


def save_or_show(filename: str) -> None:
    if SAVE_FIGURES:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        plt.savefig(OUTPUT_DIR / filename, dpi=300, bbox_inches="tight")
    if SHOW_FIGURES:
        plt.show()
    else:
        plt.close()


def concentration(sol: SimpleSolution, state: str, volume_ml: float) -> np.ndarray:
    return sol.y[IDX[state]] / volume_ml


def plot_short_distribution(sol_short: SimpleSolution, p: Dict[str, float | str]) -> None:
    # Linear-axis early view, useful for seeing smooth kinetics.
    t_min = sol_short.t * 60.0
    plt.figure(figsize=(14, 6))

    plt.subplot(1, 2, 1)
    for organ in ORGANS:
        C_v = concentration(sol_short, f"A_{organ}_v", float(p[f"V_{organ}_v"]))
        plt.plot(t_min, C_v, label=f"{organ}_vascular", color=COLORS[organ], linewidth=1.8)
    C_blood = concentration(sol_short, "A_blood", float(p["V_blood"]))
    plt.plot(t_min, C_blood, label="blood", color="black", linestyle="--", linewidth=2.0)
    plt.xlim(0, 30)
    plt.xlabel("Time (min)")
    plt.ylabel("AAV concentration (vg/mL)")
    plt.title("Vascular AAV concentrations, first 30 min")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.subplot(1, 2, 2)
    for organ in ORGANS:
        C_isf = concentration(sol_short, f"A_{organ}_isf", float(p[f"V_{organ}_isf"]))
        plt.plot(t_min, C_isf, label=f"{organ}_ISF", color=COLORS[organ], linewidth=1.8)
    plt.xlim(0, 30)
    plt.xlabel("Time (min)")
    plt.ylabel("AAV concentration (vg/mL)")
    plt.title("Interstitial AAV concentrations, first 30 min")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.tight_layout()
    save_or_show("01_short_distribution_linear_30min.png")

    # Log-axis 0 to 2 h view. Zero values become NaN, not floor = 1.
    plt.figure(figsize=(14, 6))

    plt.subplot(1, 2, 1)
    for organ in ORGANS:
        C_v = concentration(sol_short, f"A_{organ}_v", float(p[f"V_{organ}_v"]))
        plt.plot(sol_short.t, log_safe(C_v), label=f"{organ}_vascular", color=COLORS[organ], linewidth=1.8)
    C_blood = concentration(sol_short, "A_blood", float(p["V_blood"]))
    plt.plot(sol_short.t, log_safe(C_blood), label="blood", color="black", linestyle="--", linewidth=2.0)
    plt.yscale("log")
    plt.xlabel("Time (h)")
    plt.ylabel("AAV concentration (vg/mL)")
    plt.title("Vascular AAV concentrations, 0-2 h")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.subplot(1, 2, 2)
    for organ in ORGANS:
        C_isf = concentration(sol_short, f"A_{organ}_isf", float(p[f"V_{organ}_isf"]))
        plt.plot(sol_short.t, log_safe(C_isf), label=f"{organ}_ISF", color=COLORS[organ], linewidth=1.8)
    plt.yscale("log")
    plt.xlabel("Time (h)")
    plt.ylabel("AAV concentration (vg/mL)")
    plt.title("Interstitial AAV concentrations, 0-2 h")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.tight_layout()
    save_or_show("02_short_distribution_log_2h.png")


def plot_long_states(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    t_day = sol_long.t / 24.0

    plt.figure(figsize=(13, 6))
    intracellular_states = [
        ("B", "Surface-bound AAV"),
        ("EE", "Early endosome"),
        ("LE", "Late endosome"),
        ("CY", "Cytosolic capsid"),
        ("Ncap", "Nuclear capsid"),
        ("Nss", "Nuclear ssDNA"),
        ("Nds", "Nuclear dsDNA"),
        ("Epi", "Episome"),
    ]
    for state, label in intracellular_states:
        plt.plot(t_day, log_safe(sol_long.y[IDX[state]]), label=label, linewidth=1.8)
    plt.yscale("log")
    plt.xlabel("Time (day)")
    plt.ylabel("Amount / arbitrary units")
    plt.title("Liver intracellular AAV trafficking states, 0-56 days")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)
    plt.tight_layout()
    save_or_show("03_liver_intracellular_56d.png")

    plt.figure(figsize=(12, 6))
    plt.plot(t_day, log_safe(sol_long.y[IDX["Epi"]]), label="Episome (Epi)", color="red", linewidth=2.0)
    plt.plot(t_day, log_safe(sol_long.y[IDX["M"]]), label="mRNA (M)", color="blue", linewidth=2.0)
    plt.plot(t_day, log_safe(sol_long.y[IDX["P"]]), label="Protein (P)", color="green", linewidth=2.0)
    plt.yscale("log")
    plt.xlabel("Time (day)")
    plt.ylabel("Amount / arbitrary units")
    plt.title("Liver expression states, 0-56 days")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.legend()
    plt.tight_layout()
    save_or_show("04_liver_expression_56d.png")

    plt.figure(figsize=(8, 4))
    plt.plot(t_day, log_safe(sol_long.y[IDX["Ab"]]), color="magenta", linewidth=2.0)
    plt.yscale("log")
    plt.xlabel("Time (day)")
    plt.ylabel("Antibody level (a.u.)")
    plt.title("Simplified antibody kinetics, 0-56 days")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.tight_layout()
    save_or_show("05_antibody_56d.png")


def auc_trapz(y: np.ndarray, x: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def print_metrics(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    C_blood_long = concentration(sol_long, "A_blood", float(p["V_blood"]))
    print("----- Simulation settings -----")
    print(f"Administration: {p['administration']}")
    print(f"Dose: {p['dose_vg']:.3e} vg")
    print(f"Infusion duration: {float(p['T_inf_h']) * 60:.2f} min")
    print(f"Q_scale: {p['Q_scale']}")

    print("\n----- Summary metrics -----")
    print("AUC_blood:", auc_trapz(C_blood_long, sol_long.t))
    print("Cmax_blood:", np.nanmax(C_blood_long))
    for organ in ORGANS:
        peak_isf_amount = np.nanmax(sol_long.y[IDX[f"A_{organ}_isf"]])
        peak_isf_conc = np.nanmax(concentration(sol_long, f"A_{organ}_isf", float(p[f"V_{organ}_isf"])))
        print(f"Peak {organ} ISF amount:", peak_isf_amount)
        print(f"Peak {organ} ISF concentration:", peak_isf_conc)
    print("Peak liver episome:", np.nanmax(sol_long.y[IDX["Epi"]]))
    print("Peak mRNA:", np.nanmax(sol_long.y[IDX["M"]]))
    print("Peak protein:", np.nanmax(sol_long.y[IDX["P"]]))
    print("Peak antibody:", np.nanmax(sol_long.y[IDX["Ab"]]))


# ---------------------------------------------------------------------
# Main script
# ---------------------------------------------------------------------
def main() -> None:
    p = make_params()
    y0 = make_initial_condition(p)

    t_eval_short = make_short_grid()
    t_eval_long = make_long_grid()

    sol_short = solve_model(t_eval_short, y0, p, post_infusion_max_step=0.002)
    sol_long = solve_model(t_eval_long, y0, p, post_infusion_max_step=1.0)

    plot_short_distribution(sol_short, p)
    plot_long_states(sol_long, p)
    print_metrics(sol_long, p)

    if SAVE_FIGURES:
        print(f"\nFigures saved to: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
