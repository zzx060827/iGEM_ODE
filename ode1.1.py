"""
AAV PBPK + liver cellular fate + multilevel kidney proximal-tubule uptake model.

This version keeps the original liver cellular fate module and adds a separate
kidney module inspired by renal proximal-tubule biology:
blood/kidney vascular -> glomerular filtrate -> proximal tubule lumen ->
apical receptor binding/endocytosis, plus a basolateral ISF uptake route ->
endosome trafficking -> cytosol/nucleus -> episome -> mRNA/protein.

Important modeling note:
All kidney cellular parameters are phenomenological demonstration parameters.
They are intended for iGEM dry-lab model structure and sensitivity analysis,
not fitted quantitative AAV kidney PK constants.

Refined version notes:
- Default administration is a 5 min infusion so early blood exposure peaks near infusion end.
- Output folder can be changed from the command line with --output-dir.
- Early PBPK plots use corrected time windows and true log axes.
- Kidney mRNA and protein are plotted separately because their magnitudes/time scales differ.
- Scenario outputs explicitly compare filtration-limited vs basolateral-tropism renal uptake hypotheses.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple
import argparse
import csv

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


ORGANS = ["liver", "spleen", "kidney", "heart", "muscle", "lung", "rest"]

STATE_NAMES = [
    # systemic extracellular AAV
    "A_blood",
    "A_liver_v", "A_liver_isf",
    "A_spleen_v", "A_spleen_isf",
    "A_kidney_v", "A_kidney_isf",
    "A_heart_v", "A_heart_isf",
    "A_muscle_v", "A_muscle_isf",
    "A_lung_v", "A_lung_isf",
    "A_rest_v", "A_rest_isf",

    # original liver cellular fate module
    "B", "EE", "LE", "CY", "Ncap", "Nss", "Nds", "Epi", "M", "P", "Ab",

    # new multilevel kidney / proximal-tubule module
    "K_filtrate",       # filtered AAV in Bowman's space / early filtrate
    "K_pt_lumen",       # AAV in proximal tubule lumen
    "K_bound_apical",   # apical brush-border receptor-bound AAV
    "K_bound_bsl",      # basolateral/interstitial receptor-bound AAV
    "K_EE",             # proximal-tubule early endosome
    "K_REC",            # recycling endosome / recycling tubules
    "K_LE",             # late endosome / large apical vacuole-like pool
    "K_LYS",            # lysosomal AAV load
    "K_CY",             # escaped cytosolic capsids
    "K_Ncap",           # nuclear capsids
    "K_Nss",            # nuclear single-stranded vector genome
    "K_Nds",            # double-stranded vector genome
    "K_Epi",            # episomal expression-competent genome
    "K_M",              # kidney transgene mRNA
    "K_P",              # kidney transgene protein
    "K_Urine",          # cumulative urinary loss
    "K_Deg",            # cumulative intracellular degradation/loss

    # cumulative bookkeeping states for mass-balance / mechanism audit
    "Dose_in",                  # cumulative administered vector
    "Loss_blood_clear",         # central blood nonspecific clearance
    "Loss_vascular_res_clear",  # organ vascular RES/endothelial clearance
    "Loss_isf_clear",           # interstitial degradation / lymph-like loss
    "Loss_neutralized",         # antibody-neutralized vector
    "Loss_liver_cell",          # liver intracellular degradation / recycling loss
]

IDX = {name: i for i, name in enumerate(STATE_NAMES)}

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

# In this one-central-blood toy model, raw cardiac output makes vascular
# concentrations equilibrate almost immediately. We therefore keep the raw
# physiological flow fractions in the parameter table, but use Q_SCALE only as
# an effective central-organ exchange throttle. Organ permeability barriers are
# still represented explicitly by PS_* terms.
Q_SCALE = 0.05

# ------------------------------------------------------------------
# Apparent extracellular AAV decay controls
# ------------------------------------------------------------------
# These are the knobs that make the in-vivo AAV exposure become bell-shaped.
# Smaller half-life -> faster decline after the peak.
ENABLE_APPARENT_AAV_DECAY = True
BLOOD_AAV_HALF_LIFE_H = 6.0        # central blood apparent half-life
VASCULAR_AAV_HALF_LIFE_H = 9.0     # nonspecific loss from organ vascular space
ISF_AAV_HALF_LIFE_H = 15.0         # degradation / lymphatic loss from tissue ISF
LIVER_EXTRA_ISF_HALF_LIFE_H = 12.0  # liver uptake / RES-like loss from liver ISF
SPLEEN_EXTRA_ISF_HALF_LIFE_H = 10.0 # spleen uptake / RES-like loss from spleen ISF
PLOT_DECAY_WINDOW_H = 48.0

SAVE_FIGURES = True
SHOW_FIGURES = False
OUTPUT_DIR_NAME = "aav_pbpk_qsp_spatial_outputs"
OUTPUT_DIR = Path(OUTPUT_DIR_NAME)

# New in the refined version:
# - "mechanistic" keeps the bell-shaped exposure but assigns clearance to
#   interpretable RES/endothelial/ISF mechanisms.
# - "half_life_demo" recovers the older apparent half-life behavior.
CLEARANCE_MODE = "mechanistic"  # allowed: "mechanistic" or "half_life_demo"

# Default scenario panel. These scenarios are intentionally simple and can be
# calibrated later against qPCR/ddPCR biodistribution data.
RUN_DESIGN_SCENARIOS = True
RUN_SPATIAL_PK_DEMO = True


CAPSID_PRESETS = {
    "baseline_AAV": {},
    "liver_detargeted": {
        "k_res_liver": 0.012,
        "k_extra_isf_clear_liver": 0.020,
        "PS_liver": 0.14,
        "R_tot": 6e4,
    },
    "kidney_tropic": {
        "PS_kidney": 0.13,
        "Kp_kidney": 1.05,
        "Bmax_pt_bsl": 4e7,
        "k_pt_bsl_on": 1.2e-11,
        "k_kidney_escape": 0.010,
        "k_res_liver": 0.014,
    },
    "endosomal_escape_enhanced": {
        "k_escape": 0.010,
        "k_kidney_escape": 0.012,
        "k_lys": 0.075,
        "k_kidney_lys": 0.075,
    },
}

PROMOTER_PRESETS = {
    "ubiquitous": {},
    "liver_biased": {
        "k_tx": 2.8,
        "k_kidney_tx": 0.8,
    },
    "kidney_biased": {
        "k_tx": 0.9,
        "k_kidney_tx": 2.6,
        "EC50_kidney_tx": 70.0,
    },
}


# ---------------------------------------------------------------------
# Model parameters
# ---------------------------------------------------------------------
def make_params() -> Dict[str, float | str]:
    co = 25.0  # mL/h, approximate cardiac output for a 25 g mouse
    q_scale = Q_SCALE

    if CLEARANCE_MODE == "mechanistic":
        # Interpretable lumped mechanisms. Values remain demonstration-scale
        # until calibrated, but now each knob maps to a biological process.
        k_clear_blood = 0.015 if ENABLE_APPARENT_AAV_DECAY else 0.0
        k_clear_vascular = 0.010 if ENABLE_APPARENT_AAV_DECAY else 0.0
        k_clear_isf = 0.006 if ENABLE_APPARENT_AAV_DECAY else 0.0
        k_extra_liver = 0.035 if ENABLE_APPARENT_AAV_DECAY else 0.0   # Kupffer/RES-like
        k_extra_spleen = 0.045 if ENABLE_APPARENT_AAV_DECAY else 0.0  # splenic macrophage-like
    elif CLEARANCE_MODE == "half_life_demo":
        k_clear_blood = (np.log(2.0) / BLOOD_AAV_HALF_LIFE_H) if ENABLE_APPARENT_AAV_DECAY else 0.0
        k_clear_vascular = (np.log(2.0) / VASCULAR_AAV_HALF_LIFE_H) if ENABLE_APPARENT_AAV_DECAY else 0.0
        k_clear_isf = (np.log(2.0) / ISF_AAV_HALF_LIFE_H) if ENABLE_APPARENT_AAV_DECAY else 0.0
        k_extra_liver = (np.log(2.0) / LIVER_EXTRA_ISF_HALF_LIFE_H) if ENABLE_APPARENT_AAV_DECAY else 0.0
        k_extra_spleen = (np.log(2.0) / SPLEEN_EXTRA_ISF_HALF_LIFE_H) if ENABLE_APPARENT_AAV_DECAY else 0.0
    else:
        raise ValueError("CLEARANCE_MODE must be 'mechanistic' or 'half_life_demo'.")

    return {
        # Administration
        "dose_vg": DOSE_VG,
        "administration": ADMINISTRATION,
        "T_inf_h": INFUSION_DURATION_MIN / 60.0,
        "clearance_mode": CLEARANCE_MODE,

        # Central blood
        "V_blood": 1.5,      # mL
        # Kept for compatibility with the original script. The stronger decay
        # in this version mainly comes from k_clear_blood below.
        "CL_blood": 0.1,    # mL/h-like effective clearance

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
        # Raw organ flow allocation retained for transparent physiology/QSP reporting.
        "Q_raw_lung": co,
        "Q_raw_liver": 0.25 * co,
        "Q_raw_spleen": 0.06 * co,
        "Q_raw_kidney": 0.20 * co,
        "Q_raw_heart": 0.05 * co,
        "Q_raw_muscle": 0.15 * co,
        "Q_raw_rest": 0.14 * co,
        # Effective exchange flow used by this 0D central-blood demo.
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

        # Extracellular AAV clearance, 1/h. In mechanistic mode these represent
        # central nonspecific loss, endothelial uptake, tissue degradation /
        # lymph-like loss, and liver/spleen RES enrichment.
        "k_clear_blood": k_clear_blood,
        "k_clear_vascular": k_clear_vascular,
        "k_clear_isf": k_clear_isf,
        "k_extra_isf_clear_liver": k_extra_liver,
        "k_extra_isf_clear_spleen": k_extra_spleen,
        "k_extra_isf_clear_kidney": 0.0,
        "k_extra_isf_clear_heart": 0.0,
        "k_extra_isf_clear_muscle": 0.0,
        "k_extra_isf_clear_lung": 0.0,
        "k_extra_isf_clear_rest": 0.0,

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
        "k_loss_epi": 0.02,
        "k_dil": 0.0,

        # Expression module
        "k_tx": 2.0,
        "h": 1.2,
        "EC50_tx": 100.0,
        "k_tl": 5.0,
        "k_deg_m": np.log(2) / 6.0,     # mRNA 半衰期 6 h
        "k_deg_p": np.log(2) / 48.0,    # Protein 半衰期 48 h

        # Antibody module
        "k_neut": 1e-14,
        "k_ab_max": 0.05,
        "EC50_ab": 1e10,
        "k_deg_ab": 0.005,

        # -------------------------------------------------------------
        # Kidney proximal-tubule multilevel uptake module
        # -------------------------------------------------------------
        # Luminal/filtrate effective volumes, mL. These are lumped spaces,
        # not direct anatomical measurements.
        "V_kidney_filtrate": 0.02,
        "V_kidney_pt_lumen": 0.04,

        # Route 1: vascular -> glomerular filtrate -> proximal tubule lumen.
        # AAV is large, so glomerular filtration is represented as a small,
        # size-restricted apparent rate.
        "k_glom_filter": 0.004,       # 1/h, kidney vascular -> filtrate
        "k_filtrate_to_pt": 1.5,      # 1/h, filtrate -> proximal tubule lumen
        "k_urine_flow": 0.25,         # 1/h, tubular lumen -> urine sink

        # Route 2: kidney interstitial/basolateral uptake into tubular cells.
        "Bmax_pt_apical": 5e7,        # vg-equivalent apical binding capacity
        "Bmax_pt_bsl": 2e7,           # vg-equivalent basolateral capacity
        "k_pt_apical_on": 2e-11,      # mL/(vg*h), lumped megalin/cubilin-like binding
        "k_pt_apical_off": 0.05,      # 1/h
        "k_pt_apical_int": 0.25,      # 1/h
        "k_pt_bsl_on": 6e-12,         # mL/(vg*h)
        "k_pt_bsl_off": 0.04,         # 1/h
        "k_pt_bsl_int": 0.12,         # 1/h

        # Proximal-tubule endosomal trafficking.
        "k_kidney_ee_rec": 0.35,      # 1/h, early endosome -> recycling tubule
        "k_kidney_rec_return": 0.25,  # 1/h, recycling -> lumen/extracellular return
        "k_kidney_ee_le": 0.25,       # 1/h, early -> late endosome
        "k_kidney_deg_ee": 0.02,      # 1/h
        "k_kidney_escape": 0.006,     # 1/h, endosomal escape
        "k_kidney_lys": 0.10,         # 1/h, late endosome -> lysosome
        "k_kidney_le_deg": 0.01,      # 1/h
        "k_kidney_lys_deg": 0.08,     # 1/h

        # Kidney intracellular AAV fate and expression.
        "k_kidney_nuc": 0.018,
        "k_kidney_uncoat_cyto": 0.004,
        "k_kidney_uncoat_nuc": 0.018,
        "k_kidney_deg_cyto": 0.01,
        "k_kidney_deg_ncap": 0.005,
        "k_kidney_ds": 0.01,
        "k_kidney_deg_ss": 0.02,
        "k_kidney_epi": 0.01,
        "k_kidney_deg_ds": 0.005,
        "k_kidney_loss_epi": 0.015,
        "k_kidney_tx": 1.5,
        "h_kidney_tx": 1.2,
        "EC50_kidney_tx": 100.0,
        "k_kidney_tl": 4.0,
        "k_kidney_deg_m": np.log(2) / 6.0,
        "k_kidney_deg_p": np.log(2) / 48.0,
    }


def apply_design_preset(
    p: Dict[str, float | str],
    capsid: str = "baseline_AAV",
    promoter: str = "ubiquitous",
    overrides: Dict[str, float | str] | None = None,
) -> Dict[str, float | str]:
    """Return a scenario-specific copy of the base parameters.

    Capsid presets modify exposure, receptor entry, and trafficking. Promoter
    presets modify expression after episome formation. This keeps vector design
    decisions explicit instead of hiding them inside one-off parameter edits.
    """
    if capsid not in CAPSID_PRESETS:
        raise ValueError(f"Unknown capsid preset: {capsid}")
    if promoter not in PROMOTER_PRESETS:
        raise ValueError(f"Unknown promoter preset: {promoter}")

    q = dict(p)
    q.update(CAPSID_PRESETS[capsid])
    q.update(PROMOTER_PRESETS[promoter])
    if overrides:
        q.update(overrides)
    q["capsid_preset"] = capsid
    q["promoter_preset"] = promoter
    return q


def make_initial_condition(p: Dict[str, float | str]) -> np.ndarray:
    y0 = np.zeros(len(IDX), dtype=float)

    if p["administration"] == "bolus":
        y0[IDX["A_blood"]] = float(p["dose_vg"])
        y0[IDX["Dose_in"]] = float(p["dose_vg"])
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

    # Organ vascular loss and interstitial loss.
    # The first term is the original organ-specific RES / degradation term.
    # The added k_clear_* terms are phenomenological AAV decay terms that make
    # extracellular AAV decline after distribution, producing bell-shaped curves.
    J_res = (float(p[f"k_res_{organ}"]) + float(p["k_clear_vascular"])) * A_v
    J_deg_isf = (
        float(p[f"k_deg_isf_{organ}"])
        + float(p["k_clear_isf"])
        + float(p[f"k_extra_isf_clear_{organ}"])
    ) * A_isf

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
        K_filtrate, K_pt_lumen, K_bound_apical, K_bound_bsl,
        K_EE, K_REC, K_LE, K_LYS, K_CY, K_Ncap, K_Nss, K_Nds,
        K_Epi, K_M, K_P, K_Urine, K_Deg,
        Dose_in, Loss_blood_clear, Loss_vascular_res_clear,
        Loss_isf_clear, Loss_neutralized, Loss_liver_cell,
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
    J_dose_input = dose_input_rate(t, p)
    J_blood_clear = float(p["CL_blood"]) * (max(A_blood, 0.0) / float(p["V_blood"])) + float(p["k_clear_blood"]) * max(A_blood, 0.0)

    dA_blood = (
        J_dose_input
        - sum(J_blood_to_v.values())
        - J_blood_clear
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
    tx = (
        float(p["k_tx"])
        * (Epi_eff ** h)
        / (float(p["EC50_tx"]) ** h + Epi_eff ** h + 1e-30)
    )
    tx = min(tx, 2.0)
    dM = tx - p["k_deg_m"] * M          # mRNA 一阶降解
    dP = p["k_tl"] * M - p["k_deg_p"] * P  # Protein 一阶降解

    Ag = max(A_blood + sum(A_v.values()) + 0.5 * sum(A_isf.values()) + LE, 0.0)
    dAb = float(p["k_ab_max"]) * Ag / (float(p["EC50_ab"]) + Ag + 1e-30) - float(p["k_deg_ab"]) * Ab

    dDose_in = J_dose_input
    dLoss_blood_clear = J_blood_clear
    dLoss_vascular_res_clear = sum(max(v, 0.0) for v in J_res.values())
    dLoss_isf_clear = sum(max(v, 0.0) for v in J_deg_isf.values())
    dLoss_neutralized = J_neut_blood
    dLoss_liver_cell = (
        (float(p["k_rec"]) + float(p["k_deg_ee"])) * max(EE, 0.0)
        + float(p["k_lys"]) * max(LE, 0.0)
        + float(p["k_deg_cyto"]) * max(CY, 0.0)
        + float(p["k_deg_ncap"]) * max(Ncap, 0.0)
        + float(p["k_deg_ss"]) * max(Nss, 0.0)
        + float(p["k_deg_ds"]) * max(Nds, 0.0)
        + (float(p["k_loss_epi"]) + float(p["k_dil"])) * max(Epi, 0.0)
    )

    # -----------------------------------------------------------------
    # New kidney proximal-tubule module
    # -----------------------------------------------------------------
    # Route 1: small apparent filtration of AAV from kidney vascular space.
    J_glom_filter = float(p["k_glom_filter"]) * max(A_kidney_v, 0.0)
    J_filtrate_to_pt = float(p["k_filtrate_to_pt"]) * max(K_filtrate, 0.0)
    J_urine = float(p["k_urine_flow"]) * max(K_pt_lumen, 0.0)

    # Apical brush-border binding from proximal-tubule lumen.
    C_pt_lumen = max(K_pt_lumen, 0.0) / float(p["V_kidney_pt_lumen"])
    B_apical_eff = max(K_bound_apical, 0.0)
    Bmax_apical_free = max(float(p["Bmax_pt_apical"]) - B_apical_eff, 0.0)
    J_bind_apical = (
        float(p["k_pt_apical_on"]) * C_pt_lumen * Bmax_apical_free
        - float(p["k_pt_apical_off"]) * B_apical_eff
    )
    if K_bound_apical <= 0.0 and J_bind_apical < 0.0:
        J_bind_apical = 0.0

    # Basolateral/interstitial binding from kidney ISF.
    C_kidney_isf = max(A_kidney_isf, 0.0) / float(p["V_kidney_isf"])
    B_bsl_eff = max(K_bound_bsl, 0.0)
    Bmax_bsl_free = max(float(p["Bmax_pt_bsl"]) - B_bsl_eff, 0.0)
    J_bind_bsl = (
        float(p["k_pt_bsl_on"]) * C_kidney_isf * Bmax_bsl_free
        - float(p["k_pt_bsl_off"]) * B_bsl_eff
    )
    if K_bound_bsl <= 0.0 and J_bind_bsl < 0.0:
        J_bind_bsl = 0.0

    J_int_apical = float(p["k_pt_apical_int"]) * max(K_bound_apical, 0.0)
    J_int_bsl = float(p["k_pt_bsl_int"]) * max(K_bound_bsl, 0.0)

    # Subtract kidney-specific filtration and basolateral binding from systemic kidney pools.
    dA_v["kidney"] -= J_glom_filter
    dA_isf["kidney"] -= J_bind_bsl

    dK_filtrate = J_glom_filter - J_filtrate_to_pt
    dK_pt_lumen = (
        J_filtrate_to_pt
        - J_bind_apical
        - J_urine
        + float(p["k_kidney_rec_return"]) * max(K_REC, 0.0)
    )
    dK_bound_apical = J_bind_apical - J_int_apical
    dK_bound_bsl = J_bind_bsl - J_int_bsl

    dK_EE = (
        J_int_apical + J_int_bsl
        - (
            float(p["k_kidney_ee_rec"])
            + float(p["k_kidney_ee_le"])
            + float(p["k_kidney_deg_ee"])
        ) * K_EE
    )
    dK_REC = float(p["k_kidney_ee_rec"]) * K_EE - float(p["k_kidney_rec_return"]) * K_REC
    dK_LE = (
        float(p["k_kidney_ee_le"]) * K_EE
        - (
            float(p["k_kidney_escape"])
            + float(p["k_kidney_lys"])
            + float(p["k_kidney_le_deg"])
        ) * K_LE
    )
    dK_LYS = float(p["k_kidney_lys"]) * K_LE - float(p["k_kidney_lys_deg"]) * K_LYS
    dK_CY = (
        float(p["k_kidney_escape"]) * K_LE
        - (
            float(p["k_kidney_nuc"])
            + float(p["k_kidney_uncoat_cyto"])
            + float(p["k_kidney_deg_cyto"])
        ) * K_CY
    )
    dK_Ncap = (
        float(p["k_kidney_nuc"]) * K_CY
        - (
            float(p["k_kidney_uncoat_nuc"])
            + float(p["k_kidney_deg_ncap"])
        ) * K_Ncap
    )
    dK_Nss = (
        float(p["k_kidney_uncoat_cyto"]) * K_CY
        + float(p["k_kidney_uncoat_nuc"]) * K_Ncap
        - (
            float(p["k_kidney_ds"])
            + float(p["k_kidney_deg_ss"])
        ) * K_Nss
    )
    dK_Nds = float(p["k_kidney_ds"]) * K_Nss - (
        float(p["k_kidney_epi"]) + float(p["k_kidney_deg_ds"])
    ) * K_Nds
    dK_Epi = float(p["k_kidney_epi"]) * K_Nds - float(p["k_kidney_loss_epi"]) * K_Epi

    K_Epi_eff = max(K_Epi, 0.0)
    h_k = float(p["h_kidney_tx"])
    tx_k = (
        float(p["k_kidney_tx"])
        * (K_Epi_eff ** h_k)
        / (float(p["EC50_kidney_tx"]) ** h_k + K_Epi_eff ** h_k + 1e-30)
    )
    tx_k = min(tx_k, float(p["k_kidney_tx"]))
    dK_M = tx_k - float(p["k_kidney_deg_m"]) * K_M
    dK_P = float(p["k_kidney_tl"]) * K_M - float(p["k_kidney_deg_p"]) * K_P
    dK_Urine = J_urine
    dK_Deg = (
        float(p["k_kidney_deg_ee"]) * max(K_EE, 0.0)
        + float(p["k_kidney_le_deg"]) * max(K_LE, 0.0)
        + float(p["k_kidney_lys_deg"]) * max(K_LYS, 0.0)
        + float(p["k_kidney_deg_cyto"]) * max(K_CY, 0.0)
        + float(p["k_kidney_deg_ncap"]) * max(K_Ncap, 0.0)
        + float(p["k_kidney_deg_ss"]) * max(K_Nss, 0.0)
        + float(p["k_kidney_deg_ds"]) * max(K_Nds, 0.0)
        + float(p["k_kidney_loss_epi"]) * max(K_Epi, 0.0)
    )

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
        dK_filtrate, dK_pt_lumen, dK_bound_apical, dK_bound_bsl,
        dK_EE, dK_REC, dK_LE, dK_LYS, dK_CY, dK_Ncap, dK_Nss, dK_Nds,
        dK_Epi, dK_M, dK_P, dK_Urine, dK_Deg,
        dDose_in, dLoss_blood_clear, dLoss_vascular_res_clear,
        dLoss_isf_clear, dLoss_neutralized, dLoss_liver_cell,
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
    """0 to 56 days; dense enough in the first 48 h to show rise-and-decay."""
    return np.unique(np.r_[
        np.linspace(0.0, 2.0, 300, endpoint=False),          # first 2 h
        np.linspace(2.0, 48.0, 700, endpoint=False),         # first 2 days
        np.linspace(48.0, 24.0 * 7.0, 300, endpoint=False),  # day 2 to 7
        np.linspace(24.0 * 7.0, 24.0 * 56.0, 500),           # day 7 to 56
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


def total_extracellular_aav(sol: SimpleSolution) -> np.ndarray:
    """Total extracellular AAV amount: blood + organ vascular + organ ISF."""
    total = sol.y[IDX["A_blood"]].copy()
    for organ in ORGANS:
        total = total + sol.y[IDX[f"A_{organ}_v"]] + sol.y[IDX[f"A_{organ}_isf"]]
    return total


def total_liver_vector_aav(sol: SimpleSolution) -> np.ndarray:
    """AAV-like liver cellular states, excluding mRNA/protein/antibody."""
    states = ["B", "EE", "LE", "CY", "Ncap", "Nss", "Nds", "Epi"]
    total = np.zeros_like(sol.t)
    for state in states:
        total = total + sol.y[IDX[state]]
    return total


def total_kidney_vector_aav(sol: SimpleSolution) -> np.ndarray:
    """Kidney vector states, excluding mRNA/protein but including urine/deg sinks."""
    states = [
        "K_filtrate", "K_pt_lumen", "K_bound_apical", "K_bound_bsl",
        "K_EE", "K_REC", "K_LE", "K_LYS", "K_CY", "K_Ncap",
        "K_Nss", "K_Nds", "K_Epi",
    ]
    total = np.zeros_like(sol.t)
    for state in states:
        total = total + sol.y[IDX[state]]
    return total


def total_accounted_aav(sol: SimpleSolution) -> np.ndarray:
    """All tracked vector plus cumulative sinks for mass-balance auditing."""
    total = total_extracellular_aav(sol) + total_liver_vector_aav(sol) + total_kidney_vector_aav(sol)
    for state in [
        "K_Urine", "K_Deg", "Loss_blood_clear", "Loss_vascular_res_clear",
        "Loss_isf_clear", "Loss_neutralized", "Loss_liver_cell",
    ]:
        total = total + sol.y[IDX[state]]
    return total


def mass_balance_error(sol: SimpleSolution) -> np.ndarray:
    delivered = np.maximum(sol.y[IDX["Dose_in"]], 1e-30)
    return (total_accounted_aav(sol) - sol.y[IDX["Dose_in"]]) / delivered


def annotate_peak(ax, t: np.ndarray, y: np.ndarray, label: str, time_unit: str = "h") -> None:
    """Mark the peak of a curve without disturbing the data scale."""
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    valid = np.isfinite(y)
    if not np.any(valid):
        return
    idx_valid = np.where(valid)[0]
    idx_peak = idx_valid[int(np.nanargmax(y[valid]))]
    ax.scatter([t[idx_peak]], [y[idx_peak]], s=30, zorder=5)
    ax.annotate(
        f"{label} peak\n{t[idx_peak]:.2g} {time_unit}",
        xy=(t[idx_peak], y[idx_peak]),
        xytext=(6, -18),
        textcoords="offset points",
        fontsize=8,
        va="top",
    )


def plot_bell_shaped_aav_decay(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    """0-48 h extracellular AAV profile showing rise, peak, and decay."""
    t_h_all = sol_long.t
    mask = t_h_all <= PLOT_DECAY_WINDOW_H
    t_h = t_h_all[mask]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    C_blood = concentration(sol_long, "A_blood", float(p["V_blood"]))[mask]
    ax.plot(t_h, log_safe(C_blood), label="blood", color="black", linestyle="--", linewidth=2.2)

    representative_organs = ["liver", "spleen", "kidney", "muscle"]
    for organ in representative_organs:
        C_isf = concentration(sol_long, f"A_{organ}_isf", float(p[f"V_{organ}_isf"]))[mask]
        ax.plot(t_h, log_safe(C_isf), label=f"{organ}_ISF", color=COLORS[organ], linewidth=1.9)

    annotate_peak(ax, t_h, C_blood, "blood", "h")
    ax.set_yscale("log")
    ax.set_xlim(0.0, PLOT_DECAY_WINDOW_H)
    ax.set_xlabel("Time after dosing (h)")
    ax.set_ylabel("AAV concentration (vg/mL, log scale)")
    ax.set_title("Extracellular AAV exposure: rise then clearance")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    ax = axes[1]
    total = total_extracellular_aav(sol_long)[mask]
    ax.plot(t_h, log_safe(total), color="black", linewidth=2.2)
    annotate_peak(ax, t_h, total, "total", "h")
    ax.set_yscale("log")
    ax.set_xlim(0.0, PLOT_DECAY_WINDOW_H)
    ax.set_xlabel("Time after dosing (h)")
    ax.set_ylabel("Total extracellular AAV amount (vg, log scale)")
    ax.set_title("Body-level extracellular burden decays after dosing")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)

    plt.tight_layout()
    save_or_show("03_bell_shaped_aav_decay_48h.png")


def plot_short_distribution(sol_short: SimpleSolution, p: Dict[str, float | str]) -> None:
    """Early exposure plots corrected for a 5 min infusion peak and true log axes."""
    t_min = sol_short.t * 60.0
    t_inf_min = float(p["T_inf_h"]) * 60.0

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    C_blood = concentration(sol_short, "A_blood", float(p["V_blood"]))
    ax.plot(t_min, C_blood, label="blood", color="black", linestyle="--", linewidth=2.2)
    for organ in ["liver", "kidney", "spleen", "lung"]:
        C_v = concentration(sol_short, f"A_{organ}_v", float(p[f"V_{organ}_v"]))
        ax.plot(t_min, C_v, label=f"{organ}_vascular", color=COLORS[organ], linewidth=1.8)
    ax.axvline(t_inf_min, linestyle=":", linewidth=1.5, label=f"infusion end = {t_inf_min:.0f} min")
    early_mask = t_min <= 15.0
    annotate_peak(ax, t_min[early_mask], C_blood[early_mask], "blood", "min")
    ax.set_xlim(0, 15)
    ax.set_xlabel("Time after dosing (min)")
    ax.set_ylabel("AAV concentration (vg/mL)")
    ax.set_title("Early vascular exposure, 0-15 min")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    ax = axes[1]
    for organ in ["liver", "kidney", "spleen", "muscle"]:
        C_isf = concentration(sol_short, f"A_{organ}_isf", float(p[f"V_{organ}_isf"]))
        ax.plot(t_min, C_isf, label=f"{organ}_ISF", color=COLORS[organ], linewidth=1.8)
    ax.axvline(t_inf_min, linestyle=":", linewidth=1.5, label=f"infusion end = {t_inf_min:.0f} min")
    ax.set_xlim(0, 30)
    ax.set_xlabel("Time after dosing (min)")
    ax.set_ylabel("AAV concentration (vg/mL)")
    ax.set_title("Interstitial entry is delayed relative to blood")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    plt.tight_layout()
    save_or_show("01_early_distribution_linear_15to30min.png")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    for organ in ORGANS:
        C_v = concentration(sol_short, f"A_{organ}_v", float(p[f"V_{organ}_v"]))
        ax.plot(sol_short.t, log_safe(C_v), label=f"{organ}_vascular", color=COLORS[organ], linewidth=1.7)
    ax.plot(sol_short.t, log_safe(C_blood), label="blood", color="black", linestyle="--", linewidth=2.0)
    ax.axvline(float(p["T_inf_h"]), linestyle=":", linewidth=1.3)
    ax.set_yscale("log")
    ax.set_xlim(0, 2.0)
    ax.set_xlabel("Time after dosing (h)")
    ax.set_ylabel("AAV concentration (vg/mL, log scale)")
    ax.set_title("Vascular distribution, 0-2 h")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    ax = axes[1]
    for organ in ORGANS:
        C_isf = concentration(sol_short, f"A_{organ}_isf", float(p[f"V_{organ}_isf"]))
        ax.plot(sol_short.t, log_safe(C_isf), label=f"{organ}_ISF", color=COLORS[organ], linewidth=1.7)
    ax.axvline(float(p["T_inf_h"]), linestyle=":", linewidth=1.3)
    ax.set_yscale("log")
    ax.set_xlim(0, 2.0)
    ax.set_xlabel("Time after dosing (h)")
    ax.set_ylabel("AAV concentration (vg/mL, log scale)")
    ax.set_title("Permeability-limited ISF distribution, 0-2 h")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    plt.tight_layout()
    save_or_show("02_short_distribution_log_2h.png")


def plot_long_states(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    plt.figure(figsize=(15, 4))
    t_day = sol_long.t / 24.0
    # 1. Episome
    plt.subplot(1, 3, 1)
    plt.plot(t_day, sol_long.y[IDX["Epi"]], color="red", linewidth=2.0)
    plt.xlabel("Time (day)")
    plt.ylabel("Episome (a.u.)")
    plt.title("Episome (Epi)")
    plt.grid(True, linestyle="--", alpha=0.35)

    # 2. mRNA
    plt.subplot(1, 3, 2)
    plt.plot(t_day, sol_long.y[IDX["M"]], color="blue", linewidth=2.0)
    plt.xlabel("Time (day)")
    plt.ylabel("mRNA (a.u.)")
    plt.title("mRNA (M)")
    plt.grid(True, linestyle="--", alpha=0.35)

    # 3. Protein
    plt.subplot(1, 3, 3)
    plt.plot(t_day, sol_long.y[IDX["P"]], color="green", linewidth=2.0)
    plt.xlabel("Time (day)")
    plt.ylabel("Protein (a.u.)")
    plt.title("Protein (P)")
    plt.grid(True, linestyle="--", alpha=0.35)

    plt.tight_layout()
    save_or_show("04_liver_expression_56d_2.png")

    plt.figure(figsize=(8, 4))
    plt.plot(t_day, log_safe(sol_long.y[IDX["Ab"]]), color="magenta", linewidth=2.0)
    #plt.yscale("log")
    plt.xlabel("Time (day)")
    plt.ylabel("Antibody level (a.u.)")
    plt.title("Simplified antibody kinetics, 0-56 days")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.tight_layout()
    save_or_show("05_antibody_56d.png")


def auto_ylim(ax, values: Iterable[np.ndarray], pad: float = 0.08, log_scale: bool = False) -> None:
    """Set pleasant y-limits for plots with very different state magnitudes."""
    data = np.concatenate([np.asarray(v, dtype=float).ravel() for v in values])
    data = data[np.isfinite(data)]
    if log_scale:
        data = data[data > 0]
    if data.size == 0:
        return
    ymin = float(np.nanmin(data))
    ymax = float(np.nanmax(data))
    if ymax <= ymin:
        ymax = ymin + 1.0
    if log_scale:
        ax.set_ylim(max(ymin / 2.0, 1e-30), ymax * 2.0)
    else:
        span = ymax - ymin
        ax.set_ylim(max(0.0, ymin - pad * span), ymax + pad * span)


def plot_kidney_module(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    """Plot multilevel kidney proximal-tubule uptake and expression."""
    t_h = sol_long.t
    t_day = sol_long.t / 24.0

    fig, axes = plt.subplots(2, 3, figsize=(18, 9))
    axes = axes.ravel()

    ax = axes[0]
    mask_h8 = t_h <= 8.0
    kidney_states = ["K_filtrate", "K_pt_lumen", "K_bound_apical", "K_bound_bsl"]
    plotted = []
    for state in kidney_states:
        y = log_safe(sol_long.y[IDX[state]][mask_h8])
        plotted.append(y)
        ax.plot(t_h[mask_h8], y, label=state, linewidth=1.8)
    ax.set_yscale("log")
    auto_ylim(ax, plotted, log_scale=True)
    ax.set_xlabel("Time after dosing (h)")
    ax.set_ylabel("AAV amount (vg-equivalent, log scale)")
    ax.set_title("Kidney accessible pools, 0-8 h")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    ax = axes[1]
    mask_h72 = t_h <= 72.0
    intracellular_states = ["K_EE", "K_REC", "K_LE", "K_LYS", "K_CY", "K_Ncap"]
    plotted = []
    for state in intracellular_states:
        y = log_safe(sol_long.y[IDX[state]][mask_h72])
        plotted.append(y)
        ax.plot(t_h[mask_h72], y, label=state, linewidth=1.8)
    ax.set_yscale("log")
    auto_ylim(ax, plotted, log_scale=True)
    ax.set_xlabel("Time after dosing (h)")
    ax.set_ylabel("AAV amount (vg-equivalent, log scale)")
    ax.set_title("Proximal-tubule trafficking, 0-72 h")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    ax = axes[2]
    mask_d21 = t_day <= 21.0
    genome_states = ["K_Nss", "K_Nds", "K_Epi"]
    plotted = []
    for state in genome_states:
        y = log_safe(sol_long.y[IDX[state]][mask_d21])
        plotted.append(y)
        ax.plot(t_day[mask_d21], y, label=state, linewidth=1.8)
    ax.set_yscale("log")
    auto_ylim(ax, plotted, log_scale=True)
    ax.set_xlabel("Time after dosing (day)")
    ax.set_ylabel("Vector genome state (a.u., log scale)")
    ax.set_title("Kidney genome processing, 0-21 d")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    ax = axes[3]
    m = sol_long.y[IDX["K_M"]][mask_d21]
    ax.plot(t_day[mask_d21], m, label="Kidney mRNA", linewidth=2.0)
    annotate_peak(ax, t_day[mask_d21], m, "mRNA", "d")
    auto_ylim(ax, [m], log_scale=False)
    ax.set_xlabel("Time after dosing (day)")
    ax.set_ylabel("mRNA (a.u.)")
    ax.set_title("Kidney transgene mRNA")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    ax = axes[4]
    mask_d56 = t_day <= 56.0
    protein = sol_long.y[IDX["K_P"]][mask_d56]
    ax.plot(t_day[mask_d56], protein, label="Kidney protein", linewidth=2.0)
    annotate_peak(ax, t_day[mask_d56], protein, "protein", "d")
    auto_ylim(ax, [protein], log_scale=False)
    ax.set_xlabel("Time after dosing (day)")
    ax.set_ylabel("Protein output (a.u.)")
    ax.set_title("Kidney transgene protein")
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    ax = axes[5]
    urine = sol_long.y[IDX["K_Urine"]]
    deg = sol_long.y[IDX["K_Deg"]]
    kidney_vector = total_kidney_vector_aav(sol_long)
    ax.plot(t_day, log_safe(kidney_vector), label="kidney retained vector", linewidth=2.0)
    ax.plot(t_day, log_safe(urine), label="urine sink", linewidth=1.8)
    ax.plot(t_day, log_safe(deg), label="intracellular degradation", linewidth=1.8)
    ax.set_yscale("log")
    ax.set_xlim(0, 56.0)
    ax.set_xlabel("Time after dosing (day)")
    ax.set_ylabel("AAV amount (vg-equivalent, log scale)")
    ax.set_title("Kidney retention vs loss routes")
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    plt.tight_layout()
    save_or_show("06_kidney_multilevel_module.png")


def plot_liver_vs_kidney_expression(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    """Compare retained liver module with the new kidney module."""
    t_day = sol_long.t / 24.0

    plt.figure(figsize=(14, 5))

    plt.subplot(1, 2, 1)
    plt.plot(t_day, log_safe(sol_long.y[IDX["Epi"]]), label="liver Epi", linewidth=2.0)
    plt.plot(t_day, log_safe(sol_long.y[IDX["K_Epi"]]), label="kidney Epi", linewidth=2.0)
    plt.xlabel("Time (day)")
    plt.ylabel("Episome / expression-competent genome (a.u.)")
    plt.title("Liver vs kidney episome")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(t_day, sol_long.y[IDX["P"]], label="liver protein", linewidth=2.0)
    plt.plot(t_day, sol_long.y[IDX["K_P"]], label="kidney protein", linewidth=2.0)
    plt.xlabel("Time (day)")
    plt.ylabel("Protein output (a.u.)")
    plt.title("Liver vs kidney transgene expression")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend()

    plt.tight_layout()
    save_or_show("07_liver_vs_kidney_expression.png")


def plot_mass_balance(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    """Audit whether vector mass is accounted for by states and cumulative sinks."""
    t_day = sol_long.t / 24.0
    plt.figure(figsize=(14, 5))

    plt.subplot(1, 2, 1)
    plt.plot(t_day, total_extracellular_aav(sol_long), label="extracellular", linewidth=2.0)
    plt.plot(t_day, total_liver_vector_aav(sol_long), label="liver vector states", linewidth=2.0)
    plt.plot(t_day, total_kidney_vector_aav(sol_long), label="kidney vector states", linewidth=2.0)
    plt.plot(t_day, sol_long.y[IDX["K_Urine"]], label="urine sink", linewidth=1.8)
    plt.plot(t_day, sol_long.y[IDX["Loss_vascular_res_clear"]], label="vascular/RES sink", linewidth=1.8)
    plt.xlabel("Time (day)")
    plt.ylabel("AAV amount (vg-equivalent)")
    plt.title("Vector mass allocation")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.subplot(1, 2, 2)
    plt.plot(t_day, mass_balance_error(sol_long), color="black", linewidth=2.0)
    plt.axhline(0.0, color="gray", linestyle="--", linewidth=1.0)
    plt.xlabel("Time (day)")
    plt.ylabel("(accounted - delivered) / delivered")
    plt.title("Mass-balance error")
    plt.grid(True, linestyle="--", alpha=0.35)

    plt.tight_layout()
    save_or_show("08_mass_balance_audit.png")


def plot_design_scenarios(base_p: Dict[str, float | str], t_eval_long: np.ndarray) -> None:
    """Compare capsid/promoter/renal-route design hypotheses without editing core code."""
    scenarios = [
        ("baseline", "baseline_AAV", "ubiquitous", {}),
        ("liver_detargeted", "liver_detargeted", "ubiquitous", {}),
        ("kidney_tropic", "kidney_tropic", "kidney_biased", {}),
        ("escape_enhanced", "endosomal_escape_enhanced", "ubiquitous", {}),
        ("renal_filtration_limited", "baseline_AAV", "ubiquitous", {
            "k_glom_filter": 0.0015,
            "k_pt_apical_on": 1.0e-11,
            "Bmax_pt_apical": 2.5e7,
        }),
        ("renal_basolateral_tropism", "kidney_tropic", "kidney_biased", {
            "k_glom_filter": 0.0015,
            "Bmax_pt_bsl": 6.0e7,
            "k_pt_bsl_on": 1.6e-11,
            "k_pt_bsl_int": 0.16,
            "k_kidney_escape": 0.012,
        }),
    ]

    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    for label, capsid, promoter, overrides in scenarios:
        p = apply_design_preset(base_p, capsid=capsid, promoter=promoter, overrides=overrides)
        sol = solve_model(t_eval_long, make_initial_condition(p), p, post_infusion_max_step=1.0)
        t_day = sol.t / 24.0
        axes[0].plot(t_day, log_safe(sol.y[IDX["Epi"]]), label=label, linewidth=2.0)
        axes[1].plot(t_day, log_safe(sol.y[IDX["K_Epi"]]), label=label, linewidth=2.0)

        C_kidney_isf = concentration(sol, "A_kidney_isf", float(p["V_kidney_isf"]))
        C_liver_isf = concentration(sol, "A_liver_isf", float(p["V_liver_isf"]))
        rows.append({
            "scenario": label,
            "capsid": capsid,
            "promoter": promoter,
            "peak_liver_epi": float(np.nanmax(sol.y[IDX["Epi"]])),
            "peak_kidney_epi": float(np.nanmax(sol.y[IDX["K_Epi"]])),
            "peak_liver_protein": float(np.nanmax(sol.y[IDX["P"]])),
            "peak_kidney_protein": float(np.nanmax(sol.y[IDX["K_P"]])),
            "auc_liver_isf": auc_trapz(C_liver_isf, sol.t),
            "auc_kidney_isf": auc_trapz(C_kidney_isf, sol.t),
            "final_urine_loss": float(sol.y[IDX["K_Urine"], -1]),
            "final_kidney_degradation": float(sol.y[IDX["K_Deg"], -1]),
            "final_mass_balance_error": float(mass_balance_error(sol)[-1]),
        })

    axes[0].set_yscale("log")
    axes[0].set_xlabel("Time after dosing (day)")
    axes[0].set_ylabel("Liver episome (a.u., log scale)")
    axes[0].set_title("Capsid/promoter scenario: liver")
    axes[0].grid(True, which="both", linestyle="--", alpha=0.35)
    axes[0].legend(fontsize=8)

    axes[1].set_yscale("log")
    axes[1].set_xlabel("Time after dosing (day)")
    axes[1].set_ylabel("Kidney episome (a.u., log scale)")
    axes[1].set_title("Capsid/promoter/renal-route scenario: kidney")
    axes[1].grid(True, which="both", linestyle="--", alpha=0.35)
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    save_or_show("09_design_scenario_comparison.png")

    if SAVE_FIGURES:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        metrics_path = OUTPUT_DIR / "09_design_scenario_metrics.csv"
        header = list(rows[0].keys())
        with metrics_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(rows)


def plot_spatial_pk_demo(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    """Minimal 1D advection-diffusion-reaction bridge toward CFD/spatial PK."""
    n = 140
    length_cm = 1.0
    x = np.linspace(0.0, length_cm, n)
    dx = x[1] - x[0]
    dt = 0.0008
    t_end = 8.0
    steps = int(t_end / dt)

    C_blood = concentration(sol_long, "A_blood", float(p["V_blood"]))

    def run_case(flow_cm_h: float, wall_access: float) -> np.ndarray:
        C = np.zeros(n)
        B = np.zeros(n)
        I = np.zeros(n)
        E = np.zeros(n)
        D = 4.5e-4
        Bmax = 1.0
        kon = 0.9
        koff = 0.12
        kint = 0.45
        kescape = 0.045
        klys = 0.18
        kloss_epi = 0.018

        for step in range(steps):
            t = step * dt
            Cin = np.interp(t, sol_long.t, C_blood)
            C_old = C.copy()
            B_old = B.copy()
            I_old = I.copy()
            E_old = E.copy()

            C[0] = Cin / max(np.nanmax(C_blood), 1e-30)
            adv = -flow_cm_h * (C_old[1:-1] - C_old[:-2]) / dx
            diff = D * (C_old[2:] - 2.0 * C_old[1:-1] + C_old[:-2]) / (dx * dx)
            free = np.maximum(Bmax - B_old[1:-1], 0.0)
            bind = wall_access * kon * C_old[1:-1] * free
            uptake = np.maximum(bind - koff * B_old[1:-1], 0.0)
            C[1:-1] = np.maximum(C_old[1:-1] + dt * (adv + diff - uptake), 0.0)
            C[-1] = C[-2]

            free_all = np.maximum(Bmax - B_old, 0.0)
            bind_all = wall_access * kon * C_old * free_all
            internalize = kint * B_old
            escape = kescape * I_old
            B = np.maximum(B_old + dt * (bind_all - koff * B_old - internalize), 0.0)
            I = np.maximum(I_old + dt * (internalize - escape - klys * I_old), 0.0)
            E = np.maximum(E_old + dt * (escape - kloss_epi * E_old), 0.0)
        return E

    baseline = run_case(flow_cm_h=2.6, wall_access=0.55)
    enhanced = run_case(flow_cm_h=2.6, wall_access=1.20)
    slow_flow = run_case(flow_cm_h=1.1, wall_access=0.90)

    plt.figure(figsize=(12, 5))
    plt.plot(x, baseline, label="baseline IV-like", linewidth=2.0)
    plt.plot(x, enhanced, label="capsid access enhanced", linewidth=2.0)
    plt.plot(x, slow_flow, label="slow-flow local trapping", linewidth=2.0)
    plt.xlabel("Normalized vascular / tissue axis (cm)")
    plt.ylabel("Final local episome E(x), normalized a.u.")
    plt.title("1D spatial PK bridge: flow and wall uptake reshape expression")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend()
    plt.tight_layout()
    save_or_show("10_spatial_pk_1d_demo.png")


def auc_trapz(y: np.ndarray, x: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def time_of_peak(t: np.ndarray, y: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    if y.size == 0 or not np.any(np.isfinite(y)):
        return float("nan")
    return float(t[int(np.nanargmax(y))])


def collect_summary_metrics(sol_long: SimpleSolution, p: Dict[str, float | str]) -> Dict[str, float | str]:
    """Collect PI-facing metrics in a CSV-friendly dictionary."""
    C_blood_long = concentration(sol_long, "A_blood", float(p["V_blood"]))
    C_kidney_isf = concentration(sol_long, "A_kidney_isf", float(p["V_kidney_isf"]))
    C_liver_isf = concentration(sol_long, "A_liver_isf", float(p["V_liver_isf"]))
    auc_kidney_isf = auc_trapz(C_kidney_isf, sol_long.t)
    auc_liver_isf = auc_trapz(C_liver_isf, sol_long.t)
    total_aav = total_extracellular_aav(sol_long)

    return {
        "administration": str(p["administration"]),
        "dose_vg": float(p["dose_vg"]),
        "infusion_duration_min": float(p["T_inf_h"]) * 60.0,
        "q_scale_effective_exchange": float(p["Q_scale"]),
        "auc_blood_vg_h_per_ml": auc_trapz(C_blood_long, sol_long.t),
        "cmax_blood_vg_per_ml": float(np.nanmax(C_blood_long)),
        "tmax_blood_min": 60.0 * time_of_peak(sol_long.t, C_blood_long),
        "final_fraction_extracellular_remaining": float(total_aav[-1] / max(float(p["dose_vg"]), 1e-30)),
        "auc_liver_isf_vg_h_per_ml": auc_liver_isf,
        "auc_kidney_isf_vg_h_per_ml": auc_kidney_isf,
        "kidney_liver_isf_auc_ratio": float(auc_kidney_isf / max(auc_liver_isf, 1e-30)),
        "peak_liver_epi": float(np.nanmax(sol_long.y[IDX["Epi"]])),
        "t_peak_liver_epi_day": time_of_peak(sol_long.t / 24.0, sol_long.y[IDX["Epi"]]),
        "peak_liver_mrna": float(np.nanmax(sol_long.y[IDX["M"]])),
        "peak_liver_protein": float(np.nanmax(sol_long.y[IDX["P"]])),
        "peak_kidney_apical_bound": float(np.nanmax(sol_long.y[IDX["K_bound_apical"]])),
        "peak_kidney_basolateral_bound": float(np.nanmax(sol_long.y[IDX["K_bound_bsl"]])),
        "basolateral_to_apical_binding_peak_ratio": float(
            np.nanmax(sol_long.y[IDX["K_bound_bsl"]]) / max(np.nanmax(sol_long.y[IDX["K_bound_apical"]]), 1e-30)
        ),
        "peak_kidney_early_endosome": float(np.nanmax(sol_long.y[IDX["K_EE"]])),
        "peak_kidney_epi": float(np.nanmax(sol_long.y[IDX["K_Epi"]])),
        "t_peak_kidney_epi_day": time_of_peak(sol_long.t / 24.0, sol_long.y[IDX["K_Epi"]]),
        "peak_kidney_mrna": float(np.nanmax(sol_long.y[IDX["K_M"]])),
        "t_peak_kidney_mrna_day": time_of_peak(sol_long.t / 24.0, sol_long.y[IDX["K_M"]]),
        "peak_kidney_protein": float(np.nanmax(sol_long.y[IDX["K_P"]])),
        "t_peak_kidney_protein_day": time_of_peak(sol_long.t / 24.0, sol_long.y[IDX["K_P"]]),
        "final_urine_loss_vg_equiv": float(sol_long.y[IDX["K_Urine"], -1]),
        "final_kidney_intracellular_degradation_vg_equiv": float(sol_long.y[IDX["K_Deg"], -1]),
        "final_mass_balance_error": float(mass_balance_error(sol_long)[-1]),
        "peak_antibody_au": float(np.nanmax(sol_long.y[IDX["Ab"]])),
    }


def write_summary_metrics_csv(metrics: Dict[str, float | str], filename: str = "00_summary_metrics.csv") -> None:
    if not SAVE_FIGURES:
        return
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUTPUT_DIR / filename).open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, value])


def print_metrics(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    metrics = collect_summary_metrics(sol_long, p)
    C_blood_long = concentration(sol_long, "A_blood", float(p["V_blood"]))
    print("----- Simulation settings -----")
    print(f"Administration: {p['administration']}")
    print(f"Dose: {p['dose_vg']:.3e} vg")
    print(f"Infusion duration: {float(p['T_inf_h']) * 60:.2f} min")
    print(f"Q_scale: {p['Q_scale']}")
    print(f"Added blood AAV half-life: {BLOOD_AAV_HALF_LIFE_H:.2f} h")
    print(f"Added vascular AAV half-life: {VASCULAR_AAV_HALF_LIFE_H:.2f} h")
    print(f"Added ISF AAV half-life: {ISF_AAV_HALF_LIFE_H:.2f} h")

    print("\n----- Summary metrics -----")
    total_aav = total_extracellular_aav(sol_long)
    print("AUC_blood:", auc_trapz(C_blood_long, sol_long.t))
    print("Cmax_blood:", np.nanmax(C_blood_long))
    print("Tmax_blood_min:", metrics["tmax_blood_min"])
    print("Initial total extracellular AAV:", total_aav[0])
    print("Final total extracellular AAV:", total_aav[-1])
    print("Fraction extracellular AAV remaining at final time:", total_aav[-1] / max(float(p["dose_vg"]), 1e-30))
    for organ in ORGANS:
        peak_isf_amount = np.nanmax(sol_long.y[IDX[f"A_{organ}_isf"]])
        peak_isf_conc = np.nanmax(concentration(sol_long, f"A_{organ}_isf", float(p[f"V_{organ}_isf"])))
        print(f"Peak {organ} ISF amount:", peak_isf_amount)
        print(f"Peak {organ} ISF concentration:", peak_isf_conc)
    print("Peak liver episome:", np.nanmax(sol_long.y[IDX["Epi"]]))
    print("Peak liver mRNA:", np.nanmax(sol_long.y[IDX["M"]]))
    print("Peak liver protein:", np.nanmax(sol_long.y[IDX["P"]]))

    C_kidney_isf = concentration(sol_long, "A_kidney_isf", float(p["V_kidney_isf"]))
    C_liver_isf = concentration(sol_long, "A_liver_isf", float(p["V_liver_isf"]))
    auc_kidney_isf = auc_trapz(C_kidney_isf, sol_long.t)
    auc_liver_isf = auc_trapz(C_liver_isf, sol_long.t)
    print("\n----- Kidney module metrics -----")
    print("AUC_kidney_ISF:", auc_kidney_isf)
    print("AUC_liver_ISF:", auc_liver_isf)
    print("Kidney/Liver ISF AUC ratio:", auc_kidney_isf / max(auc_liver_isf, 1e-30))
    print("Peak kidney apical-bound AAV:", np.nanmax(sol_long.y[IDX["K_bound_apical"]]))
    print("Peak kidney basolateral-bound AAV:", np.nanmax(sol_long.y[IDX["K_bound_bsl"]]))
    print("Peak kidney early endosome:", np.nanmax(sol_long.y[IDX["K_EE"]]))
    print("Peak kidney episome:", np.nanmax(sol_long.y[IDX["K_Epi"]]))
    print("Peak kidney mRNA:", np.nanmax(sol_long.y[IDX["K_M"]]))
    print("Peak kidney protein:", np.nanmax(sol_long.y[IDX["K_P"]]))
    print("Cumulative urinary AAV loss:", sol_long.y[IDX["K_Urine"], -1])
    print("Cumulative kidney intracellular degradation/loss:", sol_long.y[IDX["K_Deg"], -1])
    print("Kidney entry efficiency (Peak K_Epi / AUC kidney ISF):", np.nanmax(sol_long.y[IDX["K_Epi"]]) / max(auc_kidney_isf, 1e-30))
    print("Peak antibody:", np.nanmax(sol_long.y[IDX["Ab"]]))
    print("Final mass-balance error:", metrics["final_mass_balance_error"])


def parse_cli_args() -> argparse.Namespace:
    """Small command-line interface so output folder and common settings are editable."""
    parser = argparse.ArgumentParser(
        description="AAV PBPK + liver/kidney cellular fate + spatial PK demonstration model."
    )
    parser.add_argument("--output-dir", type=str, default=None, help="Folder for figures and CSV outputs.")
    parser.add_argument("--infusion-min", type=float, default=None, help="Infusion duration in minutes; default is 5 min.")
    parser.add_argument("--dose-vg", type=float, default=None, help="Total vector dose in vg.")
    parser.add_argument("--clearance-mode", choices=["mechanistic", "half_life_demo"], default=None)
    parser.add_argument("--show", action="store_true", help="Display plots interactively after saving.")
    parser.add_argument("--no-save", action="store_true", help="Do not save figures/CSV outputs.")
    parser.add_argument("--no-scenarios", action="store_true", help="Skip capsid/promoter scenario simulations.")
    parser.add_argument("--no-spatial", action="store_true", help="Skip 1D spatial PK demo.")
    return parser.parse_args()


def apply_cli_args(args: argparse.Namespace) -> None:
    global OUTPUT_DIR, INFUSION_DURATION_MIN, DOSE_VG, SAVE_FIGURES, SHOW_FIGURES
    global RUN_DESIGN_SCENARIOS, RUN_SPATIAL_PK_DEMO, CLEARANCE_MODE

    if args.output_dir:
        OUTPUT_DIR = Path(args.output_dir)
    if args.infusion_min is not None:
        if args.infusion_min <= 0:
            raise ValueError("--infusion-min must be positive.")
        INFUSION_DURATION_MIN = float(args.infusion_min)
    if args.dose_vg is not None:
        if args.dose_vg <= 0:
            raise ValueError("--dose-vg must be positive.")
        DOSE_VG = float(args.dose_vg)
    if args.clearance_mode:
        CLEARANCE_MODE = args.clearance_mode
    if args.show:
        SHOW_FIGURES = True
    if args.no_save:
        SAVE_FIGURES = False
    if args.no_scenarios:
        RUN_DESIGN_SCENARIOS = False
    if args.no_spatial:
        RUN_SPATIAL_PK_DEMO = False


# ---------------------------------------------------------------------
# Main script
# ---------------------------------------------------------------------
def main() -> None:
    args = parse_cli_args()
    apply_cli_args(args)

    p = make_params()
    y0 = make_initial_condition(p)

    t_eval_short = make_short_grid()
    t_eval_long = make_long_grid()

    sol_short = solve_model(t_eval_short, y0, p, post_infusion_max_step=0.002)
    sol_long = solve_model(t_eval_long, y0, p, post_infusion_max_step=1.0)

    plot_short_distribution(sol_short, p)
    plot_bell_shaped_aav_decay(sol_long, p)
    plot_long_states(sol_long, p)
    plot_kidney_module(sol_long, p)
    plot_liver_vs_kidney_expression(sol_long, p)
    plot_mass_balance(sol_long, p)
    if RUN_DESIGN_SCENARIOS:
        plot_design_scenarios(p, t_eval_long)
    if RUN_SPATIAL_PK_DEMO:
        plot_spatial_pk_demo(sol_long, p)

    metrics = collect_summary_metrics(sol_long, p)
    write_summary_metrics_csv(metrics)
    print_metrics(sol_long, p)

    if SAVE_FIGURES:
        print(f"\nFigures and CSV outputs saved to: {OUTPUT_DIR.resolve()}")



if __name__ == "__main__":
    main()
