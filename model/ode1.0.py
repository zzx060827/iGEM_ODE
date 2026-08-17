"""
AAV PBPK + liver, kidney, and CNS cellular-delivery model.

This version keeps the original liver cellular fate module and adds a separate
kidney module inspired by renal proximal-tubule biology:
blood/kidney vascular -> glomerular filtrate -> proximal tubule lumen ->
apical receptor binding/endocytosis, plus a basolateral ISF uptake route ->
endosome trafficking -> cytosol/nucleus -> episome -> mRNA/protein.

Important modeling note:
All kidney and CNS cellular parameters are phenomenological demonstration
parameters. They are intended for iGEM dry-lab model structure and sensitivity
analysis, not fitted quantitative AAV organ PK constants.
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Dict, Iterable, Tuple

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

MODEL_DIR = Path(__file__).resolve().parent
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

from aav_parameter_evidence import (
    MOUSE_AAV9_BLOOD_HALF_LIFE_H,
    MOUSE_AAV9_CAPSID_HALF_LIFE_H,
    aav9_parameter_payload,
)


ORGANS = ["liver", "spleen", "kidney", "heart", "muscle", "lung", "brain", "rest"]

STATE_NAMES = [
    # systemic extracellular AAV
    "A_blood",
    "A_liver_v", "A_liver_isf",
    "A_spleen_v", "A_spleen_isf",
    "A_kidney_v", "A_kidney_isf",
    "A_heart_v", "A_heart_isf",
    "A_muscle_v", "A_muscle_isf",
    "A_lung_v", "A_lung_isf",
    "A_brain_v", "A_brain_isf",
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

    # CNS delivery: BBB handling followed by neural-cell transduction
    "C_BBB_bound",      # capsid bound to luminal BBB transport sites
    "C_BBB_endo",       # capsid internalized into brain endothelial cells
    "C_bound",          # capsid bound to CNS parenchymal-cell receptors
    "C_EE",             # CNS-cell early endosome
    "C_LE",             # CNS-cell late endosome
    "C_CY",             # escaped cytosolic capsid
    "C_Ncap",           # nuclear capsid
    "C_Nss",            # nuclear single-stranded vector genome
    "C_Nds",            # double-stranded vector genome
    "C_Epi",            # CNS expression-competent episome
    "C_M",              # CNS transgene mRNA
    "C_P",              # CNS transgene protein
    "C_Deg",            # cumulative BBB/intracellular degradation

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
    "brain": "deepskyblue",
    "rest": "brown",
}


# ---------------------------------------------------------------------
# User-facing switches
# ---------------------------------------------------------------------
DOSE_VG = 1e12

# "infusion" gives a visually smoother early curve.
# Set ADMINISTRATION = "bolus" to recover the original instantaneous dose style.
ADMINISTRATION = "infusion"  # allowed: "infusion" or "bolus"
INFUSION_DURATION_MIN = 10.0

# In this one-central-blood reduced-order model, raw cardiac output would make
# all vascular concentrations equilibrate within seconds. Keep the intended
# pulmonary exchange rate explicit and derive the dimensionless multiplier
# from physiological cardiac output inside ``make_params``.
EFFECTIVE_LUNG_EXCHANGE_ML_H = 1.25

# ------------------------------------------------------------------
# Apparent extracellular AAV decay controls
# ------------------------------------------------------------------
# These are the knobs that make the in-vivo AAV exposure become bell-shaped.
# Smaller half-life -> faster decline after the peak.
ENABLE_APPARENT_AAV_DECAY = True
BLOOD_AAV_HALF_LIFE_H = MOUSE_AAV9_BLOOD_HALF_LIFE_H
VASCULAR_AAV_HALF_LIFE_H = 9.0     # nonspecific loss from organ vascular space
ISF_AAV_HALF_LIFE_H = 15.0         # degradation / lymphatic loss from tissue ISF
LIVER_EXTRA_ISF_HALF_LIFE_H = 12.0  # liver uptake / RES-like loss from liver ISF
SPLEEN_EXTRA_ISF_HALF_LIFE_H = 10.0 # spleen uptake / RES-like loss from spleen ISF
PLOT_DECAY_WINDOW_H = 48.0

SAVE_FIGURES = True
SHOW_FIGURES = False
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "mouse_pbpk"

# New in the refined version:
# - "mechanistic" keeps the bell-shaped exposure but assigns clearance to
#   interpretable RES/endothelial/ISF mechanisms.
# - "half_life_demo" recovers the older apparent half-life behavior.
CLEARANCE_MODE = "mechanistic"  # allowed: "mechanistic" or "half_life_demo"

# Default scenario panel. These scenarios are intentionally simple and can be
# calibrated later against qPCR/ddPCR biodistribution data.
RUN_DESIGN_SCENARIOS = True
RUN_SPATIAL_PK_DEMO = True

# Mouse IV AAV9 early intact/extracellular-capsid priors. Blood is the 5.0 h
# unmodified AAV9 PET estimate from Seo et al.; the old 2.4 h number belongs to
# tetracysteine-modified AAV9-TC. Organ values are recomputed from Wang et al.
# Table S1 at import time. They are not episome or expression half-lives.
NORMAL_AAV9_CAPSID_HALF_LIFE_H = dict(MOUSE_AAV9_CAPSID_HALF_LIFE_H)
AAV9_PARAMETER_EVIDENCE = aav9_parameter_payload()


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
        "k_cns_escape": 0.015,
        "k_lys": 0.075,
        "k_kidney_lys": 0.075,
        "k_cns_lys": 0.075,
    },
    "cns_tropic": {
        "PS_brain": 0.006,
        "Kp_brain": 0.35,
        "Bmax_bbb": 1.2e8,
        "k_bbb_on": 8e-12,
        "k_bbb_trans": 0.12,
        "k_bbb_deg": 0.04,
        "Bmax_cns": 8e7,
        "k_cns_on": 1.2e-11,
        "k_cns_escape": 0.015,
        "k_res_liver": 0.014,
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
    "cns_biased": {
        "k_tx": 0.7,
        "k_kidney_tx": 0.7,
        "k_cns_tx": 2.6,
        "EC50_cns_tx": 1.2e5,
    },
}


# ---------------------------------------------------------------------
# Model parameters
# ---------------------------------------------------------------------
def make_params() -> Dict[str, float | str]:
    # Reference cardiac output for an unanesthetized 25 g mouse. The previous
    # 25 mL/h value was an effective exchange rate, not a physiological cardiac
    # output. Keep the historical 1.25 mL/h lung exchange rate by separating
    # physiology from the explicit reduced-order exchange multiplier.
    co = 840.0  # mL/h (14 mL/min; Brown et al. 1997 reference mean)
    q_scale = EFFECTIVE_LUNG_EXCHANGE_ML_H / co

    if CLEARANCE_MODE == "mechanistic":
        # Early capsid-PK calibration. Organ-specific apparent loss is split
        # between vascular/endothelial removal and tissue/ISF catabolism below.
        k_clear_blood = (
            np.log(2.0) / NORMAL_AAV9_CAPSID_HALF_LIFE_H["blood"]
            if ENABLE_APPARENT_AAV_DECAY else 0.0
        )
        k_clear_vascular = 0.0
        k_clear_isf = 0.0
        k_extra_liver = 0.0
        k_extra_spleen = 0.0
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
        # Blood clearance is represented by k_clear_blood to avoid counting
        # the calibrated AAV9 blood half-life twice.
        "CL_blood": 0.0,

        # Organ vascular and interstitial effective volumes, mL
        # Vascular and interstitial volumes use the 20 g mouse large-protein
        # PBPK system parameters in Shah & Betts, Table 6. These are
        # physiological spaces; they are not AAV-specific fitted values.
        "V_liver_v": 0.095,
        "V_liver_isf": 0.19,
        "V_spleen_v": 0.01,
        "V_spleen_isf": 0.02,
        "V_kidney_v": 0.03,
        "V_kidney_isf": 0.101,
        "V_heart_v": 0.007,
        "V_heart_isf": 0.019,
        "V_muscle_v": 0.15,
        "V_muscle_isf": 1.032,
        "V_lung_v": 0.0191,
        "V_lung_isf": 0.057,
        "V_brain_v": 0.035,
        "V_brain_isf": 0.08,
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
        "Q_brain": q_scale * 0.10 * co,
        "Q_rest": q_scale * 0.14 * co,

        # Vascular-to-interstitial permeability terms, mL/h-like
        "PS_liver": 0.20,
        "PS_spleen": 0.15,
        "PS_kidney": 0.08,
        "PS_heart": 0.05,
        "PS_muscle": 0.03,
        "PS_lung": 0.10,
        # Small passive BBB permeability. Receptor-mediated transcytosis is
        # represented separately in the CNS module below.
        "PS_brain": 0.00005,
        "PS_rest": 0.02,

        # Tissue partition coefficients
        "Kp_liver": 1.5,
        "Kp_spleen": 1.2,
        "Kp_kidney": 0.8,
        "Kp_heart": 0.6,
        "Kp_muscle": 0.5,
        "Kp_lung": 0.7,
        "Kp_brain": 0.04,
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
        "k_extra_isf_clear_brain": 0.0,
        "k_extra_isf_clear_rest": 0.0,

        # Organ-specific early capsid loss, 1/h. We assign 35% of the measured
        # apparent loss to vascular/endothelial removal and 65% to tissue/ISF
        # catabolism. This split is a modeling assumption; the summed rate is
        # ln(2)/organ half-life before exchange and cellular uptake are added.
        **{
            f"k_res_{organ}": 0.35 * np.log(2.0) / half_life
            for organ, half_life in NORMAL_AAV9_CAPSID_HALF_LIFE_H.items()
            if organ != "blood"
        },
        **{
            f"k_deg_isf_{organ}": 0.65 * np.log(2.0) / half_life
            for organ, half_life in NORMAL_AAV9_CAPSID_HALF_LIFE_H.items()
            if organ != "blood"
        },

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

        # -------------------------------------------------------------
        # CNS / blood-brain barrier module
        # -------------------------------------------------------------
        # BBB luminal binding and endothelial processing. Transcytosed capsid
        # enters A_brain_isf; recycled capsid returns to A_brain_v.
        "Bmax_bbb": 1e7,
        "k_bbb_on": 1e-13,          # mL/(vg*h)
        "k_bbb_off": 0.08,          # 1/h
        "k_bbb_int": 0.15,          # 1/h
        "k_bbb_trans": 0.003,       # 1/h, low baseline AAV9 BBB passage
        "k_bbb_recycle": 0.12,      # 1/h
        "k_bbb_deg": np.log(2) / NORMAL_AAV9_CAPSID_HALF_LIFE_H["brain"],

        # CNS parenchymal-cell uptake from brain interstitial space.
        "Bmax_cns": 3e7,
        "k_cns_on": 5e-12,          # mL/(vg*h)
        "k_cns_off": 0.05,          # 1/h
        "k_cns_int": 0.15,          # 1/h
        "k_cns_ee_le": 0.25,
        "k_cns_rec": 0.08,
        "k_cns_deg_ee": 0.02,
        "k_cns_escape": 0.008,
        "k_cns_lys": 0.10,
        "k_cns_nuc": 0.018,
        "k_cns_uncoat_cyto": 0.004,
        "k_cns_uncoat_nuc": 0.018,
        "k_cns_deg_cyto": 0.01,
        "k_cns_deg_ncap": 0.005,
        "k_cns_ds": 0.01,
        "k_cns_deg_ss": 0.02,
        "k_cns_epi": 0.01,
        "k_cns_deg_ds": 0.005,
        "k_cns_loss_epi": 0.012,
        "k_cns_tx": 1.5,
        "h_cns_tx": 1.2,
        "EC50_cns_tx": 2.0e5,
        "k_cns_tl": 4.0,
        "k_cns_deg_m": np.log(2) / 6.0,
        "k_cns_deg_p": np.log(2) / 48.0,
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
        A_brain_v, A_brain_isf,
        A_rest_v, A_rest_isf,
        B, EE, LE, CY, Ncap, Nss, Nds, Epi, M, P, Ab,
        K_filtrate, K_pt_lumen, K_bound_apical, K_bound_bsl,
        K_EE, K_REC, K_LE, K_LYS, K_CY, K_Ncap, K_Nss, K_Nds,
        K_Epi, K_M, K_P, K_Urine, K_Deg,
        C_BBB_bound, C_BBB_endo, C_bound, C_EE, C_LE, C_CY,
        C_Ncap, C_Nss, C_Nds, C_Epi, C_M, C_P, C_Deg,
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
        "brain": A_brain_v,
        "rest": A_rest_v,
    }
    A_isf = {
        "liver": A_liver_isf,
        "spleen": A_spleen_isf,
        "kidney": A_kidney_isf,
        "heart": A_heart_isf,
        "muscle": A_muscle_isf,
        "lung": A_lung_isf,
        "brain": A_brain_isf,
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
    tx = min(tx, float(p["k_tx"]))
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

    # -----------------------------------------------------------------
    # CNS delivery module: BBB handling -> brain ISF -> neural-cell fate
    # -----------------------------------------------------------------
    C_brain_v = max(A_brain_v, 0.0) / float(p["V_brain_v"])
    BBB_bound_eff = max(C_BBB_bound, 0.0)
    BBB_free = max(float(p["Bmax_bbb"]) - BBB_bound_eff, 0.0)
    J_bbb_bind = (
        float(p["k_bbb_on"]) * C_brain_v * BBB_free
        - float(p["k_bbb_off"]) * BBB_bound_eff
    )
    if C_BBB_bound <= 0.0 and J_bbb_bind < 0.0:
        J_bbb_bind = 0.0

    J_bbb_int = float(p["k_bbb_int"]) * max(C_BBB_bound, 0.0)
    J_bbb_trans = float(p["k_bbb_trans"]) * max(C_BBB_endo, 0.0)
    J_bbb_recycle = float(p["k_bbb_recycle"]) * max(C_BBB_endo, 0.0)
    J_bbb_deg = float(p["k_bbb_deg"]) * max(C_BBB_endo, 0.0)

    C_brain_isf = max(A_brain_isf, 0.0) / float(p["V_brain_isf"])
    C_bound_eff = max(C_bound, 0.0)
    C_free = max(float(p["Bmax_cns"]) - C_bound_eff, 0.0)
    J_cns_bind = (
        float(p["k_cns_on"]) * C_brain_isf * C_free
        - float(p["k_cns_off"]) * C_bound_eff
    )
    if C_bound <= 0.0 and J_cns_bind < 0.0:
        J_cns_bind = 0.0

    J_cns_int = float(p["k_cns_int"]) * max(C_bound, 0.0)
    J_cns_rec = float(p["k_cns_rec"]) * max(C_EE, 0.0)

    # Receptor recycling returns intact vector to the vascular/ISF pools.
    dA_v["brain"] += J_bbb_recycle - J_bbb_bind
    dA_isf["brain"] += J_bbb_trans + J_cns_rec - J_cns_bind

    dC_BBB_bound = J_bbb_bind - J_bbb_int
    dC_BBB_endo = J_bbb_int - J_bbb_trans - J_bbb_recycle - J_bbb_deg
    dC_bound = J_cns_bind - J_cns_int
    dC_EE = J_cns_int - (
        float(p["k_cns_ee_le"])
        + float(p["k_cns_rec"])
        + float(p["k_cns_deg_ee"])
    ) * C_EE
    dC_LE = float(p["k_cns_ee_le"]) * C_EE - (
        float(p["k_cns_escape"]) + float(p["k_cns_lys"])
    ) * C_LE
    dC_CY = float(p["k_cns_escape"]) * C_LE - (
        float(p["k_cns_nuc"])
        + float(p["k_cns_uncoat_cyto"])
        + float(p["k_cns_deg_cyto"])
    ) * C_CY
    dC_Ncap = float(p["k_cns_nuc"]) * C_CY - (
        float(p["k_cns_uncoat_nuc"]) + float(p["k_cns_deg_ncap"])
    ) * C_Ncap
    dC_Nss = (
        float(p["k_cns_uncoat_cyto"]) * C_CY
        + float(p["k_cns_uncoat_nuc"]) * C_Ncap
        - (float(p["k_cns_ds"]) + float(p["k_cns_deg_ss"])) * C_Nss
    )
    dC_Nds = float(p["k_cns_ds"]) * C_Nss - (
        float(p["k_cns_epi"]) + float(p["k_cns_deg_ds"])
    ) * C_Nds
    dC_Epi = float(p["k_cns_epi"]) * C_Nds - float(p["k_cns_loss_epi"]) * C_Epi

    C_Epi_eff = max(C_Epi, 0.0)
    h_cns = float(p["h_cns_tx"])
    tx_cns = (
        float(p["k_cns_tx"])
        * (C_Epi_eff ** h_cns)
        / (float(p["EC50_cns_tx"]) ** h_cns + C_Epi_eff ** h_cns + 1e-30)
    )
    tx_cns = min(tx_cns, float(p["k_cns_tx"]))
    dC_M = tx_cns - float(p["k_cns_deg_m"]) * C_M
    dC_P = float(p["k_cns_tl"]) * C_M - float(p["k_cns_deg_p"]) * C_P
    dC_Deg = (
        J_bbb_deg
        + float(p["k_cns_deg_ee"]) * max(C_EE, 0.0)
        + float(p["k_cns_lys"]) * max(C_LE, 0.0)
        + float(p["k_cns_deg_cyto"]) * max(C_CY, 0.0)
        + float(p["k_cns_deg_ncap"]) * max(C_Ncap, 0.0)
        + float(p["k_cns_deg_ss"]) * max(C_Nss, 0.0)
        + float(p["k_cns_deg_ds"]) * max(C_Nds, 0.0)
        + float(p["k_cns_loss_epi"]) * max(C_Epi, 0.0)
    )

    return [
        dA_blood,
        dA_v["liver"], dA_isf["liver"],
        dA_v["spleen"], dA_isf["spleen"],
        dA_v["kidney"], dA_isf["kidney"],
        dA_v["heart"], dA_isf["heart"],
        dA_v["muscle"], dA_isf["muscle"],
        dA_v["lung"], dA_isf["lung"],
        dA_v["brain"], dA_isf["brain"],
        dA_v["rest"], dA_isf["rest"],
        dB, dEE, dLE, dCY, dNcap, dNss, dNds, dEpi, dM, dP, dAb,
        dK_filtrate, dK_pt_lumen, dK_bound_apical, dK_bound_bsl,
        dK_EE, dK_REC, dK_LE, dK_LYS, dK_CY, dK_Ncap, dK_Nss, dK_Nds,
        dK_Epi, dK_M, dK_P, dK_Urine, dK_Deg,
        dC_BBB_bound, dC_BBB_endo, dC_bound, dC_EE, dC_LE, dC_CY,
        dC_Ncap, dC_Nss, dC_Nds, dC_Epi, dC_M, dC_P, dC_Deg,
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


def total_cns_vector_aav(sol: SimpleSolution) -> np.ndarray:
    """BBB and CNS vector states, excluding mRNA/protein and degradation sink."""
    states = [
        "C_BBB_bound", "C_BBB_endo", "C_bound", "C_EE", "C_LE", "C_CY",
        "C_Ncap", "C_Nss", "C_Nds", "C_Epi",
    ]
    total = np.zeros_like(sol.t)
    for state in states:
        total = total + sol.y[IDX[state]]
    return total


def total_accounted_aav(sol: SimpleSolution) -> np.ndarray:
    """All tracked vector plus cumulative sinks for mass-balance auditing."""
    total = (
        total_extracellular_aav(sol)
        + total_liver_vector_aav(sol)
        + total_kidney_vector_aav(sol)
        + total_cns_vector_aav(sol)
    )
    for state in [
        "K_Urine", "K_Deg", "C_Deg", "Loss_blood_clear", "Loss_vascular_res_clear",
        "Loss_isf_clear", "Loss_neutralized", "Loss_liver_cell",
    ]:
        total = total + sol.y[IDX[state]]
    return total


def mass_balance_error(sol: SimpleSolution) -> np.ndarray:
    delivered = np.maximum(sol.y[IDX["Dose_in"]], 1e-30)
    return (total_accounted_aav(sol) - sol.y[IDX["Dose_in"]]) / delivered


def plot_bell_shaped_aav_decay(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    """0-48 h extracellular AAV profile showing rise, peak, and decay."""
    t_h = sol_long.t
    mask = t_h <= PLOT_DECAY_WINDOW_H
    t_h = t_h[mask]

    plt.figure(figsize=(14, 6))

    # Concentration profiles: blood + representative tissue ISF compartments.
    plt.subplot(1, 2, 1)
    C_blood = concentration(sol_long, "A_blood", float(p["V_blood"]))[mask]
    plt.plot(t_h, C_blood, label="blood", color="black", linestyle="--", linewidth=2.2)

    representative_organs = ["liver", "spleen", "kidney", "muscle", "brain"]
    for organ in representative_organs:
        C_isf = concentration(sol_long, f"A_{organ}_isf", float(p[f"V_{organ}_isf"]))[mask]
        plt.plot(t_h, C_isf, label=f"{organ}_ISF", color=COLORS[organ], linewidth=1.9)

    plt.xlabel("Time (h)")
    plt.ylabel("AAV concentration (vg/mL)")
    plt.title("Extracellular AAV exposure, 0-48 h")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    # Total extracellular amount: confirms that the body-level AAV burden decays.
    plt.subplot(1, 2, 2)
    total = total_extracellular_aav(sol_long)[mask]
    plt.plot(t_h, total, color="black", linewidth=2.2)
    plt.xlabel("Time (h)")
    plt.ylabel("Total extracellular AAV amount (vg)")
    plt.title("Total extracellular AAV decays after dosing")
    plt.grid(True, linestyle="--", alpha=0.35)

    plt.tight_layout()
    save_or_show("03_bell_shaped_aav_decay_48h.png")


def plot_normal_aav9_organ_concentration_comparison(
    sol_long: SimpleSolution,
    p: Dict[str, float | str],
) -> None:
    """Compare early extracellular AAV9-like concentrations across all organs."""
    # Logarithmic time cannot include t=0; the first positive solver output is
    # retained so the infusion peak remains visible.
    mask = (sol_long.t > 0.0) & (sol_long.t <= 72.0)
    t_h = sol_long.t[mask]
    blood = concentration(sol_long, "A_blood", float(p["V_blood"]))[mask]
    organ_metrics = []

    plt.figure(figsize=(16, 10))
    plt.subplot(2, 2, 1)
    plt.plot(t_h, log_safe(blood), color="black", linestyle="--", linewidth=2.2, label="blood (2.4 h)")
    for organ in ORGANS:
        amount = sol_long.y[IDX[f"A_{organ}_v"]] + sol_long.y[IDX[f"A_{organ}_isf"]]
        isf_amount = sol_long.y[IDX[f"A_{organ}_isf"]]
        tissue_volume = float(p[f"V_{organ}_v"]) + float(p[f"V_{organ}_isf"])
        tissue_conc = amount / tissue_volume
        isf_conc = isf_amount / float(p[f"V_{organ}_isf"])
        half_life = NORMAL_AAV9_CAPSID_HALF_LIFE_H[organ]
        plt.plot(
            t_h,
            log_safe(tissue_conc[mask]),
            color=COLORS[organ],
            linewidth=1.9,
            label=f"{organ} ({half_life:g} h)",
        )

        auc_isf_amount = auc_trapz(isf_amount[mask], t_h)
        organ_metrics.append({
            "organ": organ,
            "half_life_h": half_life,
            "peak_isf_conc_vg_ml": float(np.nanmax(isf_conc[mask])),
            "peak_isf_time_h": float(t_h[int(np.nanargmax(isf_conc[mask]))]),
            "auc_isf_conc_vg_h_ml": auc_trapz(isf_conc[mask], t_h),
            "auc_isf_amount_vg_h": auc_isf_amount,
        })

    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Time (h)")
    plt.ylabel("Extracellular capsid concentration (vg/mL)")
    plt.title("Total extracellular concentration (includes vascular blood), 0-72 h")
    plt.grid(True, which="both", linestyle="--", alpha=0.3)
    plt.legend(fontsize=8, ncol=2)

    plt.subplot(2, 2, 2)
    for organ in ORGANS:
        isf_conc = concentration(sol_long, f"A_{organ}_isf", float(p[f"V_{organ}_isf"]))[mask]
        peak_idx = int(np.nanargmax(isf_conc))
        plt.plot(
            t_h,
            log_safe(isf_conc),
            color=COLORS[organ],
            linewidth=1.9,
            label=f"{organ} (Tmax {t_h[peak_idx]:.2g} h)",
        )
        plt.scatter(t_h[peak_idx], isf_conc[peak_idx], color=COLORS[organ], s=28, zorder=3)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Time (h)")
    plt.ylabel("Interstitial concentration (vg/mL)")
    plt.title("Interstitial exposure: CNS versus peripheral organs")
    plt.grid(True, which="both", linestyle="--", alpha=0.3)
    plt.legend(fontsize=8, ncol=2)

    names = [row["organ"] for row in organ_metrics]
    peak_concs = [row["peak_isf_conc_vg_ml"] for row in organ_metrics]
    colors = [COLORS[organ] for organ in names]

    plt.subplot(2, 2, 3)
    plt.bar(names, peak_concs, color=colors)
    plt.yscale("log")
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Peak concentration (vg/mL)")
    plt.title("Peak interstitial concentration (post-barrier)")
    plt.grid(True, axis="y", which="both", linestyle="--", alpha=0.3)

    total_auc_amount = sum(row["auc_isf_amount_vg_h"] for row in organ_metrics)
    shares = [100.0 * row["auc_isf_amount_vg_h"] / max(total_auc_amount, 1e-30) for row in organ_metrics]
    plt.subplot(2, 2, 4)
    bars = plt.bar(names, shares, color=colors)
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Share of 0-72 h interstitial exposure (%)")
    plt.title("Post-barrier organ exposure share (ISF amount AUC)")
    plt.grid(True, axis="y", linestyle="--", alpha=0.3)
    for bar, value in zip(bars, shares):
        label = f"{value:.2f}%" if value < 1.0 else f"{value:.1f}%"
        plt.text(bar.get_x() + bar.get_width() / 2.0, value, label, ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    save_or_show("13_normal_aav9_organ_concentration_comparison_log_axes.png")

    if SAVE_FIGURES:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        metrics_path = OUTPUT_DIR / "13_normal_aav9_organ_distribution_metrics.csv"
        header = [
            "organ", "half_life_h", "peak_isf_conc_vg_ml", "peak_isf_time_h",
            "auc_isf_conc_vg_h_ml", "auc_isf_amount_vg_h", "auc_isf_amount_share_pct",
        ]
        with metrics_path.open("w", encoding="utf-8") as f:
            f.write(",".join(header) + "\n")
            for row, share in zip(organ_metrics, shares):
                values = [row[key] for key in header[:-1]] + [share]
                f.write(",".join(str(value) for value in values) + "\n")


def plot_normal_aav9_organ_concentration_comparison_linear(
    sol_long: SimpleSolution,
    p: Dict[str, float | str],
) -> None:
    """Show the same AAV9-like comparison with ordinary linear axes."""
    mask = sol_long.t <= 72.0
    t_h = sol_long.t[mask]
    blood = concentration(sol_long, "A_blood", float(p["V_blood"]))[mask]
    organ_metrics = []

    plt.figure(figsize=(16, 10))
    plt.subplot(2, 2, 1)
    plt.plot(t_h, blood, color="black", linestyle="--", linewidth=2.2, label="blood")
    for organ in ORGANS:
        amount = sol_long.y[IDX[f"A_{organ}_v"]] + sol_long.y[IDX[f"A_{organ}_isf"]]
        isf_amount = sol_long.y[IDX[f"A_{organ}_isf"]]
        tissue_volume = float(p[f"V_{organ}_v"]) + float(p[f"V_{organ}_isf"])
        tissue_conc = amount / tissue_volume
        isf_conc = isf_amount / float(p[f"V_{organ}_isf"])
        plt.plot(t_h, tissue_conc[mask], color=COLORS[organ], linewidth=1.9, label=organ)
        organ_metrics.append({
            "organ": organ,
            "peak_isf_conc_vg_ml": float(np.nanmax(isf_conc[mask])),
            "peak_isf_time_h": float(t_h[int(np.nanargmax(isf_conc[mask]))]),
            "auc_isf_amount_vg_h": auc_trapz(isf_amount[mask], t_h),
        })

    plt.xlabel("Time (h)")
    plt.ylabel("Extracellular capsid concentration (vg/mL)")
    plt.title("Total extracellular concentration (linear axes)")
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.legend(fontsize=8, ncol=2)

    plt.subplot(2, 2, 2)
    for organ in ORGANS:
        isf_conc = concentration(sol_long, f"A_{organ}_isf", float(p[f"V_{organ}_isf"]))[mask]
        peak_idx = int(np.nanargmax(isf_conc))
        plt.plot(
            t_h,
            isf_conc,
            color=COLORS[organ],
            linewidth=1.9,
            label=f"{organ} (Tmax {t_h[peak_idx]:.2g} h)",
        )
        plt.scatter(t_h[peak_idx], isf_conc[peak_idx], color=COLORS[organ], s=28, zorder=3)
    plt.xlabel("Time (h)")
    plt.ylabel("Interstitial concentration (vg/mL)")
    plt.title("Interstitial exposure (linear axes)")
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.legend(fontsize=8, ncol=2)

    names = [row["organ"] for row in organ_metrics]
    colors = [COLORS[organ] for organ in names]
    plt.subplot(2, 2, 3)
    plt.bar(names, [row["peak_isf_conc_vg_ml"] for row in organ_metrics], color=colors)
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Peak concentration (vg/mL)")
    plt.title("Peak interstitial concentration (linear axes)")
    plt.grid(True, axis="y", linestyle="--", alpha=0.3)

    total_auc_amount = sum(row["auc_isf_amount_vg_h"] for row in organ_metrics)
    shares = [100.0 * row["auc_isf_amount_vg_h"] / max(total_auc_amount, 1e-30) for row in organ_metrics]
    plt.subplot(2, 2, 4)
    bars = plt.bar(names, shares, color=colors)
    plt.xticks(rotation=35, ha="right")
    plt.ylabel("Share of 0-72 h interstitial exposure (%)")
    plt.title("Post-barrier exposure share (linear axes)")
    plt.grid(True, axis="y", linestyle="--", alpha=0.3)
    for bar, value in zip(bars, shares):
        label = f"{value:.2f}%" if value < 1.0 else f"{value:.1f}%"
        plt.text(bar.get_x() + bar.get_width() / 2.0, value, label, ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    save_or_show("14_normal_aav9_organ_concentration_comparison_linear_axes.png")


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
    #plt.yscale("log")
    plt.xlabel("Time (h)")
    plt.ylabel("AAV concentration (vg/mL)")
    plt.title("Vascular AAV concentrations, 0-2 h")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.subplot(1, 2, 2)
    for organ in ORGANS:
        C_isf = concentration(sol_short, f"A_{organ}_isf", float(p[f"V_{organ}_isf"]))
        plt.plot(sol_short.t, log_safe(C_isf), label=f"{organ}_ISF", color=COLORS[organ], linewidth=1.8)
    #plt.yscale("log")
    plt.xlabel("Time (h)")
    plt.ylabel("AAV concentration (vg/mL)")
    plt.title("Interstitial AAV concentrations, 0-2 h")
    plt.grid(True, which="both", linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

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


def plot_liver_intracellular_uptake(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    """Diagnostic plot for the liver uptake / intracellular trafficking chain."""
    t_day = sol_long.t / 24.0
    fig = plt.figure(figsize=(16, 10))
    grid = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.32)

    # Resolve the rapid post-dose peak instead of compressing it into 56 days.
    ax_entry = fig.add_subplot(grid[0, 0])
    focus_mask = sol_long.t <= 48.0
    t_focus = sol_long.t[focus_mask]
    liver_isf_focus = log_safe(sol_long.y[IDX["A_liver_isf"]][focus_mask])
    bound_focus = log_safe(sol_long.y[IDX["B"]][focus_mask])
    isf_peak = int(np.nanargmax(liver_isf_focus))

    isf_line = ax_entry.plot(t_focus, liver_isf_focus, color="tab:blue", label="Liver ISF", linewidth=2.2)[0]
    ax_entry.scatter(t_focus[isf_peak], liver_isf_focus[isf_peak], color="tab:blue", s=28, zorder=4)
    ax_entry.annotate(
        f"ISF peak: {t_focus[isf_peak]:.2f} h",
        (t_focus[isf_peak], liver_isf_focus[isf_peak]),
        xytext=(42, -42), textcoords="offset points", fontsize=8, color="tab:blue",
        arrowprops={"arrowstyle": "->", "color": "tab:blue", "lw": 0.8},
    )
    ax_entry.set_xlim(0.0, 48.0)
    ax_entry.set_xlabel("Time after dose (h)")
    ax_entry.set_ylabel("Liver ISF AAV (vg-equivalent)", color="tab:blue")
    ax_entry.tick_params(axis="y", labelcolor="tab:blue")
    ax_entry.margins(y=0.14)
    ax_entry.grid(True, linestyle="--", alpha=0.35)

    ax_bound = ax_entry.twinx()
    bound_line = ax_bound.plot(t_focus, bound_focus, color="tab:orange", label="Surface-bound B", linewidth=2.0)[0]
    ax_bound.set_ylabel("Receptor-bound AAV (vg-equivalent)", color="tab:orange")
    ax_bound.tick_params(axis="y", labelcolor="tab:orange")
    ax_bound.margins(y=0.16)
    ax_bound.text(0.97, 0.78, "Receptor-capacity plateau", transform=ax_bound.transAxes,
                  ha="right", va="top", fontsize=8, color="tab:orange")
    ax_entry.set_title("Early liver ISF peak and receptor binding (0–48 h)")
    ax_entry.legend([isf_line, bound_line], [isf_line.get_label(), bound_line.get_label()], loc="upper right", fontsize=8)

    ax_traffic = fig.add_subplot(grid[0, 1])
    for state in ["EE", "LE", "CY", "Ncap"]:
        ax_traffic.plot(t_day, log_safe(sol_long.y[IDX[state]]), label=state, linewidth=1.9)
    ax_traffic.set_xlabel("Time (day)")
    ax_traffic.set_ylabel("AAV amount (vg-equivalent)")
    ax_traffic.set_title("Liver intracellular trafficking")
    ax_traffic.grid(True, linestyle="--", alpha=0.35)
    ax_traffic.legend(fontsize=8)

    ax_genome = fig.add_subplot(grid[1, 0])
    for state in ["Nss", "Nds", "Epi"]:
        ax_genome.plot(t_day, log_safe(sol_long.y[IDX[state]]), label=state, linewidth=1.9)
    ax_genome.set_xlabel("Time (day)")
    ax_genome.set_ylabel("Vector genome state (a.u.)")
    ax_genome.set_title("Liver nuclear genome processing")
    ax_genome.grid(True, linestyle="--", alpha=0.35)
    ax_genome.legend(fontsize=8)

    expression_grid = grid[1, 1].subgridspec(2, 1, hspace=0.38)
    ax_mrna = fig.add_subplot(expression_grid[0, 0])
    ax_mrna.plot(t_day, sol_long.y[IDX["M"]], color="tab:blue", linewidth=2.0)
    ax_mrna.set_ylabel("mRNA (a.u.)")
    ax_mrna.set_title("Liver mRNA output")
    ax_mrna.grid(True, linestyle="--", alpha=0.35)
    ax_mrna.tick_params(axis="x", labelbottom=False)

    ax_protein = fig.add_subplot(expression_grid[1, 0], sharex=ax_mrna)
    ax_protein.plot(t_day, sol_long.y[IDX["P"]], color="tab:orange", linewidth=2.0)
    ax_protein.set_xlabel("Time (day)")
    ax_protein.set_ylabel("Protein (a.u.)")
    ax_protein.set_title("Liver protein output")
    ax_protein.grid(True, linestyle="--", alpha=0.35)

    fig.subplots_adjust(left=0.07, right=0.93, top=0.95, bottom=0.08)
    save_or_show("03_liver_intracellular_56d.png")


def plot_kidney_module(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    """Plot multilevel kidney proximal-tubule uptake and expression."""
    t_day = sol_long.t / 24.0
    fig = plt.figure(figsize=(16, 10))
    grid = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.32)

    # Focus on the early renal peak and separate large luminal pools from
    # smaller receptor-bound pools so both remain visually resolvable.
    ax_entry = fig.add_subplot(grid[0, 0])
    focus_mask = sol_long.t <= 96.0
    t_focus = sol_long.t[focus_mask]
    filtrate = log_safe(sol_long.y[IDX["K_filtrate"]][focus_mask])
    lumen = log_safe(sol_long.y[IDX["K_pt_lumen"]][focus_mask])
    bound_apical = log_safe(sol_long.y[IDX["K_bound_apical"]][focus_mask])
    bound_bsl = log_safe(sol_long.y[IDX["K_bound_bsl"]][focus_mask])

    line_filtrate = ax_entry.plot(t_focus, filtrate, label="K_filtrate", color="tab:blue", linewidth=1.9)[0]
    line_lumen = ax_entry.plot(t_focus, lumen, label="K_pt_lumen", color="tab:orange", linewidth=2.1)[0]
    lumen_peak = int(np.nanargmax(lumen))
    ax_entry.scatter(t_focus[lumen_peak], lumen[lumen_peak], color="tab:orange", s=28, zorder=4)
    ax_entry.annotate(
        f"Luminal peak: {t_focus[lumen_peak]:.2f} h",
        (t_focus[lumen_peak], lumen[lumen_peak]),
        xytext=(34, -38), textcoords="offset points", fontsize=8, color="tab:orange",
        arrowprops={"arrowstyle": "->", "color": "tab:orange", "lw": 0.8},
    )
    ax_entry.set_xlim(0.0, 96.0)
    ax_entry.set_xlabel("Time after dose (h)")
    ax_entry.set_ylabel("Filtrate / tubular-lumen AAV", color="tab:blue")
    ax_entry.tick_params(axis="y", labelcolor="tab:blue")
    ax_entry.margins(y=0.14)
    ax_entry.grid(True, linestyle="--", alpha=0.35)

    ax_bound = ax_entry.twinx()
    line_apical = ax_bound.plot(t_focus, bound_apical, label="K_bound_apical", color="tab:green", linewidth=1.9)[0]
    line_bsl = ax_bound.plot(t_focus, bound_bsl, label="K_bound_bsl", color="tab:red", linewidth=1.9)[0]
    bsl_peak = int(np.nanargmax(bound_bsl))
    ax_bound.scatter(t_focus[bsl_peak], bound_bsl[bsl_peak], color="tab:red", s=28, zorder=4)
    ax_bound.annotate(
        f"Basolateral-binding peak: {t_focus[bsl_peak]:.2f} h",
        (t_focus[bsl_peak], bound_bsl[bsl_peak]),
        xytext=(38, -18), textcoords="offset points", fontsize=8, color="tab:red",
        arrowprops={"arrowstyle": "->", "color": "tab:red", "lw": 0.8},
    )
    ax_bound.set_ylabel("Receptor-bound AAV", color="tab:red")
    ax_bound.tick_params(axis="y", labelcolor="tab:red")
    ax_bound.margins(y=0.14)
    ax_entry.set_title("Early kidney exposure and receptor binding (0–96 h)")
    entry_lines = [line_filtrate, line_lumen, line_apical, line_bsl]
    ax_entry.legend(entry_lines, [line.get_label() for line in entry_lines], loc="upper right", fontsize=8)

    ax_traffic = fig.add_subplot(grid[0, 1])
    intracellular_states = ["K_EE", "K_REC", "K_LE", "K_LYS", "K_CY", "K_Ncap"]
    for state in intracellular_states:
        ax_traffic.plot(sol_long.t, log_safe(sol_long.y[IDX[state]]), label=state, linewidth=1.8)
    ax_traffic.set_xlim(0.0, 240.0)
    ax_traffic.set_xlabel("Time after dose (h)")
    ax_traffic.set_ylabel("AAV amount (vg-equivalent)")
    ax_traffic.set_title("Kidney proximal-tubule intracellular trafficking (0–10 d)")
    ax_traffic.grid(True, linestyle="--", alpha=0.35)
    ax_traffic.legend(fontsize=8)

    ax_genome = fig.add_subplot(grid[1, 0])
    genome_states = ["K_Nss", "K_Nds", "K_Epi"]
    for state in genome_states:
        ax_genome.plot(t_day, log_safe(sol_long.y[IDX[state]]), label=state, linewidth=1.8)
    ax_genome.set_xlabel("Time (day)")
    ax_genome.set_ylabel("Vector genome state (a.u.)")
    ax_genome.set_title("Kidney nuclear genome processing")
    ax_genome.grid(True, linestyle="--", alpha=0.35)
    ax_genome.legend(fontsize=8)

    expression_grid = grid[1, 1].subgridspec(2, 1, hspace=0.38)
    ax_mrna = fig.add_subplot(expression_grid[0, 0])
    ax_mrna.plot(t_day, sol_long.y[IDX["K_M"]], color="tab:blue", linewidth=2.0)
    ax_mrna.set_ylabel("Kidney mRNA (a.u.)")
    ax_mrna.set_title("Kidney mRNA output")
    ax_mrna.grid(True, linestyle="--", alpha=0.35)
    ax_mrna.tick_params(axis="x", labelbottom=False)

    ax_protein = fig.add_subplot(expression_grid[1, 0], sharex=ax_mrna)
    ax_protein.plot(t_day, sol_long.y[IDX["K_P"]], color="tab:orange", linewidth=2.0)
    ax_protein.set_xlabel("Time (day)")
    ax_protein.set_ylabel("Kidney protein (a.u.)")
    ax_protein.set_title("Kidney protein output")
    ax_protein.grid(True, linestyle="--", alpha=0.35)

    fig.subplots_adjust(left=0.07, right=0.93, top=0.95, bottom=0.08)
    save_or_show("06_kidney_multilevel_module.png")


def plot_cns_module(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    """Plot BBB transfer, CNS-cell trafficking, genome processing, and expression."""
    t_day = sol_long.t / 24.0
    fig = plt.figure(figsize=(16, 10))
    grid = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.28)

    # Vascular AAV is orders of magnitude larger than BBB/ISF pools. Separate
    # the two scales so post-BBB delivery is not flattened against zero.
    bbb_grid = grid[0, 0].subgridspec(2, 1, hspace=0.38)
    ax_vascular = fig.add_subplot(bbb_grid[0, 0])
    vascular_mask = sol_long.t <= 48.0
    vascular_t = sol_long.t[vascular_mask]
    vascular = log_safe(sol_long.y[IDX["A_brain_v"]][vascular_mask])
    vascular_peak = int(np.nanargmax(vascular))
    ax_vascular.plot(vascular_t, vascular, color="tab:blue", linewidth=2.1)
    ax_vascular.scatter(vascular_t[vascular_peak], vascular[vascular_peak], color="tab:blue", s=26, zorder=4)
    ax_vascular.annotate(
        f"Vascular peak: {vascular_t[vascular_peak]:.2f} h",
        (vascular_t[vascular_peak], vascular[vascular_peak]),
        xytext=(42, -36), textcoords="offset points", fontsize=8, color="tab:blue",
        arrowprops={"arrowstyle": "->", "color": "tab:blue", "lw": 0.8},
    )
    ax_vascular.set_xlim(0.0, 48.0)
    ax_vascular.set_ylabel("Brain vascular AAV")
    ax_vascular.set_title("Early brain vascular exposure (0–48 h)")
    ax_vascular.grid(True, linestyle="--", alpha=0.35)
    ax_vascular.tick_params(axis="x", labelbottom=False)
    ax_vascular.margins(y=0.14)

    ax_bbb = fig.add_subplot(bbb_grid[1, 0])
    for state, color in [("C_BBB_bound", "tab:orange"), ("C_BBB_endo", "tab:green"), ("A_brain_isf", "tab:red")]:
        ax_bbb.plot(sol_long.t, log_safe(sol_long.y[IDX[state]]), label=state, color=color, linewidth=1.9)
    ax_bbb.set_yscale("log")
    ax_bbb.set_xlim(0.0, 240.0)
    ax_bbb.set_xlabel("Time after dose (h)")
    ax_bbb.set_ylabel("BBB / brain ISF AAV (log)")
    ax_bbb.set_title("BBB processing and post-BBB delivery (0–10 d)")
    ax_bbb.grid(True, which="both", linestyle="--", alpha=0.35)
    ax_bbb.legend(fontsize=8, ncol=3, loc="upper right")

    ax_traffic = fig.add_subplot(grid[0, 1])
    for state in ["C_bound", "C_EE", "C_LE", "C_CY", "C_Ncap"]:
        ax_traffic.plot(sol_long.t, log_safe(sol_long.y[IDX[state]]), label=state, linewidth=1.8)
    ax_traffic.set_xlabel("Time after dose (h)")
    ax_traffic.set_ylabel("AAV amount (vg-equivalent)")
    ax_traffic.set_title("CNS-cell uptake and intracellular trafficking (0–10 d)")
    ax_traffic.set_xlim(0.0, 240.0)
    ax_traffic.grid(True, linestyle="--", alpha=0.35)
    ax_traffic.legend(fontsize=8)

    ax_genome = fig.add_subplot(grid[1, 0])
    for state in ["C_Nss", "C_Nds", "C_Epi"]:
        ax_genome.plot(t_day, log_safe(sol_long.y[IDX[state]]), label=state, linewidth=1.9)
    ax_genome.set_xlabel("Time (day)")
    ax_genome.set_ylabel("Vector genome state (a.u.)")
    ax_genome.set_title("CNS nuclear genome processing")
    ax_genome.grid(True, linestyle="--", alpha=0.35)
    ax_genome.legend(fontsize=8)

    expression_grid = grid[1, 1].subgridspec(2, 1, hspace=0.38)
    ax_mrna = fig.add_subplot(expression_grid[0, 0])
    ax_mrna.plot(t_day, sol_long.y[IDX["C_M"]], color="tab:blue", linewidth=2.0)
    ax_mrna.set_ylabel("CNS mRNA (a.u.)")
    ax_mrna.set_title("CNS mRNA output")
    ax_mrna.grid(True, linestyle="--", alpha=0.35)
    ax_mrna.tick_params(axis="x", labelbottom=False)

    ax_protein = fig.add_subplot(expression_grid[1, 0], sharex=ax_mrna)
    ax_protein.plot(t_day, sol_long.y[IDX["C_P"]], color="tab:orange", linewidth=2.0)
    ax_protein.set_xlabel("Time (day)")
    ax_protein.set_ylabel("CNS protein (a.u.)")
    ax_protein.set_title("CNS protein output")
    ax_protein.grid(True, linestyle="--", alpha=0.35)

    fig.subplots_adjust(left=0.07, right=0.96, top=0.95, bottom=0.08)
    save_or_show("11_cns_bbb_and_transduction.png")


def plot_cns_scenarios(base_p: Dict[str, float | str], t_eval_long: np.ndarray) -> None:
    """Compare baseline and CNS-directed capsid/promoter hypotheses."""
    scenarios = [
        ("baseline", "baseline_AAV", "ubiquitous"),
        ("BBB/capsid enhanced", "cns_tropic", "ubiquitous"),
        ("CNS capsid + promoter", "cns_tropic", "cns_biased"),
        ("escape + CNS promoter", "endosomal_escape_enhanced", "cns_biased"),
    ]

    rows = []
    plt.figure(figsize=(16, 5))
    for label, capsid, promoter in scenarios:
        p = apply_design_preset(base_p, capsid=capsid, promoter=promoter)
        sol = solve_model(t_eval_long, make_initial_condition(p), p, post_infusion_max_step=1.0)
        t_day = sol.t / 24.0
        C_brain_isf = concentration(sol, "A_brain_isf", float(p["V_brain_isf"]))

        early_mask = sol.t <= 48.0
        if label in {"baseline", "BBB/capsid enhanced"}:
            plt.subplot(1, 3, 1)
            plt.plot(sol.t[early_mask], log_safe(C_brain_isf[early_mask]), label=label, linewidth=1.9)
        if label in {"baseline", "BBB/capsid enhanced", "escape + CNS promoter"}:
            plt.subplot(1, 3, 2)
            epi_label = "escape enhanced" if label == "escape + CNS promoter" else label
            plt.plot(t_day, log_safe(sol.y[IDX["C_Epi"]]), label=epi_label, linewidth=1.9)
        plt.subplot(1, 3, 3)
        plt.plot(t_day, sol.y[IDX["C_P"]], label=label, linewidth=1.9)

        rows.append({
            "scenario": label,
            "capsid": capsid,
            "promoter": promoter,
            "auc_brain_isf": auc_trapz(C_brain_isf, sol.t),
            "peak_bbb_endo": float(np.nanmax(sol.y[IDX["C_BBB_endo"]])),
            "peak_cns_epi": float(np.nanmax(sol.y[IDX["C_Epi"]])),
            "peak_cns_protein": float(np.nanmax(sol.y[IDX["C_P"]])),
            "final_cns_deg": float(sol.y[IDX["C_Deg"], -1]),
            "final_mass_balance_error": float(mass_balance_error(sol)[-1]),
        })

    plt.subplot(1, 3, 1)
    plt.xlabel("Time (h)")
    plt.ylabel("Brain ISF concentration (vg/mL)")
    plt.title("BBB delivery into brain ISF")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.subplot(1, 3, 2)
    plt.xlabel("Time (day)")
    plt.ylabel("CNS episome (a.u.)")
    plt.title("Expression-competent CNS genome")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.subplot(1, 3, 3)
    plt.xlabel("Time (day)")
    plt.ylabel("CNS protein (a.u.)")
    plt.title("CNS transgene expression")
    plt.grid(True, linestyle="--", alpha=0.35)

    plt.tight_layout()
    save_or_show("12_cns_design_scenario_comparison.png")

    if SAVE_FIGURES:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        metrics_path = OUTPUT_DIR / "12_cns_design_scenario_metrics.csv"
        header = list(rows[0].keys())
        with metrics_path.open("w", encoding="utf-8") as f:
            f.write(",".join(header) + "\n")
            for row in rows:
                f.write(",".join(str(row[key]) for key in header) + "\n")


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
    plt.plot(t_day, total_cns_vector_aav(sol_long), label="CNS vector states", linewidth=2.0)
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
    """Compare capsid/promoter design hypotheses without editing core code."""
    scenarios = [
        ("baseline", "baseline_AAV", "ubiquitous"),
        ("liver_detargeted", "liver_detargeted", "ubiquitous"),
        ("kidney_tropic", "kidney_tropic", "kidney_biased"),
        ("escape_enhanced", "endosomal_escape_enhanced", "ubiquitous"),
    ]

    rows = []
    plt.figure(figsize=(15, 5))
    for label, capsid, promoter in scenarios:
        p = apply_design_preset(base_p, capsid=capsid, promoter=promoter)
        sol = solve_model(t_eval_long, make_initial_condition(p), p, post_infusion_max_step=1.0)
        t_day = sol.t / 24.0
        plt.subplot(1, 2, 1)
        plt.plot(t_day, log_safe(sol.y[IDX["Epi"]]), label=label, linewidth=2.0)
        plt.subplot(1, 2, 2)
        plt.plot(t_day, log_safe(sol.y[IDX["K_Epi"]]), label=label, linewidth=2.0)

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
            "final_mass_balance_error": float(mass_balance_error(sol)[-1]),
        })

    plt.subplot(1, 2, 1)
    plt.xlabel("Time (day)")
    plt.ylabel("Liver episome (a.u.)")
    plt.title("Capsid/promoter scenario: liver")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.subplot(1, 2, 2)
    plt.xlabel("Time (day)")
    plt.ylabel("Kidney episome (a.u.)")
    plt.title("Capsid/promoter scenario: kidney")
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend(fontsize=8)

    plt.tight_layout()
    save_or_show("09_design_scenario_comparison.png")

    if SAVE_FIGURES:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        metrics_path = OUTPUT_DIR / "09_design_scenario_metrics.csv"
        header = list(rows[0].keys())
        with metrics_path.open("w", encoding="utf-8") as f:
            f.write(",".join(header) + "\n")
            for row in rows:
                f.write(",".join(str(row[key]) for key in header) + "\n")


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


def print_metrics(sol_long: SimpleSolution, p: Dict[str, float | str]) -> None:
    C_blood_long = concentration(sol_long, "A_blood", float(p["V_blood"]))
    print("----- Simulation settings -----")
    print(f"Administration: {p['administration']}")
    print(f"Dose: {p['dose_vg']:.3e} vg")
    print(f"Infusion duration: {float(p['T_inf_h']) * 60:.2f} min")
    print(f"Q_scale: {p['Q_scale']}")
    print("Reference: normal IV AAV9-like early intact/labeled capsid PK")
    print(f"Blood capsid half-life prior: {BLOOD_AAV_HALF_LIFE_H:.2f} h")
    print("Organ capsid half-life priors (h):")
    for organ in ORGANS:
        print(f"  {organ}: {NORMAL_AAV9_CAPSID_HALF_LIFE_H[organ]:.1f}")

    print("\n----- Summary metrics -----")
    total_aav = total_extracellular_aav(sol_long)
    print("AUC_blood:", auc_trapz(C_blood_long, sol_long.t))
    print("Cmax_blood:", np.nanmax(C_blood_long))
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

    C_brain_isf = concentration(sol_long, "A_brain_isf", float(p["V_brain_isf"]))
    auc_brain_isf = auc_trapz(C_brain_isf, sol_long.t)
    print("\n----- CNS / BBB module metrics -----")
    print("AUC_brain_ISF:", auc_brain_isf)
    print("Brain/Liver ISF AUC ratio:", auc_brain_isf / max(auc_liver_isf, 1e-30))
    print("Peak BBB-bound AAV:", np.nanmax(sol_long.y[IDX["C_BBB_bound"]]))
    print("Peak BBB-endosomal AAV:", np.nanmax(sol_long.y[IDX["C_BBB_endo"]]))
    print("Peak CNS cell-bound AAV:", np.nanmax(sol_long.y[IDX["C_bound"]]))
    print("Peak CNS episome:", np.nanmax(sol_long.y[IDX["C_Epi"]]))
    print("Peak CNS mRNA:", np.nanmax(sol_long.y[IDX["C_M"]]))
    print("Peak CNS protein:", np.nanmax(sol_long.y[IDX["C_P"]]))
    print("Cumulative CNS/BBB degradation:", sol_long.y[IDX["C_Deg"], -1])
    print("CNS entry efficiency (Peak C_Epi / AUC brain ISF):", np.nanmax(sol_long.y[IDX["C_Epi"]]) / max(auc_brain_isf, 1e-30))
    print("Final mass-balance error:", mass_balance_error(sol_long)[-1])
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
    plot_bell_shaped_aav_decay(sol_long, p)
    plot_normal_aav9_organ_concentration_comparison(sol_long, p)
    plot_normal_aav9_organ_concentration_comparison_linear(sol_long, p)
    plot_long_states(sol_long, p)
    plot_liver_intracellular_uptake(sol_long, p)
    plot_kidney_module(sol_long, p)
    plot_cns_module(sol_long, p)
    plot_liver_vs_kidney_expression(sol_long, p)
    plot_mass_balance(sol_long, p)
    if RUN_DESIGN_SCENARIOS:
        plot_design_scenarios(p, t_eval_long)
        plot_cns_scenarios(p, t_eval_long)
    if RUN_SPATIAL_PK_DEMO:
        plot_spatial_pk_demo(sol_long, p)
    print_metrics(sol_long, p)

    if SAVE_FIGURES:
        print(f"\nFigures saved to: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
