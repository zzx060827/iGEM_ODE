"""Batch PBPK-to-PD export for the SINEUP Delivery Atlas.

Early biodistribution and liver/kidney/CNS intracellular transduction are read
directly from ode1.0.py. A longitudinal ODE then propagates expression-competent
episome -> SINEUP RNA -> endogenous target protein for every design point.
"""

from __future__ import annotations

import argparse
import csv
import json
import runpy
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

import human_spatial_pbpk as human_model


TROPISM_LITERATURE_PATH = Path(__file__).with_name("data") / "aav_capsid_tropism_literature.csv"

CAPSID_PRIORS = {
    "aav2": {
        "label": "AAV2", "evidence": "medium", "species": "preclinical + clinical ocular",
        "persistence_factor": 0.95,
        "tropism": {"liver": 0.30, "spleen": 0.80, "kidney": 0.45, "heart": 0.28, "muscle": 0.22, "lung": 0.80, "brain": 0.16, "rest": 0.35},
        "source": "https://pubmed.ncbi.nlm.nih.gov/18414476/",
        "additional_sources": ["https://pubmed.ncbi.nlm.nih.gov/40337478/"],
    },
    "aav5": {
        "label": "AAV5", "evidence": "medium", "species": "preclinical",
        "persistence_factor": 0.90,
        "tropism": {"liver": 0.55, "spleen": 0.85, "kidney": 0.48, "heart": 0.25, "muscle": 0.22, "lung": 2.20, "brain": 0.18, "rest": 0.45},
        "source": "https://pubmed.ncbi.nlm.nih.gov/39863928/",
        "additional_sources": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC11919325/"],
    },
    "aav8": {
        "label": "AAV8", "evidence": "strong", "species": "mouse + NHP",
        "persistence_factor": 1.00,
        "tropism": {"liver": 2.20, "spleen": 0.75, "kidney": 0.90, "heart": 0.72, "muscle": 0.95, "lung": 0.65, "brain": 0.14, "rest": 0.90},
        "source": "https://doi.org/10.1016/j.xphs.2023.10.005",
        "additional_sources": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC11919325/"],
    },
    "aav9": {
        "label": "AAV9", "evidence": "strong", "species": "mouse + NHP",
        "persistence_factor": 1.00,
        "tropism": {organ: 1.0 for organ in ["liver", "spleen", "kidney", "heart", "muscle", "lung", "brain", "rest"]},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11404148/",
        "additional_sources": [
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC7769048/",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC11919325/",
        ],
    },
    "aavrh10": {
        "label": "AAVrh.10", "evidence": "medium", "species": "NHP",
        "persistence_factor": 1.00,
        "tropism": {"liver": 0.82, "spleen": 0.78, "kidney": 0.75, "heart": 1.05, "muscle": 0.95, "lung": 0.90, "brain": 1.55, "rest": 0.95},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11919325/",
        "additional_sources": [
            "https://pubmed.ncbi.nlm.nih.gov/39863928/",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC7769048/",
        ],
    },
    "php-eb": {
        "label": "PHP.eB", "evidence": "exploratory", "species": "Ly6a-positive mouse only",
        "persistence_factor": 0.95,
        "tropism": {"liver": 0.55, "spleen": 0.72, "kidney": 0.65, "heart": 0.72, "muscle": 0.82, "lung": 0.70, "brain": 12.0, "rest": 0.90},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7193641/",
        "additional_sources": [
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC11919325/",
            "https://pubmed.ncbi.nlm.nih.gov/40337478/",
        ],
    },
    "cap-b10": {
        "label": "CAP-B10", "evidence": "exploratory", "species": "mouse",
        "persistence_factor": 0.95,
        "tropism": {"liver": 0.20, "spleen": 0.55, "kidney": 0.48, "heart": 0.58, "muscle": 0.72, "lung": 0.58, "brain": 14.0, "rest": 0.85},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11919325/",
        "additional_sources": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC9621732/"],
    },
    "lk03": {
        "label": "AAV-LK03", "evidence": "exploratory", "species": "human-hepatocyte prior",
        "persistence_factor": 1.05,
        "tropism": {"liver": 2.80, "spleen": 0.48, "kidney": 0.55, "heart": 0.38, "muscle": 0.42, "lung": 0.48, "brain": 0.10, "rest": 0.55},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11919325/",
        "additional_sources": ["https://pmc.ncbi.nlm.nih.gov/articles/PMC9621732/"],
    },
}

TARGETS = {
    "CNS": {"model_organ": "brain", "cellular_states": ("C_Epi", "C_M", "C_P"), "episome_half_life_days": 365.0, "route": "IV"},
    "Liver": {"model_organ": "liver", "cellular_states": ("Epi", "M", "P"), "episome_half_life_days": 120.0, "route": "IV"},
    "Heart": {"model_organ": "heart", "episome_half_life_days": 300.0, "route": "IV"},
    "Muscle": {"model_organ": "muscle", "episome_half_life_days": 365.0, "route": "IV"},
    "Kidney": {"model_organ": "kidney", "cellular_states": ("K_Epi", "K_M", "K_P"), "episome_half_life_days": 90.0, "route": "IV"},
    # The current PBPK has no eye compartment. Eye is an explicit local-route
    # surrogate based on the brain barrier compartment plus ocular capsid priors.
    "Eye": {"model_organ": "brain", "cellular_states": ("C_Epi", "C_M", "C_P"), "episome_half_life_days": 540.0, "route": "local ocular", "surrogate": True},
}

EYE_PRIOR = {"aav2": 12.0, "aav5": 2.4, "aav8": 6.5, "aav9": 2.0, "aavrh10": 1.8, "php-eb": 0.8, "cap-b10": 0.7, "lk03": 0.5}
OFF_TARGET_WEIGHTS = {"liver": 1.3, "spleen": 0.8, "kidney": 1.0, "heart": 1.0, "muscle": 0.5, "lung": 1.0, "brain": 1.4, "rest": 0.4}

# Reduced post-BBB spatial profiles. Layer weights correspond to
# perivascular/superficial, cortical parenchymal, and deep-nuclei exposure.
CNS_PROFILES = {
    "cortical_excitatory": {"weights": (0.60, 0.35, 0.05), "cell_access": 0.90, "persistence": 1.00, "depth_mm": 0.8},
    "cortical_inhibitory": {"weights": (0.45, 0.45, 0.10), "cell_access": 0.68, "persistence": 0.95, "depth_mm": 1.2},
    "cortical_projection": {"weights": (0.30, 0.52, 0.18), "cell_access": 0.78, "persistence": 1.00, "depth_mm": 1.8},
    "synaptic_neuron": {"weights": (0.38, 0.47, 0.15), "cell_access": 0.82, "persistence": 1.00, "depth_mm": 1.5},
    "deep_striatal": {"weights": (0.08, 0.32, 0.60), "cell_access": 0.62, "persistence": 1.00, "depth_mm": 4.5},
    "hypothalamic": {"weights": (0.05, 0.23, 0.72), "cell_access": 0.55, "persistence": 1.00, "depth_mm": 5.5},
    "broad_neuronal": {"weights": (0.34, 0.41, 0.25), "cell_access": 0.75, "persistence": 1.00, "depth_mm": 2.5},
    "neural_progenitor": {"weights": (0.48, 0.37, 0.15), "cell_access": 0.58, "persistence": 0.30, "depth_mm": 1.6},
}

HEATMAP_TIME_H = np.unique(np.r_[
    np.linspace(0.0, 2.0, 17),
    np.linspace(2.5, 12.0, 20),
    np.linspace(14.0, 24.0, 6),
    np.linspace(28.0, 72.0, 12),
])


def trapezoid(module: dict, y: np.ndarray, x: np.ndarray) -> float:
    return module["auc_trapz"](y, x)


def solve_cns_multilevel(t_pbpk: np.ndarray, brain_isf_amount: np.ndarray, dose: float, profile: dict) -> dict:
    """Solve post-BBB depth transport and neural-cell transduction for one profile."""
    weights = np.asarray(profile["weights"], dtype=float)

    def rhs(t: float, y: np.ndarray) -> list[float]:
        superficial, cortical, deep, bound, endosome, cytosol, nucleus, epi = np.maximum(y, 0.0)
        post_bbb = float(np.interp(t, t_pbpk, brain_isf_amount)) / max(dose, 1e-30)
        target_exposure = float(np.dot(weights, [superficial, cortical, deep])) * float(profile["cell_access"])
        return [
            0.55 * post_bbb - (0.18 + 0.08) * superficial + 0.025 * cortical,
            0.18 * superficial - (0.025 + 0.025 + 0.035) * cortical + 0.012 * deep,
            0.025 * cortical - (0.012 + 0.050) * deep,
            0.45 * target_exposure - (0.15 + 0.05) * bound,
            0.15 * bound - (0.008 + 0.10 + 0.02) * endosome,
            0.008 * endosome - (0.018 + 0.01) * cytosol,
            0.018 * cytosol - (0.018 + 0.005) * nucleus,
            0.018 * nucleus - 0.012 * epi,
        ]

    solution = solve_ivp(rhs, (float(t_pbpk[0]), float(t_pbpk[-1])), np.zeros(8), t_eval=t_pbpk, rtol=1e-8, atol=1e-12)
    if not solution.success:
        raise RuntimeError(f"CNS multilevel solve failed: {solution.message}")
    layers = solution.y[:3]
    target_layer_signal = weights @ layers * float(profile["cell_access"])
    mean_layer_signal = np.mean(layers, axis=0)
    target_auc = float(np.trapezoid(target_layer_signal, solution.t))
    reference_auc = float(np.trapezoid(mean_layer_signal, solution.t))
    return {
        "episome_auc_signal": float(np.trapezoid(solution.y[7], solution.t)),
        "target_layer_auc_signal": target_auc,
        "specificity_adjustment_log10": float(np.log10((target_auc + 1e-30) / (reference_auc + 1e-30))),
        "depth_mm": float(profile["depth_mm"]),
        "cell_access": float(profile["cell_access"]),
        "layer_auc_fraction_pct": float(100.0 * target_auc / max(float(np.trapezoid(np.sum(layers, axis=0), solution.t)), 1e-30)),
    }


def run_capsid(module: dict, capsid: dict) -> dict:
    p = module["make_params"]()
    for organ in module["ORGANS"]:
        multiplier = capsid["tropism"][organ]
        p[f"PS_{organ}"] = float(p[f"PS_{organ}"]) * multiplier
        p[f"Kp_{organ}"] = float(p[f"Kp_{organ}"]) * np.sqrt(multiplier)
    p["k_bbb_trans"] = float(p["k_bbb_trans"]) * np.sqrt(capsid["tropism"]["brain"])

    t_eval = np.unique(np.r_[np.linspace(0.0, 2.0, 260, endpoint=False), np.linspace(2.0, 24.0, 320, endpoint=False), np.linspace(24.0, 72.0, 260)])
    solution = module["solve_model"](t_eval, module["make_initial_condition"](p), p, post_infusion_max_step=0.5)
    dose = float(p["dose_vg"])
    organs = {}
    for organ in module["ORGANS"]:
        amount = solution.y[module["IDX"][f"A_{organ}_isf"]]
        vascular_amount = solution.y[module["IDX"][f"A_{organ}_v"]]
        concentration = amount / float(p[f"V_{organ}_isf"])
        vascular_concentration = vascular_amount / float(p[f"V_{organ}_v"])
        organs[organ] = {
            "auc_amount_vg_h": trapezoid(module, amount, solution.t),
            "auc_concentration_vg_h_ml": trapezoid(module, concentration, solution.t),
            "peak_isf_amount_vg": float(np.max(amount)),
            "peak_isf_concentration_vg_ml": float(np.max(concentration)),
            "peak_post_barrier_delivery_pct": float(100.0 * np.max(amount) / dose),
            "tmax_h": float(solution.t[int(np.argmax(amount))]),
            "isf_amount_vg": np.interp(HEATMAP_TIME_H, solution.t, amount).tolist(),
            "isf_concentration_vg_ml": np.interp(HEATMAP_TIME_H, solution.t, concentration).tolist(),
            "vascular_concentration_vg_ml": np.interp(
                HEATMAP_TIME_H, solution.t, vascular_concentration
            ).tolist(),
        }
    total_organ_auc = sum(metrics["auc_amount_vg_h"] for metrics in organs.values())
    for metrics in organs.values():
        metrics["exposure_share_pct"] = float(
            100.0 * metrics["auc_amount_vg_h"] / max(total_organ_auc, 1e-30)
        )
    brain_isf_amount = solution.y[module["IDX"]["A_brain_isf"]]
    cns_profiles = {
        profile_id: solve_cns_multilevel(solution.t, brain_isf_amount, dose, profile)
        for profile_id, profile in CNS_PROFILES.items()
    }
    cellular = {}
    for target_name, target in TARGETS.items():
        states = target.get("cellular_states")
        if not states:
            continue
        epi, mrna, protein = (solution.y[module["IDX"][state]] for state in states)
        cellular[target_name] = {
            "peak_epi": float(np.max(epi)),
            "auc_epi_vg_h": trapezoid(module, epi, solution.t),
            "epi_72h": float(epi[-1]),
            "peak_vector_rna": float(np.max(mrna)),
            "peak_vector_expression": float(np.max(protein)),
        }

    native_cellular_organs = {"liver": "Liver", "kidney": "Kidney", "brain": "CNS"}
    for organ, target_name in native_cellular_organs.items():
        organs[organ]["peak_episome"] = cellular[target_name]["peak_epi"]
        organs[organ]["episome_auc_vg_h"] = cellular[target_name]["auc_epi_vg_h"]
        organs[organ]["peak_expression"] = cellular[target_name]["peak_vector_expression"]
        organs[organ]["transduction_model"] = "native ode1.0 intracellular module"
        organs[organ]["model_status"] = "ode-derived"
    for organ in module["ORGANS"]:
        if organ not in native_cellular_organs:
            organs[organ]["peak_episome"] = None
            organs[organ]["episome_auc_vg_h"] = None
            organs[organ]["peak_expression"] = None
            organs[organ]["transduction_model"] = "PBPK ISF exposure only"
            organs[organ]["model_status"] = "exposure-only"
    return {
        "organs": organs,
        "cellular": cellular,
        "cns_profiles": cns_profiles,
        "max_mass_balance_error": float(np.max(np.abs(module["mass_balance_error"](solution)))),
    }


def build_organ_heatmap_payload(simulated: dict) -> dict:
    """Serialize capsid-by-organ ODE outputs for the React anatomical map."""
    capsids = []
    for capsid_id, capsid_prior in CAPSID_PRIORS.items():
        result = simulated[capsid_id]
        capsids.append({
            "capsid_id": capsid_id,
            "capsid": capsid_prior["label"],
            "evidence": capsid_prior["evidence"],
            "species": capsid_prior["species"],
            "source": capsid_prior["source"],
            "max_mass_balance_error": result["max_mass_balance_error"],
            "organs": result["organs"],
        })
    return {
        "time_h": HEATMAP_TIME_H.tolist(),
        "organs": ["brain", "lung", "heart", "liver", "spleen", "kidney", "muscle", "rest"],
        "capsids": capsids,
        "reference_anatomy": "human SVG projection",
        "reference_model": "adult mouse-scale PBPK",
        "interpretation": (
            "Relative organ distribution demo; the anatomy is human-shaped but the current "
            "physiology is adult mouse-scale and is not a calibrated human prediction."
        ),
    }


def human_capsid_tropism(capsid_id: str, capsid: dict) -> tuple[dict[str, float], str]:
    """Return cautious reference-human priors without mouse-only BBB inflation."""
    tropism = dict(capsid["tropism"])
    note = "Preclinical organ prior transferred to reference-human physiology; not clinically calibrated."
    if capsid_id == "php-eb":
        tropism["brain"] = 1.0
        note = (
            "Mouse LY6A-dependent CNS gain removed for human projection; brain prior reset to "
            "AAV9-like baseline because primates do not express LY6A."
        )
    elif capsid_id == "cap-b10":
        tropism["brain"] = 1.5
        note = (
            "Mouse CAP-B10 CNS gain strongly down-weighted for human projection; no calibrated "
            "human BBB transport estimate is available."
        )
    return tropism, note


def refined_tmax_h(solution, state_name: str, scale: float = 1.0) -> float:
    """Refine a sampled maximum against the dense ODE solution."""
    state_index = human_model.IDX[state_name]
    sampled = np.maximum(solution.y[state_index] / scale, 0.0)
    peak_index = int(np.argmax(sampled))
    if solution.sol is None or peak_index == 0 or peak_index == len(solution.t) - 1:
        return float(solution.t[peak_index])
    lower = float(solution.t[peak_index - 1])
    upper = float(solution.t[peak_index + 1])
    optimum = minimize_scalar(
        lambda t_h: -max(float(solution.sol(t_h)[state_index]) / scale, 0.0),
        bounds=(lower, upper),
        method="bounded",
        options={"xatol": 1e-4},
    )
    return float(optimum.x if optimum.success else solution.t[peak_index])


def build_human_spatial_payload() -> dict:
    """Solve capsid-specific reference-human regional PBPK for each route."""
    time_h = human_model.make_time_grid()
    route_payloads = []
    for route_id, route in human_model.ADMINISTRATION_ROUTES.items():
        capsids = []
        for capsid_id, capsid in CAPSID_PRIORS.items():
            tropism, translation_note = human_capsid_tropism(capsid_id, capsid)
            solution = human_model.solve_human_capsid(
                tropism,
                t_eval=time_h,
                administration=route_id,
            )
            dose = float(human_model.DOSE_VG)
            regions = {}
            for region_id, region in human_model.REGIONS.items():
                vascular = np.maximum(solution.y[human_model.IDX[f"A_{region_id}_v"]], 0.0)
                isf = np.maximum(solution.y[human_model.IDX[f"A_{region_id}_isf"]], 0.0)
                epi = np.maximum(solution.y[human_model.IDX[f"A_{region_id}_epi"]], 0.0)
                protein = np.maximum(solution.y[human_model.IDX[f"A_{region_id}_protein"]], 0.0)
                isf_concentration = isf / region.isf_ml
                vascular_concentration = vascular / region.vascular_ml
                regions[region_id] = {
                    "label": region.label,
                    "parent_organ": region.parent,
                    "route": region.route,
                    "vascular_volume_ml": region.vascular_ml,
                    "isf_volume_ml": region.isf_ml,
                    "blood_flow_ml_h": human_model.CARDIAC_OUTPUT_ML_H * region.flow_fraction,
                    "effective_exchange_flow_ml_h": (
                        human_model.CARDIAC_OUTPUT_ML_H
                        * human_model.EFFECTIVE_FLOW_SCALE
                        * region.flow_fraction
                    ),
                    "auc_isf_amount_vg_h": float(np.trapezoid(isf, solution.t)),
                    "auc_isf_concentration_vg_h_ml": float(np.trapezoid(isf_concentration, solution.t)),
                    "peak_isf_concentration_vg_ml": float(np.max(isf_concentration)),
                    "peak_post_barrier_delivery_pct": float(100.0 * np.max(isf) / dose),
                    "tmax_isf_h": refined_tmax_h(solution, f"A_{region_id}_isf", region.isf_ml),
                    "peak_episome_vg": float(np.max(epi)),
                    "peak_protein_au": float(np.max(protein)),
                    "tmax_protein_h": refined_tmax_h(solution, f"A_{region_id}_protein"),
                    "vascular_concentration_vg_ml": vascular_concentration.astype(np.float32).tolist(),
                    "isf_concentration_vg_ml": isf_concentration.astype(np.float32).tolist(),
                    "episome_vg": epi.astype(np.float32).tolist(),
                    "protein_au": protein.astype(np.float32).tolist(),
                }
            total_auc = sum(region["auc_isf_amount_vg_h"] for region in regions.values())
            for region in regions.values():
                region["exposure_share_pct"] = float(
                    100.0 * region["auc_isf_amount_vg_h"] / max(total_auc, 1e-30)
                )

            circulation = {}
            for compartment, volume_ml in human_model.CIRCULATION_VOLUMES_ML.items():
                amount = np.maximum(solution.y[human_model.IDX[f"A_{compartment}"]], 0.0)
                circulation[compartment] = {
                    "volume_ml": volume_ml,
                    "concentration_vg_ml": (amount / volume_ml).astype(np.float32).tolist(),
                }
            route_compartments = {}
            for compartment, volume_ml in human_model.ROUTE_COMPARTMENT_VOLUMES_ML.items():
                amount = np.maximum(solution.y[human_model.IDX[f"A_{compartment}"]], 0.0)
                route_compartments[compartment] = {
                    "volume_ml": volume_ml,
                    "concentration_vg_ml": (amount / volume_ml).astype(np.float32).tolist(),
                }
            for depot_id in ("im_depot", "airway_depot"):
                depot = np.maximum(solution.y[human_model.IDX[f"A_{depot_id}"]], 0.0)
                route_compartments[depot_id] = {"amount_vg": depot.astype(np.float32).tolist()}
            balance = human_model.mass_balance_error(solution)
            meaningful = solution.y[human_model.IDX["Dose_in"]] > dose * 1e-9
            capsids.append({
                "capsid_id": capsid_id,
                "capsid": capsid["label"],
                "evidence": capsid["evidence"],
                "source_species": capsid["species"],
                "source": capsid["source"],
                "additional_sources": capsid.get("additional_sources", []),
                "human_translation_note": translation_note,
                "max_mass_balance_error": float(np.max(np.abs(balance[meaningful]))),
                "regions": regions,
                "circulation": circulation,
                "route_compartments": route_compartments,
            })
        route_payloads.append({
            "route_id": route_id,
            "label": route.label,
            "label_zh": route.label_zh,
            "description": route.description,
            "description_zh": route.description_zh,
            "infusion_duration_h": route.infusion_duration_h,
            "route_class": route.route_class,
            "evidence_source": route.evidence_source,
            "capsids": capsids,
        })
    return {
        "time_h": time_h.tolist(),
        "max_time_days": float(time_h[-1] / 24.0),
        "dose_vg": float(human_model.DOSE_VG),
        "dose_vg_per_kg": float(human_model.DOSE_VG_PER_KG),
        "body_weight_kg": float(human_model.BODY_WEIGHT_KG),
        "default_route_id": "iv",
        "administration_routes": route_payloads,
        "reference_model": "70-kg reference-adult, mouse-equation-aligned multiregion PBPK",
        "physiology_status": "reference-human physiology; capsid parameters exploratory",
        "cardiac_output_ml_h": float(human_model.CARDIAC_OUTPUT_ML_H),
        "effective_flow_scale": float(human_model.EFFECTIVE_FLOW_SCALE),
        "total_modeled_blood_volume_ml": float(human_model.TOTAL_MODELED_BLOOD_VOLUME_ML),
        "reference_blood_volume_ml": float(human_model.REFERENCE_BLOOD_VOLUME_ML),
        "csf_total_volume_ml": float(human_model.CSF_TOTAL_VOLUME_ML),
        "csf_production_ml_h": float(human_model.CSF_PRODUCTION_ML_H),
        "csf_absorption_half_life_h": float(human_model.CSF_ABSORPTION_HALF_LIFE_H),
        "aav9_capsid_half_life_h": {
            "blood": float(human_model.BLOOD_CAPSID_HALF_LIFE_H),
            **{
                organ: float(value)
                for organ, value in human_model.CAPSID_HALF_LIFE_H.items()
            },
        },
        "aav9_capsid_half_life_provenance": human_model.REFERENCE_HUMAN_AAV9_PROVENANCE,
        "parameter_evidence": human_model.AAV9_PARAMETER_EVIDENCE,
        "physiology_sources": human_model.HUMAN_PHYSIOLOGY_SOURCES,
        "equation_family": "ode1.0 Q-PS-Kp-J_res-J_deg plus receptor-to-episome trafficking",
        "state_count": len(human_model.STATE_NAMES),
        "region_ids": list(human_model.REGIONS),
        "circulation_ids": list(human_model.CIRCULATION_VOLUMES_ML),
        "interpretation": (
            "Human physiology projection using the same Q-PS-Kp-clearance and intracellular "
            "equations and Q scaling as the mouse model. IV, lumbar intrathecal, intracisternal, "
            "intracerebroventricular, deltoid intramuscular and inhaled routes change the ODE "
            "input compartment and are re-solved. "
            "It is a mechanistic research demo, not a clinically validated dose prediction."
        ),
    }


def build_human_route_summary(human_payload: dict) -> list[dict]:
    """Collapse regional trajectories into route/capsid/organ design metrics."""
    summaries = []
    for route in human_payload["administration_routes"]:
        for capsid in route["capsids"]:
            by_organ = {}
            parent_ids = sorted({region["parent_organ"] for region in capsid["regions"].values()})
            for parent in parent_ids:
                regions = [
                    region for region in capsid["regions"].values()
                    if region["parent_organ"] == parent
                ]
                total_isf_ml = sum(region["isf_volume_ml"] for region in regions)
                auc_amount = sum(region["auc_isf_amount_vg_h"] for region in regions)
                by_organ[parent] = {
                    "auc_isf_concentration_vg_h_ml": float(auc_amount / max(total_isf_ml, 1e-30)),
                    "exposure_share_pct": float(sum(region["exposure_share_pct"] for region in regions)),
                    "peak_post_barrier_delivery_pct": float(sum(region["peak_post_barrier_delivery_pct"] for region in regions)),
                    "median_tmax_isf_h": float(np.median([region["tmax_isf_h"] for region in regions])),
                    "peak_protein_au": float(sum(region["peak_protein_au"] for region in regions)),
                }
            summaries.append({
                "route_id": route["route_id"],
                "route_label": route["label"],
                "route_label_zh": route["label_zh"],
                "route_class": route["route_class"],
                "evidence_source": route["evidence_source"],
                "capsid_id": capsid["capsid_id"],
                "capsid": capsid["capsid"],
                "evidence": capsid["evidence"],
                "organs": by_organ,
            })
    return summaries


def target_metrics(capsid_id: str, capsid: dict, result: dict, target_name: str, target: dict) -> dict:
    model_organ = target["model_organ"]
    organs = result["organs"]
    target_auc = organs[model_organ]["auc_concentration_vg_h_ml"]
    target_amount_auc = organs[model_organ]["auc_amount_vg_h"]
    target_peak_pct = organs[model_organ]["peak_post_barrier_delivery_pct"]
    local_factor = (
        EYE_PRIOR[capsid_id] / capsid["tropism"][model_organ]
        if target.get("surrogate") else 1.0
    )
    target_auc *= local_factor
    target_amount_auc *= local_factor
    target_peak_pct *= local_factor

    off_target = np.mean([
        metrics["auc_concentration_vg_h_ml"] * OFF_TARGET_WEIGHTS[organ]
        for organ, metrics in organs.items() if organ != model_organ
    ])
    if target.get("surrogate"):
        # Local ocular administration reduces systemic exposure in this MVP.
        off_target *= 0.08
    specificity = float(np.log10((target_auc + 1e-30) / (off_target + 1e-30)))
    total_amount_auc = sum(metrics["auc_amount_vg_h"] for metrics in organs.values())
    exposure_share = float(100.0 * target_amount_auc / max(total_amount_auc + target_amount_auc * (local_factor - 1.0), 1e-30))
    cellular = result["cellular"].get(target_name)
    if cellular is not None:
        episome_signal = cellular["auc_epi_vg_h"] * local_factor
        expression_signal = cellular["peak_vector_expression"] * local_factor
        transduction_model = "native ode1.0 intracellular module"
    else:
        # The current source model has no heart/muscle cellular states. Keep
        # these targets usable, but expose the reduced assumption explicitly.
        episome_signal = target_amount_auc
        expression_signal = target_auc
        transduction_model = "PBPK ISF-driven reduced cellular surrogate"
    return {
        "capsid_id": capsid_id,
        "capsid": capsid["label"],
        "target": target_name,
        "model_organ": model_organ,
        "route": target["route"],
        "specificity_log10": specificity,
        "target_concentration_auc_vg_h_ml": float(target_auc),
        "target_exposure_share_pct": exposure_share,
        "peak_post_barrier_delivery_pct": float(target_peak_pct),
        "tmax_h": organs[model_organ]["tmax_h"],
        "tropism_multiplier": float(capsid["tropism"][model_organ] * local_factor),
        "episome_half_life_days_prior": float(target["episome_half_life_days"]),
        "persistence_factor": float(capsid["persistence_factor"]),
        "evidence": capsid["evidence"],
        "species": capsid["species"],
        "source": capsid["source"],
        "model_status": "surrogate" if target.get("surrogate") or not target.get("cellular_states") else "ode-derived",
        "max_mass_balance_error": result["max_mass_balance_error"],
        "expression_competent_episome_auc_signal": float(episome_signal),
        "peak_vector_expression_signal": float(expression_signal),
        "transduction_model": transduction_model,
    }


def solve_sineup_pd(relative_epi: float, episome_half_life_days: float) -> dict:
    """Solve Epi -> SINEUP RNA -> endogenous protein on a 0-730 day horizon."""
    k_epi_loss = np.log(2.0) / max(episome_half_life_days, 1e-12)
    k_sineup_deg = np.log(2.0) / 0.25  # 6 h RNA half-life, matching ode1.0 mRNA logic
    k_protein_turnover = np.log(2.0) / 2.0  # 48 h protein half-life
    k_sineup_tx = 4.0
    ec50_epi = 0.30
    ec50_sineup = 0.50
    max_translation_boost = 1.0

    def rhs(_t: float, y: np.ndarray) -> list[float]:
        epi, sineup, protein = np.maximum(y, 0.0)
        tx = k_sineup_tx * epi / (ec50_epi + epi + 1e-30)
        translation_gain = max_translation_boost * sineup / (ec50_sineup + sineup + 1e-30)
        protein_setpoint = min(1.0, 0.5 * (1.0 + translation_gain))
        return [
            -k_epi_loss * epi,
            tx - k_sineup_deg * sineup,
            k_protein_turnover * (protein_setpoint - protein),
        ]

    t_eval = np.linspace(0.0, 730.0, 2921)
    solution = solve_ivp(rhs, (0.0, 730.0), [max(relative_epi, 0.0), 0.0, 0.5], t_eval=t_eval, rtol=1e-8, atol=1e-10)
    if not solution.success:
        raise RuntimeError(f"SINEUP PD solve failed: {solution.message}")
    protein = solution.y[2]
    threshold = 0.65
    above = np.flatnonzero(protein >= threshold)
    onset = float(solution.t[above[0]]) if len(above) else 0.0
    end = float(solution.t[above[-1]]) if len(above) else 0.0
    peak_index = int(np.argmax(protein))
    return {
        "predicted_protein_restoration_pct": float(100.0 * protein[peak_index]),
        "peak_restoration_day": float(solution.t[peak_index]),
        "therapeutic_onset_days": onset,
        "therapeutic_window_days": float(max(0.0, end - onset)),
        "effective_duration_days": end,
        "therapeutic_threshold_pct": float(100.0 * threshold),
        "persistence_censored_at_730d": bool(len(above) and above[-1] == len(solution.t) - 1),
        "pd_episome_loss_rate_per_day": float(k_epi_loss),
        "pd_sineup_rna_half_life_days": 0.25,
        "pd_target_protein_half_life_days": 2.0,
    }


def add_pd_ode(rows: list[dict]) -> None:
    group_keys = list(dict.fromkeys((row["target"], row.get("cns_profile")) for row in rows))
    for target_name, cns_profile in group_keys:
        target_rows = [row for row in rows if row["target"] == target_name and row.get("cns_profile") == cns_profile]
        maximum = max(row["expression_competent_episome_auc_signal"] for row in target_rows)
        for row in target_rows:
            relative_delivery = row["expression_competent_episome_auc_signal"] / max(maximum, 1e-30)
            effective_half_life = row["episome_half_life_days_prior"] * row["persistence_factor"]
            row["relative_delivery"] = float(relative_delivery)
            row.update(solve_sineup_pd(relative_delivery, effective_half_life))


def build_cns_profile_rows(base_rows: list[dict], simulated: dict) -> list[dict]:
    rows = []
    cns_base = {row["capsid_id"]: row for row in base_rows if row["target"] == "CNS"}
    for profile_id, profile in CNS_PROFILES.items():
        for capsid_id in CAPSID_PRIORS:
            row = dict(cns_base[capsid_id])
            spatial = simulated[capsid_id]["cns_profiles"][profile_id]
            row.update({
                "cns_profile": profile_id,
                "cns_depth_mm": spatial["depth_mm"],
                "cns_cell_access_factor": spatial["cell_access"],
                "cns_target_layer_auc_fraction_pct": spatial["layer_auc_fraction_pct"],
                "specificity_log10": float(row["specificity_log10"] + spatial["specificity_adjustment_log10"]),
                "expression_competent_episome_auc_signal": spatial["episome_auc_signal"],
                "episome_half_life_days_prior": float(row["episome_half_life_days_prior"] * profile["persistence"]),
                "transduction_model": "ode1.0 BBB output -> 3-depth CNS transport -> intracellular transduction ODE",
                "model_status": "ode-derived",
            })
            rows.append(row)
    add_pd_ode(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, default=Path(__file__).with_name("ode1.0.py"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    module = runpy.run_path(str(args.model))
    simulated = {capsid_id: run_capsid(module, capsid) for capsid_id, capsid in CAPSID_PRIORS.items()}
    rows = [
        target_metrics(capsid_id, capsid, simulated[capsid_id], target_name, target)
        for target_name, target in TARGETS.items()
        for capsid_id, capsid in CAPSID_PRIORS.items()
    ]
    add_pd_ode(rows)
    cns_profile_rows = build_cns_profile_rows(rows, simulated)
    human_spatial_payload = build_human_spatial_payload()
    human_route_summary = build_human_route_summary(human_spatial_payload)

    payload = {
        "schema_version": "3.2",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "reference_model": str(args.model),
        "reference_species": "adult mouse-scale PBPK; capsid priors retain source-species labels",
        "dose_vg": 1.0e12,
        "model_counts": {
            "organ_level": len(rows),
            "cns_profile": len(cns_profile_rows),
            "organ_heatmap": len(CAPSID_PRIORS) * len(module["ORGANS"]),
            "human_spatial": (
                len(human_model.ADMINISTRATION_ROUTES)
                * len(CAPSID_PRIORS)
                * len(human_model.REGIONS)
            ),
            "total": (
                len(rows) + len(cns_profile_rows) + len(CAPSID_PRIORS) * len(module["ORGANS"])
                + len(human_model.ADMINISTRATION_ROUTES)
                * len(CAPSID_PRIORS)
                * len(human_model.REGIONS)
            ),
        },
        "method": {
            "delivery_efficiency": "100 * max(target organ ISF amount) / administered dose",
            "organ_specificity": "log10(target organ ISF concentration AUC / mean toxicity-weighted off-target concentration AUC)",
            "persistence": "numerically solved Epi -> SINEUP RNA -> endogenous protein ODE; duration ends when protein falls below 65%",
            "cellular_transduction": "liver, kidney and CNS use native intracellular states from ode1.0.py; heart and muscle use an explicitly labelled PBPK-ISF-driven reduced surrogate; Eye reuses the CNS chain with a local-route prior",
            "capsid_parameterization": "PS and Kp are scaled by literature-informed capsid-organ tropism priors; BBB transcytosis is also scaled for brain",
            "capsid_prior_status": "relative tropism multipliers remain literature-informed priors; the linked head-to-head studies are recorded, but cross-study values are not treated as one calibrated numeric assay",
            "important_limit": "Eye uses a local-route barrier surrogate because the PBPK model does not yet contain an anatomical eye compartment",
            "human_spatial_model": "mouse-equation-aligned Q-PS-Kp-J_res-J_deg PBPK expanded to IV, lumbar intrathecal, intracisternal, intracerebroventricular, deltoid intramuscular and airway-depot inputs, cardiopulmonary circulation, bidirectional CSF transport, 24 human regions, receptor uptake, intracellular trafficking, episome, mRNA and protein",
            "human_translation_limit": "reference-human physiological flows and volumes are mechanistic priors; capsid-specific parameters require NHP/human fitting and uncertainty analysis",
        },
        "results": rows,
        "cns_profile_results": cns_profile_rows,
        "organ_heatmap": build_organ_heatmap_payload(simulated),
        "human_spatial_heatmap_file": "/data/human-spatial-results.json",
        "human_route_summary": human_route_summary,
        "parameter_evidence": module["AAV9_PARAMETER_EVIDENCE"],
        "capsid_tropism_literature_file": "/data/aav_capsid_tropism_literature.csv",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    args.output.with_name("aav_capsid_tropism_literature.csv").write_bytes(
        TROPISM_LITERATURE_PATH.read_bytes()
    )
    human_output = args.output.with_name("human-spatial-results.json")
    human_output.write_text(
        json.dumps(human_spatial_payload, ensure_ascii=False, separators=(",", ":")),
        encoding="utf-8",
    )

    all_rows = rows + cns_profile_rows
    fieldnames = list(dict.fromkeys(key for row in all_rows for key in row))
    csv_path = args.output.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Wrote {len(all_rows)} modeled design points to {args.output}")
    print(
        f"Wrote {len(human_model.ADMINISTRATION_ROUTES)} routes x "
        f"{len(human_model.REGIONS)}-region human trajectories to {human_output}"
    )


if __name__ == "__main__":
    main()
