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


CAPSID_PRIORS = {
    "aav2": {
        "label": "AAV2", "evidence": "medium", "species": "preclinical + clinical ocular",
        "persistence_factor": 0.95,
        "tropism": {"liver": 0.30, "spleen": 0.80, "kidney": 0.45, "heart": 0.28, "muscle": 0.22, "lung": 0.80, "brain": 0.16, "rest": 0.35},
        "source": "https://www.fda.gov/vaccines-blood-biologics/cellular-gene-therapy-products/luxturna",
    },
    "aav5": {
        "label": "AAV5", "evidence": "medium", "species": "preclinical",
        "persistence_factor": 0.90,
        "tropism": {"liver": 0.55, "spleen": 0.85, "kidney": 0.48, "heart": 0.25, "muscle": 0.22, "lung": 2.20, "brain": 0.18, "rest": 0.45},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10659018/",
    },
    "aav8": {
        "label": "AAV8", "evidence": "strong", "species": "mouse + NHP",
        "persistence_factor": 1.00,
        "tropism": {"liver": 2.20, "spleen": 0.75, "kidney": 0.90, "heart": 0.72, "muscle": 0.95, "lung": 0.65, "brain": 0.14, "rest": 0.90},
        "source": "https://www.sciencedirect.com/science/article/pii/S0022354923004148",
    },
    "aav9": {
        "label": "AAV9", "evidence": "strong", "species": "mouse + NHP",
        "persistence_factor": 1.00,
        "tropism": {organ: 1.0 for organ in ["liver", "spleen", "kidney", "heart", "muscle", "lung", "brain", "rest"]},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7769048/",
    },
    "aavrh10": {
        "label": "AAVrh.10", "evidence": "medium", "species": "NHP",
        "persistence_factor": 1.00,
        "tropism": {"liver": 0.82, "spleen": 0.78, "kidney": 0.75, "heart": 1.05, "muscle": 0.95, "lung": 0.90, "brain": 1.55, "rest": 0.95},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7769048/",
    },
    "php-eb": {
        "label": "PHP.eB", "evidence": "exploratory", "species": "Ly6a-positive mouse only",
        "persistence_factor": 0.95,
        "tropism": {"liver": 0.55, "spleen": 0.72, "kidney": 0.65, "heart": 0.72, "muscle": 0.82, "lung": 0.70, "brain": 12.0, "rest": 0.90},
        "source": "https://clover.caltech.edu/aav/faq",
    },
    "cap-b10": {
        "label": "CAP-B10", "evidence": "exploratory", "species": "mouse",
        "persistence_factor": 0.95,
        "tropism": {"liver": 0.20, "spleen": 0.55, "kidney": 0.48, "heart": 0.58, "muscle": 0.72, "lung": 0.58, "brain": 14.0, "rest": 0.85},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9621732/",
    },
    "lk03": {
        "label": "AAV-LK03", "evidence": "exploratory", "species": "human-hepatocyte prior",
        "persistence_factor": 1.05,
        "tropism": {"liver": 2.80, "spleen": 0.48, "kidney": 0.55, "heart": 0.38, "muscle": 0.42, "lung": 0.48, "brain": 0.10, "rest": 0.55},
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9621732/",
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
        concentration = amount / float(p[f"V_{organ}_isf"])
        organs[organ] = {
            "auc_amount_vg_h": trapezoid(module, amount, solution.t),
            "auc_concentration_vg_h_ml": trapezoid(module, concentration, solution.t),
            "peak_isf_amount_vg": float(np.max(amount)),
            "peak_post_barrier_delivery_pct": float(100.0 * np.max(amount) / dose),
            "tmax_h": float(solution.t[int(np.argmax(amount))]),
        }
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
    return {
        "organs": organs,
        "cellular": cellular,
        "cns_profiles": cns_profiles,
        "max_mass_balance_error": float(np.max(np.abs(module["mass_balance_error"](solution)))),
    }


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

    payload = {
        "schema_version": "2.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "reference_model": str(args.model),
        "reference_species": "adult mouse-scale PBPK; capsid priors retain source-species labels",
        "dose_vg": 1.0e12,
        "model_counts": {"organ_level": len(rows), "cns_profile": len(cns_profile_rows), "total": len(rows) + len(cns_profile_rows)},
        "method": {
            "delivery_efficiency": "100 * max(target organ ISF amount) / administered dose",
            "organ_specificity": "log10(target organ ISF concentration AUC / mean toxicity-weighted off-target concentration AUC)",
            "persistence": "numerically solved Epi -> SINEUP RNA -> endogenous protein ODE; duration ends when protein falls below 65%",
            "cellular_transduction": "liver, kidney and CNS use native intracellular states from ode1.0.py; heart and muscle use an explicitly labelled PBPK-ISF-driven reduced surrogate; Eye reuses the CNS chain with a local-route prior",
            "capsid_parameterization": "PS and Kp are scaled by literature-informed capsid-organ tropism priors; BBB transcytosis is also scaled for brain",
            "important_limit": "Eye uses a local-route barrier surrogate because the PBPK model does not yet contain an anatomical eye compartment",
        },
        "results": rows,
        "cns_profile_results": cns_profile_rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    all_rows = rows + cns_profile_rows
    fieldnames = list(dict.fromkeys(key for row in all_rows for key in row))
    csv_path = args.output.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Wrote {len(all_rows)} modeled design points to {args.output}")


if __name__ == "__main__":
    main()
