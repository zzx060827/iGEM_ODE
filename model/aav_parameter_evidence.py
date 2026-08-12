"""Traceable AAV9 capsid-PK priors used by the mouse and human projections.

The mouse organ half-lives are apparent log-linear decline constants fitted to
the mean 125I-AAV9 organ concentrations in Wang et al. (2024), Table S1. 125I
catabolites leave cells, so this signal is used as an early intact/extracellular
capsid proxy. It is not an episome or transgene-expression half-life.

The reference-human projection uses the closer-species NHP PET estimates from
Ballon et al. (2020) where available. Kidney and lung remain explicitly marked
mouse-derived provisional values because that NHP table did not report them.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np


AAV9_PK_SOURCES = {
    "wang_2024_mouse_dual_isotope": {
        "title": (
            "A novel approach to quantitate biodistribution and transduction "
            "of adeno-associated virus gene therapy using radiolabeled AAV vectors in mice"
        ),
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11404148/",
        "doi": "10.1016/j.omtm.2024.101326",
        "species": "C57BL/6 mouse",
        "route": "intravenous",
        "readout": "125I-AAV9 mean %ID/g, Table S1",
    },
    "seo_2020_mouse_pet": {
        "title": "Positron emission tomography imaging of novel AAV capsids maps rapid brain accumulation",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7193641/",
        "species": "mouse",
        "route": "intravenous",
        "readout": "64Cu-AAV blood circulation",
    },
    "ballon_2020_nhp_pet": {
        "title": "Quantitative Whole-Body Imaging of I-124-Labeled Adeno-Associated Viral Vector Biodistribution in Nonhuman Primates",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7769048/",
        "doi": "10.1089/hum.2020.116",
        "species": "nonhuman primate",
        "route": "intravenous",
        "readout": "I-124-AAV9 organ biological half-life, immune-naive animal, Table 3",
    },
}


MOUSE_AAV9_125I_TIME_H = np.array([4.0, 24.0, 72.0, 168.0])
MOUSE_AAV9_125I_MEAN_ID_PER_G = {
    "liver": np.array([26.54, 19.55, 15.76, 3.63]),
    "heart": np.array([5.76, 4.50, 0.61, 0.05]),
    "lung": np.array([1.68, 1.21, 0.07, 0.01]),
    "kidney": np.array([6.56, 3.97, 0.74, 0.64]),
    "spleen": np.array([39.10, 20.07, 13.09, 5.21]),
    "brain": np.array([0.38, 0.25, 0.02, 0.00]),
    "muscle": np.array([0.56, 0.64, 0.79, 0.02]),
    "small_intestine": np.array([1.13, 0.89, 0.19, 0.01]),
    "large_intestine": np.array([3.82, 2.78, 0.58, 0.03]),
    "stomach": np.array([0.47, 0.48, 0.09, 0.01]),
    "mammary_fat": np.array([0.90, 2.20, 0.21, 0.03]),
    "white_fat": np.array([0.57, 1.02, 0.12, 0.02]),
    "bone_marrow": np.array([9.44, 8.92, 2.28, 0.18]),
}

_REST_TISSUES = (
    "small_intestine",
    "large_intestine",
    "stomach",
    "mammary_fat",
    "white_fat",
    "bone_marrow",
)
MOUSE_AAV9_125I_MEAN_ID_PER_G["rest"] = np.exp(
    np.mean(
        np.log(np.vstack([MOUSE_AAV9_125I_MEAN_ID_PER_G[name] for name in _REST_TISSUES])),
        axis=0,
    )
)

# Brain reaches the reported zero and kidney shows a near-floor plateau at day
# 7, so their terminal point is excluded. All other fits use 4 h through 7 d.
MOUSE_AAV9_FIT_WINDOWS_H = {
    "liver": (4.0, 168.0),
    "heart": (4.0, 168.0),
    "lung": (4.0, 168.0),
    "kidney": (4.0, 72.0),
    "spleen": (4.0, 168.0),
    "brain": (4.0, 72.0),
    "muscle": (4.0, 168.0),
    "rest": (4.0, 168.0),
}


def _fit_apparent_half_life(
    times_h: np.ndarray,
    values: np.ndarray,
    window_h: Iterable[float],
) -> dict[str, float | int | list[float]]:
    lower, upper = window_h
    mask = (times_h >= lower) & (times_h <= upper) & (values > 0.0)
    fit_time = times_h[mask]
    log_values = np.log(values[mask])
    if fit_time.size < 2:
        raise ValueError("At least two positive observations are required")
    slope, intercept = np.polyfit(fit_time, log_values, 1)
    if slope >= 0.0:
        raise ValueError("The selected organ time course does not decline log-linearly")
    predicted = intercept + slope * fit_time
    residual_ss = float(np.sum((log_values - predicted) ** 2))
    total_ss = float(np.sum((log_values - np.mean(log_values)) ** 2))
    return {
        "half_life_h": float(-np.log(2.0) / slope),
        "rate_h_inv": float(-slope),
        "r_squared_log": float(1.0 - residual_ss / total_ss) if total_ss > 0.0 else 1.0,
        "n_timepoints": int(fit_time.size),
        "fit_times_h": fit_time.astype(float).tolist(),
    }


MOUSE_AAV9_ORGAN_FIT = {
    organ: _fit_apparent_half_life(
        MOUSE_AAV9_125I_TIME_H,
        MOUSE_AAV9_125I_MEAN_ID_PER_G[organ],
        MOUSE_AAV9_FIT_WINDOWS_H[organ],
    )
    for organ in MOUSE_AAV9_FIT_WINDOWS_H
}

# Seo et al. report 5.0 h for unmodified AAV9. The old 2.4 h value belongs to
# the tetracysteine-modified AAV9-TC construct, not ordinary AAV9.
MOUSE_AAV9_BLOOD_HALF_LIFE_H = 5.0
MOUSE_AAV9_CAPSID_HALF_LIFE_H = {
    "blood": MOUSE_AAV9_BLOOD_HALF_LIFE_H,
    **{organ: float(fit["half_life_h"]) for organ, fit in MOUSE_AAV9_ORGAN_FIT.items()},
}


# Immune-naive IV AAV9 values from Ballon et al., Table 3. "Source" is used as
# the circulating-source prior; "body remainder" is used for muscle/rest.
NHP_AAV9_PET_HALF_LIFE_H = {
    "blood": 1.2,
    "liver": 22.6,
    "heart": 15.6,
    "spleen": 22.9,
    "brain": 24.8,
    "body_remainder": 48.7,
}
REFERENCE_HUMAN_AAV9_CAPSID_HALF_LIFE_H = {
    "blood": NHP_AAV9_PET_HALF_LIFE_H["blood"],
    "brain": NHP_AAV9_PET_HALF_LIFE_H["brain"],
    "heart": NHP_AAV9_PET_HALF_LIFE_H["heart"],
    "liver": NHP_AAV9_PET_HALF_LIFE_H["liver"],
    "spleen": NHP_AAV9_PET_HALF_LIFE_H["spleen"],
    "kidney": MOUSE_AAV9_CAPSID_HALF_LIFE_H["kidney"],
    "muscle": NHP_AAV9_PET_HALF_LIFE_H["body_remainder"],
    "lung": MOUSE_AAV9_CAPSID_HALF_LIFE_H["lung"],
    "rest": NHP_AAV9_PET_HALF_LIFE_H["body_remainder"],
}
REFERENCE_HUMAN_AAV9_PROVENANCE = {
    organ: (
        "Ballon 2020 NHP PET"
        if organ not in {"kidney", "lung"}
        else "Wang 2024 mouse 125I fit; provisional cross-species prior"
    )
    for organ in REFERENCE_HUMAN_AAV9_CAPSID_HALF_LIFE_H
}


def aav9_parameter_payload() -> dict:
    """Return JSON-serializable data provenance for frontend/report exports."""
    return {
        "interpretation": (
            "Early apparent intact/extracellular capsid decline; not episome, "
            "vector-genome, mRNA, protein, or therapeutic-effect persistence."
        ),
        "mouse": {
            "blood_half_life_h": MOUSE_AAV9_BLOOD_HALF_LIFE_H,
            "organ_half_life_h": {
                key: float(value)
                for key, value in MOUSE_AAV9_CAPSID_HALF_LIFE_H.items()
                if key != "blood"
            },
            "organ_fit_diagnostics": MOUSE_AAV9_ORGAN_FIT,
            "time_h": MOUSE_AAV9_125I_TIME_H.tolist(),
            "organ_mean_id_per_g": {
                key: value.astype(float).tolist()
                for key, value in MOUSE_AAV9_125I_MEAN_ID_PER_G.items()
            },
            "source_ids": [
                "seo_2020_mouse_pet",
                "wang_2024_mouse_dual_isotope",
            ],
        },
        "reference_human_projection": {
            "half_life_h": {
                key: float(value)
                for key, value in REFERENCE_HUMAN_AAV9_CAPSID_HALF_LIFE_H.items()
            },
            "provenance": REFERENCE_HUMAN_AAV9_PROVENANCE,
            "source_id": "ballon_2020_nhp_pet",
            "status": "NHP-informed reference projection, not human clinical calibration",
        },
        "sources": AAV9_PK_SOURCES,
        "structural_assumption": (
            "The organ apparent loss rate is split 35% vascular/endothelial and "
            "65% ISF/catabolic because the available studies do not identify "
            "those two rates separately."
        ),
    }
