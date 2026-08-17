"""Generate an evidence-labelled AAV exposure-margin and safety-priority screen.

The output separates three questions that must not be collapsed:

1. Is the administered dose below a route-matched marketed-product context?
2. Is the modeled organ exposure below a reported signal-associated exposure?
3. Is there enough product-, capsid-, route- and payload-specific evidence to
   call the proposed construct safe?  At present the answer to (3) is no.

The calculations are useful for dose justification and toxicology planning.
They are not clinical recommendations, adverse-event probabilities or NOAELs.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODEL_DIR = Path(__file__).resolve().parent
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

import human_spatial_pbpk as human_model


HUMAN_RESULTS = ROOT / "sineup-delivery-atlas" / "public" / "data" / "human-spatial-results.json"
OUTPUT = ROOT / "sineup-delivery-atlas" / "public" / "data" / "safety-screen.json"
MARGIN_CSV = MODEL_DIR / "data" / "aav_safety_organ_margins.csv"
FRONTEND_MARGIN_CSV = ROOT / "sineup-delivery-atlas" / "public" / "data" / MARGIN_CSV.name
EVIDENCE_CSV = MODEL_DIR / "data" / "aav_safety_evidence.csv"
FRONTEND_EVIDENCE_CSV = ROOT / "sineup-delivery-atlas" / "public" / "data" / EVIDENCE_CSV.name
REPORT = ROOT / "docs" / "aav_safety_margin_assessment_2026.md"


REFERENCE_CASES: dict[str, dict[str, Any]] = {
    "zolgensma_iv": {
        "product": "ZOLGENSMA (onasemnogene abeparvovec-xioi)",
        "capsid": "AAV9",
        "route": "iv",
        "dose_value": 1.1e14,
        "dose_unit": "vg/kg",
        "total_dose_vg_for_model": human_model.BODY_WEIGHT_KG * 1.1e14,
        "evidence_kind": "marketed-product dose with serious known risks",
        "source": "https://www.fda.gov/media/126109/download",
        "interpretation": (
            "Same route and capsid as IV AAV9, but different payload, manufacturing, "
            "age and disease. The dose is not a universal safe threshold."
        ),
    },
    "zolgensma_mouse_histology_signal": {
        "product": "Onasemnogene abeparvovec mouse toxicology",
        "capsid": "AAV9",
        "route": "iv",
        "dose_value": 7.9e13,
        "dose_unit": "vg/kg",
        "total_dose_vg_for_model": human_model.BODY_WEIGHT_KG * 7.9e13,
        "evidence_kind": "lowest explicitly reported cardiac histology signal dose; not a NOAEL",
        "source": "https://www.fda.gov/media/126109/download",
        "interpretation": (
            "Neonatal-mouse product-specific signal. Directly relevant to heart and supportive "
            "for liver hazard identification; not an all-organ human threshold."
        ),
    },
    "itvisma_it": {
        "product": "ITVISMA (onasemnogene abeparvovec-brve)",
        "capsid": "AAV9",
        "route": "intrathecal",
        "dose_value": 1.2e14,
        "dose_unit": "vg/patient",
        "total_dose_vg_for_model": 1.2e14,
        "evidence_kind": "marketed-product dose with serious known risks",
        "source": "https://www.fda.gov/media/193168/download",
        "interpretation": (
            "Same route and capsid as lumbar intrathecal AAV9, but different payload, "
            "manufacturing, population and procedure. The dose is not a universal safe threshold."
        ),
    },
    "kebilidi_intraputaminal": {
        "product": "KEBILIDI (eladocagene exuparvovec-tneq)",
        "capsid": "AAV2",
        "route": "intraputaminal",
        "dose_value": 1.8e11,
        "dose_unit": "vg/patient",
        "total_dose_vg_for_model": None,
        "evidence_kind": "marketed local-CNS product dose and procedure context",
        "source": "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=9d6a6401-c6b5-4f29-af11-67707d249482",
        "interpretation": (
            "Four stereotactic putaminal infusions totaling 0.32 mL. The current model has no "
            "intraparenchymal compartment, so this case is not used as an exposure threshold."
        ),
    },
    "sineup_gdnf_mouse": {
        "product": "AAV9-miniSINEUP-GDNF proof of concept",
        "capsid": "AAV9",
        "route": "intrastriatal",
        "dose_value": 7.0e9,
        "dose_unit": "vg/mouse",
        "total_dose_vg_for_model": None,
        "evidence_kind": "mechanism and local tolerability proof of concept; not GLP toxicology",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7000958/",
        "interpretation": (
            "1 uL of 7.0e12 vg/mL AAV9; about two-fold endogenous GDNF increase for at least "
            "six months with no reported body-weight or food-intake signal. Target- and species-specific."
        ),
    },
    "itvisma_nhp_low_dose": {
        "product": "Onasemnogene abeparvovec 12-month NHP IT toxicology",
        "capsid": "AAV9",
        "route": "intrathecal",
        "dose_value": 1.2e13,
        "dose_unit": "vg/NHP",
        "total_dose_vg_for_model": None,
        "evidence_kind": "lowest NHP IT dose with microscopic neural findings; not a NOAEL",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9347375/",
        "interpretation": (
            "Minimal-to-moderate microscopic CNS/DRG/nerve findings occurred at all tested doses, "
            "including 1.2e13 vg/NHP. No clinical correlate; CSF scaling was used for human context."
        ),
    },
    "itvisma_clinical_safety": {
        "product": "ITVISMA randomized and open-label clinical safety database",
        "capsid": "AAV9",
        "route": "intrathecal",
        "dose_value": 1.2e14,
        "dose_unit": "vg/patient",
        "total_dose_vg_for_model": None,
        "evidence_kind": "clinical adverse-reaction frequencies at marketed dose",
        "source": "https://www.fda.gov/media/193168/download",
        "interpretation": (
            "At 1.2e14 vg: hepatic-enzyme increase 6/75 (8%), thrombocytopenia 3/75 (4%), "
            "sensory disturbance 2/75 (3%); two ALT values reached 20x ULN."
        ),
    },
    "spinraza_it_rna": {
        "product": "SPINRAZA (nusinersen)",
        "capsid": "none; antisense oligonucleotide",
        "route": "intrathecal",
        "dose_value": 12.0,
        "dose_unit": "mg/5 mL per low-dose administration",
        "total_dose_vg_for_model": None,
        "evidence_kind": "marketed CNS RNA-drug procedure context",
        "source": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/209531s016lbl.pdf",
        "interpretation": "Useful for lumbar-puncture procedure and monitoring context; mass cannot be converted to vg.",
    },
    "qalsody_it_rna": {
        "product": "QALSODY (tofersen)",
        "capsid": "none; antisense oligonucleotide",
        "route": "intrathecal",
        "dose_value": 100.0,
        "dose_unit": "mg/15 mL per administration",
        "total_dose_vg_for_model": None,
        "evidence_kind": "marketed CNS RNA-drug procedure context",
        "source": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/215887s000lbl.pdf",
        "interpretation": "Useful for lumbar-puncture procedure and route-specific hazards; mass cannot be converted to vg.",
    },
    "roctavian_iv": {
        "product": "ROCTAVIAN (valoctocogene roxaparvovec-rvox)",
        "capsid": "AAV5",
        "route": "iv",
        "dose_value": 6.0e13,
        "dose_unit": "vg/kg",
        "total_dose_vg_for_model": human_model.BODY_WEIGHT_KG * 6.0e13,
        "evidence_kind": "approved AAV5 systemic-dose context with required liver monitoring",
        "source": "https://www.fda.gov/media/169937/download",
        "interpretation": (
            "Same capsid and route for IV AAV5, but a different liver-directed cassette, product, "
            "manufacturing process and population. Historical approved-dose context is not a NOAEL."
        ),
    },
    "aav8_hemophilia_iv": {
        "product": "scAAV2/8-LP1-hFIXco long-term hemophilia B study",
        "capsid": "AAV8",
        "route": "iv",
        "dose_value": 2.0e12,
        "dose_unit": "vg/kg",
        "total_dose_vg_for_model": human_model.BODY_WEIGHT_KG * 2.0e12,
        "evidence_kind": "clinical AAV8 systemic dose with 13-year follow-up",
        "source": "https://doi.org/10.1056/NEJMoa2414783",
        "interpretation": (
            "Same capsid family and systemic route for IV AAV8; small hemophilia cohort and a "
            "liver-specific cassette do not establish a general-organ safety threshold."
        ),
    },
    "lk03_hemophilia_iv": {
        "product": "SPK-8011 AAV-LK03 hemophilia A phase 1-2 study",
        "capsid": "AAV-LK03",
        "route": "iv",
        "dose_value": 2.0e12,
        "dose_unit": "vg/kg",
        "total_dose_vg_for_model": human_model.BODY_WEIGHT_KG * 2.0e12,
        "evidence_kind": "clinical LK03 systemic dose context",
        "source": "https://doi.org/10.1056/NEJMoa2104205",
        "interpretation": (
            "Same engineered capsid and route, but immune events and glucocorticoid use occurred; "
            "the dose is a clinical context rather than a no-effect threshold."
        ),
    },
    "glybera_im": {
        "product": "GLYBERA (alipogene tiparvovec; historical EU authorization)",
        "capsid": "AAV1",
        "route": "intramuscular",
        "dose_value": 1.0e12,
        "dose_unit": "gc/kg",
        "total_dose_vg_for_model": human_model.BODY_WEIGHT_KG * 1.0e12,
        "evidence_kind": "human intramuscular AAV route context; authorization later withdrawn",
        "source": "https://www.ema.europa.eu/en/documents/all-authorised-presentations/glybera-epar-all-authorised-presentations_en.pdf",
        "interpretation": (
            "Route-matched but cross-capsid and historical. Genome copies and vector genomes are "
            "treated as an approximate modeling bridge, not an assay-equivalent conversion."
        ),
    },
    "aav2_cf_aerosol": {
        "product": "tgAAVCF phase 2 cystic-fibrosis aerosol study",
        "capsid": "AAV2",
        "route": "inhaled",
        "dose_value": 1.0e13,
        "dose_unit": "DNase-resistant particles/administration",
        "total_dose_vg_for_model": 1.0e13,
        "evidence_kind": "human route- and capsid-matched repeated aerosol dose context",
        "source": "https://doi.org/10.1378/chest.126.2.509",
        "interpretation": (
            "Three administrations were studied 30 days apart. Particle and vg assays, cassette, "
            "nebulizer deposition and lung disease prevent direct threshold equivalence."
        ),
    },
    "aavrh10_intracerebral": {
        "product": "AAVrh.10-SGSH-IRES-SUMF1 phase 1-2 study",
        "capsid": "AAVrh.10",
        "route": "intracerebral",
        "dose_value": 7.2e11,
        "dose_unit": "vg/patient",
        "total_dose_vg_for_model": None,
        "evidence_kind": "human same-capsid local-brain procedure context",
        "source": "https://doi.org/10.1089/hum.2013.238",
        "interpretation": (
            "Twelve direct parenchymal deposits with immunosuppression are not equivalent to IT, "
            "ICM or ICV administration; retained as qualitative AAVrh.10 evidence only."
        ),
    },
}


SCENARIO_CONTEXT = {
    "iv": {
        "aav5": ("roctavian_iv", "approved-same-capsid-route", 0.75),
        "aav8": ("aav8_hemophilia_iv", "clinical-same-capsid-route", 0.75),
        "aav9": ("zolgensma_iv", "approved-same-capsid-route", 0.75),
        "lk03": ("lk03_hemophilia_iv", "clinical-same-capsid-route", 0.75),
    },
    "intrathecal": {
        "aav9": ("itvisma_it", "approved-same-capsid-route", 0.75),
    },
    "intramuscular": {},
    "intracisternal": {
        "aav9": ("itvisma_it", "clinical-closest-csf-route", 0.10),
    },
    "intracerebroventricular": {
        "aav9": ("itvisma_it", "clinical-closest-csf-route", 0.10),
    },
    "inhaled": {
        "aav2": ("aav2_cf_aerosol", "clinical-same-capsid-route", 0.75),
    },
}

DEFAULT_ROUTE_CONTEXT = {
    "iv": ("aav8_hemophilia_iv", "clinical-same-route-cross-capsid", 0.10),
    "intrathecal": ("itvisma_it", "approved-same-route-cross-capsid", 0.10),
    "intramuscular": ("glybera_im", "historical-same-route-cross-capsid", 0.10),
    "intracisternal": ("itvisma_it", "clinical-closest-csf-route-cross-capsid", 0.03),
    "intracerebroventricular": ("itvisma_it", "clinical-closest-csf-route-cross-capsid", 0.03),
    "inhaled": ("aav2_cf_aerosol", "clinical-same-route-cross-capsid", 0.10),
}

MOUSE_RESTRICTED_CAPSIDS = {"php-eb", "cap-b10"}
EPISOME_HALF_LIFE_DAYS = {
    "brain": 365.0,
    "heart": 300.0,
    "kidney": 90.0,
    "liver": 1095.0,
    "lung": 180.0,
    "muscle": 1460.0,
    "rest": 180.0,
    "spleen": 120.0,
}


RISK_DOMAINS = {
    "liver": {
        "organs": ("liver",),
        "hazard": "Hepatotoxicity and capsid/transgene-directed immune injury",
        "monitoring": "ALT, AST, total bilirubin, albumin, PT/INR; liver histology and cytokines preclinically",
    },
    "renal_tma": {
        "organs": ("kidney", "spleen"),
        "hazard": "Complement-associated thrombotic microangiopathy and acute kidney injury",
        "monitoring": "Platelets, hemoglobin, creatinine, urinalysis, LDH/haptoglobin and complement markers",
    },
    "cardiac": {
        "organs": ("heart",),
        "hazard": "Cardiac injury signal",
        "monitoring": "Troponin and cardiac assessment when biologically or clinically indicated",
    },
    "cns_drg": {
        "organs": ("brain",),
        "hazard": "DRG or local CNS neurotoxicity, especially for CSF/local routes",
        "monitoring": "Neurological/sensory endpoints, MRI, CSF markers and DRG histopathology in relevant animals",
    },
    "systemic_immune": {
        "organs": ("liver", "spleen", "lung", "rest"),
        "hazard": "Innate/adaptive immune activation and loss of efficacy",
        "monitoring": "Pre-existing NAb, anti-capsid antibody/T cells, cytokines, complement and infection status",
    },
}


def parent_auc(capsid: dict[str, Any]) -> dict[str, float]:
    result: dict[str, float] = {}
    for region in capsid["regions"].values():
        parent = region["parent_organ"]
        result[parent] = result.get(parent, 0.0) + float(region["auc_isf_amount_vg_h"])
    return result


def parent_peak_episome(capsid: dict[str, Any]) -> dict[str, float]:
    """Aggregate regional peak episome signals into parent-organ efficacy inputs."""
    result: dict[str, float] = {}
    for region in capsid["regions"].values():
        parent = region["parent_organ"]
        result[parent] = result.get(parent, 0.0) + float(region["peak_episome_vg"])
    return result


def scenario_context(route_id: str, capsid_id: str) -> tuple[str, str, float]:
    """Choose the nearest disclosed dose context and a transparent evidence haircut."""
    reference_id, inference_class, factor = SCENARIO_CONTEXT.get(route_id, {}).get(
        capsid_id,
        DEFAULT_ROUTE_CONTEXT[route_id],
    )
    if capsid_id in MOUSE_RESTRICTED_CAPSIDS:
        factor *= 0.30
        inference_class += "-mouse-restricted-capsid"
    return reference_id, inference_class, factor


def context_grade(inference_class: str) -> str:
    if inference_class.startswith("approved-same-capsid-route"):
        return "moderate"
    if inference_class.startswith("clinical-same-capsid-route"):
        return "low"
    return "exploratory"


def solve_reference_auc(route: str, dose_vg: float) -> dict[str, float]:
    """Re-solve AAV9 at a reference dose using the active human PBPK equations."""
    parents = {region.parent for region in human_model.REGIONS.values()}
    solution = human_model.solve_human_capsid(
        {parent: 1.0 for parent in parents},
        dose_vg=dose_vg,
        administration=route,
    )
    result = {parent: 0.0 for parent in parents}
    for region_id, region in human_model.REGIONS.items():
        isf = np.maximum(solution.y[human_model.IDX[f"A_{region_id}_isf"]], 0.0)
        result[region.parent] += float(np.trapezoid(isf, solution.t))
    return result


def relative_priority(ratio: float) -> str:
    if ratio >= 2.0:
        return "higher-priority"
    if ratio >= 0.75:
        return "reference-range"
    return "lower-modeled-exposure"


def benchmark_for_route(route_id: str) -> tuple[str | None, str]:
    if route_id == "iv":
        return "zolgensma_mouse_histology_signal", "same-route AAV9 product-specific nonclinical signal context"
    if route_id == "intrathecal":
        return "itvisma_it", "same-route marketed AAV9 context"
    if route_id in {"intracisternal", "intracerebroventricular"}:
        return "itvisma_it", "closest-route lumbar intrathecal AAV9 context; route mismatch"
    return None, "no route-matched quantitative AAV comparator"


def organ_evidence_scope(route_id: str, capsid_id: str, organ: str) -> tuple[str, str]:
    if route_id == "iv" and capsid_id == "aav9":
        if organ in {"heart", "liver"}:
            return "moderate", "below product-specific mouse signal exposure; not a human NOAEL"
        if organ == "brain":
            return "low", "below NHP DRG/CNS finding-dose context; no modeled DRG compartment"
        if organ in {"kidney", "spleen"}:
            return "low", "below marketed-dose context; TMA is immune-mediated and not an organ-AUC threshold"
        return "insufficient", "no organ-specific toxic exposure threshold"
    if route_id == "intrathecal" and capsid_id == "aav9":
        return "moderate", "same route/capsid marketed-product context, but current dose exceeds it"
    if route_id in {"intracisternal", "intracerebroventricular"}:
        return "insufficient", "closest-route comparison only; delivery geometry and local Cmax differ"
    if capsid_id != "aav9":
        return "insufficient", "cross-capsid AAV9 comparison only"
    return "insufficient", "no route-matched quantitative comparator"


def margin_assessment(margin: float | None, route_id: str, capsid_id: str, organ: str) -> str:
    grade, _ = organ_evidence_scope(route_id, capsid_id, organ)
    if margin is None:
        return "insufficient-evidence"
    if route_id == "intrathecal" and capsid_id == "aav9" and margin < 1.0:
        return "exceeds-marketed-product-exposure-context"
    if route_id in {"intracisternal", "intracerebroventricular"}:
        return "closest-route-context-only"
    if capsid_id != "aav9":
        return "cross-capsid-context-only"
    if margin >= 1.0 and grade in {"moderate", "low"}:
        return "below-signal-or-clinical-context-not-proven-safe"
    if margin < 1.0:
        return "exceeds-reference-context"
    return "insufficient-evidence"


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_report(
    payload: dict[str, Any],
    reference_auc: dict[str, dict[str, float]],
    margin_rows: list[dict[str, Any]],
) -> str:
    current_dose = float(payload["dose_vg_per_kg"])
    total_dose = float(payload["dose_vg"])
    iv_rows = [
        row for row in margin_rows
        if row["route_id"] == "iv" and row["capsid_id"] == "aav9"
    ]
    it_rows = [
        row for row in margin_rows
        if row["route_id"] == "intrathecal" and row["capsid_id"] == "aav9"
    ]
    iv_contextual = sum(row["evidence_grade"] != "insufficient" for row in iv_rows)
    iv_insufficient = len(iv_rows) - iv_contextual
    iv_margin_min = min(float(row["exposure_margin_reference_over_current"]) for row in iv_rows)
    iv_margin_max = max(float(row["exposure_margin_reference_over_current"]) for row in iv_rows)
    it_margin_min = min(float(row["exposure_margin_reference_over_current"]) for row in it_rows)
    it_margin_max = max(float(row["exposure_margin_reference_over_current"]) for row in it_rows)
    ratio_it = total_dose / REFERENCE_CASES["itvisma_it"]["dose_value"]
    ratio_kebilidi = total_dose / REFERENCE_CASES["kebilidi_intraputaminal"]["dose_value"]
    return f"""# AAV–SINEUP safety-margin assessment (2026-08-13)

## Executive conclusion

The current model does **not** prove that every organ is safe.  It supports a narrower, defensible statement for the default **IV AAV9** scenario: the administered dose of `{current_dose:.2e} vg/kg` is below the marketed ZOLGENSMA dose and below the lowest explicitly reported mouse cardiac-histology signal dose.  Re-solving the same PBPK equations gives organ ISF-AUC margins of `{iv_margin_min:.2f}-{iv_margin_max:.2f}` relative to the mouse signal-dose context.  This is supportive dose justification, not a human NOAEL.

Of the eight modeled organ groups, {iv_contextual}/8 have some route/capsid-matched contextual support (heart, liver, brain/DRG context, kidney and spleen); {iv_insufficient}/8 (lung, muscle and rest-of-body) lack organ-specific toxic exposure thresholds.  Therefore the correct claim is **“lower modeled exposure than selected AAV9 signal/clinical contexts, with residual product-specific risk”**, not “all organs proven safe.”

The current lumbar intrathecal scenario is not supportable at its present dose.  The model applies the same adult total dose `{total_dose:.2e} vg` to CSF delivery, which is `{ratio_it:.1f}x` the marketed ITVISMA dose (`1.2e14 vg`).  Re-solving at the ITVISMA dose gives reference/current organ-AUC margins of only `{it_margin_min:.3f}-{it_margin_max:.3f}`; equivalently, current modeled exposure is about `{1.0/it_margin_max:.1f}-{1.0/it_margin_min:.1f}x` the marketed-product context.

## Calculation

For organ `i`, route `r` and capsid `c`:

`Margin_i = AUC_i(reference AAV9, route-matched reference dose) / AUC_i(current route, current capsid)`

- `Margin >= 1` means the current modeled ISF amount AUC is below that reference exposure.
- `Margin < 1` means the current modeled exposure exceeds it.
- A margin is **not** a probability of safety. It inherits all comparator limitations: product, payload, capsid, manufacturing, age, disease, immunity, full/empty ratio and species.

For IV AAV9, the model is nearly dose proportional over this range: the organ margins to ZOLGENSMA's `1.1e14 vg/kg` dose are approximately `2.75`, and margins to the `7.9e13 vg/kg` mouse cardiac-histology signal dose are approximately `1.98`.

## CNS AAV case studies

| Case | Delivery and dose | What it supports | What it does not support |
|---|---|---|---|
| ZOLGENSMA | Single 60-min IV AAV9, `1.1e14 vg/kg` | Systemic AAV9 clinical dose context and monitoring domains | Universal AAV safe dose; label includes fatal/nonfatal liver failure, TMA, thrombocytopenia and troponin signals |
| ITVISMA | Single lumbar IT AAV9, `1.2e14 vg` in 3 mL over about 1-2 min | Best current same-route/same-capsid CNS comparator | No-risk threshold; trial and label include liver enzyme, sensory/DRG-like, platelet, TMA and troponin signals |
| KEBILIDI | AAV2, four stereotactic putaminal infusions, total `1.8e11 vg` in 0.32 mL | Approved local-parenchymal delivery precedent | Comparator for CSF or IV exposure; current model lacks an intraputaminal compartment |
| AAV9-miniSINEUP-GDNF | Mouse striatum, `7.0e9 vg` | Direct SINEUP payload precedent: about two-fold endogenous-protein increase for >=6 months with no reported weight/food-intake signal | Human safety, GLP toxicology or safety of another binding domain/target |

The current model's `2.8e15 vg` local-CNS total dose is also about `{ratio_kebilidi:,.0f}x` the KEBILIDI total dose, but this ratio is descriptive only because route, capsid and distribution geometry are different.

## RNA-drug delivery precedents

- SPINRAZA: lumbar intrathecal antisense oligonucleotide; current 2026 label includes `12 mg/5 mL` low-dose injections and a higher-dose regimen. It provides procedure, repeat-dose and monitoring precedent, not an AAV vg conversion.
- QALSODY: lumbar intrathecal antisense oligonucleotide, `100 mg/15 mL`, three loading doses then every 28 days. Myelitis/radiculitis, papilledema/intracranial pressure and aseptic meningitis show that route-specific safety cannot be inferred from AAV biodistribution alone.

## Payload-specific conclusion

The model caps haploinsufficient protein restoration at 100% from a 50% baseline, i.e. a maximum two-fold restoration.  This is mechanism-consistent with reported SINEUP increases (often about 1.5-3-fold) and with the approximately two-fold AAV9-miniSINEUP-GDNF mouse result.  It is not proof that the current plasmid is safe because the repository does not specify one final binding-domain sequence, target transcript/isoform, promoter, cassette sequence, CpG burden, ITR integrity or human off-target screen.

The plasmid claim should therefore be: **“designed for endogenous-transcript-dependent, bounded restoration; existing SINEUP data reduce but do not eliminate overexpression concern.”**  Required construct-specific checks are target-protein upper tolerance, RNA-seq/proteomics off-targets, innate RNA sensing, promoter leakiness, replication-competent AAV, full/empty ratio and long-term integration analysis.

## Recommended safety-model architecture

1. Keep capsid particles, vector genomes, episomes, SINEUP RNA and target protein as separate states.
2. Fit organ vector exposure with PBPK; do not map AUC directly to injury without data.
3. Use a Bayesian interval-censored logistic model for binary histopathology: `logit(P(injury)) = alpha_study + beta*log10(exposure) + route + species + capsid + payload-expression`.
4. Model longitudinal ALT/AST, platelets, creatinine, complement, troponin and neurofilament as continuous submodels; connect them to time-to-event hazards for clinical events.
5. Treat DRG separately from whole brain. The current brain compartment cannot validate DRG safety.
6. Propagate uncertainty and report `P(margin>1)` plus credible intervals, not a single green/red score.
7. Validate on a held-out dose group or second study before human extrapolation.

## Current decision

- **Proceed as a research hypothesis:** IV AAV9 at `4.0e13 vg/kg`, with an explicit ~2-fold margin to the product-specific mouse histology signal context and full liver/TMA/cardiac/DRG monitoring.
- **Do not claim clinical safety:** the evidence is not a product-specific NOAEL and recent marketed-product experience confirms serious AAV liver risk.
- **Do not use the current CSF dose:** reduce and re-solve lumbar IT around `1.2e14 vg total` or lower before interpreting CNS safety; ICM/ICV need route-specific large-animal calibration.
- **Do not claim the plasmid itself safe yet:** sequence- and target-specific evidence is missing.

## Primary sources

- ZOLGENSMA US prescribing information: https://www.fda.gov/media/126109/download
- ITVISMA US prescribing information: https://www.fda.gov/media/193168/download
- KEBILIDI US prescribing information: https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=9d6a6401-c6b5-4f29-af11-67707d249482
- SPINRAZA 2026 US prescribing information: https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/209531s016lbl.pdf
- QALSODY US prescribing information: https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/215887s000lbl.pdf
- AAV9-miniSINEUP-GDNF mouse study: https://pmc.ncbi.nlm.nih.gov/articles/PMC7000958/
- Intrathecal onasemnogene NHP DRG study: https://pmc.ncbi.nlm.nih.gov/articles/PMC9347375/
- AAV8/AAV9 PBPK precedent: https://doi.org/10.1016/j.xphs.2023.10.005
- Antibody/complement-associated AAV TMA study: https://pmc.ncbi.nlm.nih.gov/articles/PMC10760971/
- FDA 2025 ELEVIDYS liver-failure safety action: https://www.fda.gov/news-events/press-announcements/fda-approves-new-safety-warning-and-revised-indication-limits-use-elevidys-following-reports-fatal
"""


def main() -> None:
    payload = json.loads(HUMAN_RESULTS.read_text(encoding="utf-8"))
    routes = payload["administration_routes"]
    iv = next(route for route in routes if route["route_id"] == "iv")
    iv_aav9 = next(capsid for capsid in iv["capsids"] if capsid["capsid_id"] == "aav9")
    current_iv_aav9_auc = parent_auc(iv_aav9)

    scenario_capsids = {
        (route["route_id"], capsid["capsid_id"]): capsid
        for route in routes
        for capsid in route["capsids"]
    }
    max_peak_episome: dict[str, float] = {}
    for capsid in scenario_capsids.values():
        for organ, peak in parent_peak_episome(capsid).items():
            max_peak_episome[organ] = max(max_peak_episome.get(organ, 0.0), peak)

    reference_auc = {
        "zolgensma_iv": solve_reference_auc(
            "iv", REFERENCE_CASES["zolgensma_iv"]["total_dose_vg_for_model"]
        ),
        "zolgensma_mouse_histology_signal": solve_reference_auc(
            "iv", REFERENCE_CASES["zolgensma_mouse_histology_signal"]["total_dose_vg_for_model"]
        ),
        "itvisma_it": solve_reference_auc(
            "intrathecal", REFERENCE_CASES["itvisma_it"]["total_dose_vg_for_model"]
        ),
    }
    reference_auc_cache: dict[tuple[str, float], dict[str, float]] = {}

    def projected_context_auc(
        route_id: str,
        capsid_id: str,
        organ_auc: dict[str, float],
        reference_id: str,
        inference_class: str,
    ) -> dict[str, float]:
        case = REFERENCE_CASES[reference_id]
        dose_vg = float(case["total_dose_vg_for_model"])
        if inference_class.startswith(("approved-same-capsid-route", "clinical-same-capsid-route")):
            scale = dose_vg / float(payload["dose_vg"])
            return {organ: value * scale for organ, value in organ_auc.items()}
        reference_route = str(case["route"])
        if reference_route not in human_model.ADMINISTRATION_ROUTES:
            reference_route = route_id
        cache_key = (reference_route, dose_vg)
        if cache_key not in reference_auc_cache:
            reference_auc_cache[cache_key] = solve_reference_auc(reference_route, dose_vg)
        return reference_auc_cache[cache_key]

    screens = []
    margin_rows: list[dict[str, Any]] = []
    for route in routes:
        route_id = route["route_id"]
        for capsid in route["capsids"]:
            capsid_id = capsid["capsid_id"]
            organ_auc = parent_auc(capsid)
            organ_epi = parent_peak_episome(capsid)
            benchmark_id, inference_class, uncertainty_factor = scenario_context(route_id, capsid_id)
            comparator_scope = inference_class.replace("-", " ")
            scenario_reference = projected_context_auc(
                route_id,
                capsid_id,
                organ_auc,
                benchmark_id,
                inference_class,
            )
            evidence_grade = context_grade(inference_class)
            organ_margins: dict[str, Any] = {}
            for organ, current in sorted(organ_auc.items()):
                reference = scenario_reference.get(organ, 0.0)
                conservative_reference = reference * uncertainty_factor
                margin = reference / max(current, 1e-30)
                conservative_margin = conservative_reference / max(current, 1e-30)
                assessment = (
                    "below-conservative-contextual-upper-bound"
                    if conservative_margin >= 1.0
                    else "exceeds-conservative-contextual-upper-bound"
                )
                evidence_note = (
                    f"{inference_class}; {uncertainty_factor:.3g} evidence haircut applied to the "
                    "PBPK-projected disclosed-dose organ AUC"
                )
                organ_margins[organ] = {
                    "current_auc_isf_amount_vg_h": current,
                    "reference_auc_isf_amount_vg_h": reference,
                    "conservative_reference_auc_isf_amount_vg_h": conservative_reference,
                    "exposure_margin_reference_over_current": margin,
                    "conservative_margin_over_current": conservative_margin,
                    "assessment": assessment,
                    "evidence_grade": evidence_grade,
                    "evidence_note": evidence_note,
                }
                margin_rows.append({
                    "route_id": route_id,
                    "route": route["label"],
                    "capsid_id": capsid_id,
                    "capsid": capsid["capsid"],
                    "organ": organ,
                    "current_dose_vg": payload["dose_vg"],
                    "current_dose_vg_per_kg": payload["dose_vg_per_kg"],
                    "current_auc_isf_amount_vg_h": current,
                    "reference_id": benchmark_id,
                    "reference_product": REFERENCE_CASES[benchmark_id]["product"],
                    "reference_auc_isf_amount_vg_h": reference,
                    "conservative_reference_auc_isf_amount_vg_h": conservative_reference,
                    "exposure_margin_reference_over_current": margin,
                    "conservative_margin_over_current": conservative_margin,
                    "assessment": assessment,
                    "evidence_grade": evidence_grade,
                    "inference_class": inference_class,
                    "uncertainty_factor": uncertainty_factor,
                    "comparator_scope": comparator_scope,
                    "evidence_note": evidence_note,
                    "source_url": REFERENCE_CASES[benchmark_id]["source"],
                })

            conservative_upper_dose = float(payload["dose_vg"]) * min(
                value["conservative_margin_over_current"] for value in organ_margins.values()
            )
            efficacy_targets = {
                organ: {
                    "anchor_peak_episome_vg": peak,
                    "anchor_relative_episome": peak / max(max_peak_episome.get(organ, 0.0), 1e-30),
                    "episome_half_life_days": EPISOME_HALF_LIFE_DAYS[organ],
                    "baseline_protein_pct": 50.0,
                    "therapeutic_threshold_pct": 65.0,
                    "maximum_modeled_protein_pct": 100.0,
                    "sineup_activity_factor": 2.0,
                    "interpretation": (
                        "Dose-responsive SINEUP efficacy proxy normalized across route/capsid scenarios "
                        "for this target organ. A two-fold activity factor is an explicit sensitivity "
                        "assumption aligned to the approximate miniSINEUP-GDNF proof-of-concept gain; "
                        "it is not a calibrated human response."
                    ),
                }
                for organ, peak in organ_epi.items()
            }

            domains = {}
            for domain_id, domain in RISK_DOMAINS.items():
                current = sum(organ_auc.get(organ, 0.0) for organ in domain["organs"])
                reference = sum(current_iv_aav9_auc.get(organ, 0.0) for organ in domain["organs"])
                ratio = current / max(reference, 1e-30)
                domains[domain_id] = {
                    **domain,
                    "modeled_exposure_relative_to_current_iv_aav9": ratio,
                    "screening_priority": relative_priority(ratio),
                }
            if route["route_class"] == "csf":
                domains["cns_drg"]["screening_priority"] = "route-specific-priority"
                domains["cns_drg"]["route_note"] = (
                    "CSF delivery makes DRG/CNS assessment mandatory; whole-brain AUC cannot clear DRG risk."
                )

            conclusion = (
                "below-conservative-contextual-upper-bound-not-proven-safe"
                if conservative_upper_dose >= float(payload["dose_vg"])
                else "dose-reduction-required-for-conservative-contextual-bound"
            )
            screens.append({
                "route_id": route_id,
                "route": route["label"],
                "capsid_id": capsid_id,
                "capsid": capsid["capsid"],
                "dose_vg": payload["dose_vg"],
                "dose_vg_per_kg": payload["dose_vg_per_kg"],
                "primary_benchmark_id": benchmark_id,
                "comparator_scope": comparator_scope,
                "inference_class": inference_class,
                "uncertainty_factor": uncertainty_factor,
                "evidence_grade": evidence_grade,
                "disclosed_context_dose_vg": REFERENCE_CASES[benchmark_id]["total_dose_vg_for_model"],
                "conservative_contextual_upper_dose_vg": conservative_upper_dose,
                "dose_slider_min_vg": max(1.0e9, conservative_upper_dose / 1000.0),
                "dose_slider_max_vg": max(float(payload["dose_vg"]), conservative_upper_dose * 3.0),
                "organ_exposure_margins": organ_margins,
                "efficacy_targets": efficacy_targets,
                "risk_domains": domains,
                "conclusion": conclusion,
                "interpretation": (
                    "Research-use dose-window and assay-priority screen. The conservative bound is a "
                    "haircut applied to a PBPK-projected disclosed-dose context; it is not a NOAEL, "
                    "clinical recommendation or adverse-event probability."
                ),
            })

    evidence_rows = []
    for reference_id, case in REFERENCE_CASES.items():
        evidence_rows.append({
            "reference_id": reference_id,
            "product_or_study": case["product"],
            "capsid_or_modality": case["capsid"],
            "route": case["route"],
            "dose_value": case["dose_value"],
            "dose_unit": case["dose_unit"],
            "evidence_kind": case["evidence_kind"],
            "model_use": "quantitative reference" if case["total_dose_vg_for_model"] is not None else "case/procedure context only",
            "interpretation_limit": case["interpretation"],
            "source_url": case["source"],
        })

    output = {
        "schema_version": "3.0",
        "generated_on": "2026-08-17",
        "current_scenario": {
            "body_weight_kg": payload["body_weight_kg"],
            "dose_vg": payload["dose_vg"],
            "dose_vg_per_kg": payload["dose_vg_per_kg"],
            "important_route_warning": (
                "The same total dose is currently applied to every route. This is not acceptable for "
                "interpreting local-CNS safety without route-specific dose calibration."
            ),
        },
        "reference_cases": REFERENCE_CASES,
        "important_limits": [
            "A marketed dose is not a universal safe threshold.",
            "The 7.9e13 vg/kg mouse signal dose is not a NOAEL and is not an all-organ human threshold.",
            "The model has no calibrated exposure-to-ALT, TMA, troponin or DRG injury function.",
            "The model has no DRG compartment and no immune/complement state.",
            "Vector genome, capsid particle, empty/full ratio, promoter, payload, impurities, age and disease change risk.",
            "The final plasmid sequence, target isoform, promoter and human off-target screen are not present in the repository.",
            "Evidence haircuts (0.75, 0.10 or 0.03, with an additional mouse-restricted capsid penalty) are transparent modeling policy factors, not literature-derived probabilities.",
            "The two-fold SINEUP activity factor is a sensitivity assumption informed by the approximate miniSINEUP-GDNF proof-of-concept gain; it does not calibrate a human effective dose.",
            "A below-bound result means below a conservative contextual exposure projection, not clinically proven safe.",
        ],
        "dose_window_method": {
            "pk_scaling": "Organ AUC is scaled approximately linearly from the 2.8e15-vg anchor trajectory.",
            "reference_projection": "The active PBPK model projects organ AUC at a disclosed human or route-nearest dose context.",
            "evidence_haircut": "A transparent factor reduces that contextual AUC before it is used as an exploratory upper bound.",
            "efficacy": "A three-state Epi-SINEUP RNA-protein ODE is re-solved from the dose-scaled normalized peak episome signal.",
            "decision_rule": "A model-feasible window requires protein >=65% and all modeled organ AUCs below their conservative contextual upper bounds.",
        },
        "screens": screens,
    }
    OUTPUT.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    write_csv(MARGIN_CSV, list(margin_rows[0]), margin_rows)
    write_csv(FRONTEND_MARGIN_CSV, list(margin_rows[0]), margin_rows)
    write_csv(EVIDENCE_CSV, list(evidence_rows[0]), evidence_rows)
    write_csv(FRONTEND_EVIDENCE_CSV, list(evidence_rows[0]), evidence_rows)
    REPORT.write_text(build_report(payload, reference_auc, margin_rows), encoding="utf-8")
    print(f"Wrote {len(screens)} safety screens to {OUTPUT}")
    print(f"Wrote {len(margin_rows)} organ-margin rows to {MARGIN_CSV}")
    print(f"Wrote {len(evidence_rows)} evidence rows to {EVIDENCE_CSV}")
    print(f"Wrote assessment report to {REPORT}")


if __name__ == "__main__":
    main()
