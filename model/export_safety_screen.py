"""Generate an evidence-labelled AAV exposure and safety-priority screen.

This module deliberately does not classify a candidate as clinically safe.
It compares modeled organ exposure against the project's IV-AAV9 reference and
links each signal to monitoring endpoints reported in regulatory sources.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HUMAN_RESULTS = ROOT / "sineup-delivery-atlas" / "public" / "data" / "human-spatial-results.json"
OUTPUT = ROOT / "sineup-delivery-atlas" / "public" / "data" / "safety-screen.json"

REGULATORY_CONTEXT = {
    "zolgensma_label": {
        "url": "https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=68cd4f06-70e1-40d8-bedb-609ec0afa471",
        "context": "Product-specific AAV9 IV dose is 1.1e14 vg/kg. Serious liver injury, thrombocytopenia, TMA and elevated troponin are monitored. This is context, not a universal threshold.",
    },
    "fda_clinical_development": {
        "url": "https://www.fda.gov/media/167536/download",
        "context": "FDA examples link IV AAV to hepatotoxicity/TMA and intrathecal or intraparenchymal AAV to DRG or local neurotoxicity.",
    },
    "fda_neurodegenerative_guidance": {
        "url": "https://www.fda.gov/media/144886/download",
        "context": "Product purity, potency, delivery-device compatibility, unintended immune responses and unwanted expression require product-specific assessment.",
    },
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
        "monitoring": "Platelet count, hemoglobin, creatinine, urinalysis, LDH/haptoglobin and complement activation markers",
    },
    "cardiac": {
        "organs": ("heart",),
        "hazard": "Cardiac injury signal",
        "monitoring": "Troponin and cardiac assessment when biologically or clinically indicated",
    },
    "cns_drg": {
        "organs": ("brain",),
        "hazard": "DRG or local CNS neurotoxicity, especially for CSF/local routes",
        "monitoring": "Neurological examination, sensory endpoints, MRI, CSF markers and DRG histopathology in relevant animals",
    },
    "systemic_immune": {
        "organs": ("liver", "spleen", "lung", "rest"),
        "hazard": "Innate/adaptive immune activation and loss of efficacy",
        "monitoring": "Pre-existing NAb, anti-capsid antibody/T cells, cytokines and complement; infection status",
    },
}


def parent_auc(capsid: dict) -> dict[str, float]:
    result: dict[str, float] = {}
    for region in capsid["regions"].values():
        parent = region["parent_organ"]
        result[parent] = result.get(parent, 0.0) + float(region["auc_isf_amount_vg_h"])
    return result


def level(ratio: float) -> str:
    if ratio >= 2.0:
        return "higher-priority"
    if ratio >= 0.75:
        return "reference-range"
    return "lower-modeled-exposure"


def main() -> None:
    payload = json.loads(HUMAN_RESULTS.read_text(encoding="utf-8"))
    routes = payload["administration_routes"]
    iv = next(route for route in routes if route["route_id"] == "iv")
    reference_capsid = next(capsid for capsid in iv["capsids"] if capsid["capsid_id"] == "aav9")
    reference_auc = parent_auc(reference_capsid)
    screens = []
    for route in routes:
        for capsid in route["capsids"]:
            organ_auc = parent_auc(capsid)
            domains = {}
            for domain_id, domain in RISK_DOMAINS.items():
                current = sum(organ_auc.get(organ, 0.0) for organ in domain["organs"])
                reference = sum(reference_auc.get(organ, 0.0) for organ in domain["organs"])
                ratio = current / max(reference, 1e-30)
                domains[domain_id] = {
                    **domain,
                    "modeled_exposure_relative_to_iv_aav9": ratio,
                    "screening_priority": level(ratio),
                }
            if route["route_class"] == "csf":
                domains["cns_drg"]["screening_priority"] = "route-specific-priority"
                domains["cns_drg"]["route_note"] = "CSF delivery makes DRG/CNS safety assessment mandatory regardless of the whole-brain exposure ratio."
            screens.append({
                "route_id": route["route_id"],
                "route": route["label"],
                "capsid_id": capsid["capsid_id"],
                "capsid": capsid["capsid"],
                "dose_vg_per_kg": payload["dose_vg_per_kg"],
                "dose_context_ratio_to_zolgensma_label": payload["dose_vg_per_kg"] / 1.1e14,
                "risk_domains": domains,
                "interpretation": "Exposure-priority screen only; it does not establish clinical safety or a no-adverse-effect level.",
            })
    output = {
        "schema_version": "1.0",
        "reference": "IV AAV9 in the same 70 kg exploratory PBPK model",
        "regulatory_context": REGULATORY_CONTEXT,
        "important_limits": [
            "The model has no calibrated exposure-to-ALT, TMA, troponin or DRG injury function.",
            "Vector genome, capsid particle, empty/full ratio, promoter, payload, impurities, patient age and disease can change risk.",
            "Relative exposure is suitable for prioritizing assays, not approving a dose.",
        ],
        "screens": screens,
    }
    OUTPUT.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote {len(screens)} safety screens to {OUTPUT}")


if __name__ == "__main__":
    main()
