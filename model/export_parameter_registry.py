"""Export every active model parameter with provenance and calibration status.

The register separates a documented biological mechanism from a directly
measured numerical value. A parameter may therefore cite a mechanistic review
and still be labelled ``assumed`` when its present value has not been fitted.
"""

from __future__ import annotations

import csv
import json
import re
import runpy
from pathlib import Path
from typing import Any

import human_spatial_pbpk as human
import export_delivery_design_space as atlas


MODEL_DIR = Path(__file__).resolve().parent
OUTPUT_CSV = MODEL_DIR / "data" / "model_parameter_registry.csv"
OUTPUT_JSON = MODEL_DIR / "data" / "model_parameter_registry.json"
OUTPUT_TEX = MODEL_DIR.parents[0] / "docs" / "latex" / "generated_parameter_table.tex"

SOURCES = {
    "mouse_pk": "Wang et al. 2024, doi:10.1016/j.omtm.2024.101326; Seo et al. 2020, PMCID:PMC7193641",
    "nhp_pk": "Ballon et al. 2020, doi:10.1089/hum.2020.116",
    "pbpk": "Liu et al. 2024, doi:10.1016/j.xphs.2023.10.005",
    "physiology": "ICRP Publication 89 (2002), reference adult anatomy and physiology",
    "csf": "Damkier et al. 2013, doi:10.1152/physrev.00004.2013",
    "aavr": "Pillay et al. 2016, doi:10.1038/nature16465",
    "trafficking": "Nonnenmacher and Weber 2012, doi:10.1038/gt.2012.6; Riyad and Weber 2021, doi:10.1038/s41434-021-00243-z",
    "kidney": "Christensen and Birn 2002, doi:10.1038/nrm778 (mechanistic proximal-tubule endocytosis only)",
    "sineup": "Zucchelli et al. 2015, PMID:26259533 (SINEUP mechanism only)",
    "model": "Current model structural prior; requires project-specific calibration",
}

FIELDS = (
    "scope", "parameter", "description", "value", "unit", "evidence_type",
    "species", "source", "rationale", "confidence", "calibration_priority",
    "code_location",
)


def row(**kwargs: Any) -> dict[str, Any]:
    result = {field: "" for field in FIELDS}
    result.update(kwargs)
    return result


def value_text(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def tex_escape(value: Any) -> str:
    text = value_text(value)
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
        "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{",
        "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(character, character) for character in text)


def tex_identifier(value: Any) -> str:
    """Escape identifiers while allowing TeX to wrap long scopes and names."""
    text = value_text(value)
    rendered: list[str] = []
    cursor = 0
    for match in re.finditer(r"[A-Za-z0-9]{9,}", text):
        rendered.append(tex_escape(text[cursor:match.start()]))
        token = match.group(0)
        rendered.append(r"\allowbreak{}".join(tex_escape(token[i:i + 7]) for i in range(0, len(token), 7)))
        cursor = match.end()
    rendered.append(tex_escape(text[cursor:]))
    return (
        "".join(rendered)
        .replace(r"\_", r"\_\allowbreak{}")
        .replace(":", r":\allowbreak{}")
        .replace("/", r"/\allowbreak{}")
    )


def tex_table_value(value: Any) -> str:
    """Use report-level precision without changing machine-readable exports."""
    text = value_text(value)
    try:
        numeric = float(text)
    except (TypeError, ValueError):
        return tex_identifier(text)
    return tex_escape(f"{numeric:.3g}")


def mouse_metadata(name: str, value: Any) -> dict[str, str]:
    if name == "dose_vg":
        return dict(description="Administered vector-genome dose", unit="vg", evidence_type="assumed", species="mouse-scale demo", source=SOURCES["model"], rationale="Scenario input, not a recommended dose", confidence="low", calibration_priority="high")
    if name.startswith("V_"):
        return dict(description="Effective compartment volume", unit="mL", evidence_type="scaled", species="adult mouse", source=SOURCES["physiology"], rationale="Organ-scale values assembled for mass-conserving PBPK geometry", confidence="medium-low", calibration_priority="medium")
    if name.startswith("Q_") or name in {"CO", "Q_scale"}:
        unit = "mL/h" if name != "Q_scale" else "dimensionless"
        return dict(description="Effective perfusion or exchange-flow parameter", unit=unit, evidence_type="scaled" if name != "Q_scale" else "assumed", species="adult mouse", source=SOURCES["pbpk"], rationale="Physiological flow is multiplied by an explicit effective-exchange scale", confidence="low", calibration_priority="high")
    if name.startswith("PS_"):
        return dict(description="Permeability-surface area product", unit="mL/h", evidence_type="assumed", species="adult mouse", source=SOURCES["pbpk"], rationale="Mechanistic PBPK parameter; present value chosen for organ ordering and requires fit", confidence="low", calibration_priority="high")
    if name.startswith("Kp_"):
        return dict(description="Tissue-to-plasma partition coefficient", unit="dimensionless", evidence_type="assumed", species="adult mouse", source=SOURCES["pbpk"], rationale="Exploratory capsid partition prior", confidence="low", calibration_priority="high")
    if name == "k_clear_blood":
        return dict(description="Apparent circulating capsid loss rate", unit="1/h", evidence_type="direct/derived", species="mouse", source=SOURCES["mouse_pk"], rationale="ln(2) divided by reported AAV9 circulation half-life", confidence="medium", calibration_priority="medium")
    if name.startswith("k_res_") or name.startswith("k_deg_isf_"):
        return dict(description="Organ capsid-loss rate component", unit="1/h", evidence_type="fitted/partitioned", species="mouse", source=SOURCES["mouse_pk"], rationale="Organ log-linear half-life fit split 35% vascular and 65% ISF; split itself is structural", confidence="medium-low", calibration_priority="high")
    if name.startswith("k_deg_m") or name.startswith("k_deg_p") or name.endswith("_deg_m") or name.endswith("_deg_p"):
        return dict(description="RNA or protein first-order turnover", unit="1/h", evidence_type="assumed/derived", species="generic mammalian cell", source=SOURCES["model"], rationale="Calculated from 6 h RNA or 48 h protein half-life; target-specific values are needed", confidence="low", calibration_priority="high")
    if name.startswith("k_bbb") or name.startswith("Bmax_bbb"):
        return dict(description="BBB binding, trafficking or capacity parameter", unit="mixed; see equation", evidence_type="assumed", species="mouse-scale", source=SOURCES["trafficking"], rationale="BBB mechanism is represented explicitly; numerical rate is not fitted", confidence="low", calibration_priority="high")
    if name.startswith("k_cns") or name.startswith("Bmax_cns"):
        return dict(description="CNS-cell uptake, trafficking or expression parameter", unit="mixed; see equation", evidence_type="assumed", species="mouse-scale", source=SOURCES["trafficking"], rationale="Transferred intracellular chain with CNS-specific exploratory values", confidence="low", calibration_priority="high")
    if name.startswith("k_kidney") or name.startswith("k_pt_") or name.startswith("Bmax_pt") or name.startswith("k_glom") or name.startswith("k_filtrate") or name.startswith("k_urine"):
        return dict(description="Kidney filtration, proximal-tubule uptake or trafficking parameter", unit="mixed; see equation", evidence_type="assumed", species="mouse-scale", source=SOURCES["kidney"], rationale="Dual-entry mechanism is biologically motivated; intact-AAV rates remain provisional", confidence="low", calibration_priority="high")
    if name.startswith("k_") or name.startswith("R_") or name.startswith("Bmax") or name.startswith("EC50") or name in {"h", "h_kidney_tx", "h_cns_tx"}:
        return dict(description="Cellular uptake, intracellular trafficking, immune or expression parameter", unit="mixed; see equation", evidence_type="assumed", species="generic mammalian cell", source=f"{SOURCES['aavr']}; {SOURCES['trafficking']}", rationale="Mechanistic step is literature-supported; present numerical value is an exploratory prior", confidence="low", calibration_priority="high")
    if name in {"administration", "T_inf_h", "clearance_mode", "enable_apparent_decay"}:
        return dict(description="Administration or model-control setting", unit="h or categorical", evidence_type="design", species="model", source=SOURCES["model"], rationale="Defines the simulated scenario", confidence="not applicable", calibration_priority="low")
    return dict(description="Model parameter", unit="see code", evidence_type="assumed", species="model", source=SOURCES["model"], rationale="Retained for full auditability", confidence="low", calibration_priority="medium")


def mouse_rows() -> list[dict[str, Any]]:
    module = runpy.run_path(str(MODEL_DIR / "ode1.0.py"))
    params = module["make_params"]()
    rows = []
    for name, value in sorted(params.items()):
        metadata = mouse_metadata(name, value)
        rows.append(row(scope="mouse_pbpk", parameter=name, value=value_text(value), code_location="model/ode1.0.py:make_params", **metadata))
    return rows


def human_rows() -> list[dict[str, Any]]:
    rows = [
        row(scope="human_global", parameter="BODY_WEIGHT_KG", description="Reference adult body weight", value=value_text(human.BODY_WEIGHT_KG), unit="kg", evidence_type="design reference", species="human", source=SOURCES["physiology"], rationale="Defines a transparent reference adult, not an individual patient", confidence="medium", calibration_priority="low", code_location="model/human_spatial_pbpk.py"),
        row(scope="human_global", parameter="DOSE_VG_PER_KG", description="Exploratory administered dose", value=value_text(human.DOSE_VG_PER_KG), unit="vg/kg", evidence_type="assumed", species="human projection", source="Clinical context only: ZOLGENSMA US label 1.1e14 vg/kg; not a safety equivalence", rationale="Scenario input selected below the cited product dose; product, payload and population differ", confidence="low", calibration_priority="high", code_location="model/human_spatial_pbpk.py"),
        row(scope="human_global", parameter="CARDIAC_OUTPUT_ML_H", description="Reference adult cardiac output", value=value_text(human.CARDIAC_OUTPUT_ML_H), unit="mL/h", evidence_type="direct/scaled", species="human", source=SOURCES["physiology"], rationale="5.5 L/min reference cardiac output", confidence="medium", calibration_priority="low", code_location="model/human_spatial_pbpk.py"),
        row(scope="human_global", parameter="EFFECTIVE_FLOW_SCALE", description="Effective exchange-flow multiplier", value=value_text(human.EFFECTIVE_FLOW_SCALE), unit="dimensionless", evidence_type="assumed", species="human projection", source=SOURCES["model"], rationale="Maintains continuity with mouse equation family; must not be called physiological flow", confidence="low", calibration_priority="high", code_location="model/human_spatial_pbpk.py"),
        row(scope="human_global", parameter="CSF_TOTAL_VOLUME_ML", description="Adult total CSF volume", value=value_text(human.CSF_TOTAL_VOLUME_ML), unit="mL", evidence_type="direct", species="human", source=SOURCES["csf"], rationale="Adult physiological reference", confidence="high", calibration_priority="low", code_location="model/human_spatial_pbpk.py"),
        row(scope="human_global", parameter="CSF_PRODUCTION_ML_H", description="Adult CSF production", value=value_text(human.CSF_PRODUCTION_ML_H), unit="mL/h", evidence_type="direct/derived", species="human", source=SOURCES["csf"], rationale="500 mL/day converted to hourly rate", confidence="medium-high", calibration_priority="low", code_location="model/human_spatial_pbpk.py"),
        row(scope="human_global", parameter="CSF_ABSORPTION_HALF_LIFE_H", description="Equivalent first-order CSF turnover half-life", value=value_text(human.CSF_ABSORPTION_HALF_LIFE_H), unit="h", evidence_type="derived", species="human", source=SOURCES["csf"], rationale="ln(2)/(production/volume)", confidence="medium", calibration_priority="medium", code_location="model/human_spatial_pbpk.py"),
    ]
    for organ, value in sorted({"blood": human.BLOOD_CAPSID_HALF_LIFE_H, **human.CAPSID_HALF_LIFE_H}.items()):
        provenance = human.REFERENCE_HUMAN_AAV9_PROVENANCE[organ]
        direct_nhp = provenance == "Ballon 2020 NHP PET"
        rows.append(row(scope="human_capsid_pk", parameter=f"AAV9_HALF_LIFE_{organ}", description=f"Apparent early AAV9 capsid half-life: {organ}", value=value_text(value), unit="h", evidence_type="direct NHP PET" if direct_nhp else "scaled provisional", species="NHP" if direct_nhp else "mouse-to-human projection", source=SOURCES["nhp_pk"] if direct_nhp else SOURCES["mouse_pk"], rationale=provenance, confidence="medium" if direct_nhp else "low", calibration_priority="high", code_location="model/aav_parameter_evidence.py"))
    for region_id, region in human.REGIONS.items():
        fields = {
            "flow_fraction": (region.flow_fraction, "fraction of cardiac output", "scaled"),
            "vascular_ml": (region.vascular_ml, "mL", "scaled"),
            "isf_ml": (region.isf_ml, "mL", "scaled"),
            "ps_ml_h": (region.ps_ml_h, "mL/h", "assumed"),
            "kp": (region.kp, "dimensionless", "assumed"),
            "internalization_half_life_h": (region.internalization_half_life_h, "h", "assumed"),
            "episome_half_life_days": (region.episome_half_life_days, "day", "assumed"),
        }
        for field, (value, unit, evidence) in fields.items():
            source = SOURCES["physiology"] if evidence == "scaled" else SOURCES["model"]
            rows.append(row(scope=f"human_region:{region_id}", parameter=field, description=f"{region.label}: {field}", value=value_text(value), unit=unit, evidence_type=evidence, species="reference human", source=source, rationale="Region-resolved reference geometry" if evidence == "scaled" else "Exploratory regional transport/transduction prior", confidence="medium-low" if evidence == "scaled" else "low", calibration_priority="medium" if evidence == "scaled" else "high", code_location="model/human_spatial_pbpk.py:REGIONS"))
    for route_id, route in human.ADMINISTRATION_ROUTES.items():
        rows.append(row(scope=f"administration:{route_id}", parameter="infusion_duration_h", description=route.label, value=value_text(route.infusion_duration_h), unit="h", evidence_type="design/assumed", species="reference human", source=route.evidence_source, rationale=route.description, confidence="low-medium", calibration_priority="high", code_location="model/human_spatial_pbpk.py:ADMINISTRATION_ROUTES"))
    return rows


def design_rows() -> list[dict[str, Any]]:
    rows = []
    for capsid_id, capsid in atlas.CAPSID_PRIORS.items():
        for organ, value in capsid["tropism"].items():
            rows.append(row(scope=f"capsid:{capsid_id}", parameter=f"tropism_{organ}", description=f"Relative {capsid['label']} prior for {organ}", value=value_text(value), unit="relative to AAV9", evidence_type="literature-constrained prior", species=capsid["species"], source="; ".join([capsid["source"], *capsid.get("additional_sources", [])]), rationale="Applied to PS and sqrt-scaled Kp; not a cross-study calibrated effect size", confidence=capsid["evidence"], calibration_priority="high", code_location="model/export_delivery_design_space.py:CAPSID_PRIORS"))
    for profile_id, profile in atlas.CNS_PROFILES.items():
        for name, value in profile.items():
            rows.append(row(scope=f"cns_profile:{profile_id}", parameter=name, description=f"Reduced-order CNS depth parameter: {name}", value=value_text(value), unit="dimensionless or mm", evidence_type="assumed", species="generic CNS", source=SOURCES["model"], rationale="Produces disease-specific depth contrasts; not an anatomical diffusion fit", confidence="low", calibration_priority="high", code_location="model/export_delivery_design_space.py:CNS_PROFILES"))
    pd = {
        "sineup_rna_half_life_days": 0.25,
        "target_protein_half_life_days": 2.0,
        "k_sineup_tx": 4.0,
        "ec50_epi": 0.30,
        "ec50_sineup": 0.50,
        "max_translation_boost": 1.0,
        "therapeutic_threshold": 0.65,
    }
    for name, value in pd.items():
        rows.append(row(scope="sineup_pd", parameter=name, description="SINEUP pharmacodynamic parameter", value=value_text(value), unit="day, fraction or a.u.", evidence_type="assumed", species="generic target", source=SOURCES["sineup"], rationale="Mechanism is supported; numerical value is a transparent design prior", confidence="low", calibration_priority="high", code_location="model/export_delivery_design_space.py:solve_sineup_pd"))
    return rows


def main() -> None:
    rows = mouse_rows() + human_rows() + design_rows()
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    OUTPUT_JSON.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    lines = [
        r"\begingroup\scriptsize\sloppy",
        r"\setlength{\tabcolsep}{2.5pt}",
        r"\begin{longtable}{>{\raggedright\arraybackslash}p{0.10\textwidth}>{\raggedright\arraybackslash}p{0.13\textwidth}>{\raggedright\arraybackslash}p{0.15\textwidth}>{\raggedright\arraybackslash}p{0.08\textwidth}>{\raggedright\arraybackslash}p{0.11\textwidth}>{\raggedright\arraybackslash}p{0.11\textwidth}>{\raggedright\arraybackslash}p{0.22\textwidth}}",
        r"\caption{完整模型参数注册表。evidence type 区分 direct、fitted、derived、scaled 与 assumed；assumed 不表示机制没有依据，而表示当前数值尚未由本项目数据校准。}\label{tab:all-parameters}\\",
        r"\toprule Scope & Parameter & Description & Value & Unit & Evidence & Reference / rationale \\",
        r"\midrule\endfirsthead",
        r"\toprule Scope & Parameter & Description & Value & Unit & Evidence & Reference / rationale \\",
        r"\midrule\endhead",
        r"\midrule\multicolumn{7}{r}{Continued on next page}\\\endfoot",
        r"\bottomrule\endlastfoot",
    ]
    for item in rows:
        reference = f"{item['source']}. {item['rationale']}"
        cells = [item["scope"], item["parameter"], item["description"], item["value"], item["unit"], item["evidence_type"], reference]
        rendered = [tex_identifier(cells[0]), tex_identifier(cells[1])]
        rendered.extend([tex_identifier(cells[2]), tex_table_value(cells[3])])
        rendered.extend(tex_identifier(cell) for cell in cells[4:])
        lines.append(" & ".join(rendered) + r" \\")
    lines.extend([r"\end{longtable}", r"\endgroup", ""])
    OUTPUT_TEX.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_TEX.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {len(rows)} parameters to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
