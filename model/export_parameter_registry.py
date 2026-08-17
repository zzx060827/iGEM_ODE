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
    "physiology": "ICRP Publication 89 (2002); Jones and Rowland-Yeo 2013; Shah and Betts 2012/2013",
    "mouse_physiology": "Brown et al. 1997; Davies and Morris 1993; Shah and Betts 2012/2013",
    "csf": "Damkier et al. 2013, doi:10.1152/physrev.00004.2013",
    "aavr": "Pillay et al. 2016, doi:10.1038/nature16465",
    "trafficking": "Nonnenmacher and Weber 2012, doi:10.1038/gt.2012.6; Riyad and Weber 2021, doi:10.1038/s41434-021-00243-z",
    "kidney": "Christensen and Birn 2002, doi:10.1038/nrm778 (mechanistic proximal-tubule endocytosis only)",
    "sineup": "Zucchelli et al. 2015, PMID:26259533 (SINEUP mechanism only)",
    "model": "Current model structural prior; requires project-specific calibration",
    "cell_entry": "Bartlett et al. 2000, doi:10.1128/JVI.74.6.2777-2785.2000",
    "human_liver_epi": "Fong et al. 2022, PMCID:PMC9018415",
    "human_muscle_epi": "Mueller et al. 2017, PMCID:PMC5374867",
}

SOURCE_URLS = {
    "mouse_pk": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11404148/; https://pmc.ncbi.nlm.nih.gov/articles/PMC7193641/",
    "nhp_pk": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7769048/",
    "pbpk": "https://doi.org/10.1016/j.xphs.2023.10.005",
    "physiology": "https://www.icrp.org/publication.asp?id=icrp%20publication%2089; https://pmc.ncbi.nlm.nih.gov/articles/PMC6890583/",
    "mouse_physiology": "https://pubmed.ncbi.nlm.nih.gov/8378254/; https://pmc.ncbi.nlm.nih.gov/articles/PMC3727051/",
    "csf": "https://doi.org/10.1152/physrev.00004.2013",
    "aavr": "https://doi.org/10.1038/nature16465",
    "trafficking": "https://doi.org/10.1038/gt.2012.6; https://doi.org/10.1038/s41434-021-00243-z",
    "kidney": "https://doi.org/10.1038/nrm778",
    "sineup": "https://pubmed.ncbi.nlm.nih.gov/26259533/",
    "cell_entry": "https://doi.org/10.1128/JVI.74.6.2777-2785.2000",
    "human_liver_epi": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9018415/",
    "human_muscle_epi": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5374867/",
    "zolgensma": "https://www.fda.gov/media/126109/download",
}

BASE_FIELDS = (
    "scope", "parameter", "description", "value", "unit", "evidence_type",
    "species", "source", "rationale", "confidence", "calibration_priority",
    "code_location", "evidence_scope", "literature_value",
    "comparison_to_literature", "source_url",
)

ZH_FIELDS = (
    "scope_zh", "description_zh", "unit_zh", "evidence_type_zh",
    "species_zh", "rationale_zh", "confidence_zh",
    "calibration_priority_zh", "evidence_scope_zh",
    "literature_value_zh", "comparison_to_literature_zh",
)

FIELDS = BASE_FIELDS + ZH_FIELDS

ORGAN_ZH = {
    "blood": "血液", "brain": "脑", "heart": "心脏", "kidney": "肾脏",
    "liver": "肝脏", "lung": "肺", "muscle": "肌肉", "rest": "其他组织",
    "spleen": "脾脏",
}

REGION_ZH = {
    "brain_frontal": "额叶皮层", "brain_parietal": "顶叶皮层",
    "brain_temporal": "颞叶皮层", "brain_occipital": "枕叶皮层",
    "brain_deep_gray": "深部灰质核团", "brain_cerebellum": "小脑",
    "brainstem_spinal": "脑干与脊髓", "heart": "心脏", "liver": "肝脏",
    "spleen": "脾脏", "kidney_left_cortex": "左肾皮质",
    "kidney_left_medulla": "左肾髓质", "kidney_right_cortex": "右肾皮质",
    "kidney_right_medulla": "右肾髓质", "muscle_injected_arm": "注射侧上臂肌肉",
    "muscle_contralateral_arm": "对侧上臂肌肉", "muscle_trunk": "躯干肌肉",
    "muscle_legs": "下肢肌肉", "gut": "胃肠道", "skin_adipose": "皮肤与脂肪",
    "bone_marrow": "骨与骨髓", "rest": "其他组织", "lung_left": "左肺",
    "lung_right": "右肺",
}

FIELD_ZH = {
    "flow_fraction": "血流分数", "vascular_ml": "血管容积",
    "isf_ml": "组织间液容积", "ps_ml_h": "通透性-表面积乘积",
    "kp": "组织-血浆分配系数", "internalization_half_life_h": "内化半衰期",
    "episome_half_life_days": "游离体半衰期", "weights": "三级深度权重",
    "cell_access": "细胞可及性", "persistence": "持久性修正",
    "depth_mm": "目标深度",
}

DESCRIPTION_ZH = {
    "Model parameter": "模型参数",
    "Administered vector-genome dose": "给药载体基因组剂量",
    "Effective compartment volume": "有效隔室容积",
    "Effective perfusion or exchange-flow parameter": "有效灌流或交换血流参数",
    "Permeability-surface area product": "通透性-表面积乘积",
    "Tissue-to-plasma partition coefficient": "组织-血浆分配系数",
    "Apparent circulating capsid loss rate": "循环衣壳表观损失速率",
    "Organ capsid-loss rate component": "器官衣壳损失速率分量",
    "RNA or protein first-order turnover": "RNA 或蛋白一阶周转速率",
    "BBB binding, trafficking or capacity parameter": "血脑屏障结合、转运或容量参数",
    "CNS-cell uptake, trafficking or expression parameter": "CNS 细胞摄取、转运或表达参数",
    "Kidney filtration, proximal-tubule uptake or trafficking parameter": "肾滤过、近端小管摄取或转运参数",
    "Cellular uptake, intracellular trafficking, immune or expression parameter": "细胞摄取、胞内转运、免疫或表达参数",
    "Administration or model-control setting": "给药或模型控制设置",
    "Reference adult body weight": "参考成人体重",
    "Exploratory administered dose": "探索性给药剂量",
    "Reference adult cardiac output": "参考成人心输出量",
    "Reference cardiac output": "参考心输出量",
    "Effective exchange-flow multiplier": "有效交换血流缩放系数",
    "Adult total CSF volume": "成人脑脊液总容积",
    "Adult CSF production": "成人脑脊液生成速率",
    "Equivalent first-order CSF turnover half-life": "等效一阶脑脊液周转半衰期",
    "SINEUP pharmacodynamic parameter": "SINEUP 药效动力学参数",
}

UNIT_ZH = {
    "mixed; see equation": "混合单位，见对应方程", "mL": "mL",
    "relative to AAV9": "相对于 AAV9", "mL/h": "mL/h", "h": "h",
    "dimensionless": "无量纲", "dimensionless or mm": "无量纲或 mm",
    "fraction of cardiac output": "心输出量分数", "day": "d", "1/h": "h^-1",
    "day, fraction or a.u.": "d、比例或任意单位", "h or categorical": "h 或类别变量",
    "see code": "见代码", "vg": "vg", "kg": "kg", "vg/kg": "vg/kg",
}

EVIDENCE_ZH = {
    "assumed": "假设值", "scaled": "缩放值",
    "literature-constrained prior": "文献约束先验", "fitted/partitioned": "拟合后分配",
    "direct NHP PET": "NHP PET 直接估计", "assumed/derived": "假设/推导值",
    "design/assumed": "设计输入/假设值", "design": "设计输入",
    "direct/derived": "直接证据/推导值", "scaled provisional": "暂定缩放值",
    "design reference": "设计参考值", "direct/scaled": "直接证据/缩放值",
    "direct": "直接文献值", "derived": "推导值",
}

SPECIES_ZH = {
    "reference human": "参考成人", "mouse-scale": "小鼠尺度",
    "generic mammalian cell": "通用哺乳动物细胞", "adult mouse": "成年小鼠",
    "generic CNS": "通用中枢神经系统", "mouse": "小鼠",
    "mouse + NHP": "小鼠与非人灵长类", "NHP": "非人灵长类",
    "preclinical + clinical ocular": "临床前与临床眼部证据", "preclinical": "临床前",
    "Ly6a-positive mouse only": "仅限 Ly6a 阳性小鼠",
    "human-hepatocyte prior": "人肝细胞先验", "generic target": "通用靶标",
    "human": "人", "model": "模型", "human projection": "人体投影",
    "mouse-to-human projection": "小鼠至人体投影", "mouse-scale demo": "小鼠尺度演示",
}

RATIONALE_ZH = {
    "Exploratory regional transport/transduction prior": "区域转运与转导的探索性先验",
    "Region-resolved reference geometry": "分区参考解剖几何",
    "Applied to PS and sqrt-scaled Kp; not a cross-study calibrated effect size": "用于缩放 PS，并以平方根缩放 Kp；不是跨研究校准的效应量",
    "Mechanistic step is literature-supported; present numerical value is an exploratory prior": "机制步骤有文献支持，但当前数值仍是探索性先验",
    "Produces disease-specific depth contrasts; not an anatomical diffusion fit": "用于产生疾病相关的深度差异，并非脑解剖扩散拟合",
    "Dual-entry mechanism is biologically motivated; intact-AAV rates remain provisional": "双入口机制具有生物学依据，但完整 AAV 的速率仍为暂定值",
    "Transferred intracellular chain with CNS-specific exploratory values": "沿用胞内转运链，并采用 CNS 特异的探索性数值",
    "Organ-scale values assembled for mass-conserving PBPK geometry": "为满足质量守恒的 PBPK 几何而汇总的器官尺度数值",
    "Organ log-linear half-life fit split 35% vascular and 65% ISF; split itself is structural": "器官对数线性半衰期拟合后按血管 35%、ISF 65% 分配；该比例本身是结构假设",
    "Physiological flow is multiplied by an explicit effective-exchange scale": "生理血流乘以显式的有效交换缩放系数",
    "Exploratory capsid partition prior": "探索性衣壳分配先验",
    "Mechanistic PBPK parameter; present value chosen for organ ordering and requires fit": "机制性 PBPK 参数；当前数值用于形成器官排序，仍需拟合",
    "BBB mechanism is represented explicitly; numerical rate is not fitted": "显式表示 BBB 机制，但数值速率尚未拟合",
    "Mechanism is supported; numerical value is a transparent design prior": "机制有文献支持，数值是公开透明的设计先验",
    "Calculated from 6 h RNA or 48 h protein half-life; target-specific values are needed": "由 6 h RNA 或 48 h 蛋白半衰期推导；需要靶标特异数值",
    "Defines the simulated scenario": "定义模拟场景",
    "Wang 2024 mouse 125I fit; provisional cross-species prior": "Wang 2024 小鼠 125I 数据拟合；暂作跨物种先验",
    "Retained for full auditability": "为完整审计而保留",
    "Scenario input, not a recommended dose": "场景输入，并非推荐剂量",
    "ln(2) divided by reported AAV9 circulation half-life": "由 ln(2) 除以文献报道的 AAV9 循环半衰期得到",
    "Defines a transparent reference adult, not an individual patient": "定义透明的参考成人，并不代表具体患者",
    "Scenario input selected below the cited product dose; product, payload and population differ": "场景输入低于所引用产品剂量，但产品、载荷和人群均不同",
    "5.5 L/min reference cardiac output": "采用 5.5 L/min 的参考心输出量",
    "Maintains continuity with mouse equation family; must not be called physiological flow": "保持与小鼠方程族一致，但不能称为真实生理血流",
    "Adult physiological reference": "成人生理参考值",
    "500 mL/day converted to hourly rate": "由 500 mL/d 换算为小时速率",
    "ln(2)/(production/volume)": "按 ln(2)/(生成速率/容积) 推导",
    "Ballon 2020 NHP PET": "来自 Ballon 2020 的 NHP PET 估计",
    "Left arm vein to right heart, lung, left heart, then systemic organs": "由左臂静脉进入右心、肺、左心，再到全身器官",
    "Lumbar CSF to spinal/cranial CSF and CNS ISF, with venous CSF drainage": "由腰段 CSF 进入脊髓与颅内 CSF/CNS ISF，并经静脉回流",
    "Injected deltoid depot to local muscle ISF with slower systemic escape": "三角肌注射仓进入局部肌肉 ISF，并缓慢泄漏至全身",
    "Cisterna magna to cranial CSF, favoring cerebellar and brainstem surfaces": "由枕大池进入颅内 CSF，偏向小脑与脑干表面",
    "Ventricular CSF input with greater access to periventricular and deep-gray regions": "由脑室 CSF 输入，更易到达脑室周围和深部灰质",
    "Airway depot to lung ISF with limited leakage into pulmonary venous blood": "气道给药仓进入肺 ISF，并有限进入肺静脉血",
}

CONFIDENCE_ZH = {
    "low": "低", "medium-low": "中低", "medium": "中",
    "exploratory": "探索性", "strong": "强", "low-medium": "低至中",
    "not applicable": "不适用", "high": "高", "medium-high": "中高",
}

PRIORITY_ZH = {"high": "高", "medium": "中", "low": "低"}

EVIDENCE_SCOPE_ZH = {
    "quantitative-direct": "直接定量证据",
    "quantitative-derived": "定量推导/拟合证据",
    "physiology-scaled": "生理学缩放证据",
    "comparative-prior": "同条件比较形成的先验",
    "mechanistic-only": "仅支持机制，不能直接约束数值",
    "clinical-context": "临床背景，不构成等效校准",
    "structural/design": "模型结构或设计输入",
}


def row(**kwargs: Any) -> dict[str, Any]:
    result = {field: "" for field in FIELDS}
    result.update(kwargs)
    return result


def scope_zh(scope: str) -> str:
    if scope == "mouse_pbpk":
        return "小鼠 PBPK"
    if scope == "human_global":
        return "人体全局参数"
    if scope == "human_capsid_pk":
        return "人体衣壳 PK"
    if scope == "sineup_pd":
        return "SINEUP 药效动力学"
    if scope.startswith("human_region:"):
        region_id = scope.split(":", 1)[1]
        return f"人体分区：{REGION_ZH.get(region_id, region_id)}"
    if scope.startswith("administration:"):
        route_id = scope.split(":", 1)[1]
        route = human.ADMINISTRATION_ROUTES.get(route_id)
        return f"给药途径：{route.label_zh if route else route_id}"
    if scope.startswith("capsid:"):
        return f"衣壳先验：{scope.split(':', 1)[1]}"
    if scope.startswith("cns_profile:"):
        return f"CNS 深度配置：{scope.split(':', 1)[1]}"
    return scope


def description_zh(item: dict[str, Any]) -> str:
    description = str(item["description"])
    if description in DESCRIPTION_ZH:
        return DESCRIPTION_ZH[description]
    match = re.fullmatch(r"Apparent early AAV9 capsid half-life: (.+)", description)
    if match:
        return f"AAV9 早期表观衣壳半衰期：{ORGAN_ZH.get(match.group(1), match.group(1))}"
    match = re.fullmatch(r"(.+): (flow_fraction|vascular_ml|isf_ml|ps_ml_h|kp|internalization_half_life_h|episome_half_life_days)", description)
    if match and str(item["scope"]).startswith("human_region:"):
        region_id = str(item["scope"]).split(":", 1)[1]
        return f"{REGION_ZH.get(region_id, match.group(1))}：{FIELD_ZH[match.group(2)]}"
    match = re.fullmatch(r"Relative (.+) prior for (.+)", description)
    if match:
        return f"{match.group(1)} 对{ORGAN_ZH.get(match.group(2), match.group(2))}的相对先验"
    match = re.fullmatch(r"Reduced-order CNS depth parameter: (.+)", description)
    if match:
        return f"降阶 CNS 深度参数：{FIELD_ZH.get(match.group(1), match.group(1))}"
    if str(item["scope"]).startswith("administration:"):
        route_id = str(item["scope"]).split(":", 1)[1]
        route = human.ADMINISTRATION_ROUTES.get(route_id)
        if route:
            return route.label_zh
    return description


def add_chinese_fields(item: dict[str, Any]) -> dict[str, Any]:
    """Append Chinese metadata without translating identifiers or citations."""
    item.update({
        "scope_zh": scope_zh(str(item["scope"])),
        "description_zh": description_zh(item),
        "unit_zh": UNIT_ZH.get(str(item["unit"]), str(item["unit"])),
        "evidence_type_zh": EVIDENCE_ZH.get(str(item["evidence_type"]), str(item["evidence_type"])),
        "species_zh": SPECIES_ZH.get(str(item["species"]), str(item["species"])),
        "rationale_zh": RATIONALE_ZH.get(str(item["rationale"]), str(item["rationale"])),
        "confidence_zh": CONFIDENCE_ZH.get(str(item["confidence"]), str(item["confidence"])),
        "calibration_priority_zh": PRIORITY_ZH.get(
            str(item["calibration_priority"]), str(item["calibration_priority"]),
        ),
        "evidence_scope_zh": EVIDENCE_SCOPE_ZH.get(
            str(item["evidence_scope"]), str(item["evidence_scope"]),
        ),
        "literature_value_zh": str(item["literature_value"]),
        "comparison_to_literature_zh": str(item["comparison_to_literature"]),
    })
    return item


def add_evidence_audit(item: dict[str, Any]) -> dict[str, Any]:
    """Classify whether a citation constrains a value or only its mechanism."""
    if not item["evidence_scope"]:
        evidence_type = str(item["evidence_type"])
        if evidence_type in {"direct", "direct NHP PET"}:
            item["evidence_scope"] = "quantitative-direct"
        elif evidence_type in {"direct/derived", "direct/scaled", "derived", "fitted/partitioned"}:
            item["evidence_scope"] = "quantitative-derived"
        elif evidence_type == "scaled":
            item["evidence_scope"] = "physiology-scaled"
        elif evidence_type == "literature-constrained prior":
            item["evidence_scope"] = "comparative-prior"
        elif evidence_type in {"design", "design/assumed", "design reference"}:
            item["evidence_scope"] = "structural/design"
        elif "Current model structural prior" in str(item["source"]):
            item["evidence_scope"] = "structural/design"
        else:
            item["evidence_scope"] = "mechanistic-only"
    if not item["source_url"]:
        source = str(item["source"])
        for key, citation in SOURCES.items():
            if citation and citation in source and key in SOURCE_URLS:
                item["source_url"] = SOURCE_URLS[key]
                break
        if not item["source_url"]:
            urls = re.findall(r"https?://[^;\s]+", source)
            item["source_url"] = "; ".join(urls)
    return item


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
        return dict(description="Administered vector-genome dose", unit="vg", evidence_type="assumed", species="mouse-scale demo", source=SOURCES["model"], rationale="Scenario input, not a recommended dose", confidence="low", calibration_priority="high", evidence_scope="structural/design")
    if name.startswith("V_"):
        return dict(description="Effective compartment volume", unit="mL", evidence_type="scaled", species="adult mouse", source=SOURCES["mouse_physiology"], rationale="Organ-scale values assembled for mass-conserving PBPK geometry", confidence="medium", calibration_priority="medium", evidence_scope="physiology-scaled", source_url=SOURCE_URLS["mouse_physiology"])
    if name.startswith("Q_") or name in {"CO", "Q_scale"}:
        unit = "mL/h" if name != "Q_scale" else "dimensionless"
        if name == "CO":
            return dict(description="Reference cardiac output", unit=unit, evidence_type="direct/scaled", species="adult mouse", source=SOURCES["mouse_physiology"], rationale="14 mL/min reference for an unanesthetized 25 g mouse", confidence="medium", calibration_priority="low", evidence_scope="physiology-scaled", literature_value="12-16 mL/min; reference mean 13.98 mL/min", comparison_to_literature="updated from 25 to 840 mL/h; previous value was an effective exchange rate", source_url=SOURCE_URLS["mouse_physiology"])
        return dict(description="Effective perfusion or exchange-flow parameter", unit=unit, evidence_type="scaled" if name != "Q_scale" else "assumed", species="adult mouse", source=SOURCES["pbpk"], rationale="Physiological flow is multiplied by an explicit effective-exchange scale", confidence="low", calibration_priority="high", evidence_scope="structural/design" if name == "Q_scale" else "physiology-scaled", source_url=SOURCE_URLS["pbpk"])
    if name.startswith("PS_"):
        return dict(description="Permeability-surface area product", unit="mL/h", evidence_type="assumed", species="adult mouse", source=SOURCES["pbpk"], rationale="Mechanistic PBPK parameter; present value chosen for organ ordering and requires fit", confidence="low", calibration_priority="high", evidence_scope="mechanistic-only", source_url=SOURCE_URLS["pbpk"])
    if name.startswith("Kp_"):
        organ = name.removeprefix("Kp_")
        observed = {"liver": "~50", "heart": "~10", "muscle": "~10", "brain": "~2", "lung": "~2", "kidney": "~2", "spleen": "~2", "rest": "not directly comparable"}.get(organ, "")
        return dict(description="Tissue-to-plasma partition coefficient", unit="dimensionless", evidence_type="assumed", species="adult mouse", source=SOURCES["pbpk"], rationale="Exploratory ISF exchange prior; Liu whole-tissue T/B ratios are calibration targets, not direct Kp replacements", confidence="low", calibration_priority="high", evidence_scope="mechanistic-only", literature_value=f"Liu 2024 whole-tissue T/B {observed}", comparison_to_literature="not directly comparable: observed whole tissue includes vascular, cellular and nuclear vector", source_url=SOURCE_URLS["pbpk"])
    if name == "k_clear_blood":
        return dict(description="Apparent circulating capsid loss rate", unit="1/h", evidence_type="direct/derived", species="mouse", source=SOURCES["mouse_pk"], rationale="ln(2) divided by reported AAV9 circulation half-life", confidence="medium", calibration_priority="medium", evidence_scope="quantitative-derived", literature_value="AAV9 blood half-life 5.0 h", comparison_to_literature="matches", source_url=SOURCE_URLS["mouse_pk"])
    if name.startswith("k_res_") or name.startswith("k_deg_isf_"):
        return dict(description="Organ capsid-loss rate component", unit="1/h", evidence_type="fitted/partitioned", species="mouse", source=SOURCES["mouse_pk"], rationale="Organ log-linear half-life fit split 35% vascular and 65% ISF; split itself is structural", confidence="medium-low", calibration_priority="high", evidence_scope="quantitative-derived", comparison_to_literature="total fitted organ loss matches; 35/65 decomposition remains unidentifiable", source_url=SOURCE_URLS["mouse_pk"])
    if name.startswith("k_deg_m") or name.startswith("k_deg_p") or name.endswith("_deg_m") or name.endswith("_deg_p"):
        return dict(description="RNA or protein first-order turnover", unit="1/h", evidence_type="assumed/derived", species="generic mammalian cell", source=SOURCES["model"], rationale="Calculated from 6 h RNA or 48 h protein half-life; target-specific values are needed", confidence="low", calibration_priority="high", evidence_scope="structural/design")
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
        if name.startswith("k_int") or name.endswith("_int"):
            metadata.update(
                source=f"{metadata['source']}; {SOURCES['cell_entry']}",
                source_url=f"{metadata.get('source_url', '')}; {SOURCE_URLS['cell_entry']}".strip("; "),
                literature_value="AAV2-HeLa cell-surface internalization t1/2 <10 min",
                comparison_to_literature="model rates are much slower; retained because serotype, cell type and in-vivo accessibility are not comparable",
                evidence_scope="mechanistic-only",
            )
        rows.append(row(scope="mouse_pbpk", parameter=name, value=value_text(value), code_location="model/ode1.0.py:make_params", **metadata))
    return rows


def human_rows() -> list[dict[str, Any]]:
    rows = [
        row(scope="human_global", parameter="BODY_WEIGHT_KG", description="Reference adult body weight", value=value_text(human.BODY_WEIGHT_KG), unit="kg", evidence_type="design reference", species="human", source=SOURCES["physiology"], rationale="Defines a transparent reference adult, not an individual patient", confidence="medium", calibration_priority="low", code_location="model/human_spatial_pbpk.py", source_url=SOURCE_URLS["physiology"]),
        row(scope="human_global", parameter="DOSE_VG_PER_KG", description="Exploratory administered dose", value=value_text(human.DOSE_VG_PER_KG), unit="vg/kg", evidence_type="assumed", species="human projection", source="Clinical context only: ZOLGENSMA US label 1.1e14 vg/kg; not a safety equivalence", rationale="Scenario input selected below the cited product dose; product, payload and population differ", confidence="low", calibration_priority="high", code_location="model/human_spatial_pbpk.py", evidence_scope="clinical-context", literature_value="ZOLGENSMA label: 1.1e14 vg/kg", comparison_to_literature="model scenario is 0.364x label dose; product, payload, age and indication differ", source_url=SOURCE_URLS["zolgensma"]),
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
            extra: dict[str, Any] = {"source_url": SOURCE_URLS["physiology"]} if evidence == "scaled" else {}
            if field == "episome_half_life_days" and region.parent == "liver":
                source = SOURCES["human_liver_epi"]
                extra = dict(
                    evidence_scope="comparative-prior",
                    literature_value="transcriptionally competent liver episomes detected 2.6-4.1 y post-dose",
                    comparison_to_literature="updated from 120 to 1095 d; conservative effective persistence prior, not a fitted half-life",
                    source_url=SOURCE_URLS["human_liver_epi"],
                )
                evidence = "literature-constrained prior"
            elif field == "episome_half_life_days" and region.parent == "muscle":
                source = SOURCES["human_muscle_epi"]
                extra = dict(
                    evidence_scope="comparative-prior",
                    literature_value="circular muscle episomes and ~0.5 vg/diploid genome persisted at least 4 y",
                    comparison_to_literature="updated from 365 to 1460 d; lower-bound effective persistence prior",
                    source_url=SOURCE_URLS["human_muscle_epi"],
                )
                evidence = "literature-constrained prior"
            rows.append(row(scope=f"human_region:{region_id}", parameter=field, description=f"{region.label}: {field}", value=value_text(value), unit=unit, evidence_type=evidence, species="reference human", source=source, rationale="Region-resolved reference geometry" if evidence == "scaled" else "Exploratory regional transport/transduction prior", confidence="medium-low" if evidence == "scaled" else "low", calibration_priority="medium" if evidence == "scaled" else "high", code_location="model/human_spatial_pbpk.py:REGIONS", **extra))
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
    rows = [add_chinese_fields(add_evidence_audit(item)) for item in rows]
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    # utf-8-sig writes a BOM so Chinese Excel/Windows reliably recognizes the
    # CSV as UTF-8. JSON and TeX remain BOM-free UTF-8.
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    OUTPUT_JSON.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    lines = [
        r"\begingroup\scriptsize\sloppy",
        r"\setlength{\tabcolsep}{2.5pt}",
        r"\begin{longtable}{>{\raggedright\arraybackslash}p{0.10\textwidth}>{\raggedright\arraybackslash}p{0.13\textwidth}>{\raggedright\arraybackslash}p{0.15\textwidth}>{\raggedright\arraybackslash}p{0.08\textwidth}>{\raggedright\arraybackslash}p{0.11\textwidth}>{\raggedright\arraybackslash}p{0.11\textwidth}>{\raggedright\arraybackslash}p{0.22\textwidth}}",
        r"\caption{完整中英双语模型参数注册表。evidence type 区分 direct、fitted、derived、scaled 与 assumed；assumed 不表示机制没有依据，而表示当前数值尚未由本项目数据校准。}\label{tab:all-parameters}\\",
        r"\toprule Scope / 范围 & Parameter & Description / 描述 & Value & Unit & Evidence / 证据 & Reference / rationale / 参考与理由 \\",
        r"\midrule\endfirsthead",
        r"\toprule Scope / 范围 & Parameter & Description / 描述 & Value & Unit & Evidence / 证据 & Reference / rationale / 参考与理由 \\",
        r"\midrule\endhead",
        r"\midrule\multicolumn{7}{r}{Continued on next page}\\\endfoot",
        r"\bottomrule\endlastfoot",
    ]
    for item in rows:
        reference = f"{item['source']}. {item['rationale']} / {item['rationale_zh']}"
        cells = [
            f"{item['scope']} / {item['scope_zh']}", item["parameter"],
            f"{item['description']} / {item['description_zh']}", item["value"],
            f"{item['unit']} / {item['unit_zh']}",
            f"{item['evidence_type']} / {item['evidence_type_zh']}", reference,
        ]
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
