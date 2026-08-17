"use client";

import { Activity, BookOpen, CircleAlert, Gauge, ShieldCheck, SlidersHorizontal } from "lucide-react";
import { useMemo, useState } from "react";
import safetyPayload from "../public/data/safety-screen.json";
import type { Language } from "./disease-data";

type EvidenceGrade = "moderate" | "low" | "exploratory";
type OrganMargin = {
  current_auc_isf_amount_vg_h: number;
  reference_auc_isf_amount_vg_h: number;
  conservative_reference_auc_isf_amount_vg_h: number;
  exposure_margin_reference_over_current: number;
  conservative_margin_over_current: number;
  assessment: string;
  evidence_grade: EvidenceGrade;
  evidence_note: string;
};
type EfficacyTarget = {
  anchor_peak_episome_vg: number;
  anchor_relative_episome: number;
  episome_half_life_days: number;
  baseline_protein_pct: number;
  therapeutic_threshold_pct: number;
  maximum_modeled_protein_pct: number;
  sineup_activity_factor: number;
  interpretation: string;
};
type SafetyScreen = {
  route_id: string;
  route: string;
  capsid_id: string;
  capsid: string;
  dose_vg: number;
  dose_vg_per_kg: number;
  primary_benchmark_id: string;
  comparator_scope: string;
  inference_class: string;
  uncertainty_factor: number;
  evidence_grade: EvidenceGrade;
  disclosed_context_dose_vg: number;
  conservative_contextual_upper_dose_vg: number;
  dose_slider_min_vg: number;
  dose_slider_max_vg: number;
  organ_exposure_margins: Record<string, OrganMargin>;
  efficacy_targets: Record<string, EfficacyTarget>;
  conclusion: string;
  interpretation: string;
};
type ReferenceCase = {
  product: string;
  capsid: string;
  route: string;
  dose_value: number;
  dose_unit: string;
  evidence_kind: string;
  source: string;
  interpretation: string;
};
type Payload = {
  current_scenario: { body_weight_kg: number; dose_vg: number; dose_vg_per_kg: number };
  reference_cases: Record<string, ReferenceCase>;
  important_limits: string[];
  screens: SafetyScreen[];
};

const data = safetyPayload as unknown as Payload;
const routeOrder = ["iv", "intrathecal", "intramuscular", "intracisternal", "intracerebroventricular", "inhaled"];
const capsidOrder = ["aav2", "aav5", "aav8", "aav9", "aavrh10", "php-eb", "cap-b10", "lk03"];
const routeZh: Record<string, string> = {
  iv: "外周静脉注射",
  intrathecal: "腰椎鞘内注射",
  intramuscular: "肌内注射",
  intracisternal: "枕大池注射",
  intracerebroventricular: "脑室内注射",
  inhaled: "吸入 / 气道给药",
};
const organName: Record<string, { zh: string; en: string }> = {
  brain: { zh: "脑 / CNS", en: "Brain / CNS" },
  heart: { zh: "心脏", en: "Heart" },
  kidney: { zh: "肾脏", en: "Kidney" },
  liver: { zh: "肝脏", en: "Liver" },
  lung: { zh: "肺", en: "Lung" },
  muscle: { zh: "肌肉", en: "Muscle" },
  rest: { zh: "其余组织", en: "Rest of body" },
  spleen: { zh: "脾脏", en: "Spleen" },
};

const copy = {
  zh: {
    eyebrow: "剂量—暴露—药效联动",
    title: "AAV 安全—有效窗口",
    description: "调节总 vg 剂量，实时更新器官 AUC、保守参考裕度与 SINEUP 蛋白恢复时间。绿色只表示低于模型参考上限，不代表临床安全。",
    route: "给药途径", capsid: "衣壳", target: "治疗靶器官", dose: "总注射剂量", perKg: "折算剂量",
    contextual: "低于保守参考暴露", effective: "达到模型疗效阈值", feasible: "模型内安全—有效重叠",
    unsafe: "超过保守参考暴露", ineffective: "未达到 65% 蛋白阈值", noWindow: "当前证据与参数下无重叠区间",
    upper: "保守暴露上限", lower: "最低模型有效剂量", window: "安全—有效候选区间",
    current: "当前选择", peak: "峰值蛋白恢复", duration: "有效持续时间", margin: "最小器官裕度",
    proteinCurve: "SINEUP 蛋白恢复轨迹", threshold: "65% 治疗阈值", days: "天",
    organTable: "器官暴露审计", organ: "器官", projected: "当前剂量预测 AUC", bound: "保守参考 AUC 上限", ratio: "裕度",
    below: "低于上限", above: "超过上限", evidence: "比较器与外推证据", disclosed: "披露剂量背景", factor: "证据折扣",
    inference: "外推类型", source: "查看一手来源", limits: "解释边界", anchor: "模型锚定剂量",
    referenceWarning: "该区间是研究用优先级筛查，不是 NOAEL、临床处方或不良事件概率。",
  },
  en: {
    eyebrow: "Dose–exposure–response coupling",
    title: "AAV safety–efficacy window",
    description: "Adjust total vector genomes to update organ AUC, conservative contextual margins, and SINEUP protein duration. Green means below a modeled reference bound, not clinically proven safe.",
    route: "Administration", capsid: "Capsid", target: "Target organ", dose: "Total injected dose", perKg: "Dose equivalent",
    contextual: "Below conservative exposure context", effective: "Meets modeled efficacy threshold", feasible: "Modeled safety–efficacy overlap",
    unsafe: "Above conservative exposure context", ineffective: "Does not reach 65% protein", noWindow: "No overlapping interval under current evidence and assumptions",
    upper: "Conservative exposure upper bound", lower: "Minimum modeled effective dose", window: "Candidate safety–efficacy interval",
    current: "Selected dose", peak: "Peak protein restoration", duration: "Effective duration", margin: "Minimum organ margin",
    proteinCurve: "SINEUP protein restoration trajectory", threshold: "65% therapeutic threshold", days: "days",
    organTable: "Organ exposure audit", organ: "Organ", projected: "Projected AUC at selected dose", bound: "Conservative reference AUC bound", ratio: "Margin",
    below: "Below bound", above: "Above bound", evidence: "Comparator and inference evidence", disclosed: "Disclosed dose context", factor: "Evidence haircut",
    inference: "Inference class", source: "Open primary source", limits: "Interpretation boundary", anchor: "Model anchor dose",
    referenceWarning: "This is a research prioritization window, not a NOAEL, clinical prescription, or adverse-event probability.",
  },
};

type PdResult = { peakPct: number; onsetDay: number | null; effectiveDuration: number; endDay: number | null; series: Array<[number, number]> };

function solveSineup(relativeEpi: number, halfLifeDays: number): PdResult {
  const dt = 0.25;
  const steps = Math.round(730 / dt);
  const kEpi = Math.log(2) / Math.max(halfLifeDays, 1e-12);
  const kS = Math.log(2) / 0.25;
  const kP = Math.log(2) / 2;
  const derivative = (e: number, s: number, p: number) => {
    const tx = 4 * Math.max(e, 0) / (0.30 + Math.max(e, 0) + 1e-30);
    const gain = Math.max(s, 0) / (0.50 + Math.max(s, 0) + 1e-30);
    const setpoint = Math.min(1, 0.5 * (1 + gain));
    return [-kEpi * Math.max(e, 0), tx - kS * Math.max(s, 0), kP * (setpoint - Math.max(p, 0))] as const;
  };
  let e = Math.max(relativeEpi, 0), s = 0, p = 0.5;
  let peak = p, onset: number | null = null, end: number | null = null;
  const series: Array<[number, number]> = [[0, p * 100]];
  for (let index = 1; index <= steps; index += 1) {
    const k1 = derivative(e, s, p);
    const k2 = derivative(e + dt * k1[0] / 2, s + dt * k1[1] / 2, p + dt * k1[2] / 2);
    const k3 = derivative(e + dt * k2[0] / 2, s + dt * k2[1] / 2, p + dt * k2[2] / 2);
    const k4 = derivative(e + dt * k3[0], s + dt * k3[1], p + dt * k3[2]);
    e = Math.max(0, e + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6);
    s = Math.max(0, s + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6);
    p = Math.max(0, p + dt * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6);
    const day = index * dt;
    peak = Math.max(peak, p);
    if (p >= 0.65) { if (onset === null) onset = day; end = day; }
    if (index % 24 === 0 || index === steps) series.push([day, p * 100]);
  }
  return { peakPct: peak * 100, onsetDay: onset, effectiveDuration: onset !== null && end !== null ? end - onset : 0, endDay: end, series };
}

function scientific(value: number, digits = 2) {
  if (!Number.isFinite(value) || value === 0) return "0";
  return value.toExponential(digits).replace("e+", "e");
}

function findEffectiveDose(screen: SafetyScreen, target: EfficacyTarget) {
  const effective = (dose: number) => solveSineup(target.anchor_relative_episome * target.sineup_activity_factor * dose / screen.dose_vg, target.episome_half_life_days).peakPct >= target.therapeutic_threshold_pct;
  const lowLimit = Math.min(screen.dose_slider_min_vg, 1e9);
  const highLimit = Math.max(screen.dose_slider_max_vg, screen.dose_vg * 10);
  if (!effective(highLimit)) return null;
  let low = Math.log10(lowLimit), high = Math.log10(highLimit);
  for (let index = 0; index < 30; index += 1) {
    const mid = (low + high) / 2;
    if (effective(10 ** mid)) high = mid; else low = mid;
  }
  return 10 ** high;
}

function suggestedDoseFor(screen: SafetyScreen, target: EfficacyTarget) {
  const effectiveLower = findEffectiveDose(screen, target);
  const safeUpper = screen.conservative_contextual_upper_dose_vg;
  return effectiveLower !== null && effectiveLower <= safeUpper
    ? Math.sqrt(effectiveLower * safeUpper)
    : Math.min(safeUpper, screen.disclosed_context_dose_vg);
}

function curvePath(series: Array<[number, number]>) {
  return series.map(([day, value], index) => `${index ? "L" : "M"}${8 + day / 730 * 484},${96 - (value - 45) / 60 * 78}`).join(" ");
}

export function SafetyDashboard({ language }: { language: Language }) {
  const t = copy[language];
  const [routeId, setRouteId] = useState("intrathecal");
  const [capsidId, setCapsidId] = useState("aav9");
  const [targetOrgan, setTargetOrgan] = useState("brain");
  const selectedScreen = data.screens.find((item) => item.route_id === routeId && item.capsid_id === capsidId) ?? data.screens[0];
  const target = selectedScreen.efficacy_targets[targetOrgan] ?? Object.values(selectedScreen.efficacy_targets)[0];
  const effectiveLower = useMemo(() => findEffectiveDose(selectedScreen, target), [selectedScreen, target]);
  const safeUpper = selectedScreen.conservative_contextual_upper_dose_vg;
  const hasWindow = effectiveLower !== null && effectiveLower <= safeUpper;
  const suggestedDose = suggestedDoseFor(selectedScreen, target);
  const [doseLog, setDoseLog] = useState(Math.log10(Math.max(suggestedDose, 1e9)));

  const dose = 10 ** doseLog;
  const doseRatio = dose / selectedScreen.dose_vg;
  const pd = useMemo(() => solveSineup(target.anchor_relative_episome * target.sineup_activity_factor * doseRatio, target.episome_half_life_days), [target, doseRatio]);
  const organs = Object.entries(selectedScreen.organ_exposure_margins).map(([organ, metric]) => ({
    organ,
    projected: metric.current_auc_isf_amount_vg_h * doseRatio,
    bound: metric.conservative_reference_auc_isf_amount_vg_h,
    margin: metric.conservative_margin_over_current / Math.max(doseRatio, 1e-30),
  })).sort((left, right) => left.margin - right.margin);
  const minimumMargin = organs[0]?.margin ?? 0;
  const belowContext = minimumMargin >= 1;
  const effective = pd.peakPct >= target.therapeutic_threshold_pct;
  const feasible = belowContext && effective;
  const reference = data.reference_cases[selectedScreen.primary_benchmark_id];
  const minLog = Math.log10(Math.max(selectedScreen.dose_slider_min_vg, 1e9));
  const maxLog = Math.log10(Math.max(selectedScreen.dose_slider_max_vg, safeUpper * 3, effectiveLower ?? 0));
  const windowLeft = effectiveLower === null ? 100 : Math.max(0, Math.min(100, (Math.log10(effectiveLower) - minLog) / Math.max(maxLog - minLog, 1e-9) * 100));
  const windowRight = Math.max(0, Math.min(100, (Math.log10(safeUpper) - minLog) / Math.max(maxLog - minLog, 1e-9) * 100));

  const selectRoute = (next: string) => {
    const nextTarget = next === "intramuscular" ? "muscle" : next === "inhaled" ? "lung" : next === "iv" ? "liver" : "brain";
    const nextScreen = data.screens.find((item) => item.route_id === next && item.capsid_id === capsidId) ?? selectedScreen;
    const nextTargetInput = nextScreen.efficacy_targets[nextTarget] ?? Object.values(nextScreen.efficacy_targets)[0];
    setRouteId(next); setTargetOrgan(nextTarget); setDoseLog(Math.log10(Math.max(suggestedDoseFor(nextScreen, nextTargetInput), 1e9)));
  };

  const selectCapsid = (next: string) => {
    const nextScreen = data.screens.find((item) => item.route_id === routeId && item.capsid_id === next) ?? selectedScreen;
    const nextTargetInput = nextScreen.efficacy_targets[targetOrgan] ?? Object.values(nextScreen.efficacy_targets)[0];
    setCapsidId(next); setDoseLog(Math.log10(Math.max(suggestedDoseFor(nextScreen, nextTargetInput), 1e9)));
  };

  const selectTarget = (next: string) => {
    const nextTargetInput = selectedScreen.efficacy_targets[next] ?? target;
    setTargetOrgan(next); setDoseLog(Math.log10(Math.max(suggestedDoseFor(selectedScreen, nextTargetInput), 1e9)));
  };

  return <section className="safety-page">
    <div className="safety-heading">
      <div><span className="eyebrow"><ShieldCheck size={15} />{t.eyebrow}</span><h1>{t.title}</h1><p>{t.description}</p></div>
      <div className={`safety-verdict ${feasible ? "feasible" : belowContext ? "ineffective" : "above"}`}>
        <ShieldCheck size={18} /><div><strong>{feasible ? t.feasible : belowContext ? t.ineffective : t.unsafe}</strong><span>{hasWindow ? `${scientific(effectiveLower as number)}–${scientific(safeUpper)} vg` : t.noWindow}</span></div>
      </div>
    </div>

    <div className="safety-controls">
      <label><span>{t.route}</span><select value={routeId} onChange={(event) => selectRoute(event.target.value)}>{routeOrder.map((id) => <option value={id} key={id}>{language === "zh" ? routeZh[id] : data.screens.find((item) => item.route_id === id)?.route}</option>)}</select></label>
      <label><span>{t.capsid}</span><select value={capsidId} onChange={(event) => selectCapsid(event.target.value)}>{capsidOrder.map((id) => <option value={id} key={id}>{data.screens.find((item) => item.capsid_id === id)?.capsid}</option>)}</select></label>
      <label><span>{t.target}</span><select value={targetOrgan} onChange={(event) => selectTarget(event.target.value)}>{Object.keys(selectedScreen.efficacy_targets).map((organ) => <option value={organ} key={organ}>{organName[organ]?.[language] ?? organ}</option>)}</select></label>
      <div className="safety-evidence-badge"><span>{t.evidence}</span><strong className={`grade-${selectedScreen.evidence_grade}`}>{selectedScreen.evidence_grade}</strong></div>
    </div>

    <div className="dose-console">
      <div className="dose-slider-card">
        <div className="dose-label"><div><SlidersHorizontal size={17} /><span>{t.dose}</span></div><strong>{scientific(dose, 3)} vg</strong><small>{t.perKg}: {scientific(dose / data.current_scenario.body_weight_kg, 2)} vg/kg</small></div>
        <input aria-label={t.dose} type="range" min={minLog} max={maxLog} step={0.005} value={doseLog} onChange={(event) => setDoseLog(Number(event.target.value))} />
        <div className="dose-scale"><span>{scientific(10 ** minLog)}</span><span>{t.anchor}: {scientific(selectedScreen.dose_vg)}</span><span>{scientific(10 ** maxLog)}</span></div>
        <div className="window-track" aria-label={t.window}><i className="safety-span" style={{ left: 0, width: `${windowRight}%` }} /><i className="efficacy-span" style={{ left: `${windowLeft}%`, width: `${Math.max(0, 100 - windowLeft)}%` }} />{hasWindow && <i className="overlap-span" style={{ left: `${windowLeft}%`, width: `${Math.max(0, windowRight - windowLeft)}%` }} />}<b style={{ left: `${Math.max(0, Math.min(100, (doseLog - minLog) / Math.max(maxLog - minLog, 1e-9) * 100))}%` }} /></div>
        <div className="window-legend"><span><i className="legend-safe" />{t.upper}: {scientific(safeUpper)}</span><span><i className="legend-effective" />{t.lower}: {effectiveLower ? scientific(effectiveLower) : "—"}</span><strong>{hasWindow ? t.window : t.noWindow}</strong></div>
      </div>
      <div className="dose-metrics">
        <div><span>{t.peak}</span><strong>{pd.peakPct.toFixed(1)}%</strong><small>{effective ? t.effective : t.ineffective}</small></div>
        <div><span>{t.duration}</span><strong>{pd.effectiveDuration.toFixed(0)} {t.days}</strong><small>{pd.onsetDay !== null ? `onset ${pd.onsetDay.toFixed(1)} d` : "—"}</small></div>
        <div><span>{t.margin}</span><strong>{minimumMargin.toFixed(2)}×</strong><small>{belowContext ? t.contextual : t.unsafe}</small></div>
      </div>
    </div>

    <div className="safety-grid">
      <section className="protein-panel">
        <div className="safety-panel-title"><Activity size={16} /><strong>{t.proteinCurve}</strong><span>{organName[targetOrgan]?.[language]}</span></div>
        <svg viewBox="0 0 508 112" role="img" aria-label={t.proteinCurve}><path d="M8 96H492" className="safety-chart-axis" /><line x1="8" x2="492" y1="70" y2="70" className="safety-threshold" /><text x="12" y="66">{t.threshold}</text><path d={curvePath(pd.series)} className="safety-protein-line" /></svg>
        <div className="curve-axis-labels"><span>0</span><span>365 d</span><span>730 d</span></div>
        <div className="pd-summary"><span>Epi t½ {target.episome_half_life_days.toFixed(0)} d</span><span>RNA t½ 0.25 d</span><span>Protein t½ 2 d</span><span>SINEUP activity {target.sineup_activity_factor.toFixed(1)}×</span><span>{target.interpretation}</span></div>
      </section>

      <section className="evidence-panel">
        <div className="safety-panel-title"><BookOpen size={16} /><strong>{t.evidence}</strong></div>
        <h2>{reference.product}</h2><p>{reference.evidence_kind}</p>
        <dl><div><dt>{t.disclosed}</dt><dd>{scientific(reference.dose_value)} {reference.dose_unit}</dd></div><div><dt>{t.factor}</dt><dd>{selectedScreen.uncertainty_factor.toFixed(3)}×</dd></div><div><dt>{t.inference}</dt><dd>{selectedScreen.inference_class}</dd></div></dl>
        <p className="evidence-limit">{reference.interpretation}</p><a href={reference.source} target="_blank" rel="noreferrer">{t.source}</a>
      </section>
    </div>

    <section className="organ-audit">
      <div className="safety-panel-title"><Gauge size={16} /><strong>{t.organTable}</strong><span>AUC (vg·h)</span></div>
      <div className="organ-audit-table" role="table">
        <div className="audit-head" role="row"><span>{t.organ}</span><span>{t.projected}</span><span>{t.bound}</span><span>{t.ratio}</span><span /></div>
        {organs.map((item) => <div className="audit-row" role="row" key={item.organ}><strong>{organName[item.organ]?.[language] ?? item.organ}</strong><code>{scientific(item.projected)}</code><code>{scientific(item.bound)}</code><code>{item.margin.toFixed(2)}×</code><span className={item.margin >= 1 ? "audit-below" : "audit-above"}>{item.margin >= 1 ? t.below : t.above}</span></div>)}
      </div>
    </section>

    <div className="safety-limits"><CircleAlert size={17} /><div><strong>{t.limits}</strong><p>{t.referenceWarning}</p><ul>{data.important_limits.slice(0, 6).map((limit) => <li key={limit}>{limit}</li>)}</ul></div></div>
  </section>;
}
