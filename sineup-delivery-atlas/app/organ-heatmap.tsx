"use client";

import {
  Activity,
  BarChart3,
  Clock3,
  ExternalLink,
  Info,
  Pause,
  Play,
} from "lucide-react";
import { useEffect, useMemo, useState } from "react";
import modelPayload from "../public/data/model-results.json";
import type { Language } from "./disease-data";

type OrganId = "brain" | "lung" | "heart" | "liver" | "spleen" | "kidney" | "muscle" | "rest";
type MetricId = "time" | "auc" | "share" | "peak" | "expression";
type HumanMetricId = "vascular" | "isf" | "episome" | "protein" | "auc" | "share" | "peak" | "tmax";
type ColorMappingMode = "absolute" | "relative";
type OrganMetrics = {
  auc_amount_vg_h: number;
  auc_concentration_vg_h_ml: number;
  peak_isf_amount_vg: number;
  peak_isf_concentration_vg_ml: number;
  peak_post_barrier_delivery_pct: number;
  tmax_h: number;
  exposure_share_pct: number;
  isf_amount_vg: number[];
  isf_concentration_vg_ml: number[];
  vascular_concentration_vg_ml: number[];
  peak_episome: number | null;
  episome_auc_vg_h: number | null;
  peak_expression: number | null;
  transduction_model: string;
  model_status: "ode-derived" | "exposure-only";
};
type CapsidHeatmap = {
  capsid_id: string;
  capsid: string;
  evidence: "strong" | "medium" | "exploratory";
  species: string;
  source: string;
  max_mass_balance_error: number;
  organs: Record<OrganId, OrganMetrics>;
};
type HeatmapPayload = {
  time_h: number[];
  organs: OrganId[];
  capsids: CapsidHeatmap[];
  reference_anatomy: string;
  reference_model: string;
  interpretation: string;
};

type HumanRegionMetrics = {
  label: string;
  parent_organ: string;
  route: "systemic" | "pulmonary";
  vascular_volume_ml: number;
  isf_volume_ml: number;
  blood_flow_ml_h: number;
  effective_exchange_flow_ml_h: number;
  auc_isf_amount_vg_h: number;
  auc_isf_concentration_vg_h_ml: number;
  peak_isf_concentration_vg_ml: number;
  peak_post_barrier_delivery_pct: number;
  tmax_isf_h: number;
  peak_episome_vg: number;
  peak_protein_au: number;
  tmax_protein_h: number;
  exposure_share_pct: number;
  vascular_concentration_vg_ml: number[];
  isf_concentration_vg_ml: number[];
  episome_vg: number[];
  protein_au: number[];
};
type HumanCapsid = {
  capsid_id: string;
  capsid: string;
  evidence: "strong" | "medium" | "exploratory";
  source_species: string;
  source: string;
  human_translation_note: string;
  max_mass_balance_error: number;
  regions: Record<string, HumanRegionMetrics>;
  circulation: Record<string, { volume_ml: number; concentration_vg_ml: number[] }>;
  route_compartments: Record<string, { volume_ml?: number; concentration_vg_ml?: number[]; amount_vg?: number[] }>;
};
type HumanAdministrationRoute = {
  route_id: "iv" | "intrathecal" | "intramuscular" | "intracisternal" | "intracerebroventricular" | "inhaled";
  label: string;
  label_zh: string;
  description: string;
  description_zh: string;
  infusion_duration_h: number;
  route_class: "systemic" | "csf" | "local_depot";
  evidence_source: string;
  capsids: HumanCapsid[];
};
type HumanSpatialPayload = {
  time_h: number[];
  max_time_days: number;
  dose_vg: number;
  dose_vg_per_kg: number;
  body_weight_kg: number;
  default_route_id: HumanAdministrationRoute["route_id"];
  administration_routes: HumanAdministrationRoute[];
  reference_model: string;
  physiology_status: string;
  cardiac_output_ml_h: number;
  effective_flow_scale: number;
  equation_family: string;
  state_count: number;
  region_ids: string[];
  circulation_ids: string[];
  interpretation: string;
};

const heatmap = (modelPayload as unknown as { organ_heatmap: HeatmapPayload }).organ_heatmap;

const copy = {
  zh: {
    eyebrow: "全身 ODE 投射",
    title: "AAV 衣壳器官递送热图",
    description: "每种衣壳重新求解同一套 PBPK。颜色来自器官 ISF 暴露、剂量份额或原生胞内转导状态。",
    capsid: "衣壳",
    metric: "显示指标",
    mapping: "颜色映射",
    absolute: "绝对值",
    relative: "相对值",
    administration: "给药途径",
    time: "时点 ISF 浓度",
    auc: "ISF 浓度 AUC",
    share: "器官暴露份额",
    peak: "峰值后屏障递送率",
    expression: "峰值表达输出",
    timeline: "给药后时间",
    play: "播放时间序列",
    pause: "暂停时间序列",
    humanProjection: "人体解剖投射",
    mouseModel: "成年小鼠尺度 PBPK",
    selectedOrgan: "当前器官",
    currentValue: "热图数值",
    peakConcentration: "峰值 ISF 浓度",
    tmax: "达到峰值",
    exposureShare: "器官暴露份额",
    postBarrier: "峰值后屏障递送",
    transduction: "转导状态",
    odeDerived: "原生胞内 ODE",
    exposureOnly: "仅 PBPK 暴露",
    rank: "器官排序",
    curve: "0–72 h ISF 浓度",
    evidence: "衣壳证据",
    source: "查看衣壳来源",
    logScale: "当前时间窗固定的跨衣壳对数色标",
    linearScale: "当前时间窗固定的跨衣壳线性色标",
    relativeScale: "相对色标：当前衣壳的最高区域定义为 100%",
    tmaxScale: "同一路径下的连续 ODE 达峰时间",
    earlyPeak: "早达峰",
    latePeak: "晚达峰",
    unavailable: "该器官尚无独立胞内表达模块",
    disclaimer: "这是 mouse-scale PBPK 的相对分布演示，人体轮廓仅用于解剖定位，不代表已校准的人体递送预测。",
    anatomyCredit: "解剖底图：Reactome，CC BY 4.0",
    humanMode: "人体多区域",
    mouseMode: "小鼠器官级",
    humanLoading: "正在载入人体多区域 ODE 数据…",
    vascular: "血管内浓度",
    isf: "区域 ISF 浓度",
    episome: "区域 episome",
    protein: "区域蛋白表达",
    timeScale: "时间尺度",
    logarithmicTime: "对数",
    linearTime: "线性",
    injectionRoute: "左臂静脉 → 右心 → 肺 → 左心 → 全身",
    routeDoseNote: "各途径使用相同总 vg 剂量，差异来自 ODE 输入位置和转运路径。",
    humanDisclaimer: "人体与小鼠使用同一套 Q–PS–Kp–清除–胞内转运方程，Q scale 已统一为 0.05；人体额外保留串联心肺循环、24 个区域及 CSF/IM depot。衣壳区域参数仍主要来自临床前数据，因此不是临床剂量预测。",
    humanCurve: "完整 0–730 天区域轨迹",
    translationNote: "跨物种说明",
    flow: "区域血流",
    physiology: "70 kg · 4×10¹³ vg/kg",
    earlyWindow: "早期递送 · 0–72 h",
    longWindow: "长期表达 · 1–730 d",
    earlyWindowHint: "血管与 ISF 中的完整 AAV 颗粒",
    longWindowHint: "入核 episome 与持续蛋白上调",
    persistenceClarification: "AAV 颗粒的清除不等于表达终止；长期窗口显示的是已经进入细胞核的 episome 和下游蛋白。",
  },
  en: {
    eyebrow: "Whole-body ODE projection",
    title: "AAV capsid organ-delivery heatmap",
    description: "Each capsid re-solves the same PBPK system. Color encodes organ ISF exposure, dose share, or native intracellular transduction.",
    capsid: "Capsid",
    metric: "Metric",
    mapping: "Color mapping",
    absolute: "Absolute",
    relative: "Relative",
    administration: "Administration",
    time: "ISF concentration at time",
    auc: "ISF concentration AUC",
    share: "Organ exposure share",
    peak: "Peak post-barrier delivery",
    expression: "Peak expression output",
    timeline: "Time after dose",
    play: "Play time course",
    pause: "Pause time course",
    humanProjection: "Human anatomy projection",
    mouseModel: "Adult mouse-scale PBPK",
    selectedOrgan: "Selected organ",
    currentValue: "Mapped value",
    peakConcentration: "Peak ISF concentration",
    tmax: "Time to peak",
    exposureShare: "Organ exposure share",
    postBarrier: "Peak post-barrier delivery",
    transduction: "Transduction status",
    odeDerived: "Native intracellular ODE",
    exposureOnly: "PBPK exposure only",
    rank: "Organ ranking",
    curve: "0–72 h ISF concentration",
    evidence: "Capsid evidence",
    source: "Open capsid source",
    logScale: "Fixed cross-capsid scale for this time window",
    linearScale: "Fixed cross-capsid scale for this time window",
    relativeScale: "Relative scale: highest region for this capsid = 100%",
    tmaxScale: "Continuous ODE peak time within this route",
    earlyPeak: "Earlier peak",
    latePeak: "Later peak",
    unavailable: "No organ-specific intracellular expression module yet",
    disclaimer: "This is a relative mouse-scale PBPK demonstration. The human outline provides anatomical orientation and is not a calibrated human prediction.",
    anatomyCredit: "Anatomy reference: Reactome, CC BY 4.0",
    humanMode: "Human multiregion",
    mouseMode: "Mouse organ-level",
    humanLoading: "Loading human multiregion ODE data…",
    vascular: "Vascular concentration",
    isf: "Regional ISF concentration",
    episome: "Regional episome",
    protein: "Regional protein output",
    timeScale: "Time scale",
    logarithmicTime: "Logarithmic",
    linearTime: "Linear",
    injectionRoute: "Left arm vein → right heart → lung → left heart → body",
    routeDoseNote: "Routes use the same total vg dose; differences arise from the ODE input location and transport path.",
    humanDisclaimer: "Human and mouse views use the same Q–PS–Kp, clearance, and intracellular-trafficking equations with Q scale aligned at 0.05. Human anatomy adds serial cardiopulmonary transit, 24 regions, and CSF/IM depots. Capsid priors remain predominantly preclinical, so this is not a clinical dose prediction.",
    humanCurve: "Complete 0–730 day regional trajectory",
    translationNote: "Cross-species note",
    flow: "Regional blood flow",
    physiology: "70 kg · 4×10¹³ vg/kg",
    earlyWindow: "Early delivery · 0–72 h",
    longWindow: "Long-term expression · 1–730 d",
    earlyWindowHint: "Intact AAV particles in vascular and ISF spaces",
    longWindowHint: "Nuclear episome and sustained protein restoration",
    persistenceClarification: "Clearance of intact AAV particles does not mean expression has ended. The long-term window follows nuclear episome and downstream protein.",
  },
};

const organName: Record<OrganId, Record<Language, string>> = {
  brain: { zh: "脑 / CNS", en: "Brain / CNS" },
  lung: { zh: "肺", en: "Lung" },
  heart: { zh: "心脏", en: "Heart" },
  liver: { zh: "肝脏", en: "Liver" },
  spleen: { zh: "脾脏", en: "Spleen" },
  kidney: { zh: "肾脏", en: "Kidney" },
  muscle: { zh: "骨骼肌", en: "Skeletal muscle" },
  rest: { zh: "其余组织", en: "Rest of body" },
};

const metricOrder: MetricId[] = ["time", "auc", "share", "peak", "expression"];
const colorStops = ["#edf1ef", "#9bc8b8", "#2f8b73", "#e4b743", "#c94f3d"];

function hexToRgb(hex: string) {
  const value = Number.parseInt(hex.slice(1), 16);
  return [(value >> 16) & 255, (value >> 8) & 255, value & 255];
}

function colorScale(value: number) {
  const clamped = Math.min(1, Math.max(0, value));
  const position = clamped * (colorStops.length - 1);
  const index = Math.min(colorStops.length - 2, Math.floor(position));
  const fraction = position - index;
  const start = hexToRgb(colorStops[index]);
  const end = hexToRgb(colorStops[index + 1]);
  const rgb = start.map((channel, i) => Math.round(channel + (end[i] - channel) * fraction));
  return `rgb(${rgb.join(",")})`;
}

function colorScaleBounds(values: number[], logarithmic: boolean) {
  const positive = values.filter((value) => Number.isFinite(value) && value > 0);
  if (!positive.length) {
    return { min: 0, max: 1, rawMin: 0, rawMax: 1, logarithmic };
  }
  const rawMin = Math.min(...positive);
  const rawMax = Math.max(...positive);
  return logarithmic
    ? {
        min: Math.log10(rawMin),
        max: Math.log10(rawMax),
        rawMin,
        rawMax,
        logarithmic,
      }
    : { min: 0, max: rawMax, rawMin: 0, rawMax, logarithmic };
}

function superscriptExponent(exponent: number) {
  const characters: Record<string, string> = {
    "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
    "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹",
  };
  return String(exponent).split("").map((character) => characters[character] ?? character).join("");
}

function compact(value: number, digits = 3) {
  if (!Number.isFinite(value)) return "—";
  if (value === 0) return (0).toFixed(digits);
  const magnitude = Math.abs(value);
  if (magnitude >= 1e4 || magnitude < 1e-3) {
    const exponent = Math.floor(Math.log10(magnitude));
    const coefficient = value / 10 ** exponent;
    return `${coefficient.toFixed(digits)} × 10${superscriptExponent(exponent)}`;
  }
  return value.toLocaleString("en-US", {
    minimumFractionDigits: digits,
    maximumFractionDigits: digits,
  });
}

function metricValue(metrics: OrganMetrics, metric: MetricId, timeIndex: number) {
  if (metric === "time") return metrics.isf_concentration_vg_ml[timeIndex];
  if (metric === "auc") return metrics.auc_concentration_vg_h_ml;
  if (metric === "share") return metrics.exposure_share_pct;
  if (metric === "peak") return metrics.peak_post_barrier_delivery_pct;
  return metrics.peak_expression;
}

function metricUnit(metric: MetricId) {
  if (metric === "time") return "vg-eq/mL";
  if (metric === "auc") return "vg-eq·h/mL";
  if (metric === "share" || metric === "peak") return "%";
  return "a.u.";
}

function sparklinePath(values: number[]) {
  const positive = values.filter((value) => value > 0);
  const minLog = Math.log10(Math.min(...positive, 1));
  const maxLog = Math.log10(Math.max(...positive, 1));
  return values.map((value, index) => {
    const x = 8 + (index / Math.max(values.length - 1, 1)) * 284;
    const transformed = value > 0 ? Math.log10(value) : minLog;
    const y = 74 - ((transformed - minLog) / Math.max(maxLog - minLog, 1e-9)) * 58;
    return `${index === 0 ? "M" : "L"}${x.toFixed(2)},${y.toFixed(2)}`;
  }).join(" ");
}

const reactomeAnatomyHref = "/reactome-male-body-organs.svg";
const reactomeAnatomyTransform = "translate(110 20) scale(5.14706 5.14444)";

function ReactomeAnatomyMasks({ prefix }: { prefix: string }) {
  const mask = (id: string, fragment: string) => (
    <mask id={`${prefix}-${id}`} key={id} maskUnits="userSpaceOnUse" x="0" y="0" width="600" height="1000" style={{ maskType: "alpha" }}>
      <use href={`${reactomeAnatomyHref}#${fragment}`} transform={reactomeAnatomyTransform} />
    </mask>
  );

  return (
    <defs>
      {mask("body-mask", "BG")}
      {mask("brain-mask", "R-ICO-013680")}
      {mask("lung-mask", "R-ICO-013935")}
      {mask("liver-mask", "R-ICO-012959")}
      {mask("kidney-mask", "R-ICO-012931")}
      <mask id={`${prefix}-gut-mask`} maskUnits="userSpaceOnUse" x="0" y="0" width="600" height="1000" style={{ maskType: "alpha" }}>
        <use href={`${reactomeAnatomyHref}#R-ICO-013406`} transform={reactomeAnatomyTransform} />
        <use href={`${reactomeAnatomyHref}#R-ICO-012904`} transform={reactomeAnatomyTransform} />
        <use href={`${reactomeAnatomyHref}#R-ICO-013151`} transform={reactomeAnatomyTransform} />
      </mask>
      <filter id={`${prefix}-selected-glow`} x="-30%" y="-30%" width="160%" height="160%">
        <feDropShadow dx="0" dy="0" stdDeviation="4" floodColor="#12211c" floodOpacity="0.9" />
      </filter>
    </defs>
  );
}

function AnatomyReferenceLayer({ opacity = 0.3 }: { opacity?: number }) {
  return (
    <image
      href={reactomeAnatomyHref}
      x="110"
      y="20"
      width="350"
      height="926"
      opacity={opacity}
      pointerEvents="none"
      preserveAspectRatio="xMidYMid meet"
    />
  );
}

function AnatomicalMap({
  values,
  selected,
  onSelect,
  language,
}: {
  values: Record<OrganId, number | null>;
  selected: OrganId;
  onSelect: (organ: OrganId) => void;
  language: Language;
}) {
  const tintFor = (organ: OrganId, opacity = 0.84) => ({
    fill: values[organ] === null ? "#d7dcda" : colorScale(values[organ] ?? 0),
    opacity,
    filter: selected === organ ? "url(#mouse-selected-glow)" : undefined,
  });
  const activate = (event: React.KeyboardEvent<SVGGElement>, organ: OrganId) => {
    if (event.key === "Enter" || event.key === " ") onSelect(organ);
  };
  const groupProps = (organ: OrganId) => ({
    role: "button" as const,
    tabIndex: 0,
    "aria-label": organName[organ][language],
    onClick: () => onSelect(organ),
    onKeyDown: (event: React.KeyboardEvent<SVGGElement>) => activate(event, organ),
    className: `anatomy-organ anatomy-hit-target ${selected === organ ? "selected" : ""}`,
  });

  return (
    <svg className="anatomy-map" viewBox="0 0 600 1000" aria-label="Interactive AAV organ distribution map">
      <ReactomeAnatomyMasks prefix="mouse" />
      <AnatomyReferenceLayer opacity={0.28} />

      <g aria-hidden="true" pointerEvents="none">
        <rect x="0" y="0" width="600" height="1000" mask="url(#mouse-body-mask)" style={tintFor("rest", 0.28)} />
        <g mask="url(#mouse-body-mask)" style={tintFor("muscle", 0.74)}>
          <path d="M147 170 Q188 145 219 185 L202 520 Q173 535 128 510 L147 170 Z M423 170 Q382 145 351 185 L368 520 Q397 535 442 510 L423 170 Z" />
          <path d="M205 525 Q232 505 270 548 L268 950 L205 950 Z M365 525 Q338 505 300 548 L302 950 L365 950 Z" />
        </g>
        <rect x="0" y="0" width="600" height="1000" mask="url(#mouse-brain-mask)" style={tintFor("brain")} />
        <rect x="0" y="0" width="600" height="1000" mask="url(#mouse-lung-mask)" style={tintFor("lung")} />
        <path d="M286 273 C273 255 249 269 252 294 C255 317 285 338 287 339 C291 336 319 317 321 293 C323 268 299 256 286 273 Z" style={tintFor("heart")} />
        <rect x="0" y="0" width="600" height="1000" mask="url(#mouse-liver-mask)" style={tintFor("liver")} />
        <path d="M352 344 C365 337 377 352 373 374 C369 392 357 402 348 391 C341 379 342 353 352 344 Z" style={tintFor("spleen")} />
        <rect x="0" y="0" width="600" height="1000" mask="url(#mouse-kidney-mask)" style={tintFor("kidney")} />
      </g>

      <g {...groupProps("rest")}>
        <circle cx="285" cy="70" r="58" />
        <path d="M238 150 Q285 126 332 150 L382 205 L365 520 L337 590 L335 930 L292 930 L282 610 L270 930 L227 930 L228 590 L202 520 L188 205 Z" />
        <path d="M194 206 L132 520 L91 510 L153 188 Z M376 206 L438 520 L479 510 L417 188 Z" />
      </g>
      <g {...groupProps("muscle")}>
        <path d="M160 205 Q136 262 119 360 L93 496 Q104 515 126 512 L169 383 Q185 300 201 224 Z" />
        <path d="M410 205 Q434 262 451 360 L477 496 Q466 515 444 512 L401 383 Q385 300 369 224 Z" />
        <path d="M225 575 Q211 700 219 924 L260 924 L277 612 L270 570 Z" />
        <path d="M345 575 Q359 700 351 924 L310 924 L293 612 L300 570 Z" />
        <title>{organName.muscle[language]}</title>
      </g>
      <g {...groupProps("brain")}>
        <path d="M247 24 C253 12 313 12 320 31 L318 88 Q283 101 247 89 Z" />
        <title>{organName.brain[language]}</title>
      </g>
      <g {...groupProps("lung")}>
        <path d="M273 203 C247 194 220 220 213 257 C207 294 216 335 245 346 C264 348 273 329 274 301 Z" />
        <path d="M298 203 C324 194 351 220 358 257 C364 294 355 335 326 346 C307 348 298 329 297 301 Z" />
        <title>{organName.lung[language]}</title>
      </g>
      <g {...groupProps("heart")}>
        <path d="M286 273 C273 255 249 269 252 294 C255 317 285 338 287 339 C291 336 319 317 321 293 C323 268 299 256 286 273 Z" />
        <title>{organName.heart[language]}</title>
      </g>
      <g {...groupProps("liver")}>
        <path d="M228 305 C252 289 328 291 345 311 L336 367 Q282 389 222 369 Z" />
        <title>{organName.liver[language]}</title>
      </g>
      <g {...groupProps("spleen")}>
        <path d="M352 344 C365 337 377 352 373 374 C369 392 357 402 348 391 C341 379 342 353 352 344 Z" />
        <title>{organName.spleen[language]}</title>
      </g>
      <g {...groupProps("kidney")}>
        <path d="M248 337 C235 339 232 366 242 383 C251 394 263 379 264 358 C264 345 258 337 248 337 Z M322 337 C335 339 338 366 328 383 C319 394 307 379 306 358 C306 345 312 337 322 337 Z" />
        <title>{organName.kidney[language]}</title>
      </g>
    </svg>
  );
}

function MouseOrganHeatmap({ language }: { language: Language }) {
  const t = copy[language];
  const [capsidId, setCapsidId] = useState("aav9");
  const [metric, setMetric] = useState<MetricId>("time");
  const [mappingMode, setMappingMode] = useState<ColorMappingMode>("relative");
  const [selectedOrgan, setSelectedOrgan] = useState<OrganId>("liver");
  const initialTime = heatmap.time_h.reduce((best, value, index) =>
    Math.abs(value - 8) < Math.abs(heatmap.time_h[best] - 8) ? index : best, 0);
  const [timeIndex, setTimeIndex] = useState(initialTime);
  const [playing, setPlaying] = useState(false);
  const selectedCapsid = heatmap.capsids.find((capsid) => capsid.capsid_id === capsidId) ?? heatmap.capsids[0];

  useEffect(() => {
    if (!playing || metric !== "time") return;
    const timer = window.setInterval(() => {
      setTimeIndex((index) => index >= heatmap.time_h.length - 1 ? 0 : index + 1);
    }, 420);
    return () => window.clearInterval(timer);
  }, [playing, metric]);

  const scale = useMemo(() => {
    const logarithmic = metric === "time" || metric === "auc" || metric === "peak" || metric === "expression";
    const values = metric === "time"
      ? heatmap.capsids.flatMap((capsid) =>
          heatmap.organs.flatMap((organ) => capsid.organs[organ].isf_concentration_vg_ml),
        )
      : heatmap.capsids.flatMap((capsid) =>
          heatmap.organs.map((organ) => metricValue(capsid.organs[organ], metric, timeIndex) ?? 0),
        );
    return colorScaleBounds(values, logarithmic);
  }, [metric, timeIndex]);

  const rawValues = Object.fromEntries(heatmap.organs.map((organ) => [
    organ,
    metricValue(selectedCapsid.organs[organ], metric, timeIndex),
  ])) as Record<OrganId, number | null>;
  const relativeMaximum = Math.max(
    ...Object.values(rawValues).filter((value): value is number => value !== null && value > 0),
    1e-30,
  );
  const normalized = Object.fromEntries(heatmap.organs.map((organ) => {
    const value = rawValues[organ];
    if (value === null || value <= 0) return [organ, null];
    if (mappingMode === "relative") return [organ, Math.min(1, value / relativeMaximum)];
    const transformed = scale.logarithmic ? Math.log10(value) : value;
    return [organ, Math.min(1, Math.max(0, (transformed - scale.min) / Math.max(scale.max - scale.min, 1e-12)))];
  })) as Record<OrganId, number | null>;
  const ranking = [...heatmap.organs].sort((a, b) => (rawValues[b] ?? -1) - (rawValues[a] ?? -1));
  const selectedMetrics = selectedCapsid.organs[selectedOrgan];
  const selectedValue = rawValues[selectedOrgan];
  const curve = sparklinePath(selectedMetrics.isf_concentration_vg_ml);
  const currentTimeX = 8 + (timeIndex / Math.max(heatmap.time_h.length - 1, 1)) * 284;

  return (
    <section className="heatmap-page">
      <div className="heatmap-heading">
        <div>
          <span className="eyebrow"><Activity size={15} /> {t.eyebrow}</span>
          <h1>{t.title}</h1>
          <p>{t.description}</p>
        </div>
        <div className="model-scope-badges">
          <span>{t.humanProjection}</span>
          <strong>{t.mouseModel}</strong>
        </div>
      </div>

      <div className="heatmap-controls">
        <label>
          <span>{t.capsid}</span>
          <select value={capsidId} onChange={(event) => setCapsidId(event.target.value)}>
            {heatmap.capsids.map((capsid) => <option key={capsid.capsid_id} value={capsid.capsid_id}>{capsid.capsid}</option>)}
          </select>
        </label>
        <div className="metric-control">
          <span>{t.metric}</span>
          <div className="metric-segments">
            {metricOrder.map((id) => <button key={id} type="button" className={metric === id ? "active" : ""} onClick={() => { setMetric(id); if (id !== "time") setPlaying(false); }}>{t[id]}</button>)}
          </div>
        </div>
        <div className="metric-control mapping-control">
          <span>{t.mapping}</span>
          <div className="metric-segments">
            {(["absolute", "relative"] as ColorMappingMode[]).map((mode) => <button key={mode} type="button" className={mappingMode === mode ? "active" : ""} onClick={() => setMappingMode(mode)}>{t[mode]}</button>)}
          </div>
        </div>
      </div>

      <div className="heatmap-workspace">
        <aside className="organ-ranking">
          <div className="map-panel-title"><BarChart3 size={16} /><strong>{t.rank}</strong><span>{metricUnit(metric)}</span></div>
          <div className="ranking-list">
            {ranking.map((organ, index) => {
              const value = rawValues[organ];
              return <button type="button" key={organ} className={selectedOrgan === organ ? "selected" : ""} onClick={() => setSelectedOrgan(organ)}>
                <span className="rank-number">{String(index + 1).padStart(2, "0")}</span>
                <span className="rank-organ"><strong>{organName[organ][language]}</strong><i style={{ background: normalized[organ] === null ? "#d7dcda" : colorScale(normalized[organ] ?? 0) }} /></span>
                <span className="rank-value">{value === null ? "—" : compact(value)}</span>
              </button>;
            })}
          </div>
          <div className="scale-legend">
            <div className="scale-bar" />
            <div><span>{mappingMode === "relative" ? "0%" : compact(scale.rawMin)}</span><span>{mappingMode === "relative" ? "100%" : compact(scale.rawMax)}</span></div>
            <small>{mappingMode === "relative" ? t.relativeScale : scale.logarithmic ? t.logScale : t.linearScale}</small>
          </div>
        </aside>

        <div className="anatomy-panel">
          <div className="anatomy-status">
            <span className={`evidence-pill evidence-${selectedCapsid.evidence}`}>{selectedCapsid.evidence}</span>
            <strong>{selectedCapsid.capsid}</strong>
            <small>{selectedCapsid.species}</small>
          </div>
          <AnatomicalMap values={normalized} selected={selectedOrgan} onSelect={setSelectedOrgan} language={language} />
          <a className="anatomy-credit" href="https://reactome.org/content/detail/R-ICO-013956" target="_blank" rel="noreferrer">{t.anatomyCredit}<ExternalLink size={11} /></a>
        </div>

        <aside className="organ-inspector">
          <div className="map-panel-title"><Info size={16} /><strong>{t.selectedOrgan}</strong></div>
          <h2>{organName[selectedOrgan][language]}</h2>
          <div className="mapped-value">
            <span>{t.currentValue}</span>
            <strong>{selectedValue === null ? "—" : compact(selectedValue)}</strong>
            <small>{metricUnit(metric)}</small>
          </div>
          {selectedValue === null && <p className="unavailable-note">{t.unavailable}</p>}
          <dl className="organ-metrics">
            <div><dt>{t.peakConcentration}</dt><dd>{compact(selectedMetrics.peak_isf_concentration_vg_ml)} vg-eq/mL</dd></div>
            <div><dt>{t.tmax}</dt><dd>{compact(selectedMetrics.tmax_h)} h</dd></div>
            <div><dt>{t.exposureShare}</dt><dd>{compact(selectedMetrics.exposure_share_pct)}%</dd></div>
            <div><dt>{t.postBarrier}</dt><dd>{compact(selectedMetrics.peak_post_barrier_delivery_pct, 3)}%</dd></div>
            <div><dt>{t.transduction}</dt><dd className={selectedMetrics.model_status === "ode-derived" ? "status-derived" : "status-exposure"}>{selectedMetrics.model_status === "ode-derived" ? t.odeDerived : t.exposureOnly}</dd></div>
          </dl>
          <div className="curve-card">
            <span>{t.curve}</span>
            <svg viewBox="0 0 300 84" role="img" aria-label={`${organName[selectedOrgan][language]} concentration curve`}>
              <path d="M8 74 H292" className="curve-axis" />
              <path d={curve} className="curve-line" />
              {metric === "time" && <line x1={currentTimeX} x2={currentTimeX} y1="10" y2="74" className="curve-cursor" />}
            </svg>
          </div>
          <a className="source-link" href={selectedCapsid.source} target="_blank" rel="noreferrer">{t.source}<ExternalLink size={14} /></a>
        </aside>
      </div>

      <div className="heatmap-timeline">
        <button type="button" className="timeline-play" onClick={() => { setMetric("time"); setPlaying((value) => !value); }} aria-label={playing ? t.pause : t.play} title={playing ? t.pause : t.play}>
          {playing ? <Pause size={16} /> : <Play size={16} />}
        </button>
        <Clock3 size={16} />
        <span>{t.timeline}</span>
        <input type="range" min="0" max={heatmap.time_h.length - 1} value={timeIndex} onChange={(event) => { setTimeIndex(Number(event.target.value)); setMetric("time"); setPlaying(false); }} />
        <strong>{compact(heatmap.time_h[timeIndex], 1)} h</strong>
      </div>

      <div className="heatmap-disclaimer"><Info size={15} /><p>{t.disclaimer}</p></div>
    </section>
  );
}

const humanRegionNameZh: Record<string, string> = {
  brain_frontal: "额叶皮层",
  brain_parietal: "顶叶皮层",
  brain_temporal: "颞叶皮层",
  brain_occipital: "枕叶皮层",
  brain_deep_gray: "深部灰质核团",
  brain_cerebellum: "小脑",
  brainstem_spinal: "脑干与脊髓",
  heart: "心肌",
  liver: "肝脏",
  spleen: "脾脏",
  kidney_left_cortex: "左肾皮质",
  kidney_left_medulla: "左肾髓质",
  kidney_right_cortex: "右肾皮质",
  kidney_right_medulla: "右肾髓质",
  muscle_injected_arm: "注射侧上肢肌肉",
  muscle_contralateral_arm: "对侧上肢肌肉",
  muscle_trunk: "躯干肌肉",
  muscle_legs: "下肢肌肉",
  gut: "胃肠道",
  skin_adipose: "皮肤与脂肪",
  bone_marrow: "骨与骨髓",
  rest: "其他组织",
  lung_left: "左肺",
  lung_right: "右肺",
};

function humanRegionName(regionId: string, region: HumanRegionMetrics, language: Language) {
  return language === "zh" ? (humanRegionNameZh[regionId] ?? region.label) : region.label;
}

function interpolateAt(times: number[], values: number[], target: number) {
  if (target <= times[0]) return values[0];
  if (target >= times[times.length - 1]) return values[values.length - 1];
  let low = 0;
  let high = times.length - 1;
  while (high - low > 1) {
    const middle = Math.floor((low + high) / 2);
    if (times[middle] <= target) low = middle;
    else high = middle;
  }
  const fraction = (target - times[low]) / Math.max(times[high] - times[low], 1e-30);
  return values[low] + (values[high] - values[low]) * fraction;
}

function humanMetricValue(metrics: HumanRegionMetrics, metric: HumanMetricId, timeH: number, times: number[]) {
  if (metric === "vascular") return interpolateAt(times, metrics.vascular_concentration_vg_ml, timeH);
  if (metric === "isf") return interpolateAt(times, metrics.isf_concentration_vg_ml, timeH);
  if (metric === "episome") return interpolateAt(times, metrics.episome_vg, timeH);
  if (metric === "protein") return interpolateAt(times, metrics.protein_au, timeH);
  if (metric === "auc") return metrics.auc_isf_concentration_vg_h_ml;
  if (metric === "share") return metrics.exposure_share_pct;
  if (metric === "tmax") return metrics.tmax_isf_h;
  return metrics.peak_post_barrier_delivery_pct;
}

function humanMetricUnit(metric: HumanMetricId) {
  if (metric === "vascular" || metric === "isf") return "vg-eq/mL";
  if (metric === "episome") return "vg-eq";
  if (metric === "protein") return "a.u.";
  if (metric === "auc") return "vg-eq·h/mL";
  if (metric === "tmax") return "h";
  return "%";
}

function formatModelTime(hours: number, language: Language) {
  if (hours < 1 / 60) return `${compact(hours * 3600, 1)} s`;
  if (hours < 1) return `${compact(hours * 60, 1)} min`;
  if (hours < 48) return `${compact(hours, 2)} h`;
  if (hours < 24 * 365) return `${compact(hours / 24, 1)} ${language === "zh" ? "天" : "d"}`;
  return `${compact(hours / (24 * 365), 2)} ${language === "zh" ? "年" : "y"}`;
}

function sparklinePathLinearTime(times: number[], values: number[], startH: number, endH: number) {
  const points = times.map((hours, index) => ({ hours, value: values[index] }))
    .filter(({ hours }) => hours >= startH && hours <= endH);
  const maxValue = Math.max(...points.map(({ value }) => value), 1e-30);
  return points.map(({ hours, value }, index) => {
    const x = 8 + ((hours - startH) / Math.max(endH - startH, 1e-30)) * 284;
    const y = 74 - (Math.max(value, 0) / maxValue) * 58;
    return `${index === 0 ? "M" : "L"}${x.toFixed(2)},${y.toFixed(2)}`;
  }).join(" ");
}

const humanRegionShape = {
  contralateralArm: "M154 174 C176 160 196 170 207 194 C202 253 192 316 174 379 L142 486 C139 499 143 511 137 525 C131 538 117 546 106 542 C97 539 94 529 98 520 C101 513 100 505 96 499 L120 374 C133 303 139 235 154 174 Z",
  injectedArm: "M416 174 C394 160 374 170 363 194 C368 253 378 316 396 379 L428 486 C431 499 427 511 433 525 C439 538 453 546 464 542 C473 539 476 529 472 520 C469 513 470 505 474 499 L450 374 C437 303 431 235 416 174 Z",
  legs: "M200 530 C221 523 252 527 278 545 C275 610 270 673 263 742 C258 792 260 840 267 890 C272 921 264 945 249 960 C235 972 216 970 207 958 C201 949 203 939 208 930 C202 846 194 765 188 686 C184 620 186 567 200 530 Z M370 530 C349 523 318 527 292 545 C295 610 300 673 307 742 C312 792 310 840 303 890 C298 921 306 945 321 960 C335 972 354 970 363 958 C369 949 367 939 362 930 C368 846 376 765 382 686 C386 620 384 567 370 530 Z",
  boneMarrow: "M282 155 H288 L291 338 Q285 345 279 338 Z M236 565 Q240 559 244 565 L244 900 Q240 911 236 900 Z M326 565 Q330 559 334 565 L334 900 Q330 911 326 900 Z",
} as const;

function HumanRegionalMap({
  values,
  circulation,
  selected,
  onSelect,
  language,
  playing,
  administrationId,
}: {
  values: Record<string, number>;
  circulation: Record<string, number>;
  selected: string;
  onSelect: (regionId: string) => void;
  language: Language;
  playing: boolean;
  administrationId: HumanAdministrationRoute["route_id"];
}) {
  const tintFor = (regionId: string, opacity = 0.86) => ({
    fill: colorScale(values[regionId] ?? 0),
    opacity,
    filter: selected === regionId ? "url(#human-selected-glow)" : undefined,
  });
  const propsFor = (regionId: string) => ({
    role: "button" as const,
    tabIndex: 0,
    className: `anatomy-organ anatomy-hit-target ${selected === regionId ? "selected" : ""}`,
    "aria-label": humanRegionNameZh[regionId] ?? regionId,
    onClick: () => onSelect(regionId),
    onKeyDown: (event: React.KeyboardEvent<SVGGElement>) => {
      if (event.key === "Enter" || event.key === " ") onSelect(regionId);
    },
  });
  const vessel = (id: string) => ({ stroke: colorScale(circulation[id] ?? 0) });

  return (
    <svg className="anatomy-map human-regional-map" viewBox="0 0 600 1000" aria-label="Human multiregion AAV PBPK map">
      <ReactomeAnatomyMasks prefix="human" />
      <AnatomyReferenceLayer opacity={0.2} />

      <g className="anatomy-tint-layer" aria-hidden="true" pointerEvents="none">
        <rect x="0" y="0" width="600" height="1000" mask="url(#human-body-mask)" style={tintFor("rest", 0.2)} />
        <rect x="0" y="0" width="600" height="1000" mask="url(#human-body-mask)" style={tintFor("skin_adipose", 0.2)} />

        <g mask="url(#human-body-mask)">
          <path d={humanRegionShape.contralateralArm} style={tintFor("muscle_contralateral_arm", 0.76)} />
          <path d={humanRegionShape.injectedArm} style={tintFor("muscle_injected_arm", 0.76)} />
          <path d="M200 145 Q285 124 370 145 L365 570 Q285 600 205 570 Z" style={tintFor("muscle_trunk", 0.52)} />
          <path d={humanRegionShape.legs} style={tintFor("muscle_legs", 0.76)} />
        </g>

        <g mask="url(#human-brain-mask)">
          <path d="M238 18 H283 V59 H238 Z" style={tintFor("brain_frontal")} />
          <path d="M283 18 H330 V57 L283 59 Z" style={tintFor("brain_parietal")} />
          <path d="M238 59 H283 V98 H238 Z" style={tintFor("brain_temporal")} />
          <path d="M283 57 H330 V98 H283 Z" style={tintFor("brain_occipital")} />
          <ellipse cx="285" cy="59" rx="10" ry="9" style={tintFor("brain_deep_gray")} />
          <path d="M278 77 Q296 69 318 78 L318 98 H278 Z" style={tintFor("brain_cerebellum")} />
        </g>
        <path d="M282 81 Q288 76 294 82 L298 190 L291 250 L284 190 Z" style={tintFor("brainstem_spinal")} />

        <rect x="285" y="0" width="315" height="1000" mask="url(#human-lung-mask)" style={tintFor("lung_left")} />
        <rect x="0" y="0" width="285" height="1000" mask="url(#human-lung-mask)" style={tintFor("lung_right")} />
        <path d="M286 273 C273 255 249 269 252 294 C255 317 285 338 287 339 C291 336 319 317 321 293 C323 268 299 256 286 273 Z" style={tintFor("heart")} />
        <rect x="0" y="0" width="600" height="1000" mask="url(#human-liver-mask)" style={tintFor("liver")} />
        <path d="M352 344 C365 337 377 352 373 374 C369 392 357 402 348 391 C341 379 342 353 352 344 Z" style={tintFor("spleen")} />

        <rect x="285" y="0" width="315" height="1000" mask="url(#human-kidney-mask)" style={tintFor("kidney_left_cortex")} />
        <ellipse cx="311" cy="361" rx="6" ry="12" mask="url(#human-kidney-mask)" style={tintFor("kidney_left_medulla")} />
        <rect x="0" y="0" width="285" height="1000" mask="url(#human-kidney-mask)" style={tintFor("kidney_right_cortex")} />
        <ellipse cx="259" cy="361" rx="6" ry="12" mask="url(#human-kidney-mask)" style={tintFor("kidney_right_medulla")} />
        <rect x="0" y="0" width="600" height="1000" mask="url(#human-gut-mask)" style={tintFor("gut")} />
        <path d={humanRegionShape.boneMarrow} style={tintFor("bone_marrow", 0.5)} />
      </g>

      <g className="anatomy-interaction-layer">
        <g {...propsFor("rest")}><circle cx="285" cy="70" r="58" /><path d="M238 150 Q285 126 332 150 L382 205 L365 520 L337 590 L335 930 L292 930 L282 610 L270 930 L227 930 L228 590 L202 520 L188 205 Z" /><path d="M194 206 L132 520 L91 510 L153 188 Z M376 206 L438 520 L479 510 L417 188 Z" /></g>
        <g {...propsFor("skin_adipose")}><path d="M209 178 Q285 135 361 178 L389 330 L365 540 Q339 580 317 590 L306 548 Q335 515 340 465 L349 210 Q285 178 221 210 L230 465 Q235 515 264 548 L253 590 Q231 580 205 540 L181 330 Z" fillRule="evenodd" /></g>
        <g {...propsFor("muscle_contralateral_arm")}><path d={humanRegionShape.contralateralArm} /></g>
        <g {...propsFor("muscle_injected_arm")}><path d={humanRegionShape.injectedArm} /></g>
        <g {...propsFor("muscle_trunk")}><path d="M214 190 Q285 158 356 190 L349 360 Q329 350 311 362 L285 405 L259 362 Q241 350 221 360 Z" /></g>
        <g {...propsFor("muscle_legs")}><path d={humanRegionShape.legs} /></g>
        <g {...propsFor("bone_marrow")}><path d={humanRegionShape.boneMarrow} /></g>

        <g {...propsFor("brain_frontal")}><path d="M247 24 H283 V59 H247 Z" /></g>
        <g {...propsFor("brain_parietal")}><path d="M283 24 H319 V57 L283 59 Z" /></g>
        <g {...propsFor("brain_temporal")}><path d="M247 59 H283 V91 H247 Z" /></g>
        <g {...propsFor("brain_occipital")}><path d="M283 57 H319 V91 H283 Z" /></g>
        <g {...propsFor("brain_deep_gray")}><ellipse cx="285" cy="59" rx="10" ry="9" /></g>
        <g {...propsFor("brain_cerebellum")}><path d="M278 77 Q296 69 318 78 L316 93 H278 Z" /></g>
        <g {...propsFor("brainstem_spinal")}><path d="M282 81 Q288 76 294 82 L298 190 L291 250 L284 190 Z" /></g>

        <g {...propsFor("lung_left")}><path d="M297 203 C323 194 351 220 358 257 C364 294 355 335 326 346 C307 348 298 329 297 301 Z" /></g>
        <g {...propsFor("lung_right")}><path d="M274 203 C248 194 220 220 213 257 C207 294 216 335 245 346 C264 348 273 329 274 301 Z" /></g>
        <g {...propsFor("heart")}><path d="M286 273 C273 255 249 269 252 294 C255 317 285 338 287 339 C291 336 319 317 321 293 C323 268 299 256 286 273 Z" /></g>
        <g {...propsFor("liver")}><path d="M228 305 C252 289 328 291 345 311 L336 367 Q282 389 222 369 Z" /></g>
        <g {...propsFor("spleen")}><path d="M352 344 C365 337 377 352 373 374 C369 392 357 402 348 391 C341 379 342 353 352 344 Z" /></g>
        <g {...propsFor("gut")}><path d="M227 354 Q285 338 343 356 L338 468 Q286 486 232 466 Z" /></g>
        <g {...propsFor("kidney_left_cortex")}><path d="M322 337 C335 339 338 366 328 383 C319 394 307 379 306 358 C306 345 312 337 322 337 Z" /></g>
        <g {...propsFor("kidney_left_medulla")}><ellipse cx="311" cy="361" rx="6" ry="12" /></g>
        <g {...propsFor("kidney_right_cortex")}><path d="M248 337 C235 339 232 366 242 383 C251 394 263 379 264 358 C264 345 258 337 248 337 Z" /></g>
        <g {...propsFor("kidney_right_medulla")}><ellipse cx="259" cy="361" rx="6" ry="12" /></g>
      </g>

      {administrationId === "iv" && <>
        <g className={`vascular-route ${playing ? "playing" : ""}`} pointerEvents="none">
          <path d="M472 485 C448 441 426 349 398 274 C365 244 330 249 304 276" style={vessel("arm_vein")} />
          <path d="M304 276 C294 274 289 282 286 297" style={vessel("right_heart")} />
          <path d="M286 297 C270 270 245 250 231 245 M286 297 C307 267 333 249 349 245" style={vessel("pulmonary_artery")} />
          <path d="M231 258 C253 276 270 287 287 305 M349 258 C326 278 307 290 287 305" style={vessel("pulmonary_vein")} />
          <path d="M287 305 C300 300 307 310 303 324" style={vessel("left_heart")} />
          <path d="M303 324 C319 354 317 432 304 548 M304 548 C294 570 263 588 249 622 C241 701 242 802 240 895 M304 548 C316 570 337 588 337 622 C345 701 337 802 335 895" style={vessel("arterial")} />
          <path d="M236 895 C234 802 235 702 246 620 C250 580 266 560 276 548 M334 895 C336 802 335 702 324 620 C320 580 286 560 276 548 M276 548 C258 432 260 355 276 322" style={vessel("venous")} />
          <circle className="injection-site" cx="472" cy="485" r="5" style={{ fill: colorScale(circulation.arm_vein ?? 0) }} />
        </g>
        <text x="482" y="480" className="injection-label">{language === "zh" ? "左臂 IV" : "Left-arm IV"}</text>
      </>}

      {administrationId === "intrathecal" && <>
        <g className={`vascular-route csf-route ${playing ? "playing" : ""}`} pointerEvents="none">
          <path d="M286 535 C282 470 288 405 286 338 C284 270 291 205 288 130" style={vessel("csf_lumbar")} />
          <path d="M288 130 C286 105 286 91 286 81 M286 81 C263 94 246 82 244 58 M286 81 C310 94 326 82 327 58" style={vessel("csf_cranial")} />
          <path d="M310 92 C337 126 326 205 304 276 M286 535 C301 482 299 414 276 322" style={vessel("venous")} />
          <circle className="injection-site" cx="286" cy="535" r="5" style={{ fill: colorScale(circulation.csf_lumbar ?? 0) }} />
        </g>
        <text x="300" y="548" className="injection-label">{language === "zh" ? "腰椎 IT" : "Lumbar IT"}</text>
      </>}

      {administrationId === "intramuscular" && <>
        <g className={`vascular-route im-route ${playing ? "playing" : ""}`} pointerEvents="none">
          <path d="M407 188 C416 211 420 240 416 271" style={vessel("im_depot")} />
          <path d="M407 188 C386 220 346 251 304 276" style={vessel("arm_vein")} />
          <path d="M304 276 C294 274 289 282 286 297" style={vessel("right_heart")} />
          <path d="M286 297 C270 270 245 250 231 245 M286 297 C307 267 333 249 349 245" style={vessel("pulmonary_artery")} />
          <path d="M231 258 C253 276 270 287 287 305 M349 258 C326 278 307 290 287 305" style={vessel("pulmonary_vein")} />
          <path d="M287 305 C300 300 307 310 303 324 C319 354 317 432 304 548" style={vessel("arterial")} />
          <circle className="injection-site" cx="407" cy="188" r="6" style={{ fill: colorScale(circulation.im_depot ?? 0) }} />
        </g>
        <text x="422" y="185" className="injection-label">{language === "zh" ? "三角肌 IM" : "Deltoid IM"}</text>
      </>}

      {administrationId === "intracisternal" && <>
        <g className={`vascular-route csf-route ${playing ? "playing" : ""}`} pointerEvents="none">
          <path d="M314 91 C308 102 298 116 288 130 C286 105 286 91 286 81 M286 81 C263 94 246 82 244 58 M286 81 C310 94 326 82 327 58" style={vessel("csf_cranial")} />
          <path d="M288 130 C291 205 284 270 286 338 C288 405 282 470 286 535" style={vessel("csf_lumbar")} />
          <path d="M310 92 C337 126 326 205 304 276" style={vessel("venous")} />
          <circle className="injection-site" cx="314" cy="91" r="5" style={{ fill: colorScale(circulation.csf_cranial ?? 0) }} />
        </g>
        <text x="328" y="88" className="injection-label">{language === "zh" ? "枕大池 ICM" : "Cisterna magna"}</text>
      </>}

      {administrationId === "intracerebroventricular" && <>
        <g className={`vascular-route csf-route ${playing ? "playing" : ""}`} pointerEvents="none">
          <path d="M285 59 C286 70 286 76 286 81 M286 81 C263 94 246 82 244 58 M286 81 C310 94 326 82 327 58" style={vessel("csf_cranial")} />
          <path d="M288 130 C291 205 284 270 286 338 C288 405 282 470 286 535" style={vessel("csf_lumbar")} />
          <path d="M310 92 C337 126 326 205 304 276" style={vessel("venous")} />
          <circle className="injection-site" cx="285" cy="59" r="5" style={{ fill: colorScale(circulation.csf_cranial ?? 0) }} />
        </g>
        <text x="300" y="45" className="injection-label">{language === "zh" ? "脑室 ICV" : "Ventricular ICV"}</text>
      </>}

      {administrationId === "inhaled" && <>
        <g className={`vascular-route airway-route ${playing ? "playing" : ""}`} pointerEvents="none">
          <path d="M286 137 L286 202 M286 202 C267 218 246 232 231 245 M286 202 C306 218 333 232 349 245" style={vessel("airway_depot")} />
          <path d="M231 258 C253 276 270 287 287 305 M349 258 C326 278 307 290 287 305" style={vessel("pulmonary_vein")} />
          <circle className="injection-site" cx="286" cy="137" r="5" style={{ fill: colorScale(circulation.airway_depot ?? 0) }} />
        </g>
        <text x="302" y="141" className="injection-label">{language === "zh" ? "吸入" : "Inhaled"}</text>
      </>}
    </svg>
  );
}

function HumanSpatialHeatmap({ data, language }: { data: HumanSpatialPayload; language: Language }) {
  const t = copy[language];
  const [timeWindow, setTimeWindow] = useState<"early" | "long">("early");
  const metrics: HumanMetricId[] = timeWindow === "early"
    ? ["vascular", "isf", "tmax", "auc", "share", "peak"]
    : ["episome", "protein"];
  const [administrationId, setAdministrationId] = useState<HumanAdministrationRoute["route_id"]>(data.default_route_id);
  const [capsidId, setCapsidId] = useState("aav9");
  const [metric, setMetric] = useState<HumanMetricId>("isf");
  const [mappingMode, setMappingMode] = useState<ColorMappingMode>("relative");
  const [selectedRegion, setSelectedRegion] = useState("muscle_injected_arm");
  const [timeH, setTimeH] = useState(6);
  const [playing, setPlaying] = useState(false);
  const selectedRoute = data.administration_routes.find((route) => route.route_id === administrationId)
    ?? data.administration_routes[0];
  const selectedCapsid = selectedRoute.capsids.find((capsid) => capsid.capsid_id === capsidId)
    ?? selectedRoute.capsids[0];
  const selectedMetrics = selectedCapsid.regions[selectedRegion];
  const windowStartH = timeWindow === "early" ? 0 : 24;
  const windowEndH = timeWindow === "early" ? 72 : 730 * 24;

  useEffect(() => {
    if (!playing) return;
    const timer = window.setInterval(() => {
      const stepH = timeWindow === "early" ? 0.4 : 4 * 24;
      setTimeH((hours) => hours >= windowEndH ? windowStartH : Math.min(hours + stepH, windowEndH));
    }, 180);
    return () => window.clearInterval(timer);
  }, [playing, timeWindow, windowEndH, windowStartH]);

  const scale = useMemo(() => {
    const logarithmic = metric !== "share" && metric !== "tmax";
    const seriesKey = metric === "vascular" ? "vascular_concentration_vg_ml"
      : metric === "isf" ? "isf_concentration_vg_ml"
        : metric === "episome" ? "episome_vg"
          : metric === "protein" ? "protein_au" : null;
    const values = seriesKey
      ? selectedRoute.capsids.flatMap((capsid) =>
          data.region_ids.flatMap((regionId) =>
            capsid.regions[regionId][seriesKey].filter(
              (_, index) => data.time_h[index] >= windowStartH && data.time_h[index] <= windowEndH,
            ),
          ),
        )
      : selectedRoute.capsids.flatMap((capsid) =>
          data.region_ids.map((regionId) =>
            humanMetricValue(capsid.regions[regionId], metric, timeH, data.time_h),
          ),
        );
    return colorScaleBounds(values, logarithmic);
  }, [
    data.region_ids,
    data.time_h,
    metric,
    selectedRoute.capsids,
    timeH,
    windowEndH,
    windowStartH,
  ]);

  const rawValues = Object.fromEntries(data.region_ids.map((regionId) => [
    regionId,
    humanMetricValue(selectedCapsid.regions[regionId], metric, timeH, data.time_h),
  ])) as Record<string, number>;
  const relativeMaximum = Math.max(
    ...Object.values(rawValues).filter((value) => Number.isFinite(value) && value > 0),
    1e-30,
  );
  const normalized = Object.fromEntries(data.region_ids.map((regionId) => {
    const value = rawValues[regionId];
    if (!(value > 0)) return [regionId, 0];
    if (mappingMode === "relative") return [regionId, Math.min(1, value / relativeMaximum)];
    const transformed = scale.logarithmic ? Math.log10(value) : value;
    return [regionId, Math.min(1, Math.max(0, (transformed - scale.min) / Math.max(scale.max - scale.min, 1e-12)))];
  })) as Record<string, number>;
  const routeSeries = {
    ...Object.fromEntries(data.circulation_ids.map((id) => [id, selectedCapsid.circulation[id].concentration_vg_ml])),
    ...Object.fromEntries(Object.entries(selectedCapsid.route_compartments).map(([id, compartment]) => [
      id,
      compartment.concentration_vg_ml ?? compartment.amount_vg ?? [],
    ])),
  } as Record<string, number[]>;
  const routeNormalized = Object.fromEntries(Object.entries(routeSeries).map(([id, series]) => {
    const current = series.length ? interpolateAt(data.time_h, series, timeH) : 0;
    const peak = Math.max(...series, 1e-30);
    return [id, Math.sqrt(Math.max(current, 0) / peak)];
  })) as Record<string, number>;
  const ranking = [...data.region_ids].sort((a, b) => rawValues[b] - rawValues[a]);
  const fullSeries = metric === "vascular" ? selectedMetrics.vascular_concentration_vg_ml
    : metric === "episome" ? selectedMetrics.episome_vg
      : metric === "protein" ? selectedMetrics.protein_au
        : selectedMetrics.isf_concentration_vg_ml;
  const curve = sparklinePathLinearTime(data.time_h, fullSeries, windowStartH, windowEndH);
  const currentTimeX = 8 + ((timeH - windowStartH) / Math.max(windowEndH - windowStartH, 1e-30)) * 284;
  const presets = timeWindow === "early"
    ? [5 / 60, 1, 6, 24, 72]
    : [24, 7 * 24, 30 * 24, 180 * 24, 365 * 24, 730 * 24];

  const selectWindow = (next: "early" | "long") => {
    setTimeWindow(next);
    setMetric(next === "early" ? "isf" : "protein");
    setTimeH(next === "early" ? 6 : 30 * 24);
    setPlaying(false);
  };

  const selectAdministration = (next: HumanAdministrationRoute["route_id"]) => {
    setAdministrationId(next);
    setTimeWindow("early");
    setMetric("isf");
    const csfRoute = next === "intrathecal" || next === "intracisternal" || next === "intracerebroventricular";
    setTimeH(next === "intramuscular" ? 6 : 1);
    setSelectedRegion(
      next === "intrathecal" ? "brainstem_spinal"
        : next === "intracisternal" ? "brain_cerebellum"
          : next === "intracerebroventricular" ? "brain_deep_gray"
            : next === "intramuscular" ? "muscle_injected_arm"
              : next === "inhaled" ? "lung_right" : "liver",
    );
    if (csfRoute) setCapsidId("aav9");
    setPlaying(false);
  };

  return (
    <section className="heatmap-page human-heatmap-page">
      <div className="heatmap-heading">
        <div>
          <span className="eyebrow"><Activity size={15} /> {t.eyebrow}</span>
          <h1>{language === "zh" ? "人体多区域 AAV 时空递送" : "Human multiregion AAV delivery"}</h1>
          <p>{language === "zh" ? selectedRoute.description_zh : selectedRoute.description}</p>
        </div>
        <div className="model-scope-badges"><span>{t.humanProjection}</span><strong>{t.physiology}</strong></div>
      </div>

      <div className="heatmap-controls">
        <div className="window-control"><span>{t.timeline}</span><div className="time-window-toggle"><button type="button" className={timeWindow === "early" ? "active" : ""} onClick={() => selectWindow("early")}>{t.earlyWindow}</button><button type="button" className={timeWindow === "long" ? "active" : ""} onClick={() => selectWindow("long")}>{t.longWindow}</button></div></div>
        <label><span>{t.administration}</span><select value={administrationId} onChange={(event) => selectAdministration(event.target.value as HumanAdministrationRoute["route_id"])}>{data.administration_routes.map((route) => <option key={route.route_id} value={route.route_id}>{language === "zh" ? route.label_zh : route.label}</option>)}</select></label>
        <label><span>{t.capsid}</span><select value={capsidId} onChange={(event) => setCapsidId(event.target.value)}>{selectedRoute.capsids.map((capsid) => <option key={capsid.capsid_id} value={capsid.capsid_id}>{capsid.capsid}</option>)}</select></label>
        <div className="metric-control"><span>{t.metric}</span><div className="metric-segments human-metric-segments">{metrics.map((id) => <button key={id} type="button" className={metric === id ? "active" : ""} onClick={() => setMetric(id)}>{t[id]}</button>)}</div></div>
        <div className="metric-control mapping-control"><span>{t.mapping}</span><div className="metric-segments">{(["absolute", "relative"] as ColorMappingMode[]).map((mode) => <button key={mode} type="button" className={mappingMode === mode ? "active" : ""} onClick={() => setMappingMode(mode)}>{t[mode]}</button>)}</div></div>
      </div>

      <div className="heatmap-workspace human-workspace">
        <aside className="organ-ranking">
          <div className="map-panel-title"><BarChart3 size={16} /><strong>{t.rank}</strong><span>{humanMetricUnit(metric)}</span></div>
          <div className="ranking-list regional-ranking-list">{ranking.map((regionId, index) => <button type="button" key={regionId} className={selectedRegion === regionId ? "selected" : ""} onClick={() => setSelectedRegion(regionId)}><span className="rank-number">{String(index + 1).padStart(2, "0")}</span><span className="rank-organ"><strong>{humanRegionName(regionId, selectedCapsid.regions[regionId], language)}</strong><i style={{ background: colorScale(normalized[regionId]) }} /></span><span className="rank-value">{compact(rawValues[regionId])}</span></button>)}</div>
          <div className="scale-legend"><div className="scale-bar" /><div><span>{mappingMode === "relative" ? "0%" : metric === "tmax" ? t.earlyPeak : compact(scale.rawMin)}</span><span>{mappingMode === "relative" ? "100%" : metric === "tmax" ? t.latePeak : compact(scale.rawMax)}</span></div><small>{mappingMode === "relative" ? t.relativeScale : scale.logarithmic ? t.logScale : metric === "tmax" ? t.tmaxScale : t.linearScale}</small></div>
        </aside>

        <div className="anatomy-panel">
          <div className="anatomy-status"><span className={`evidence-pill evidence-${selectedCapsid.evidence}`}>{selectedCapsid.evidence}</span><strong>{selectedCapsid.capsid}</strong><small>{formatModelTime(timeH, language)}</small></div>
          <HumanRegionalMap values={normalized} circulation={routeNormalized} selected={selectedRegion} onSelect={setSelectedRegion} language={language} playing={playing} administrationId={administrationId} />
          <div className="route-caption"><Activity size={13} />{language === "zh" ? selectedRoute.description_zh : selectedRoute.description}</div>
          <a className="anatomy-credit" href="https://reactome.org/content/detail/R-ICO-013956" target="_blank" rel="noreferrer">{t.anatomyCredit}<ExternalLink size={11} /></a>
        </div>

        <aside className="organ-inspector">
          <div className="map-panel-title"><Info size={16} /><strong>{t.selectedOrgan}</strong></div>
          <h2>{humanRegionName(selectedRegion, selectedMetrics, language)}</h2>
          <div className="mapped-value"><span>{t.currentValue}</span><strong>{compact(rawValues[selectedRegion])}</strong><small>{humanMetricUnit(metric)}</small></div>
          <dl className="organ-metrics">
            <div><dt>{t.peakConcentration}</dt><dd>{compact(selectedMetrics.peak_isf_concentration_vg_ml)} vg-eq/mL</dd></div>
            <div><dt>{t.tmax}</dt><dd>{formatModelTime(selectedMetrics.tmax_isf_h, language)}</dd></div>
            <div><dt>{t.flow}</dt><dd>{compact(selectedMetrics.blood_flow_ml_h / 1000)} L/h</dd></div>
            <div><dt>{t.exposureShare}</dt><dd>{compact(selectedMetrics.exposure_share_pct)}%</dd></div>
            <div><dt>{t.postBarrier}</dt><dd>{compact(selectedMetrics.peak_post_barrier_delivery_pct, 3)}%</dd></div>
          </dl>
          <div className="curve-card"><span>{timeWindow === "early" ? t.earlyWindow : t.longWindow}</span><svg viewBox="0 0 300 84" role="img" aria-label={`${selectedMetrics.label} trajectory`}><path d="M8 74 H292" className="curve-axis" /><path d={curve} className="curve-line" /><line x1={currentTimeX} x2={currentTimeX} y1="10" y2="74" className="curve-cursor" /></svg></div>
          <div className="translation-note"><strong>{t.translationNote}</strong><p>{selectedCapsid.human_translation_note} {t.routeDoseNote}</p></div>
          <a className="source-link" href={selectedCapsid.source} target="_blank" rel="noreferrer">{t.source}<ExternalLink size={14} /></a>
        </aside>
      </div>

      <div className="human-time-controls">
        <div className="heatmap-timeline">
          <button type="button" className="timeline-play" onClick={() => setPlaying((value) => !value)} aria-label={playing ? t.pause : t.play}>{playing ? <Pause size={16} /> : <Play size={16} />}</button>
          <Clock3 size={16} /><span>{t.timeline}</span>
          <input type="range" min={timeWindow === "early" ? 0 : 1} max={timeWindow === "early" ? 72 : 730} step={timeWindow === "early" ? 0.1 : 1} value={timeWindow === "early" ? timeH : timeH / 24} onChange={(event) => { const value = Number(event.target.value); setTimeH(timeWindow === "early" ? value : value * 24); setPlaying(false); }} />
          <strong>{formatModelTime(timeH, language)}</strong>
        </div>
        <div className="time-toolbar"><span>{timeWindow === "early" ? t.earlyWindowHint : t.longWindowHint}</span><div className="time-presets">{presets.map((hours) => <button type="button" key={hours} onClick={() => { setTimeH(hours); setPlaying(false); }}>{formatModelTime(hours, language)}</button>)}</div></div>
      </div>
      <div className="heatmap-disclaimer"><Info size={15} /><p>{t.persistenceClarification} {t.humanDisclaimer}</p></div>
    </section>
  );
}

export function OrganHeatmap({ language }: { language: Language }) {
  const t = copy[language];
  const [mode, setMode] = useState<"human" | "mouse">("human");
  const [humanData, setHumanData] = useState<HumanSpatialPayload | null>(null);
  const [loadError, setLoadError] = useState(false);

  useEffect(() => {
    let active = true;
    fetch("/data/human-spatial-results.json")
      .then((response) => {
        if (!response.ok) throw new Error(`HTTP ${response.status}`);
        return response.json() as Promise<HumanSpatialPayload>;
      })
      .then((payload) => { if (active) setHumanData(payload); })
      .catch(() => { if (active) setLoadError(true); });
    return () => { active = false; };
  }, []);

  return (
    <div className="model-mode-shell">
      <div className="model-mode-switch" aria-label="PBPK model species">
        <button type="button" className={mode === "human" ? "active" : ""} onClick={() => setMode("human")}>{t.humanMode}</button>
        <button type="button" className={mode === "mouse" ? "active" : ""} onClick={() => setMode("mouse")}>{t.mouseMode}</button>
      </div>
      {mode === "mouse" ? <MouseOrganHeatmap language={language} /> : humanData ? <HumanSpatialHeatmap data={humanData} language={language} /> : <div className="human-loading"><Activity size={18} />{loadError ? (language === "zh" ? "人体模型数据载入失败" : "Failed to load human model data") : t.humanLoading}</div>}
    </div>
  );
}
