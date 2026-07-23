"use client";

import {
  Activity,
  BookOpen,
  Check,
  ChevronDown,
  ChevronRight,
  CircleAlert,
  Database,
  Dna,
  ExternalLink,
  Eye,
  EyeOff,
  FlaskConical,
  Info,
  Languages,
  Search,
  Sigma,
  Target,
} from "lucide-react";
import { useMemo, useRef, useState, type WheelEvent } from "react";
import modelPayload from "../public/data/model-results.json";
import expressionPayload from "../public/data/gene-expression.json";
import { diseases, type DiseaseRecord, type Language, type Organ } from "./disease-data";
import { OrganHeatmap } from "./organ-heatmap";

type Evidence = "strong" | "medium" | "exploratory";
type ModelPoint = {
  capsid_id: string;
  capsid: string;
  target: Organ;
  model_organ: string;
  route: string;
  specificity_log10: number;
  target_concentration_auc_vg_h_ml: number;
  target_exposure_share_pct: number;
  peak_post_barrier_delivery_pct: number;
  tmax_h: number;
  tropism_multiplier: number;
  episome_half_life_days_prior: number;
  persistence_factor: number;
  evidence: Evidence;
  species: string;
  source: string;
  model_status: "ode-derived" | "surrogate";
  max_mass_balance_error: number;
  relative_delivery: number;
  predicted_protein_restoration_pct: number;
  effective_duration_days: number;
  therapeutic_threshold_pct: number;
  transduction_model: string;
  therapeutic_onset_days: number;
  therapeutic_window_days: number;
  peak_restoration_day: number;
  persistence_censored_at_730d: boolean;
  pd_sineup_rna_half_life_days: number;
  pd_target_protein_half_life_days: number;
  cns_profile?: string;
  cns_depth_mm?: number;
  cns_cell_access_factor?: number;
  cns_target_layer_auc_fraction_pct?: number;
};

type ExpressionEvidence = {
  gencode_id: string;
  organ_median_tpm: Record<string, number | null>;
  top_modeled_organ: string | null;
  top_modeled_organ_tpm: number | null;
  top_gtex_tissue: string | null;
  top_gtex_tissue_tpm: number | null;
  tissue_tau: number | null;
  hpa: null | {
    tissue_specificity: string | null;
    tissue_distribution: string | null;
    protein_tissue_specificity: string | null;
    entry_url: string;
  };
};

type RouteOrganMetric = {
  auc_isf_concentration_vg_h_ml: number;
  exposure_share_pct: number;
  peak_post_barrier_delivery_pct: number;
  median_tmax_isf_h: number;
  peak_protein_au: number;
};

type HumanRouteSummary = {
  route_id: string;
  route_label: string;
  route_label_zh: string;
  route_class: string;
  evidence_source: string;
  capsid_id: string;
  capsid: string;
  evidence: Evidence;
  organs: Record<string, RouteOrganMetric>;
};

type RegimenAgent = HumanRouteSummary & { coveredGenes: string[] };
type RegimenRecommendation = {
  agents: RegimenAgent[];
  score: number;
  modeledGenes: string[];
  unmodeledGenes: string[];
};

const modelResults = modelPayload.results as ModelPoint[];
const cnsProfileResults = modelPayload.cns_profile_results as ModelPoint[];
const modelOutputCount = modelPayload.model_counts.total;
const expressionByGene = (expressionPayload as { genes: Record<string, ExpressionEvidence> }).genes;
const humanRouteSummary = ((modelPayload as unknown as { human_route_summary?: HumanRouteSummary[] }).human_route_summary ?? []);
const capsidOrder = ["aav2", "aav5", "aav8", "aav9", "aavrh10", "php-eb", "cap-b10", "lk03"];
const capsidLabels = Object.fromEntries(
  modelResults.filter((point) => point.target === "CNS").map((point) => [point.capsid_id, point.capsid]),
);

const copy = {
  zh: {
    subtitle: "疾病驱动的 AAV 空间递送设计",
    status: `PBPK–空间 CNS–SINEUP ODE · ${modelOutputCount} 个模型输出`,
    designNav: "疾病设计空间",
    heatmapNav: "器官热图",
    disease: "疾病",
    gene: "缺失基因",
    location: "递送位置",
    library: "疾病库",
    search: "搜索疾病、基因或位点",
    noResult: "没有匹配结果",
    geneCount: "个靶基因",
    eligibility: "疾病资格首先检查单倍剂量不足证据、正常等位基因转录本和可干预治疗窗口。",
    designSpace: "ODE 设计空间",
    candidate: "递送候选",
    chartHelp: "每种衣壳均重新求解 PBPK、细胞转导和 SINEUP-PD；点击点查看状态、参数和来源。",
    strong: "较强",
    medium: "中等",
    exploratory: "探索性",
    capsids: "衣壳显示",
    showAll: "全部显示",
    yAxis: "单次给药有效持续时间（天）",
    xAxis: "器官特异性 log10（目标 ISF 浓度 AUC / 加权脱靶 AUC）",
    preferred: "优选：精准且持久",
    preciseShort: "精准但短效",
    durableOff: "持久但脱靶负担高",
    lowUtility: "低特异性 / 短效",
    modelDerived: "ODE 计算",
    surrogate: "代理模型",
    specificity: "器官特异性",
    duration: "有效持续时间",
    delivery: "峰值后屏障递送率",
    restoration: "预测蛋白恢复",
    route: "给药路径",
    species: "证据物种",
    targetCell: "靶细胞",
    exposure: "目标器官 AUC 份额",
    tmax: "ISF Tmax",
    halfLife: "episome 半衰期先验",
    onset: "达到 65% 阈值",
    peakDay: "蛋白峰值时间",
    transduction: "细胞转导模块",
    targetDepth: "CNS 目标深度",
    layerExposure: "目标层暴露份额",
    nativeChain: "原 ode1.0 细胞内状态",
    multilevelChain: "BBB 后三级深度 + 细胞内 ODE",
    reducedChain: "ISF 驱动的降阶代理",
    source: "查看主要证据",
    hi: "单倍剂量证据",
    sineup: "SINEUP 前提",
    modelProof: "数学模型证据",
    earlyPk: "早期递送 ODE",
    earlyPkText: "先求解血液–血管–ISF PBPK；CNS 继续求解浅层、皮层和深部三级空间转运及细胞内转导，再数值求解 Epi–SINEUP RNA–靶蛋白 ODE。",
    efficiencyFormula: "递送效率",
    specificityFormula: "特异性",
    persistenceFormula: "持久性",
    balance: "最大质量守恒误差",
    generated: "模型生成时间",
    disclaimer: "研究级设计工具，不是临床剂量建议。Eye 当前使用局部给药屏障代理；骨髓与周围神经疾病用于显示下一步需要增加的 PBPK 器官室。",
    expressionTitle: "正常组织表达证据",
    targetTpm: "靶器官 GTEx 中位 TPM",
    topTissue: "GTEx 最高组织",
    tissueTau: "组织特异性 τ",
    hpaClass: "HPA 组织分类",
    expressionNote: "成人 bulk RNA 只表示正常转录本存在先验，不代表靶细胞表达、发育期表达或治疗有效性。",
    combinationTitle: "ODE 方案组合筛选",
    singlePlan: "优先单方案",
    dualPlan: "探索性双方案",
    coverageScore: "加权覆盖分数",
    coveredGenes: "覆盖基因",
    unmodeledGenes: "尚无对应器官 ODE",
    combinationNote: "双方案只有在校正侵入性、第二次给药负担和衣壳免疫风险后仍明显改善覆盖时才会出现。当前每个方案按完整参考剂量独立求解，尚未做总剂量拆分或免疫相互作用 ODE，因此不是临床联合用药建议。",
    dataSource: "数据库来源",
  },
  en: {
    subtitle: "Disease-guided AAV spatial delivery design",
    status: `PBPK–spatial CNS–SINEUP ODE · ${modelOutputCount} model outputs`,
    designNav: "Disease design space",
    heatmapNav: "Organ heatmap",
    disease: "Disease",
    gene: "Deleted gene",
    location: "Delivery site",
    library: "Disease library",
    search: "Search disease, gene or locus",
    noResult: "No matching records",
    geneCount: "target genes",
    eligibility: "Eligibility starts with dosage evidence, a remaining normal transcript, and an actionable treatment window.",
    designSpace: "ODE design space",
    candidate: "delivery candidates",
    chartHelp: "Every capsid re-solves PBPK, cellular transduction and SINEUP-PD. Select a point to inspect states, parameters and provenance.",
    strong: "Strong",
    medium: "Moderate",
    exploratory: "Exploratory",
    capsids: "Capsid visibility",
    showAll: "Show all",
    yAxis: "Effective duration after one dose (days)",
    xAxis: "Organ specificity log10 (target ISF concentration AUC / weighted off-target AUC)",
    preferred: "Preferred: precise + durable",
    preciseShort: "Precise but short-lived",
    durableOff: "Durable with off-target burden",
    lowUtility: "Low specificity / short-lived",
    modelDerived: "ODE-derived",
    surrogate: "Surrogate model",
    specificity: "Organ specificity",
    duration: "Effective duration",
    delivery: "Peak post-barrier delivery",
    restoration: "Predicted protein restoration",
    route: "Route",
    species: "Evidence species",
    targetCell: "Target cells",
    exposure: "Target-organ AUC share",
    tmax: "ISF Tmax",
    halfLife: "Episome half-life prior",
    onset: "Time to 65% threshold",
    peakDay: "Protein peak time",
    transduction: "Cellular transduction",
    targetDepth: "CNS target depth",
    layerExposure: "Target-layer exposure share",
    nativeChain: "Native ode1.0 intracellular states",
    multilevelChain: "Post-BBB three-depth + intracellular ODE",
    reducedChain: "ISF-driven reduced surrogate",
    source: "Open primary evidence",
    hi: "Dosage evidence",
    sineup: "SINEUP prerequisite",
    modelProof: "Mathematical model evidence",
    earlyPk: "Early-delivery ODE",
    earlyPkText: "The blood–vascular–ISF PBPK is solved first. CNS then resolves superficial, cortical and deep transport plus intracellular transduction before integrating the Epi–SINEUP RNA–target protein ODE.",
    efficiencyFormula: "Delivery efficiency",
    specificityFormula: "Specificity",
    persistenceFormula: "Persistence",
    balance: "Maximum mass-balance error",
    generated: "Model generated",
    disclaimer: "Research design tool, not clinical dosing guidance. Eye currently uses a local-route barrier surrogate; marrow and peripheral-nerve records expose PBPK compartments still to be added.",
    expressionTitle: "Normal-tissue expression evidence",
    targetTpm: "Target-organ GTEx median TPM",
    topTissue: "Top GTEx tissue",
    tissueTau: "Tissue specificity τ",
    hpaClass: "HPA tissue class",
    expressionNote: "Adult bulk RNA is only a prior for normal-transcript presence; it is not cell-type expression, developmental expression, or therapeutic efficacy.",
    combinationTitle: "ODE regimen-combination screen",
    singlePlan: "Single plan preferred",
    dualPlan: "Exploratory dual plan",
    coverageScore: "Weighted coverage score",
    coveredGenes: "Covered genes",
    unmodeledGenes: "No matching organ ODE",
    combinationNote: "A dual plan appears only when coverage improves after penalties for invasiveness, a second administration and capsid immune risk. Each agent is currently solved independently at the full reference dose; total-dose splitting and immune-interaction ODEs are not yet modeled, so this is not clinical combination-dosing guidance.",
    dataSource: "Database source",
  },
};

const organName: Record<Organ, { zh: string; en: string }> = {
  CNS: { zh: "中枢神经系统", en: "CNS" },
  Liver: { zh: "肝脏", en: "Liver" },
  Eye: { zh: "眼部", en: "Eye" },
  Heart: { zh: "心脏", en: "Heart" },
  Muscle: { zh: "肌肉 / 外周组织", en: "Muscle / peripheral tissue" },
  Kidney: { zh: "肾脏", en: "Kidney" },
};

const cnsProfileName: Record<string, { zh: string; en: string }> = {
  cortical_excitatory: { zh: "浅层皮层兴奋性神经元", en: "Superficial cortical excitatory neurons" },
  cortical_inhibitory: { zh: "皮层抑制性中间神经元", en: "Cortical inhibitory interneurons" },
  cortical_projection: { zh: "中深层皮层投射神经元", en: "Mid/deep cortical projection neurons" },
  synaptic_neuron: { zh: "皮层突触神经元", en: "Cortical synaptic neurons" },
  deep_striatal: { zh: "深部纹状体神经元", en: "Deep striatal neurons" },
  hypothalamic: { zh: "深部下丘脑神经元", en: "Deep hypothalamic neurons" },
  broad_neuronal: { zh: "跨层广泛神经元", en: "Broad neuronal populations across layers" },
  neural_progenitor: { zh: "神经祖细胞区", en: "Neural progenitor zones" },
};

const mechanismName = {
  zh: { haploinsufficiency: "单倍剂量不足", "whole-gene deletion": "完整基因缺失", "contiguous deletion": "连续基因缺失" },
  en: { haploinsufficiency: "Haploinsufficiency", "whole-gene deletion": "Whole-gene deletion", "contiguous deletion": "Contiguous deletion" },
};

const evidenceClass: Record<Evidence, string> = { strong: "evidence-strong", medium: "evidence-medium", exploratory: "evidence-limited" };
const organToHumanParent: Record<Organ, string | null> = {
  CNS: "brain", Liver: "liver", Eye: null, Heart: "heart", Muscle: "muscle", Kidney: "kidney",
};
const routeEligibility: Record<Organ, string[]> = {
  CNS: ["iv", "intrathecal", "intracisternal", "intracerebroventricular"],
  Liver: ["iv"],
  Eye: [],
  Heart: ["iv"],
  Muscle: ["iv", "intramuscular"],
  Kidney: ["iv"],
};
const routeBurden: Record<string, number> = {
  iv: 0.01, intramuscular: 0.04, inhaled: 0.03, intrathecal: 0.08,
  intracisternal: 0.12, intracerebroventricular: 0.15,
};
const xTicks = [-3, -2, -1, 0, 1];
const yTicks = [0, 90, 180, 365, 730];
const xMin = -3;
const xMax = 1.5;
const yMax = 730;
const specificityThresholdPct = ((0 - xMin) / (xMax - xMin)) * 100;
const durationThresholdPct = (180 / yMax) * 100;

function fmt(value: number, digits = 2) {
  return value.toLocaleString("en-US", { maximumFractionDigits: digits });
}

function fmtScientific(value: number) {
  if (value === 0) return "0.000";
  const superscript: Record<string, string> = {
    "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
    "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹",
  };
  const exponent = Math.floor(Math.log10(Math.abs(value)));
  const power = String(exponent).split("").map((character) => superscript[character] ?? character).join("");
  return `${(value / 10 ** exponent).toFixed(3)} × 10${power}`;
}

const generatedAtUtc = `${modelPayload.generated_at.slice(0, 19).replace("T", " ")} UTC`;

function expressionTargetScore(record: DiseaseRecord) {
  const expression = expressionByGene[record.gene];
  const parent = record.organ === "CNS" ? "CNS" : record.organ;
  const target = expression?.organ_median_tpm[parent];
  const top = expression?.top_modeled_organ_tpm;
  if (!(target !== null && target !== undefined && top && top > 0)) return 0.55;
  return Math.min(1, Math.max(0.15, Math.log1p(target) / Math.log1p(top)));
}

function recommendRegimen(records: DiseaseRecord[]): RegimenRecommendation {
  const modeled = records.filter((record) => organToHumanParent[record.organ] !== null);
  const unmodeled = records.filter((record) => organToHumanParent[record.organ] === null);
  if (!humanRouteSummary.length || !modeled.length) {
    return { agents: [], score: 0, modeledGenes: modeled.map((record) => record.gene), unmodeledGenes: unmodeled.map((record) => record.gene) };
  }

  const organScale = new Map<string, { min: number; max: number }>();
  for (const parent of new Set(modeled.map((record) => organToHumanParent[record.organ] as string))) {
    const values = humanRouteSummary.map((option) => {
      const metric = option.organs[parent];
      return metric ? Math.log10(metric.auc_isf_concentration_vg_h_ml + 1) + 0.35 * Math.log10(metric.peak_protein_au + 1) : 0;
    });
    organScale.set(parent, { min: Math.min(...values), max: Math.max(...values) });
  }

  const geneScore = (option: HumanRouteSummary, record: DiseaseRecord) => {
    if (!routeEligibility[record.organ].includes(option.route_id)) return 0;
    const parent = organToHumanParent[record.organ];
    if (!parent || !option.organs[parent]) return 0;
    const metric = option.organs[parent];
    const raw = Math.log10(metric.auc_isf_concentration_vg_h_ml + 1) + 0.35 * Math.log10(metric.peak_protein_au + 1);
    const scale = organScale.get(parent)!;
    const exposure = (raw - scale.min) / Math.max(scale.max - scale.min, 1e-12);
    const specificity = Math.min(1, metric.exposure_share_pct / 45);
    return (0.72 * exposure + 0.28 * specificity) * (0.65 + 0.35 * expressionTargetScore(record));
  };

  const candidates = humanRouteSummary.filter((option) => modeled.some((record) => routeEligibility[record.organ].includes(option.route_id)));
  const scoreAgents = (agents: HumanRouteSummary[]) => {
    const coverage = modeled.map((record) => Math.max(...agents.map((agent) => geneScore(agent, record))));
    const immunePenalty = agents.length > 1 ? 0.08 + (agents[0].capsid_id !== agents[1].capsid_id ? 0.04 : 0.02) : 0;
    const burden = agents.reduce((sum, agent) => sum + (routeBurden[agent.route_id] ?? 0.08), 0);
    return coverage.reduce((sum, value) => sum + value, 0) / coverage.length - immunePenalty - burden;
  };

  let bestAgents: HumanRouteSummary[] = [];
  let bestScore = -Infinity;
  for (const candidate of candidates) {
    const score = scoreAgents([candidate]);
    if (score > bestScore) { bestAgents = [candidate]; bestScore = score; }
  }
  const bestSingleScore = bestScore;
  for (let left = 0; left < candidates.length; left += 1) {
    for (let right = left + 1; right < candidates.length; right += 1) {
      const pair = [candidates[left], candidates[right]];
      const score = scoreAgents(pair);
      if (score > bestScore && score > bestSingleScore + 0.10) { bestAgents = pair; bestScore = score; }
    }
  }
  const agents = bestAgents.map((agent) => ({
    ...agent,
    coveredGenes: modeled.filter((record) => geneScore(agent, record) === Math.max(...bestAgents.map((item) => geneScore(item, record)))).map((record) => record.gene),
  }));
  return {
    agents,
    score: Math.max(0, Math.min(1, bestScore)),
    modeledGenes: modeled.map((record) => record.gene),
    unmodeledGenes: unmodeled.map((record) => record.gene),
  };
}

export function DesignSpaceApp() {
  const diseaseListRef = useRef<HTMLDivElement>(null);
  const [language, setLanguage] = useState<Language>("zh");
  const [activeView, setActiveView] = useState<"design" | "heatmap">("heatmap");
  const t = copy[language];
  const [query, setQuery] = useState("");
  const [diseaseId, setDiseaseId] = useState("whs-nsd2");
  const [expandedDisease, setExpandedDisease] = useState("Wolf-Hirschhorn syndrome");
  const [capsidId, setCapsidId] = useState("aav9");
  const [visibleCapsids, setVisibleCapsids] = useState<string[]>(capsidOrder);
  const disease = diseases.find((item) => item.id === diseaseId) ?? diseases[0];
  const organPoints = disease.organ === "CNS" && disease.cnsProfile
    ? cnsProfileResults.filter((point) => point.cns_profile === disease.cnsProfile)
    : modelResults.filter((point) => point.target === disease.organ);
  const visiblePoints = organPoints.filter((point) => visibleCapsids.includes(point.capsid_id));
  const selectedPoint = organPoints.find((point) => point.capsid_id === capsidId) ?? organPoints[0];
  const selectedDiseaseGenes = useMemo(
    () => diseases.filter((item) => item.name.en === disease.name.en),
    [disease.name.en],
  );
  const expression = expressionByGene[disease.gene];
  const expressionOrgan = disease.organ === "CNS" ? "CNS" : disease.organ;
  const regimen = useMemo(() => recommendRegimen(selectedDiseaseGenes), [selectedDiseaseGenes]);

  const diseaseGroups = useMemo(() => {
    const normalized = query.trim().toLowerCase();
    const filtered = normalized ? diseases.filter((item) =>
      `${item.name.zh} ${item.name.en} ${item.gene} ${item.locus} ${item.phenotype.zh} ${item.phenotype.en}`.toLowerCase().includes(normalized),
    ) : diseases;
    const groups = new Map<string, DiseaseRecord[]>();
    for (const item of filtered) {
      const key = item.name.en;
      groups.set(key, [...(groups.get(key) ?? []), item]);
    }
    return Array.from(groups.entries()).map(([key, genes]) => ({ key, genes }));
  }, [query]);

  function selectDisease(item: DiseaseRecord) {
    setDiseaseId(item.id);
    setExpandedDisease(item.name.en);
    setCapsidId("aav9");
  }

  function toggleCapsid(id: string) {
    if (visibleCapsids.includes(id)) {
      if (visibleCapsids.length === 1) return;
      const next = visibleCapsids.filter((item) => item !== id);
      setVisibleCapsids(next);
      if (capsidId === id) setCapsidId(next[0]);
    } else {
      setVisibleCapsids([...visibleCapsids, id]);
    }
  }

  function scrollDiseaseLibrary(event: WheelEvent<HTMLElement>) {
    const list = diseaseListRef.current;
    if (!list || list.contains(event.target as Node) || list.scrollHeight <= list.clientHeight) return;

    const maxScrollTop = list.scrollHeight - list.clientHeight;
    const nextScrollTop = Math.max(0, Math.min(maxScrollTop, list.scrollTop + event.deltaY));
    if (nextScrollTop === list.scrollTop) return;

    list.scrollTop = nextScrollTop;
    event.preventDefault();
  }

  return (
    <main className="app-shell">
      <header className="topbar">
        <div className="brand-lockup">
          <div className="brand-mark" aria-hidden="true"><Dna size={20} /></div>
          <div><strong>SINEUP Delivery Atlas</strong><span>{t.subtitle}</span></div>
        </div>
        <div className="topbar-actions">
          <nav className="view-switch" aria-label="Primary view">
            <button type="button" className={activeView === "design" ? "active" : ""} onClick={() => setActiveView("design")}><Target size={14} />{t.designNav}</button>
            <button type="button" className={activeView === "heatmap" ? "active" : ""} onClick={() => setActiveView("heatmap")}><Activity size={14} />{t.heatmapNav}</button>
          </nav>
          <div className="status-line"><span className="status-dot" />{t.status}</div>
          <div className="language-switch" aria-label="Language">
            <Languages size={15} aria-hidden="true" />
            <button type="button" className={language === "zh" ? "active" : ""} onClick={() => setLanguage("zh")}>中文</button>
            <button type="button" className={language === "en" ? "active" : ""} onClick={() => setLanguage("en")}>EN</button>
          </div>
        </div>
      </header>

      {activeView === "heatmap" ? <OrganHeatmap language={language} /> : <>
      <section className="context-strip" aria-label="Design context">
        <div><span>{t.disease}</span><strong>{disease.name[language]}</strong></div>
        <ChevronRight size={16} aria-hidden="true" />
        <div><span>{t.gene}</span><strong>{disease.gene}</strong></div>
        <ChevronRight size={16} aria-hidden="true" />
        <div><span>{t.location}</span><strong>{organName[disease.organ][language]} · {disease.cnsProfile ? cnsProfileName[disease.cnsProfile][language] : disease.targetCell[language]}</strong></div>
      </section>

      <div className="workspace">
        <aside className="disease-panel" aria-label={t.library} onWheel={scrollDiseaseLibrary}>
          <div className="panel-heading"><div><Database size={17} /><strong>{t.library}</strong></div><span>{diseaseGroups.length}</span></div>
          <label className="search-field">
            <Search size={17} aria-hidden="true" />
            <span className="sr-only">{t.search}</span>
            <input value={query} onChange={(event) => setQuery(event.target.value)} placeholder={t.search} />
          </label>
          <div className="disease-list" ref={diseaseListRef} tabIndex={0} aria-label={t.library}>
            {diseaseGroups.map((group) => {
              const expanded = Boolean(query.trim()) || expandedDisease === group.key;
              const active = group.genes.some((item) => item.id === disease.id);
              return <div className={`disease-group ${active ? "active" : ""}`} key={group.key}>
                <button type="button" className="disease-group-row" aria-expanded={expanded} onClick={() => setExpandedDisease(expanded && !query.trim() ? "" : group.key)}>
                  <span className="disease-copy"><strong>{group.genes[0].name[language]}</strong><small>{group.genes.length} {t.geneCount}</small></span>
                  {expanded ? <ChevronDown size={16} aria-hidden="true" /> : <ChevronRight size={16} aria-hidden="true" />}
                </button>
                {expanded && <div className="gene-sublist">
                  {group.genes.map((item) => <button type="button" key={item.id} className={item.id === disease.id ? "disease-row selected" : "disease-row"} onClick={() => selectDisease(item)}>
                    <span className="gene-symbol">{item.gene}</span>
                    <span className="disease-copy"><strong>{item.targetCell[language]}</strong><small>{item.locus} · {mechanismName[language][item.mechanism]}</small></span>
                    <ChevronRight size={15} aria-hidden="true" />
                  </button>)}
                </div>}
              </div>;
            })}
            {diseaseGroups.length === 0 && <p className="empty-state">{t.noResult}</p>}
          </div>
          <div className="data-note"><Info size={16} /><p>{t.eligibility}</p></div>
        </aside>

        <section className="analysis-panel">
          <div className="analysis-header">
            <div>
              <span className="eyebrow"><Target size={15} /> {t.designSpace}</span>
              <h1>{disease.gene} · {organName[disease.organ][language]} {t.candidate}</h1>
              <p>{t.chartHelp}</p>
            </div>
            <div className="evidence-legend" aria-label="Evidence legend">
              <span><i className="evidence-strong" />{t.strong}</span><span><i className="evidence-medium" />{t.medium}</span><span><i className="evidence-limited" />{t.exploratory}</span>
            </div>
          </div>

          <div className="capsid-toolbar">
            <span>{t.capsids}</span>
            <div className="capsid-toggles">
              {capsidOrder.map((id) => {
                const active = visibleCapsids.includes(id);
                const evidence = organPoints.find((point) => point.capsid_id === id)?.evidence ?? "exploratory";
                return <button type="button" key={id} className={active ? "active" : ""} onClick={() => toggleCapsid(id)} aria-pressed={active}>
                  {active ? <Eye size={14} /> : <EyeOff size={14} />}<i className={evidenceClass[evidence]} />{capsidLabels[id]}
                </button>;
              })}
              <button type="button" className="show-all" onClick={() => setVisibleCapsids(capsidOrder)}><Check size={14} />{t.showAll}</button>
            </div>
          </div>

          <div className="chart-and-detail">
            <div className="chart-wrap">
              <div className="chart-y-title">{t.yAxis}</div>
              <div className="chart-area" data-testid="design-chart">
                <div className="utility-zone durable-off" style={{ left: 0, right: `${100 - specificityThresholdPct}%`, top: 0, bottom: `${durationThresholdPct}%` }}><span>{t.durableOff}</span></div>
                <div className="utility-zone preferred" style={{ left: `${specificityThresholdPct}%`, right: 0, top: 0, bottom: `${durationThresholdPct}%` }}><span>{t.preferred}</span></div>
                <div className="utility-zone low-utility" style={{ left: 0, right: `${100 - specificityThresholdPct}%`, top: `${100 - durationThresholdPct}%`, bottom: 0 }}><span>{t.lowUtility}</span></div>
                <div className="utility-zone precise-short" style={{ left: `${specificityThresholdPct}%`, right: 0, top: `${100 - durationThresholdPct}%`, bottom: 0 }}><span>{t.preciseShort}</span></div>
                {yTicks.map((tick) => <div className="grid-line horizontal" key={`y-${tick}`} style={{ bottom: `${(tick / yMax) * 100}%` }}><span>{tick}</span></div>)}
                {xTicks.map((tick) => <div className="grid-line vertical" key={`x-${tick}`} style={{ left: `${((tick - xMin) / (xMax - xMin)) * 100}%` }}><span>{tick.toFixed(1)}</span></div>)}
                {visiblePoints.map((point) => {
                  const left = Math.min(100, Math.max(0, ((point.specificity_log10 - xMin) / (xMax - xMin)) * 100));
                  const bottom = Math.min(100, Math.max(0, (point.effective_duration_days / yMax) * 100));
                  return <button type="button" key={point.capsid_id} className={`plot-point ${evidenceClass[point.evidence]} ${selectedPoint.capsid_id === point.capsid_id ? "selected" : ""}`} style={{ left: `${left}%`, bottom: `${bottom}%` }} onClick={() => setCapsidId(point.capsid_id)} aria-label={`${point.capsid}, ${t.specificity} ${fmt(point.specificity_log10)}, ${t.duration} ${fmt(point.effective_duration_days)} days`} title={point.capsid}><span>{point.capsid}</span></button>;
                })}
              </div>
              <div className="chart-x-title">{t.xAxis}</div>
            </div>

            <aside className="detail-panel" aria-label="Modeled capsid details">
              <div className="detail-title-row">
                <div><span className={`evidence-pill ${selectedPoint.model_status === "ode-derived" ? "model-derived" : "model-surrogate"}`}>{selectedPoint.model_status === "ode-derived" ? t.modelDerived : t.surrogate}</span><h2>{selectedPoint.capsid}</h2></div>
                <Activity size={22} aria-hidden="true" />
              </div>
              <dl className="metric-list">
                <div><dt>{t.specificity}</dt><dd>{fmt(selectedPoint.specificity_log10)}</dd></div>
                <div><dt>{t.duration}</dt><dd>{fmt(selectedPoint.effective_duration_days, 0)} d</dd></div>
                <div><dt>{t.delivery}</dt><dd>{fmt(selectedPoint.peak_post_barrier_delivery_pct, 4)}%</dd></div>
                <div><dt>{t.restoration}</dt><dd>{fmt(selectedPoint.predicted_protein_restoration_pct, 1)}%</dd></div>
              </dl>
              <div className="parameter-block">
                <div><span>{t.route}</span><strong>{selectedPoint.route}</strong></div>
                <div><span>{t.species}</span><strong>{selectedPoint.species}</strong></div>
                <div><span>{t.targetCell}</span><strong>{disease.targetCell[language]}</strong></div>
                <div><span>{t.exposure}</span><strong>{fmt(selectedPoint.target_exposure_share_pct, 2)}%</strong></div>
                <div><span>{t.tmax}</span><strong>{fmt(selectedPoint.tmax_h, 2)} h</strong></div>
                <div><span>{t.halfLife}</span><strong>{fmt(selectedPoint.episome_half_life_days_prior, 0)} d</strong></div>
                <div><span>{t.onset}</span><strong>{selectedPoint.therapeutic_onset_days > 0 ? `${fmt(selectedPoint.therapeutic_onset_days, 1)} d` : "—"}</strong></div>
                <div><span>{t.peakDay}</span><strong>{fmt(selectedPoint.peak_restoration_day, 1)} d</strong></div>
                {selectedPoint.cns_depth_mm !== undefined && <div><span>{t.targetDepth}</span><strong>{fmt(selectedPoint.cns_depth_mm, 1)} mm</strong></div>}
                {selectedPoint.cns_target_layer_auc_fraction_pct !== undefined && <div><span>{t.layerExposure}</span><strong>{fmt(selectedPoint.cns_target_layer_auc_fraction_pct, 1)}%</strong></div>}
                <div><span>{t.transduction}</span><strong>{selectedPoint.cns_profile ? t.multilevelChain : selectedPoint.model_status === "ode-derived" ? t.nativeChain : t.reducedChain}</strong></div>
              </div>
              <a className="source-link" href={selectedPoint.source} target="_blank" rel="noreferrer"><BookOpen size={16} />{t.source}<ExternalLink size={14} /></a>
            </aside>
          </div>

          <div className="disease-evidence-band">
            <div className="evidence-item"><span><Dna size={16} />{t.hi}</span><strong>{disease.hiEvidence[language]}</strong></div>
            <div className="evidence-item"><span><FlaskConical size={16} />{t.sineup}</span><strong>{disease.residual[language]}</strong></div>
            <div className="caution-item"><CircleAlert size={17} /><p>{disease.caution[language]}</p></div>
          </div>

          <div className="evidence-design-grid">
            <section className="expression-card" aria-label={t.expressionTitle}>
              <div className="evidence-card-title"><Database size={17} /><strong>{disease.gene} · {t.expressionTitle}</strong></div>
              <dl>
                <div><dt>{t.targetTpm}</dt><dd>{expression?.organ_median_tpm[expressionOrgan] !== null && expression?.organ_median_tpm[expressionOrgan] !== undefined ? fmt(expression.organ_median_tpm[expressionOrgan] as number, 2) : "—"}</dd></div>
                <div><dt>{t.topTissue}</dt><dd>{expression?.top_gtex_tissue?.replaceAll("_", " ") ?? "—"} {expression?.top_gtex_tissue_tpm ? `· ${fmt(expression.top_gtex_tissue_tpm, 2)} TPM` : ""}</dd></div>
                <div><dt>{t.tissueTau}</dt><dd>{expression?.tissue_tau !== null && expression?.tissue_tau !== undefined ? fmt(expression.tissue_tau, 3) : "—"}</dd></div>
                <div><dt>{t.hpaClass}</dt><dd>{expression?.hpa?.tissue_specificity ?? "—"}</dd></div>
              </dl>
              <p>{t.expressionNote}</p>
              <div className="database-links"><a href="https://gtexportal.org/api/v2/docs" target="_blank" rel="noreferrer">GTEx v10 <ExternalLink size={12} /></a>{expression?.hpa && <a href={expression.hpa.entry_url} target="_blank" rel="noreferrer">Human Protein Atlas <ExternalLink size={12} /></a>}</div>
            </section>

            <section className="combination-card" aria-label={t.combinationTitle}>
              <div className="evidence-card-title"><FlaskConical size={17} /><strong>{t.combinationTitle}</strong><span>{regimen.agents.length > 1 ? t.dualPlan : t.singlePlan}</span></div>
              <div className="regimen-score"><span>{t.coverageScore}</span><strong>{Math.round(regimen.score * 100)} / 100</strong></div>
              <div className="regimen-agents">
                {regimen.agents.map((agent, index) => <div className="regimen-agent" key={`${agent.route_id}-${agent.capsid_id}`}>
                  <span>{String.fromCharCode(65 + index)}</span><div><strong>{agent.capsid} · {language === "zh" ? agent.route_label_zh : agent.route_label}</strong><small>{t.coveredGenes}: {agent.coveredGenes.join(", ") || "—"}</small></div><a href={agent.evidence_source} target="_blank" rel="noreferrer" aria-label={t.source}><ExternalLink size={14} /></a>
                </div>)}
                {!regimen.agents.length && <div className="regimen-empty">{t.unmodeledGenes}: {regimen.unmodeledGenes.join(", ")}</div>}
              </div>
              {regimen.unmodeledGenes.length > 0 && <div className="unmodeled-line"><CircleAlert size={14} />{t.unmodeledGenes}: {regimen.unmodeledGenes.join(", ")}</div>}
              <p>{t.combinationNote}</p>
            </section>
          </div>

          <section className="model-proof" aria-label={t.modelProof}>
            <div className="model-proof-title"><Sigma size={18} /><div><strong>{t.modelProof}</strong><span>{t.earlyPkText}</span></div></div>
            <div className="formula-grid">
              <div><span>{t.efficiencyFormula}</span><code>100 × max(A_target,ISF) / Dose</code><strong>{fmt(selectedPoint.peak_post_barrier_delivery_pct, 4)}%</strong></div>
              <div><span>{t.specificityFormula}</span><code>log10(AUC_target / AUC_off-target)</code><strong>{fmt(selectedPoint.specificity_log10)}</strong></div>
              <div><span>{t.persistenceFormula}</span><code>dE/dt → dSINEUP/dt → dP/dt</code><strong>{fmt(selectedPoint.effective_duration_days, 0)} d</strong></div>
            </div>
            <div className="model-meta"><span>{t.balance}: <strong>{fmtScientific(selectedPoint.max_mass_balance_error)}</strong></span><span>{t.generated}: <strong>{generatedAtUtc}</strong></span></div>
          </section>

          <footer className="method-footer"><CircleAlert size={15} />{t.disclaimer}</footer>
        </section>
      </div>
      </>}
    </main>
  );
}
