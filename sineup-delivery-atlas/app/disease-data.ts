export type Language = "zh" | "en";
export type Organ = "CNS" | "Liver" | "Eye" | "Heart" | "Muscle" | "Kidney";
export type LocalizedText = { zh: string; en: string };
export type CnsProfile = "cortical_excitatory" | "cortical_inhibitory" | "cortical_projection" | "synaptic_neuron" | "deep_striatal" | "hypothalamic" | "broad_neuronal" | "neural_progenitor";

export type DiseaseRecord = {
  id: string;
  name: LocalizedText;
  gene: string;
  locus: string;
  organ: Organ;
  cnsProfile?: CnsProfile;
  targetCell: LocalizedText;
  mechanism: "haploinsufficiency" | "whole-gene deletion" | "contiguous deletion";
  hiEvidence: LocalizedText;
  phenotype: LocalizedText;
  residual: LocalizedText;
  caution: LocalizedText;
  source: string;
};

const sufficient = {
  zh: "ClinGen：单倍剂量不足证据充分",
  en: "ClinGen: sufficient evidence for haploinsufficiency",
};

const clingen = (hgnc: string) => `https://search.clinicalgenome.org/kb/gene-dosage/${hgnc}`;
const digeorgeRegion = "https://search.clinicalgenome.org/kb/gene-dosage/region/ISCA-37446";
const williamsReview = "https://www.ncbi.nlm.nih.gov/books/NBK1249/";
const wagrReview = "https://www.ncbi.nlm.nih.gov/books/NBK55674/";
const millerDiekerReview = "https://www.ncbi.nlm.nih.gov/books/NBK5189/";
const oneP36Review = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4730066/";

export const diseases: DiseaseRecord[] = [
  {
    id: "whs-nsd2", name: { zh: "Wolf-Hirschhorn 综合征", en: "Wolf-Hirschhorn syndrome" }, gene: "NSD2", locus: "4p16.3", organ: "CNS", cnsProfile: "neural_progenitor",
    targetCell: { zh: "神经祖细胞 / 神经元", en: "Neural progenitors / neurons" }, mechanism: "contiguous deletion", hiEvidence: sufficient,
    phenotype: { zh: "神经发育、颅面与生长表型", en: "Neurodevelopmental, craniofacial and growth phenotypes" },
    residual: { zh: "另一条同源染色体通常保留正常 NSD2 转录本", en: "The homologous chromosome usually retains a normal NSD2 transcript" },
    caution: { zh: "发育期表型可能无法由出生后蛋白恢复逆转。", en: "Developmental phenotypes may not be reversible by postnatal protein restoration." }, source: clingen("HGNC:12766"),
  },
  {
    id: "whs-letm1", name: { zh: "Wolf-Hirschhorn 综合征", en: "Wolf-Hirschhorn syndrome" }, gene: "LETM1", locus: "4p16.3", organ: "CNS", cnsProfile: "broad_neuronal",
    targetCell: { zh: "广泛 CNS / 神经元", en: "Broad CNS / neurons" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "WHS 癫痫候选基因；ClinGen HI 证据不足", en: "WHS seizure candidate; no ClinGen HI evidence" },
    phenotype: { zh: "癫痫与线粒体稳态候选表型", en: "Candidate seizure and mitochondrial-homeostasis phenotypes" },
    residual: { zh: "需确认正常等位基因转录本和目标 isoform", en: "Normal-allele transcript and target isoform must be confirmed" },
    caution: { zh: "属于候选贡献基因，不应解释为 WHS 癫痫的单一原因。", en: "A candidate contributor, not a proven single cause of WHS seizures." }, source: clingen("HGNC:6556"),
  },
  {
    id: "whs-msx1", name: { zh: "Wolf-Hirschhorn 综合征", en: "Wolf-Hirschhorn syndrome" }, gene: "MSX1", locus: "4p16.2", organ: "Muscle",
    targetCell: { zh: "颅颌面间充质 / 成骨相关细胞", en: "Craniofacial mesenchyme / osteogenic cells" }, mechanism: "contiguous deletion", hiEvidence: sufficient,
    phenotype: { zh: "牙齿与颅颌面发育表型", en: "Dental and craniofacial developmental phenotypes" },
    residual: { zh: "正常等位基因通常仍可转录", en: "The normal allele is generally still transcribed" },
    caution: { zh: "主要作用窗口在发育期，成年系统给药价值有限。", en: "The key window is developmental; adult systemic dosing may have limited value." }, source: clingen("HGNC:7391"),
  },
  {
    id: "whs-fgfr3", name: { zh: "Wolf-Hirschhorn 综合征", en: "Wolf-Hirschhorn syndrome" }, gene: "FGFR3", locus: "4p16.3", organ: "Muscle",
    targetCell: { zh: "软骨细胞 / 骨骼生长板", en: "Chondrocytes / skeletal growth plate" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "WHS 区域候选；ClinGen HI 证据不足", en: "WHS-region candidate; no ClinGen HI evidence" },
    phenotype: { zh: "骨骼生长候选表型", en: "Candidate skeletal-growth phenotype" },
    residual: { zh: "需测定正常等位基因表达及剂量反应", en: "Normal-allele expression and dose response require measurement" },
    caution: { zh: "FGFR3 剂量调节复杂，过度增强也可能有害。", en: "FGFR3 dosage is complex and over-restoration may also be harmful." }, source: clingen("HGNC:3690"),
  },
  {
    id: "whs-cplx1", name: { zh: "Wolf-Hirschhorn 综合征", en: "Wolf-Hirschhorn syndrome" }, gene: "CPLX1", locus: "4p16.3", organ: "CNS", cnsProfile: "synaptic_neuron",
    targetCell: { zh: "突触神经元", en: "Synaptic neurons" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "WHS 神经表型候选，尚无 ClinGen HI 结论", en: "WHS neurological candidate without a ClinGen HI assertion" },
    phenotype: { zh: "突触释放与癫痫候选表型", en: "Candidate synaptic-release and seizure phenotypes" },
    residual: { zh: "需先验证患者细胞中的剩余 mRNA", en: "Residual mRNA must first be verified in patient cells" },
    caution: { zh: "目前应作为研究候选点，而非治疗优先级结论。", en: "Currently a research candidate, not a treatment-priority conclusion." }, source: "https://pmc.ncbi.nlm.nih.gov/articles/PMC3953918/",
  },
  {
    id: "chd8", name: { zh: "CHD8 相关神经发育障碍", en: "CHD8-related neurodevelopmental disorder" }, gene: "CHD8", locus: "14q11.2", organ: "CNS", cnsProfile: "cortical_excitatory",
    targetCell: { zh: "皮层神经元", en: "Cortical neurons" }, mechanism: "haploinsufficiency", hiEvidence: sufficient,
    phenotype: { zh: "神经发育与自闭症相关表型", en: "Neurodevelopmental and autism-associated phenotypes" },
    residual: { zh: "靶 isoform 与起始密码子需逐构建设计", en: "Target isoform and start codon require construct-specific design" },
    caution: { zh: "已有 SINEUP 概念验证，但人体治疗窗口未知。", en: "SINEUP proof-of-concept exists, but the human treatment window is unknown." }, source: clingen("HGNC:20153"),
  },
  {
    id: "scn1a", name: { zh: "Dravet 综合征", en: "Dravet syndrome" }, gene: "SCN1A", locus: "2q24.3", organ: "CNS", cnsProfile: "cortical_inhibitory",
    targetCell: { zh: "GABA 能抑制性中间神经元", en: "GABAergic inhibitory interneurons" }, mechanism: "haploinsufficiency", hiEvidence: sufficient,
    phenotype: { zh: "癫痫与神经发育表型", en: "Epilepsy and neurodevelopmental phenotypes" }, residual: { zh: "仅适用于仍产生功能性正常转录本的基因型", en: "Applicable only when a functional normal transcript remains" },
    caution: { zh: "必须加入细胞类型特异性，单纯 CNS AUC 不足以预测疗效。", en: "Cell-type specificity is essential; CNS AUC alone cannot predict efficacy." }, source: clingen("HGNC:10585"),
  },
  {
    id: "shank3", name: { zh: "Phelan-McDermid 综合征", en: "Phelan-McDermid syndrome" }, gene: "SHANK3", locus: "22q13.33", organ: "CNS", cnsProfile: "synaptic_neuron",
    targetCell: { zh: "兴奋性突触神经元", en: "Excitatory synaptic neurons" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "发育迟缓、语言与自闭症相关表型", en: "Developmental delay, speech and autism-associated phenotypes" }, residual: { zh: "完整基因缺失时另一条染色体仍提供靶 mRNA", en: "A whole-gene deletion leaves target mRNA from the homologous chromosome" },
    caution: { zh: "需要严格控制翻译增益以避免突触蛋白过表达。", en: "Translation gain must be controlled to avoid synaptic protein overexpression." }, source: clingen("HGNC:14294"),
  },
  {
    id: "ehmt1", name: { zh: "Kleefstra 综合征", en: "Kleefstra syndrome" }, gene: "EHMT1", locus: "9q34.3", organ: "CNS", cnsProfile: "broad_neuronal",
    targetCell: { zh: "广泛神经元", en: "Broad neuronal populations" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "智力、语言与行为表型", en: "Intellectual, speech and behavioral phenotypes" }, residual: { zh: "需确认患者正常等位基因表达", en: "Normal-allele expression must be confirmed" },
    caution: { zh: "染色质调控蛋白的剂量窗口可能较窄。", en: "Chromatin regulators may have a narrow dosage window." }, source: clingen("HGNC:24650"),
  },
  {
    id: "rai1", name: { zh: "Smith-Magenis 综合征", en: "Smith-Magenis syndrome" }, gene: "RAI1", locus: "17p11.2", organ: "CNS", cnsProfile: "hypothalamic",
    targetCell: { zh: "下丘脑 / 广泛神经元", en: "Hypothalamic / broad neuronal populations" }, mechanism: "contiguous deletion", hiEvidence: sufficient,
    phenotype: { zh: "睡眠、行为与神经发育表型", en: "Sleep, behavioral and neurodevelopmental phenotypes" }, residual: { zh: "正常等位基因通常保留", en: "The normal allele is usually retained" },
    caution: { zh: "治疗终点应按昼夜节律和行为分别定义。", en: "Therapeutic endpoints should separate circadian and behavioral outcomes." }, source: clingen("HGNC:9834"),
  },
  {
    id: "nsd1", name: { zh: "Sotos 综合征", en: "Sotos syndrome" }, gene: "NSD1", locus: "5q35.3", organ: "CNS", cnsProfile: "neural_progenitor",
    targetCell: { zh: "神经祖细胞 / 神经元", en: "Neural progenitors / neurons" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "过度生长与神经发育表型", en: "Overgrowth and neurodevelopmental phenotypes" }, residual: { zh: "需测定剩余等位基因表达", en: "Residual-allele expression requires measurement" },
    caution: { zh: "出生后蛋白恢复未必修正既有生长轨迹。", en: "Postnatal restoration may not correct an established growth trajectory." }, source: clingen("HGNC:14234"),
  },
  {
    id: "satb2", name: { zh: "SATB2 相关综合征", en: "SATB2-associated syndrome" }, gene: "SATB2", locus: "2q33.1", organ: "CNS", cnsProfile: "cortical_projection",
    targetCell: { zh: "皮层投射神经元", en: "Cortical projection neurons" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "语言、认知与颅面表型", en: "Speech, cognitive and craniofacial phenotypes" }, residual: { zh: "正常转录本通常保留约一半剂量", en: "Normal transcript is usually retained at approximately half dosage" },
    caution: { zh: "语言回路治疗窗口和细胞覆盖仍需研究。", en: "The treatment window and cellular coverage of speech circuits remain uncertain." }, source: clingen("HGNC:21637"),
  },
  {
    id: "mef2c", name: { zh: "MEF2C 单倍剂量不足综合征", en: "MEF2C haploinsufficiency syndrome" }, gene: "MEF2C", locus: "5q14.3", organ: "CNS", cnsProfile: "cortical_excitatory",
    targetCell: { zh: "皮层神经元", en: "Cortical neurons" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "严重语言、运动与癫痫表型", en: "Severe speech, motor and epilepsy phenotypes" }, residual: { zh: "需确认剩余转录本与细胞类型", en: "Residual transcript and cell type must be confirmed" },
    caution: { zh: "转录因子过量表达存在安全窗口问题。", en: "Transcription-factor overexpression creates a safety-window concern." }, source: clingen("HGNC:6996"),
  },
  {
    id: "foxp1", name: { zh: "FOXP1 综合征", en: "FOXP1 syndrome" }, gene: "FOXP1", locus: "3p13", organ: "CNS", cnsProfile: "deep_striatal",
    targetCell: { zh: "皮层及纹状体神经元", en: "Cortical and striatal neurons" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "语言与神经发育表型", en: "Speech and neurodevelopmental phenotypes" }, residual: { zh: "正常等位基因转录本通常存在", en: "A normal-allele transcript is generally present" },
    caution: { zh: "需避免影响 FOXP 家族转录网络的剂量平衡。", en: "Dosage balance across the FOXP transcriptional network must be preserved." }, source: clingen("HGNC:3823"),
  },
  {
    id: "tcf4", name: { zh: "Pitt-Hopkins 综合征", en: "Pitt-Hopkins syndrome" }, gene: "TCF4", locus: "18q21.2", organ: "CNS", cnsProfile: "broad_neuronal",
    targetCell: { zh: "广泛神经元", en: "Broad neuronal populations" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "认知、呼吸与癫痫表型", en: "Cognitive, breathing and epilepsy phenotypes" }, residual: { zh: "复杂 isoform 结构要求转录本特异设计", en: "Complex isoform structure requires transcript-specific design" },
    caution: { zh: "不同 TCF4 isoform 的治疗增益可能不同。", en: "Therapeutic gain may differ across TCF4 isoforms." }, source: clingen("HGNC:11634"),
  },
  {
    id: "syngap1", name: { zh: "SYNGAP1 相关智力障碍", en: "SYNGAP1-related intellectual disability" }, gene: "SYNGAP1", locus: "6p21.32", organ: "CNS", cnsProfile: "synaptic_neuron",
    targetCell: { zh: "兴奋性突触神经元", en: "Excitatory synaptic neurons" }, mechanism: "haploinsufficiency", hiEvidence: sufficient,
    phenotype: { zh: "癫痫、认知与行为表型", en: "Epilepsy, cognitive and behavioral phenotypes" }, residual: { zh: "需选择正确脑区与发育阶段 isoform", en: "The correct region- and stage-specific isoform is required" },
    caution: { zh: "突触剂量恢复应以功能阈值而非最大表达为目标。", en: "Synaptic dosage restoration should target a functional threshold, not maximal expression." }, source: clingen("HGNC:11497"),
  },
  {
    id: "dyrk1a", name: { zh: "DYRK1A 综合征", en: "DYRK1A syndrome" }, gene: "DYRK1A", locus: "21q22.13", organ: "CNS", cnsProfile: "neural_progenitor",
    targetCell: { zh: "神经祖细胞 / 神经元", en: "Neural progenitors / neurons" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "小头畸形、认知与运动表型", en: "Microcephaly, cognitive and motor phenotypes" }, residual: { zh: "正常等位基因通常保留", en: "The normal allele is generally retained" },
    caution: { zh: "激酶剂量过高同样可能致病。", en: "Excess kinase dosage may also be pathogenic." }, source: clingen("HGNC:3091"),
  },
  {
    id: "nrnx1", name: { zh: "NRXN1 缺失相关障碍", en: "NRXN1 deletion-associated disorder" }, gene: "NRXN1", locus: "2p16.3", organ: "CNS", cnsProfile: "synaptic_neuron",
    targetCell: { zh: "广泛突触神经元", en: "Broad synaptic neuronal populations" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "神经发育、语言及精神行为表型", en: "Neurodevelopmental, speech and psychiatric phenotypes" }, residual: { zh: "外显子缺失与 isoform 情况必须个体化", en: "Exon deletion and isoform context must be individualized" },
    caution: { zh: "NRXN1 转录本复杂，通用 SINEUP 可能不适用所有缺失。", en: "NRXN1 transcript complexity limits a universal SINEUP design." }, source: clingen("HGNC:8008"),
  },
  {
    id: "pafah1b1", name: { zh: "PAFAH1B1 相关无脑回畸形", en: "PAFAH1B1-related lissencephaly" }, gene: "PAFAH1B1", locus: "17p13.3", organ: "CNS", cnsProfile: "neural_progenitor",
    targetCell: { zh: "神经迁移相关祖细胞", en: "Neuronal migration progenitors" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "皮层迁移与癫痫表型", en: "Cortical migration and epilepsy phenotypes" }, residual: { zh: "正常等位基因通常存在", en: "A normal allele is generally present" },
    caution: { zh: "关键病理发生于胎儿神经迁移期，出生后可逆性很低。", en: "Core pathology occurs during fetal neuronal migration and is unlikely to be reversible postnatally." }, source: clingen("HGNC:8574"),
  },
  {
    id: "ankrd11", name: { zh: "KBG 综合征", en: "KBG syndrome" }, gene: "ANKRD11", locus: "16q24.3", organ: "CNS", cnsProfile: "broad_neuronal",
    targetCell: { zh: "神经元 / 成骨相关细胞", en: "Neurons / osteogenic cells" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "认知、癫痫与骨骼表型", en: "Cognitive, epilepsy and skeletal phenotypes" }, residual: { zh: "需确认无显性负效应转录本", en: "Dominant-negative transcripts must be excluded" },
    caution: { zh: "并非所有 ANKRD11 变异都适合翻译增强。", en: "Not all ANKRD11 variants are suitable for translation enhancement." }, source: clingen("HGNC:21316"),
  },
  {
    id: "pax6", name: { zh: "PAX6 相关无虹膜症", en: "PAX6-related aniridia" }, gene: "PAX6", locus: "11p13", organ: "Eye",
    targetCell: { zh: "角膜缘 / 视网膜相关细胞", en: "Limbal / retinal-associated cells" }, mechanism: "haploinsufficiency", hiEvidence: sufficient,
    phenotype: { zh: "无虹膜与眼表表型", en: "Aniridia and ocular-surface phenotypes" }, residual: { zh: "靶组织需表达正常 PAX6 转录本", en: "Target tissue must express the normal PAX6 transcript" },
    caution: { zh: "眼部细胞类型和局部给药路径必须分开建模。", en: "Ocular cell type and local administration route require separate models." }, source: clingen("HGNC:8620"),
  },
  {
    id: "tbx5", name: { zh: "Holt-Oram 综合征", en: "Holt-Oram syndrome" }, gene: "TBX5", locus: "12q24.21", organ: "Heart",
    targetCell: { zh: "心肌细胞 / 传导系统", en: "Cardiomyocytes / conduction system" }, mechanism: "haploinsufficiency", hiEvidence: sufficient,
    phenotype: { zh: "心脏结构与传导表型", en: "Cardiac structural and conduction phenotypes" }, residual: { zh: "正常等位基因表达量需实验测定", en: "Normal-allele expression requires experimental measurement" },
    caution: { zh: "结构缺陷多在发育期形成，成年表达仅可能影响部分功能终点。", en: "Structural defects form developmentally; adult expression may affect only selected functional endpoints." }, source: clingen("HGNC:11604"),
  },
  {
    id: "eln", name: { zh: "Williams-Beuren 综合征", en: "Williams-Beuren syndrome" }, gene: "ELN", locus: "7q11.23", organ: "Heart",
    targetCell: { zh: "血管平滑肌 / 弹性组织", en: "Vascular smooth muscle / elastic tissue" }, mechanism: "contiguous deletion", hiEvidence: sufficient,
    phenotype: { zh: "主动脉瓣上狭窄与血管表型", en: "Supravalvular aortic stenosis and vascular phenotypes" }, residual: { zh: "另一条染色体通常保留 ELN", en: "The homologous chromosome generally retains ELN" },
    caution: { zh: "细胞外基质形成有发育窗口，心脏器官 AUC 不能替代血管壁分布。", en: "Extracellular-matrix formation has a developmental window; cardiac AUC does not substitute for vessel-wall distribution." }, source: clingen("HGNC:3327"),
  },
  {
    id: "gata3", name: { zh: "HDR 综合征", en: "HDR syndrome" }, gene: "GATA3", locus: "10p14", organ: "Kidney",
    targetCell: { zh: "肾脏远端上皮 / 发育相关细胞", en: "Distal renal epithelium / developmental cells" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "甲状旁腺、耳聋与肾脏表型", en: "Hypoparathyroidism, deafness and renal phenotypes" }, residual: { zh: "不同器官正常等位基因表达需分别测量", en: "Normal-allele expression requires organ-specific measurement" },
    caution: { zh: "这是多器官疾病，单一肾脏靶向不能覆盖全部表型。", en: "This is a multi-organ disease; kidney targeting cannot cover all phenotypes." }, source: clingen("HGNC:4172"),
  },
  {
    id: "ldlr", name: { zh: "家族性高胆固醇血症", en: "Familial hypercholesterolemia" }, gene: "LDLR", locus: "19p13.2", organ: "Liver",
    targetCell: { zh: "肝细胞", en: "Hepatocytes" }, mechanism: "haploinsufficiency", hiEvidence: sufficient,
    phenotype: { zh: "LDL 清除不足", en: "Reduced LDL clearance" }, residual: { zh: "仅适用于仍产生功能性正常转录本的基因型", en: "Applicable only when a functional normal transcript remains" },
    caution: { zh: "显性负效应和无义介导降解基因型应单独排除。", en: "Dominant-negative and nonsense-mediated-decay genotypes require separate exclusion." }, source: clingen("HGNC:6547"),
  },
  {
    id: "rps19", name: { zh: "Diamond-Blackfan 贫血", en: "Diamond-Blackfan anemia" }, gene: "RPS19", locus: "19q13.2", organ: "Liver",
    targetCell: { zh: "造血祖细胞（肝脏仅作系统递送代理）", en: "Hematopoietic progenitors (liver is only a systemic-delivery proxy)" }, mechanism: "haploinsufficiency", hiEvidence: sufficient,
    phenotype: { zh: "红系造血不足", en: "Erythroid failure" }, residual: { zh: "需靶向骨髓造血祖细胞中的正常转录本", en: "The normal transcript must be targeted in marrow progenitors" },
    caution: { zh: "当前 PBPK 没有骨髓室，因此本条仅用于指出模型缺口。", en: "The current PBPK lacks a marrow compartment; this record exposes a model gap." }, source: clingen("HGNC:10402"),
  },
  {
    id: "pmp22", name: { zh: "遗传性压迫易感性神经病", en: "Hereditary neuropathy with liability to pressure palsies" }, gene: "PMP22", locus: "17p12", organ: "Muscle",
    targetCell: { zh: "周围神经 Schwann 细胞", en: "Peripheral-nerve Schwann cells" }, mechanism: "whole-gene deletion", hiEvidence: sufficient,
    phenotype: { zh: "周围神经髓鞘与压迫易感", en: "Peripheral myelin and pressure-palsy susceptibility" }, residual: { zh: "另一等位基因通常保留正常 PMP22", en: "The other allele generally retains normal PMP22" },
    caution: { zh: "肌肉室只是周围组织代理，必须增加周围神经室才能用于决策。", en: "Muscle is only a peripheral-tissue proxy; a peripheral-nerve compartment is required for decisions." }, source: clingen("HGNC:9118"),
  },
  {
    id: "dgs-tbx1", name: { zh: "22q11.2 缺失综合征", en: "22q11.2 deletion syndrome" }, gene: "TBX1", locus: "22q11.21", organ: "Heart",
    targetCell: { zh: "心脏流出道 / 咽弓发育细胞", en: "Cardiac outflow tract / pharyngeal-arch cells" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "22q11.2 区域核心候选基因", en: "Core candidate within the 22q11.2 region" }, phenotype: { zh: "先天性心脏与咽弓表型", en: "Congenital cardiac and pharyngeal-arch phenotypes" },
    residual: { zh: "正常等位基因通常保留，但主要作用窗口在胚胎期", en: "A normal allele usually remains, but the main action window is embryonic" },
    caution: { zh: "不能用出生后心脏表达代表对先天结构缺陷的逆转。", en: "Postnatal cardiac expression should not be interpreted as reversal of congenital anatomy." }, source: digeorgeRegion,
  },
  {
    id: "dgs-dgcr8", name: { zh: "22q11.2 缺失综合征", en: "22q11.2 deletion syndrome" }, gene: "DGCR8", locus: "22q11.21", organ: "CNS", cnsProfile: "broad_neuronal",
    targetCell: { zh: "广泛神经元 / miRNA 加工细胞", en: "Broad neurons / miRNA-processing cells" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "神经表型候选贡献基因", en: "Candidate contributor to neurological phenotypes" }, phenotype: { zh: "认知与精神行为候选表型", en: "Candidate cognitive and neuropsychiatric phenotypes" },
    residual: { zh: "需验证患者神经细胞中的剩余正常转录本", en: "Residual normal transcript must be verified in patient neural cells" },
    caution: { zh: "miRNA 加工蛋白的剂量窗口较窄，需设置表达上限。", en: "A miRNA-processing protein may have a narrow dosage window; an expression ceiling is required." }, source: digeorgeRegion,
  },
  {
    id: "dgs-crkl", name: { zh: "22q11.2 缺失综合征", en: "22q11.2 deletion syndrome" }, gene: "CRKL", locus: "22q11.21", organ: "Heart",
    targetCell: { zh: "心血管发育相关细胞", en: "Cardiovascular developmental cells" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "区域候选修饰基因", en: "Candidate regional modifier" }, phenotype: { zh: "心血管与肾脏表型候选", en: "Candidate cardiovascular and renal phenotypes" },
    residual: { zh: "需先证明目标组织中的剂量不足", en: "Dosage deficiency in the target tissue must first be demonstrated" },
    caution: { zh: "目前不应作为独立治疗靶点结论。", en: "This should not yet be treated as an independent therapeutic target." }, source: digeorgeRegion,
  },
  {
    id: "dgs-comt", name: { zh: "22q11.2 缺失综合征", en: "22q11.2 deletion syndrome" }, gene: "COMT", locus: "22q11.21", organ: "CNS", cnsProfile: "cortical_projection",
    targetCell: { zh: "前额叶与儿茶酚胺代谢细胞", en: "Prefrontal and catecholamine-metabolizing cells" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "精神行为表型候选修饰基因", en: "Candidate modifier of neuropsychiatric phenotypes" }, phenotype: { zh: "执行功能与精神行为候选表型", en: "Candidate executive-function and neuropsychiatric phenotypes" },
    residual: { zh: "必须结合 COMT 单倍型与基线酶活", en: "COMT haplotype and baseline enzyme activity must be considered" },
    caution: { zh: "表达越高不等于越好，儿茶酚胺代谢存在双向风险。", en: "More expression is not necessarily better; catecholamine metabolism has bidirectional risk." }, source: digeorgeRegion,
  },
  {
    id: "dgs-prodh", name: { zh: "22q11.2 缺失综合征", en: "22q11.2 deletion syndrome" }, gene: "PRODH", locus: "22q11.21", organ: "CNS", cnsProfile: "broad_neuronal",
    targetCell: { zh: "神经元 / 脯氨酸代谢细胞", en: "Neurons / proline-metabolizing cells" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "代谢与神经表型候选贡献基因", en: "Candidate metabolic and neurological contributor" }, phenotype: { zh: "高脯氨酸与神经行为候选表型", en: "Candidate hyperprolinemia and neurobehavioral phenotypes" },
    residual: { zh: "需联合测量脯氨酸与剩余酶活", en: "Proline and residual enzyme activity should be measured together" },
    caution: { zh: "仅凭脑分布不能判断代谢获益。", en: "Brain distribution alone cannot establish metabolic benefit." }, source: digeorgeRegion,
  },
  {
    id: "wbs-gtf2i", name: { zh: "Williams-Beuren 综合征", en: "Williams-Beuren syndrome" }, gene: "GTF2I", locus: "7q11.23", organ: "CNS", cnsProfile: "broad_neuronal",
    targetCell: { zh: "广泛神经元", en: "Broad neuronal populations" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "WBS 神经行为核心候选基因", en: "Core neurobehavioral candidate in WBS" }, phenotype: { zh: "社交、认知与神经行为候选表型", en: "Candidate social, cognitive and behavioral phenotypes" },
    residual: { zh: "另一条染色体通常保留正常转录本", en: "The homologous chromosome generally retains a normal transcript" }, caution: { zh: "复杂行为表型不能由单一基因递送分数解释。", en: "Complex behavioral phenotypes cannot be explained by one delivery score." }, source: williamsReview,
  },
  {
    id: "wbs-gtf2ird1", name: { zh: "Williams-Beuren 综合征", en: "Williams-Beuren syndrome" }, gene: "GTF2IRD1", locus: "7q11.23", organ: "CNS", cnsProfile: "neural_progenitor",
    targetCell: { zh: "神经发育相关细胞", en: "Neurodevelopmental cells" }, mechanism: "contiguous deletion", hiEvidence: { zh: "WBS 候选贡献基因", en: "Candidate contributor in WBS" },
    phenotype: { zh: "颅面与神经行为候选表型", en: "Candidate craniofacial and neurobehavioral phenotypes" }, residual: { zh: "需确认主要脑内 isoform", en: "The dominant brain isoform must be confirmed" },
    caution: { zh: "治疗优先级低于证据更明确的功能终点。", en: "Priority is lower than targets with clearer functional endpoints." }, source: williamsReview,
  },
  {
    id: "wbs-baz1b", name: { zh: "Williams-Beuren 综合征", en: "Williams-Beuren syndrome" }, gene: "BAZ1B", locus: "7q11.23", organ: "CNS", cnsProfile: "neural_progenitor",
    targetCell: { zh: "神经嵴 / 神经祖细胞", en: "Neural crest / neural progenitors" }, mechanism: "contiguous deletion", hiEvidence: { zh: "神经嵴与颅面候选贡献基因", en: "Candidate neural-crest and craniofacial contributor" },
    phenotype: { zh: "颅面与神经发育候选表型", en: "Candidate craniofacial and neurodevelopmental phenotypes" }, residual: { zh: "正常等位基因通常保留", en: "A normal allele generally remains" },
    caution: { zh: "关键窗口偏发育期，成人 CNS 暴露并不代表可逆。", en: "The key window is developmental; adult CNS exposure does not imply reversibility." }, source: williamsReview,
  },
  {
    id: "wbs-limk1", name: { zh: "Williams-Beuren 综合征", en: "Williams-Beuren syndrome" }, gene: "LIMK1", locus: "7q11.23", organ: "CNS", cnsProfile: "synaptic_neuron",
    targetCell: { zh: "突触神经元", en: "Synaptic neurons" }, mechanism: "contiguous deletion", hiEvidence: { zh: "视觉空间表型候选基因", en: "Candidate visuospatial-phenotype gene" },
    phenotype: { zh: "视觉空间认知候选表型", en: "Candidate visuospatial-cognition phenotype" }, residual: { zh: "需验证脑区特异的正常转录本", en: "Region-specific normal transcript requires verification" },
    caution: { zh: "目前关联不足以支持单独剂量推荐。", en: "Current association is insufficient for a stand-alone dose recommendation." }, source: williamsReview,
  },
  {
    id: "wagr-pax6", name: { zh: "WAGR 综合征", en: "WAGR syndrome" }, gene: "PAX6", locus: "11p13", organ: "Eye",
    targetCell: { zh: "角膜缘 / 视网膜相关细胞", en: "Limbal / retinal-associated cells" }, mechanism: "contiguous deletion", hiEvidence: sufficient,
    phenotype: { zh: "无虹膜与眼表表型", en: "Aniridia and ocular-surface phenotypes" }, residual: { zh: "另一等位基因通常保留 PAX6 转录本", en: "The other allele generally retains a PAX6 transcript" },
    caution: { zh: "需眼部局部 ODE；当前人体热图不应用于眼内剂量。", en: "A local ocular ODE is required; the current human heatmap should not guide intraocular dose." }, source: wagrReview,
  },
  {
    id: "wagr-wt1", name: { zh: "WAGR 综合征", en: "WAGR syndrome" }, gene: "WT1", locus: "11p13", organ: "Kidney",
    targetCell: { zh: "肾小球足细胞 / 肾脏发育细胞", en: "Glomerular podocytes / renal developmental cells" }, mechanism: "contiguous deletion", hiEvidence: sufficient,
    phenotype: { zh: "Wilms 瘤易感与泌尿生殖表型", en: "Wilms-tumor susceptibility and genitourinary phenotypes" }, residual: { zh: "必须区分肿瘤抑制与剂量恢复目标", en: "Tumor suppression must be separated from dosage-restoration goals" },
    caution: { zh: "肿瘤易感综合征不能依据递送效率自动推荐增强表达。", en: "A tumor-predisposition syndrome must not receive an automatic expression-enhancement recommendation." }, source: wagrReview,
  },
  {
    id: "mds-pafah1b1", name: { zh: "Miller-Dieker 综合征", en: "Miller-Dieker syndrome" }, gene: "PAFAH1B1", locus: "17p13.3", organ: "CNS", cnsProfile: "neural_progenitor",
    targetCell: { zh: "胎儿神经迁移祖细胞", en: "Fetal neuronal-migration progenitors" }, mechanism: "contiguous deletion", hiEvidence: sufficient,
    phenotype: { zh: "无脑回与严重癫痫表型", en: "Lissencephaly and severe epilepsy phenotypes" }, residual: { zh: "正常等位基因通常保留", en: "A normal allele generally remains" },
    caution: { zh: "主要病理在胎儿期形成，出生后递送不应被标成可逆治疗。", en: "Core pathology forms prenatally; postnatal delivery should not be labelled restorative." }, source: millerDiekerReview,
  },
  {
    id: "mds-ywhae", name: { zh: "Miller-Dieker 综合征", en: "Miller-Dieker syndrome" }, gene: "YWHAE", locus: "17p13.3", organ: "CNS", cnsProfile: "broad_neuronal",
    targetCell: { zh: "神经祖细胞 / 广泛神经元", en: "Neural progenitors / broad neurons" }, mechanism: "contiguous deletion",
    hiEvidence: { zh: "17p13.3 神经表型候选修饰基因", en: "Candidate neurological modifier in 17p13.3" }, phenotype: { zh: "脑发育严重度候选修饰表型", en: "Candidate modifier of brain-development severity" },
    residual: { zh: "需证明患者来源细胞的剂量-表型关系", en: "A dosage-phenotype relationship must be shown in patient-derived cells" }, caution: { zh: "不可替代 PAFAH1B1 的核心病因解释。", en: "It does not replace PAFAH1B1 as the core causal explanation." }, source: millerDiekerReview,
  },
  {
    id: "1p36-rere", name: { zh: "1p36 缺失综合征", en: "1p36 deletion syndrome" }, gene: "RERE", locus: "1p36.23", organ: "CNS", cnsProfile: "neural_progenitor",
    targetCell: { zh: "神经祖细胞 / 广泛神经元", en: "Neural progenitors / broad neurons" }, mechanism: "contiguous deletion", hiEvidence: sufficient,
    phenotype: { zh: "神经发育、眼与心脏表型", en: "Neurodevelopmental, ocular and cardiac phenotypes" }, residual: { zh: "正常等位基因通常保留", en: "A normal allele generally remains" }, caution: { zh: "多器官发育表型需要分器官终点。", en: "Multiorgan developmental phenotypes require organ-specific endpoints." }, source: oneP36Review,
  },
  {
    id: "1p36-prdm16", name: { zh: "1p36 缺失综合征", en: "1p36 deletion syndrome" }, gene: "PRDM16", locus: "1p36.32", organ: "Heart",
    targetCell: { zh: "心肌细胞", en: "Cardiomyocytes" }, mechanism: "contiguous deletion", hiEvidence: { zh: "心肌病候选贡献基因", en: "Candidate cardiomyopathy contributor" },
    phenotype: { zh: "扩张型 / 非致密化心肌病候选表型", en: "Candidate dilated / noncompaction cardiomyopathy phenotypes" }, residual: { zh: "需患者特异的心肌表达与功能验证", en: "Patient-specific myocardial expression and function require validation" },
    caution: { zh: "不能由区域缺失自动推定为主要病因。", en: "Regional deletion alone does not establish it as the primary cause." }, source: oneP36Review,
  },
  {
    id: "1p36-kcnab2", name: { zh: "1p36 缺失综合征", en: "1p36 deletion syndrome" }, gene: "KCNAB2", locus: "1p36.32", organ: "CNS", cnsProfile: "cortical_inhibitory",
    targetCell: { zh: "皮层神经元", en: "Cortical neurons" }, mechanism: "contiguous deletion", hiEvidence: { zh: "癫痫候选贡献基因", en: "Candidate seizure contributor" },
    phenotype: { zh: "癫痫与神经兴奋性候选表型", en: "Candidate seizure and excitability phenotypes" }, residual: { zh: "需先测量正常等位基因与通道功能", en: "Normal-allele expression and channel function should be measured first" },
    caution: { zh: "离子通道剂量恢复必须避免过度兴奋或抑制。", en: "Ion-channel restoration must avoid excessive excitation or inhibition." }, source: oneP36Review,
  },
  {
    id: "1p36-ski", name: { zh: "1p36 缺失综合征", en: "1p36 deletion syndrome" }, gene: "SKI", locus: "1p36.33", organ: "CNS", cnsProfile: "neural_progenitor",
    targetCell: { zh: "神经发育相关细胞", en: "Neurodevelopmental cells" }, mechanism: "contiguous deletion", hiEvidence: { zh: "颅面 / 神经发育候选贡献基因", en: "Candidate craniofacial / neurodevelopmental contributor" },
    phenotype: { zh: "颅面与神经发育候选表型", en: "Candidate craniofacial and neurodevelopmental phenotypes" }, residual: { zh: "需排除不同变异机制导致的反向剂量效应", en: "Opposing dosage effects from other variant mechanisms must be excluded" },
    caution: { zh: "SKI 通路剂量关系复杂，不自动进入组合推荐。", en: "SKI pathway dosage is complex and is not automatically included in combination recommendations." }, source: oneP36Review,
  },
];
