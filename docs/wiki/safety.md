# AAV safety research: current evidence and model scope

This is an early hazard-identification framework. It does not demonstrate that
the proposed AAV-SINEUP product is safe, and it must not be used to select a
clinical dose.

## Why safety is product-specific

Risk depends on total capsid and vector-genome dose, full/empty ratio, serotype,
promoter, payload, impurities, route, infusion procedure, age, disease,
pre-existing immunity and concomitant infection or treatment. A dose used by
one licensed product is therefore context, not a universal threshold.

The current 70 kg simulation uses `4.0e13 vg/kg`. For comparison only, the US
ZOLGENSMA label specifies `1.1e14 vg/kg` as a single 60-minute IV infusion in a
defined paediatric SMA population, with corticosteroid prophylaxis and intensive
monitoring. The difference does not prove a safety margin for our construct.

### 2026 route-matched re-screen

The active safety exporter now re-solves the same PBPK equations at selected
reference doses instead of comparing dose numbers alone. For IV AAV9, the
modeled organ ISF-AUC is approximately 2.75-fold below the ZOLGENSMA dose
context and 1.98-fold below the lowest explicitly reported mouse cardiac
histology signal-dose context (`7.9e13 vg/kg`). These are contextual exposure
margins, not a human NOAEL.

The same total `2.8e15 vg` was previously applied to every administration
route. For lumbar intrathecal AAV9 this is 23.3-fold the marketed ITVISMA fixed
dose of `1.2e14 vg`; the re-solved organ AUCs are approximately 23.3-23.7-fold
the ITVISMA-dose model context. The current CSF dose is therefore not supported
for safety interpretation and must be reduced and re-solved.

The full calculation, organ evidence grades and primary sources are in
[`aav_safety_margin_assessment_2026.md`](../aav_safety_margin_assessment_2026.md).

## Clinically important AAV risk domains

| Risk | Evidence context | What should be monitored or modelled |
|---|---|---|
| Hepatotoxicity | A prominent risk after systemic AAV; serious liver injury and failure are included in product warnings | Liver capsid/transgene exposure, ALT, AST, bilirubin, albumin, PT/INR, anti-capsid T cells |
| Thrombocytopenia and TMA | Reported after systemic AAV9; complement activation may accompany microvascular injury | Platelets, hemoglobin, creatinine, urinalysis, LDH/haptoglobin and complement markers |
| Innate/adaptive immunity | Pre-existing NAb can reduce efficacy; capsid and transgene responses can alter safety and persistence | NAb titre, antibody/T-cell response, cytokines, complement, infection status |
| Cardiac signal | Troponin elevation is monitored for onasemnogene abeparvovec | Heart exposure, troponin and cardiac assessment |
| DRG/CNS neurotoxicity | DRG pathology has occurred in nonclinical/clinical CNS-directed AAV programmes; local brain MRI findings are route/procedure dependent | CSF/DRG exposure, neurological/sensory endpoints, MRI, CSF biomarkers and relevant-animal histopathology |
| Unwanted expression | Broad tropism or strong promoter can create on-target/off-tissue effects; excess restoration may also be harmful | Organ-specific promoter activity, target protein dose-response and off-target tissues |
| Integration/tumorigenicity | AAV is mainly episomal but integration risk is not zero | Long-term follow-up and product-specific integration assessment |
| CMC and procedure | Empty particles, residual host-cell material, endotoxin and delivery devices can change risk | Identity, purity, potency, full:empty ratio, endotoxin, sterility and device compatibility |

Primary regulatory context: the FDA's
[AAV clinical-development toxicity examples](https://www.fda.gov/media/167536/download),
[guidance for human gene therapy in neurodegenerative diseases](https://www.fda.gov/media/144886/download),
and the current US
[ZOLGENSMA prescribing information](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=68cd4f06-70e1-40d8-bedb-609ec0afa471).

## What the model now does

`model/export_safety_screen.py` compares each route/capsid organ AUC with IV
AAV9 in the same reference-adult model. It maps relative liver, kidney/spleen,
heart, CNS and systemic exposure to assay priorities. CSF routes automatically
trigger a route-specific CNS/DRG flag.

The screen is useful for questions such as “which organ should receive extra
toxicology measurements?” It does not contain a calibrated relationship from
AUC to ALT elevation, TMA, troponin or neuronal injury, so it cannot answer
“what dose is safe?”.

## Proposed next safety model

The next defensible layer is an exposure-response model fitted to a single
well-characterised vector programme:

1. retain capsid, vector genome, transgene and immune analytes separately;
2. use a saturable or sigmoid dose-exposure model where supported;
3. link liver AUC to longitudinal ALT/AST and immune markers;
4. link platelet/creatinine/complement trajectories to a TMA hazard endpoint;
5. treat DRG/CNS injury as a route- and species-conditioned probability;
6. propagate parameter uncertainty and report probability intervals;
7. validate on held-out dose groups or a second study.

Recent mouse dose-ranging work indicates non-proportional tissue exposure,
sigmoidal transgene expression and greater immune/hepatic signals at high dose.
That supports a nonlinear framework, but mouse AAV8 data cannot directly set a
human AAV9-SINEUP toxicity threshold.

## Risk-reduction decisions available now

- minimise total capsid needed for the target protein window;
- compare local/CSF routes with systemic exposure rather than assuming “local”
  means no peripheral distribution;
- use tissue/cell-selective regulatory elements and evaluate excess expression;
- screen pre-existing immunity and define route-appropriate monitoring;
- characterise full/empty ratio, potency and impurities before interpreting
  dose-response data;
- avoid claims about repeat administration until immune and product-specific
  evidence exists.

---

# 中文版本：AAV 安全性研究、当前证据与模型边界

本页是早期危害识别框架，不能证明拟议的 AAV-SINEUP 产品安全，也不能用于确定临床剂量。

## 为什么安全性必须针对具体产品判断

风险取决于总衣壳量与载体基因组剂量、full/empty ratio、血清型、启动子、载荷、杂质、给药途径、输注过程、年龄、疾病、既往免疫和合并感染/治疗。一个获批产品使用的剂量只能提供背景，不能成为所有 AAV 产品通用的安全阈值。

当前 70 kg 模拟使用 `4.0e13 vg/kg`。仅作比较，美国 ZOLGENSMA 标签规定在特定儿科 SMA 人群中按 `1.1e14 vg/kg` 单次 60 min 静脉输注，并配合糖皮质激素预处理和密切监测。两者载荷、产品工艺与人群不同，当前模型剂量较低并不能证明存在安全裕度。

### 2026 同途径重筛查

当前安全导出器已不再只比较名义剂量，而是用相同 PBPK 方程在参照剂量下重新求解。对于 IV AAV9，模型器官 ISF-AUC 相对 ZOLGENSMA 剂量背景约低 2.75 倍，相对标签中最低明确报告的小鼠心脏组织学信号剂量（`7.9e13 vg/kg`）约低 1.98 倍。这些是有边界的暴露背景裕度，不是人体 NOAEL。

此前模型把同一个 `2.8e15 vg` 总量用于所有给药途径。对于腰椎鞘内 AAV9，这相当于已上市 ITVISMA 固定剂量 `1.2e14 vg` 的 23.3 倍；参照剂量重求解后的各器官 AUC 仅为当前结果的约 4.2%，即当前结果约高 23.3–23.7 倍。因此当前 CSF 剂量不支持安全性解释，必须先降剂量并重新求解。

完整计算、器官证据等级和一手来源见 [`aav_safety_margin_assessment_2026.md`](../aav_safety_margin_assessment_2026.md)。

## 临床上重要的 AAV 风险领域

| 风险 | 证据背景 | 应监测或建模的内容 |
|---|---|---|
| 肝毒性 | 系统性 AAV 的突出风险，产品警示中包含严重肝损伤和肝衰竭 | 肝衣壳/转基因暴露、ALT、AST、胆红素、白蛋白、PT/INR、抗衣壳 T 细胞 |
| 血小板减少与 TMA | 系统性 AAV9 后有报告；补体激活可能伴随微血管损伤 | 血小板、血红蛋白、肌酐、尿检、LDH/结合珠蛋白和补体标志物 |
| 先天/适应性免疫 | 既往 NAb 可降低疗效；衣壳与转基因免疫反应可改变安全性和持久性 | NAb 滴度、抗体/T 细胞、细胞因子、补体、感染状态 |
| 心脏信号 | onasemnogene abeparvovec 需要监测肌钙蛋白升高 | 心脏暴露、肌钙蛋白与心脏评估 |
| DRG/CNS 神经毒性 | CNS 定向 AAV 项目中出现过 DRG 病理；局部脑部 MRI 改变还与途径和操作有关 | CSF/DRG 暴露、神经/感觉终点、MRI、CSF 标志物和相关动物组织病理 |
| 非预期表达 | 广泛嗜性或强启动子可导致同一靶标在错误组织表达，过度恢复也可能有害 | 器官特异启动子活性、靶蛋白剂量反应和脱靶组织 |
| 整合与肿瘤发生 | AAV 主要以 episome 存在，但整合风险并非零 | 长期随访与产品特异的整合评估 |
| CMC 与给药过程 | 空衣壳、宿主细胞残留、内毒素和输送装置都会改变风险 | 身份、纯度、效价、full:empty ratio、内毒素、无菌和装置相容性 |

主要监管背景包括 FDA 的 [AAV 临床开发毒性实例](https://www.fda.gov/media/167536/download)、[神经退行性疾病基因治疗指南](https://www.fda.gov/media/144886/download)和美国 [ZOLGENSMA 处方信息](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=68cd4f06-70e1-40d8-bedb-609ec0afa471)。

## 当前模型已经能做什么

`model/export_safety_screen.py` 将每种给药途径/衣壳的器官 AUC 与同一参考成人模型中的 IV AAV9 比较，并把肝、肾/脾、心脏、CNS 与全身暴露的相对变化映射到检测优先级。所有 CSF 给药途径都会自动触发 CNS/DRG 特异风险标志。

这个筛查可以回答“下一轮毒理实验应在哪个器官增加测量”，但当前没有从 AUC 到 ALT、TMA、肌钙蛋白或神经损伤的校准关系，因此不能回答“多大剂量安全”。

## 下一步可辩护的安全模型

应在单一、表征充分的载体项目内拟合暴露—反应层：

1. 分别保留衣壳、载体基因组、转基因产物和免疫分析物；
2. 有数据支持时使用饱和或 S 型剂量—暴露关系；
3. 将肝 AUC 与纵向 ALT/AST 和免疫标志物连接；
4. 将血小板、肌酐与补体轨迹连接到 TMA 风险终点；
5. 将 DRG/CNS 损伤建模为依赖给药途径和物种的概率；
6. 传播参数不确定性并报告概率区间；
7. 使用保留剂量组或第二项研究进行外部验证。

近期小鼠剂量研究支持组织暴露非比例、转基因表达呈 S 型，以及高剂量下免疫/肝脏信号更强，因此未来应采用非线性框架。但小鼠 AAV8 数据不能直接定义人体 AAV9-SINEUP 的毒性阈值。

## 当前即可采用的风险降低决策

- 在达到目标蛋白治疗窗口的前提下，尽量降低总衣壳量；
- 比较局部/CSF 给药与系统暴露，不能假设“局部给药”等于“没有外周分布”；
- 使用组织或细胞选择性调控元件，并检测过度表达；
- 筛查既往免疫并定义与途径相适应的监测方案；
- 在解释剂量反应前表征 full/empty ratio、效价和杂质；
- 在缺乏免疫学与产品特异证据时，不宣称可以重复给药。
