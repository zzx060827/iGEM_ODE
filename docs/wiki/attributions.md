# Attributions: what we built, what we reused, and what remains to verify

This page is the working record behind the official iGEM Attributions Form. It
should be updated whenever a person, dataset, software package or external
service changes the project. The final form, not this draft alone, is the
competition deliverable.

## Work demonstrably completed in this repository

The current tracked implementation was completed primarily by **Zixuan Zhu**.
The work includes:

- formulation and implementation of the mouse-scale PBPK/ODE system;
- liver, proximal-tubule and BBB/CNS intracellular modules;
- the 70 kg reference-adult regional model and route-specific inputs;
- capsid batch simulation, SINEUP-PD linkage and JSON/CSV export;
- the disease/gene library, 2D design space and anatomical React heat map;
- parameter provenance, mass-balance checks, documentation and local demo;
- iterative visualization changes made in response to internal feedback.

This statement describes repository authorship. It does not imply that the
underlying biological facts, datasets, software libraries or anatomical assets
were created by the team.

## Prior scientific work used

| External work | How it was used | What we changed or added |
|---|---|---|
| Liu et al. (2024) AAV whole-body PBPK | Mechanistic precedent for organ distribution, receptor uptake, intracellular processing and transgene output | Added kidney dual-entry, BBB/CNS, multiple administration routes, SINEUP-PD and a disease-facing frontend |
| Wang et al. (2024) radiolabelled AAV9 mouse data | Fitted early apparent organ capsid half-lives | Kept raw time points, fit windows, log-fit diagnostics and the capsid-versus-episome caveat |
| Ballon et al. (2020) NHP PET | NHP-informed early AAV9 organ priors and CSF/systemic distribution context | Used only as a labelled reference-human projection; kidney/lung gaps remain explicit |
| Zincarelli et al. (2008), Walkey et al. (2025), Yang et al. (2025), Abele et al. (2025) | Head-to-head capsid and route evidence | Built a machine-readable catalog and cautious relative priors rather than pooling incompatible assays |
| GTEx and Human Protein Atlas | Healthy-tissue gene-expression priors | Aggregated tissues into model organs and exposed provenance/limitations in the interface |
| ClinGen | Haploinsufficiency and disease-gene evidence links | Organised disease entries as expandable disease-to-gene records |
| Reactome male-body SVG | Low-opacity anatomical reference | Overlaid independently calculated model regions; source SVG geometry is unchanged |
| DBCLS human anatomy SVG | Earlier anatomical reference retained for comparison | No source-geometry change |
| [Heidelberg 2025 SPARC and PHOENICS Builder](https://2025.igem.wiki/heidelberg/model) | Design precedent for turning a model-generated behaviour space into an interactive candidate-selection tool: users specify desired properties and the builder retrieves nearby precomputed circuit architectures | Adapted the general model-to-database-to-interface workflow to a different biological question. Our implementation uses prior-conditioned PBPK/ODE outputs for AAV capsids, routes, organs, CNS depth and SINEUP persistence; it does not reuse Heidelberg's code, parameters, circuit database or biological predictions |

Full citations are in the scientific report and source URLs are stored beside
the parameter/data records.

The Heidelberg attribution refers specifically to an interface and workflow
precedent. Their SPARC system models ligand sensing, receptor dimerisation and
phosphorylation circuits, whereas our atlas models AAV biodistribution,
intracellular trafficking and SINEUP-linked expression. The common idea is to
precompute a mechanistic design space and make it explorable; the equations,
data and intended decisions are distinct.

## Software, services and tools

| Tool | Role | Attribution note |
|---|---|---|
| Python, NumPy, SciPy, Matplotlib | ODE solution, integration, fitting and plots | Open-source scientific software; versions are constrained in `requirements.txt` |
| TypeScript, React 19, vinext/Vite and Next-compatible app APIs | Interactive frontend | The application is React-based; vinext/Vite provides the build and runtime layer, while the source retains a Next-compatible `app/` structure. Dependencies are listed in `sineup-delivery-atlas/package.json` |
| Codex / OpenAI tools | Code review, implementation assistance, literature-search assistance, document editing and local testing | AI assistance was supervised by the student author; references and numerical claims require human verification |
| Git and GitHub | Version control and public distribution | Repository history is the audit trail for code authorship |
| LaTeX / Tectonic, TeX Live and Overleaf | Scientific report typesetting and collaborative editing | TeX source and generated PDF are both versioned; Overleaf use should also be declared in the official form |

No AlphaFold server or Figma artifact is evidenced in this repository at the
time of writing. Add them to the official form if they were used elsewhere in
the team project.

## Visual asset licences

- Reactome `Male body with organs`, stable identifier `R-ICO-013956`, curated
  by Marija Orlic-Milacic and designed by Cristoffer Sevilla, CC BY 4.0.
- DBCLS `202403 human anatomy organs.svg`, CC BY 4.0.

The full URLs and change descriptions are in
`sineup-delivery-atlas/public/ASSET_ATTRIBUTION.md`.

## Items the full team must verify before submission

The repository cannot determine these contributions. Replace each placeholder
with names, dates and a concrete description in the official Attributions Form.

| Area | Person(s) | Exact contribution | Evidence/status |
|---|---|---|---|
| Project conception and supervision | **TEAM TO VERIFY** | Who framed the AAV-SINEUP therapeutic question and approved scope? | pending |
| Wet-lab design and experiments | **TEAM TO VERIFY** | Constructs, protocols, measurements, analysis and negative results | pending |
| RNA-binder model v2/v3 | **TEAM TO VERIFY** | Architecture, training data, code, compute and interpretation | outside this repository |
| Advisor/expert feedback | **TEAM TO VERIFY** | What decision changed after each consultation? | pending meeting records |
| Wiki integration and visual design | **TEAM TO VERIFY** | Layout, illustrations, copy editing and deployment | pending |
| Institutional facilities and funding | **TEAM TO VERIFY** | Laboratory, computing, reagents, grants and sponsorship | pending |

## Engineering contributors currently recorded by the team

The following entries mirror the team-supplied Chinese record. Names, dates and
supporting evidence should still be checked before transfer to the official
Attributions Form.

| Area | Person(s) | Exact contribution | Evidence/status |
|---|---|---|---|
| Literature review | Ziheng Wei; Zixuan Zhu | Surveyed scientific literature relevant to AAV delivery, pharmacokinetics and model design | completed; bibliography and parameter records available |
| ODE modeling concept and design | Ziheng Wei; Zixuan Zhu | Translated the wet-lab AAV-SINEUP therapeutic question into a computational delivery model. In view of [iGEM animal-use approval requirements](https://responsibility.igem.org/safety-policies/animal-use), the project timeline and available resources, the team used appropriately cited public animal data and in-silico simulation instead of conducting new animal experiments | substantially complete |
| Two-dimensional visualization | Ziheng Wei; Zixuan Zhu | Jointly conceived the model-facing visualization and developed an interactive React-based interface for disease, capsid, route and organ-level result exploration | substantially complete; implementation available in the repository |
| Wet-lab design and validation | **TEAM TO VERIFY** | Constructs, protocols, measurements, analysis and negative results | pending |
| Advisor/expert discussion | Dr. Steven A. Benner | Reviewed the project direction and encouraged broader interdisciplinary reasoning | meeting record and resulting design decisions to be documented |
| Institutional facilities and computing support | Prof. Liqin Zhang's laboratory | Provided laboratory space and computational resources | completed; exact facilities and funding to be confirmed |

## Statement of intellectual honesty

Model outputs are the team's calculations, but most parameters are not the
team's measurements. Every output therefore carries a model/evidence label.
Reference-human results are projections rather than clinical predictions, and
the safety screen prioritises measurements rather than declaring a dose safe.

---

# 我们完成了什么、复用了什么、还有什么需要确认

本页是 iGEM 官方 Attributions Form 的工作底稿。任何成员分工、数据集、软件包或外部服务发生变化时，都会继续同步更新。

## 本仓库中可以直接核验的工作

当前仓库内的工程实现主要由朱子轩（Zixuan Zhu）完成，包括：

- 小鼠尺度 PBPK/ODE 系统的建模与代码实现；
- 肝脏、肾近端小管和 BBB/CNS 胞内模块；
- 70 kg 参考成人区域模型与给药途径入口；
- 多衣壳批量模拟、SINEUP-PD 连接与 JSON/CSV 导出；
- 疾病/基因库、二维设计空间和 React 人体热图；
- 参数来源、质量守恒检查、文档和本地演示；
- 根据组内反馈完成的多轮可视化迭代。

这段说明描述的是仓库作者身份，并不表示模型所依赖的生物学知识、公开数据集、开源软件或人体解剖素材由本团队原创。

## 使用的既有科学工作

| 外部工作                                                                                    | 在项目中的用途                                                                                    | 我们新增或改变的内容                                                                                                                                                               |
| --------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Liu 等（2024），AAV8/AAV9 whole-body PBPK**                                               | 提供 AAV 全身 PBPK 的机制先例，包括器官分布、组织摄取、胞内处理和转基因输出；其 AAV8/AAV9 三周组织数据和全组织/血液比被记录为后续联合校准目标         | 在其基础上扩展肾脏双入口、BBB/CNS 空间结构、多给药途径、AAV 衣壳差异、SINEUP-PD 和面向疾病的交互前端。文献中的 whole-tissue qPCR ratio **没有直接替换模型 `Kp`**，因为全组织载体 DNA 与模型中的 ISF 分配系数并不是同一个物理量                         |
| **Wang 等（2024），放射标记 AAV9 小鼠数据**                                                         | 用于估计早期器官 AAV9 衣壳信号的表观半衰期，并为小鼠器官分布提供时间序列约束                                                  | 保存原始时间点、拟合窗口和对数拟合诊断；明确区分放射标记衣壳信号、vector genome 与长期 episome，避免把早期 capsid clearance 解释为转基因持久性                                                                              |
| **Ballon 等（2020），NHP 全身 AAV PET**                                                       | 提供非人灵长类早期 AAV9 全身分布和器官 PK 背景，并为参考人体模型的早期分布提供跨物种先验                                          | 仅作为明确标记的 reference-human projection 使用，没有把 NHP 数据描述成人体直接测量；肾、肺等缺乏充分定量数据的区域继续保留不确定性标签                                                                                     |
| **Zincarelli 等（2008），AAV1–AAV9 同条件比较**                                                  | 提供经典的多衣壳头对头组织嗜性证据；用于支持 AAV9 广泛分布以及 AAV4、AAV6 等衣壳的相对器官偏好                                    | 将其整理为机器可读的 capsid–organ 文献记录，并转换为**相对先验而非绝对转导效率**；不跨检测方法直接合并数值，也不把单一动物实验解释为人体 tropism                                                                                    |
| **Walkey 等（2025），10 种 AAV 衣壳 × 22 种组织比较**                                               | 提供更大规模的多衣壳、多组织数据，并同时包含性别差异、vector DNA 与功能性 transduction 信息                                 | 将 sex、tissue、measurement type 和 functional/DNA readout 分开记录，使前端与参数注册表能够区分“载体到达组织”和“产生功能表达”；没有将 DNA abundance 等同于功能转导                                                     |
| **Yang 等（2025），21 种天然/工程 AAV 衣壳跨小鼠与猕猴比较**                                               | 为 AAV9、AAVrh10、LK03、PHP.eB、CAP-B10 等候选衣壳提供跨品系、跨物种的相对性能证据                                   | 将 species、mouse strain、capsid engineering status 和 assay context 纳入文献目录，用于约束相对衣壳先验；跨物种差异被显式保留，而不是将小鼠或猕猴结果直接作为人体参数                                                        |
| **Abele 等（2025），34 种 AAV 衣壳及 IV/IP 给药比较**                                               | 提供大规模衣壳筛选证据，并证明同一衣壳的器官偏好会随给药途径改变                                                           | 在模型和数据库中把 **capsid × route** 作为联合条件处理，而不是为每个衣壳赋予一个固定的全局 tropism score；为未来不同给药途径重新校准相对先验提供结构                                                                              |
| **Bartlett 等（2000），AAV2 细胞内运输研究**                                                       | 提供 AAV 进入细胞后快速内化和向核周区域运输的实验时间尺度背景，包括 `<10 min` 内化及约 `2 h` 核周定位的观察                          | 将这些数据加入胞内运输参数的 evidence-comparison 字段，用于检查模型时间尺度是否合理；由于衣壳为 AAV2、细胞系和实验条件与体内 AAV9 不同，**没有直接覆盖 AAV9 的器官特异性 internalization 参数**                                            |
| **人体 AAV 肝活检长期 episome 研究**                                                             | 提供人体肝组织中 AAV vector genome 长期存在及转录活性的直接证据；给药后约 2.6–4.1 年仍可检测到具有转录能力的 episomal genomes      | 将 reference-human liver episome persistence prior 从 `120 d` 更新为 `1095 d`。该值被明确标记为保守的长期持久性先验，**不是由该研究精确拟合得到的 episome half-life**                                          |
| **人体骨骼肌 AAV 长期 persistence 研究**                                                         | 提供人体骨骼肌中环状、可转录 AAV genomes 可持续至少约 4 年的证据                                                   | 将 reference-human skeletal-muscle episome persistence prior 从 `365 d` 更新为 `1460 d`，同时保留其为 evidence-informed prior，而不是精确动力学半衰期                                            |
| **FDA ZOLGENSMA prescribing information**                                               | 提供已上市 systemic AAV9 的临床剂量背景（`1.1×10^14 vg/kg`）、组织 biodistribution 信息以及肝脏、血小板/TMA、心脏等安全监测背景 | 仅作为**临床情景与安全比较背景**，没有把上市产品剂量或肝脏 vector DNA 数值直接用于校准 70 kg reference-adult PBPK 模型；模型输出仍被标记为研究级 projection                                                                |
| **Human AAV9 post-mortem / biodistribution data**                                       | 为人体接受系统性 AAV9 后不同组织中的 vector distribution 提供直接人体背景证据                                       | 用于检查参考人体模型预测的器官排序与数量级是否明显违背人体观测，而不将尸检终点的 tissue vector DNA 直接等同于动态 ISF concentration、`Kp` 或 episome half-life；正式 Attribution 中应补上实际采用论文的题目和链接                            |
| **AAV9-miniSINEUP-GDNF 小鼠研究**                                                           | 提供与本项目载荷机制直接相关的 AAV9–miniSINEUP 先例，包括体内 SINEUP 递送以及约两倍内源蛋白提升和长期表达背景                        | 将其作为 SINEUP-PD 幅度和机制合理性的外部依据，但没有把 GDNF 的具体效应量直接视为其他 binding domain、其他靶基因或人体中的固定 SINEUP efficacy                                                                          |
| **[Heidelberg 2025：SPARC 与 PHOENICS Builder](https://2025.igem.wiki/heidelberg/model)** | 作为“模型生成设计空间，再由交互式软件筛选候选方案”的 iGEM 设计先例；使用者给定目标性质后，从预计算模型结果中寻找适合的候选设计                        | 将通用的 **model → database → interactive design tool** 思路应用到不同生物学问题：本项目展示由文献先验约束、经 PBPK/ODE 求解的 AAV 衣壳、给药途径、器官/CNS 空间分布及 SINEUP 持久性结果；**未复用 Heidelberg 的代码、参数、回路数据库或生物学预测** |
| **GTEx 与 Human Protein Atlas**                                                          | 提供健康人体不同组织的 mRNA / protein expression 背景，用于判断目标基因是否存在可供 SINEUP 利用的正常转录本，以及确定潜在治疗靶器官        | 将数据库组织映射到模型器官，在疾病设计空间中展示表达证据、来源与限制；这些健康组织数据没有被描述为患者特异性表达数据                                                                                                               |
| **ClinGen**                                                                             | 提供 haploinsufficiency、gene–disease validity 等疾病—基因证据                                       | 将其整理成可展开的 disease–gene records，并与靶器官、转录本表达和 SINEUP 使用前提连接；没有重新生成或修改 ClinGen 的原始临床遗传学判断                                                                                   |
| **Reactome `Male body with organs` SVG**                                                | 作为人体解剖可视化的低透明度底图                                                                           | 在来源 SVG 上叠加由模型独立计算得到的空间区域、颜色和交互数据；没有改变来源 SVG 的器官几何，并保留 CC BY 4.0 attribution                                                                                             |
| **DBCLS human anatomy SVG**                                                             | 作为早期人体解剖可视化和布局参考                                                                           | 保留于早期原型用于追溯和比较；未修改来源几何，当前主要模型结果并不来自该 SVG                                                                                                                                 |



完整引文位于科学报告中，关键 URL 也保存在参数与数据记录旁。

这里对海德堡工作的引用仅表示界面与工作流层面的借鉴。其 SPARC 系统描述配体识别、受体二聚化和磷酸化回路；我们的模型描述 AAV 体内分布、胞内转运和 SINEUP 相关表达。两者共有的思想是先计算机制性设计空间，再让使用者交互探索，但方程、数据和所支持的设计决策均不相同。

## 软件、服务与工具

| 工具 | 作用 | Attribution 说明 |
|---|---|---|
| Python、NumPy、SciPy、Matplotlib | ODE 求解、积分、拟合和绘图 | 开源科学软件，版本范围写入 `requirements.txt` |
| TypeScript、React 19、vinext/Vite 与 Next 兼容的 app API | 交互式前端 | 网页确实基于 React；vinext/Vite 提供构建与运行层，源码保留 Next 兼容的 `app/` 目录结构。开源依赖列于 `sineup-delivery-atlas/package.json` |
| Codex / OpenAI 工具 | 代码审阅与实现辅助、文献检索辅助 | AI 输出由作者监督；参考文献和数值结论由人核验 |
| Git 与 GitHub | 版本管理和公开发布 | 仓库历史构成代码作者与改动的审计记录 |
| LaTeX / Tectonic 和 Overleaf | 科学报告排版 | TeX 源文件与生成的 PDF 均纳入版本控制 |


## 可视化素材

- Reactome `Male body with organs`，稳定编号 `R-ICO-013956`，策展人 Marija Orlic-Milacic，设计者 Cristoffer Sevilla，许可 CC BY 4.0。
- DBCLS `202403 human anatomy organs.svg`，许可 CC BY 4.0。

完整 URL 与修改说明位于 `sineup-delivery-atlas/public/ASSET_ATTRIBUTION.md`。

## 工程

| 领域 | 人员 | 具体贡献 | 当前状态 |
|---|---|---|---|
| 文献调研 | 魏子恒、朱子轩 | 调研 AAV 递送、药代动力学和模型设计相关文献 | 已完成 |
| ODE 建模思路与设计 | 魏子恒、朱子轩 | 根据湿实验提出的 AAV-SINEUP 治疗问题，考虑到 [iGEM 动物实验的事前审批要求](https://responsibility.igem.org/safety-policies/animal-use)、项目周期和可用资源，使用规范引用的公开动物数据与计算模拟替代新增动物实验，并据此构思递送模型 | 较完整 |
| 二维可视化 | 魏子恒、朱子轩 | 共同构思面向模型结果的可视化，并使用 React 构建可交互前端，用于探索疾病、衣壳、给药途径和器官结果 | 较完整；实现已纳入仓库 |
| 湿实验设计与验证 | **待确认** | 待补充 | 待完成 |
| 与导师/专家交流 | Dr. Steven A. Benner | 向教授汇报项目方向并获得建议，促进团队从更广的跨学科角度审视问题 | 待整理会议记录及由交流引起的具体设计变化 |
| 机构设施与计算支持 | 张力勤老师课题组 | 提供实验室与计算资源 | 已完成|

## 学术诚信声明

模型输出是本团队的计算结果，但多数参数不是本团队的实测数据，因此每项输出都附有模型和证据标签。参考人体结果是外推，不是临床预测；安全筛查用于确定需要优先测量的风险器官，不用于宣称某剂量安全。
