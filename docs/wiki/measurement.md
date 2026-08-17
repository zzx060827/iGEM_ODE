# Measurement: how the model becomes testable

Measurement is the bridge between this model and a defensible biological
claim. The project does not yet possess a complete wet-lab calibration dataset,
so this page specifies the minimum measurement system needed to close the loop.

## Measurands must not be collapsed

| Biological quantity | Recommended assay | Suggested unit | Model state informed |
|---|---|---|---|
| Administered full-vector dose | ddPCR plus capsid/full-empty characterization | vg, capsid particles, full:empty ratio | `Dose_in` and input uncertainty |
| Circulating vector genome | matrix-qualified qPCR/ddPCR | vg/mL blood or plasma | blood clearance and organ input |
| Tissue vector genome | qPCR/ddPCR normalised to tissue mass and cell number | vg/g, vg/diploid genome | organ vascular/ISF plus internalised vector |
| Intact or labelled capsid | ELISA, imaging or validated label | capsid/mL, %ID/g | early capsid PK, not episome persistence |
| Cell-type uptake | co-localised imaging or sorted-cell ddPCR | fraction positive, vg/cell | receptor capacity and cell-access fraction |
| Nuclear/episomal vector | nuclear fraction, nuclease-resistant/circular genome assay | episome/cell | `Nss`, `Nds`, `Epi` |
| SINEUP RNA and target mRNA | RT-qPCR with construct-specific primers | copies/cell or relative expression | RNA production and degradation |
| Target protein | quantitative western, ELISA or functional assay | concentration or calibrated relative units | translation gain and protein turnover |
| Immune/safety biomarkers | NAb, anti-capsid T cells, cytokines, complement, ALT/AST, platelets, creatinine, troponin | assay-specific SI/clinical units | safety-priority layer |

Vector genome, capsid, episome, RNA and protein answer different questions and
should not share one “half-life”.

## Minimum time course

A practical first pass is 0.5 h, 2 h, 8 h, 24 h, 72 h and 7 d for vector/capsid,
with 7 d and 28 d (or a biologically appropriate later time) for episome, RNA
and protein. At least three biological replicates per condition are preferable;
technical replicates quantify assay precision but do not replace biological
replication.

## Controls and calibration

- untreated matrix and no-template controls;
- known vector-spike recovery and dilution linearity;
- a reference capsid/route shared across experiments;
- full versus empty capsid and genome integrity where feasible;
- tissue processing recovery controls;
- standard curves and lower limit of quantification;
- blinded image segmentation and pre-defined region boundaries;
- raw observations, replicate structure and uncertainty reported alongside the
  fitted mean.

## Proposed fitting and validation split

Use early blood and tissue observations to estimate transport/loss parameters;
use cell-type and subcellular data to estimate uptake/escape; use later
episome/RNA/protein observations to estimate expression and persistence. Hold
back at least one time point, dose or capsid for prediction rather than fitting.

The first quantitative success criterion is not a perfect curve. It is that a
model calibrated on one subset predicts the held-out condition within a
pre-declared error band and ranks the tested designs correctly with uncertainty.

## Measurement contribution that is feasible this season

A compact but strong deliverable is an assay-to-model protocol plus one pilot
dataset comparing AAV9 and one target-biased capsid in three organs. Publishing
the raw table, units, recovery controls, fitting script and failed measurements
would make the work reproducible and useful to future teams even before full
therapeutic validation.

---

# 中文版本：如何使模型成为可检验的科学工具

Measurement 是模型与可靠生物学结论之间的桥梁。项目目前还没有完整的湿实验校准数据，因此本页定义的是闭合模型—实验循环所需的最低测量体系，而不是宣称这些测量已经全部完成。

## 不同被测量必须分开

| 生物学量 | 推荐检测方法 | 建议单位 | 对应模型状态或参数 |
|---|---|---|---|
| 给入的完整载体剂量 | ddPCR，并结合衣壳总数和 full/empty 表征 | vg、capsid particles、full:empty ratio | `Dose_in` 与输入不确定性 |
| 循环载体基因组 | 经基质验证的 qPCR/ddPCR | vg/mL 血液或血浆 | 血液清除和器官输入 |
| 组织载体基因组 | qPCR/ddPCR，并按组织质量和细胞数归一化 | vg/g、vg/diploid genome | 器官血管/ISF 与内化载体的总和 |
| 完整或带标签衣壳 | ELISA、成像或经过验证的标记方法 | capsid/mL、%ID/g | 早期衣壳 PK，而不是 episome 持久性 |
| 细胞类型摄取 | 共定位成像或分选细胞 ddPCR | 阳性比例、vg/cell | 受体容量和细胞可及性 |
| 核内/episomal 载体 | 核分离、核酸酶耐受或环状基因组检测 | episome/cell | `Nss`、`Nds`、`Epi` |
| SINEUP RNA 与靶 mRNA | 使用构建体特异引物的 RT-qPCR | copies/cell 或相对表达 | RNA 生成与降解 |
| 靶蛋白 | 定量 western、ELISA 或功能检测 | 浓度或校准相对单位 | 翻译增益与蛋白周转 |
| 免疫/安全标志物 | NAb、抗衣壳 T 细胞、细胞因子、补体、ALT/AST、血小板、肌酐、肌钙蛋白 | 检测特异的 SI/临床单位 | 安全优先级层 |

载体基因组、衣壳、episome、RNA 和蛋白回答的是不同问题，不能共用同一个“半衰期”。组织 qPCR 也可能同时包含血管残留、ISF 和已内化基因组，因此需要灌流、分级或成像数据辅助解释。

## 最低时间序列

初步实验可在 0.5 h、2 h、8 h、24 h、72 h 和 7 d 测量载体/衣壳，并在 7 d、28 d 或疾病适合的更晚时点测量 episome、RNA 和蛋白。每组最好至少有 3 个生物学重复；技术重复只能说明检测精密度，不能替代生物学重复。

## 对照与校准

- 未处理基质和无模板对照；
- 已知载体加标回收率与稀释线性；
- 所有实验共享一个参考衣壳/给药途径；
- 在可行时测定 full/empty ratio 和基因组完整性；
- 组织处理回收率对照；
- 标准曲线与定量下限；
- 盲法图像分割和预先定义的解剖区域；
- 同时报告原始观测、重复结构、不确定性和拟合均值。

## 建议的拟合与验证划分

用早期血液和组织数据估计运输与损失参数；用细胞类型和亚细胞数据估计摄取与逃逸；用后期 episome/RNA/蛋白数据估计表达和持久性。至少保留一个时间点、剂量或衣壳作为预测集，不能全部用于拟合。

第一个量化成功标准不应是“曲线看起来完美”，而应事先声明：模型只用一部分数据校准后，能够在给定误差带内预测保留条件，并且在考虑不确定性后正确排序受测方案。

## 本赛季可完成的 Measurement Contribution

一个规模小但有说服力的交付物，是完成“检测—模型状态”实验方案，并生成一份 AAV9 与一个靶器官偏好衣壳在三个器官中的先导数据。公开原始表格、单位、回收率对照、拟合脚本和失败测量，即使还没有完成治疗验证，也能够使工作可复现并帮助未来团队。
