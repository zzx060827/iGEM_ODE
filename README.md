# SINEUP Delivery Atlas

[中文](#中文说明) | [English](#english-guide)

## 中文说明

SINEUP Delivery Atlas 是一个面向单倍剂量不足疾病的 AAV–SINEUP 个性化递送设计平台。平台将疾病与靶基因信息、AAV 衣壳特性、给药途径、PBPK 空间分布和 SINEUP 药效模型整合在同一界面中，用于比较候选方案并帮助设计后续实验。

> 本平台是研究级设计与假设生成工具，输出用于候选方案比较，不构成临床剂量或治疗建议。

## 快速开始

### 在线使用

打开平台主页后，可直接在顶部切换两个功能页面：

- **疾病设计空间**：根据疾病和缺失基因比较不同 AAV 衣壳的靶器官特异性、递送持续时间及预测蛋白恢复。
- **器官热图**：查看不同衣壳、给药途径和时间点下的全身或局部组织暴露。

### 本地运行

运行环境需要 Node.js 22.13.0 或更高版本。

```bash
cd sineup-delivery-atlas
npm install
npm run dev
```

根据终端提示，在浏览器中打开本地地址。若需检查生产构建，可运行：

```bash
npm run build
```

## 如何使用

### 1. 从疾病和基因开始

进入 **疾病设计空间**，在左侧疾病库中展开疾病，或搜索疾病名称、基因符号和染色体位点。选择目标基因后，页面会自动加载对应的靶器官、靶细胞、单倍剂量不足证据和 SINEUP 使用前提。

连续基因缺失可能包含多个候选基因。此时应分别检查每个基因是否仍存在可供 SINEUP 结合的正常转录本，以及提高蛋白表达是否具有可干预的治疗窗口。

### 2. 比较 AAV 衣壳

中央设计空间中的每个点代表一种衣壳方案：

- 横轴越靠右，表示目标器官相对脱靶器官的暴露越集中；
- 纵轴越靠上，表示模型预测的有效持续时间越长；
- 点的证据颜色表示其文献支持程度；
- 可使用“衣壳显示”按钮隐藏或恢复候选衣壳。

点击衣壳点后，右侧会显示给药途径、峰值递送、目标器官暴露份额、预测蛋白恢复、达峰时间、episome 持久性先验及主要证据来源。优先关注同时具有较高特异性、较长持续时间和较强证据支持的方案，而不是只比较单一指标。

### 3. 检查表达和组合方案

设计空间下方提供两类补充信息：

- **正常组织表达证据**：来自 GTEx 和 Human Protein Atlas，用于判断正常组织中是否存在目标转录本表达先验；这些数据不能替代患者或目标细胞数据。
- **方案组合筛选**：当单一衣壳和给药途径难以覆盖多个靶基因时，平台会显示探索性双方案。组合分数是模型排序指标，不是临床成功概率。

### 4. 查看器官热图

切换到 **器官热图** 后，可选择人体多区域模型或小鼠器官模型，并依次设置：

1. 给药途径；
2. AAV 衣壳；
3. 显示指标；
4. 时间点或早期/长期时间窗口。

点击人体图中的器官、区域或左侧排序列表，可查看该区域的峰值浓度、达峰时间、暴露份额、屏障后递送比例和完整时间曲线。播放按钮和时间滑块可用于观察递送信号随时间在循环系统及不同组织间的变化。

不同指标的含义如下：

| 指标 | 用途 |
|---|---|
| 血管/ISF 浓度 | 查看当前时间点的循环或组织间液暴露 |
| AUC | 比较整个观察期的累计暴露 |
| 暴露份额 | 判断总暴露在各器官间的分配 |
| Tmax | 比较不同区域的到达和达峰速度 |
| Episome | 查看进入细胞核后载体持久性的模型结果 |
| Protein | 查看 SINEUP–蛋白药效链的相对输出 |

## 如何解读结果

平台输出基于公开文献、GTEx/HPA/ClinGen 数据和机制模型。文献实测值、跨物种先验及工程假设的证据强度并不相同，因此结果适合用于：

- 筛选值得优先验证的衣壳和给药途径；
- 比较靶器官暴露与脱靶负担；
- 识别对结论影响最大的未知参数；
- 决定后续 qPCR/ddPCR、蛋白检测和空间取样的器官与时间点。

预测蛋白恢复率、65%阈值、人体衣壳倍率和联合方案评分仍需要目标基因、患者来源细胞及动物实验数据校准。使用结果时应同时查看证据来源和模型状态，避免将代理模型结果解释为临床疗效。

## 复现模型

```bash
python -m pip install -r requirements.txt
python model/export_parameter_registry.py
python model/export_delivery_design_space.py \
  --output sineup-delivery-atlas/public/data/model-results.json
python model/export_safety_screen.py
```

上述脚本依次导出参数证据表、前端所需的 ODE/PBPK 轨迹及研究级安全筛查结果。完整参数表见
[`model/data/model_parameter_registry.csv`](model/data/model_parameter_registry.csv)，其中每项参数均标注数值、单位、物种、证据类型、来源、选择理由、置信度和代码位置。

## 项目结构

| 路径 | 内容 |
|---|---|
| `model/` | 当前使用的 ODE/PBPK、参数证据、批量导出器和安全筛查 |
| `model/data/` | 衣壳文献数据与机器可读的完整参数注册表 |
| `sineup-delivery-atlas/` | React/vinext 交互式疾病设计空间与人体热图 |
| `docs/latex/`、`docs/pdf/` | 主报告、技术附录及其编译 PDF |
| `docs/wiki/` | Attribution、Engineering、Contribution、Measurement 和 Safety 草稿 |
| `docs/presentation/` | 简约中文汇报讲稿 |
| `results/` | 经过整理的模型图和表格输出 |
| `prototypes/` | 早期独立可视化原型 |
| `archive/legacy/` | 为追溯保留、但不再作为运行入口的旧脚本 |
| `releases/` | 可转移的演示压缩包 |

参数证据分为 `direct`（直接文献值）、`fitted`（由数据拟合）、`derived`（由已知量推导）、`scaled`（跨解剖或物种缩放）和 `assumed`（待实验校准的显式假设）。主报告、模型假设、局限性和湿实验闭环见
[`docs/pdf/aav_sineup_spatial_pk_project_report.pdf`](docs/pdf/aav_sineup_spatial_pk_project_report.pdf) 与
[`docs/pdf/ode_report.pdf`](docs/pdf/ode_report.pdf)。

## 安全边界

平台中的安全分数用于比较模型场景，不是临床安全剂量或不良事件概率。任何给药方案仍需结合组织学、临床化学、免疫原性、载体脱落和生物分布实验验证；详见 [`docs/wiki/safety.md`](docs/wiki/safety.md)。

---

## English Guide

SINEUP Delivery Atlas is a personalized AAV–SINEUP delivery design platform for haploinsufficiency disorders. It brings disease and target-gene information, AAV capsid properties, administration routes, PBPK spatial distribution, and SINEUP pharmacodynamic modeling into one interface to compare candidate strategies and guide follow-up experiments.

> This is a research tool for design exploration and hypothesis generation. Its outputs compare candidate strategies and are not clinical dosing or treatment recommendations.

## Quick Start

### Use the web interface

The top navigation provides two main views:

- **Disease Design Space** compares capsids by target-organ specificity, effective duration, and predicted protein restoration for a selected disease and gene.
- **Organ Heatmap** visualizes systemic or local tissue exposure across capsids, administration routes, metrics, and time points.

### Run locally

Node.js 22.13.0 or later is required.

```bash
cd sineup-delivery-atlas
npm install
npm run dev
```

Open the local address shown in the terminal. To verify a production build, run:

```bash
npm run build
```

## How to Use

### 1. Start with a disease and gene

Open **Disease Design Space**. Expand a record in the disease library, or search by disease name, gene symbol, or genomic locus. Selecting a target gene loads its target organ and cell type, haploinsufficiency evidence, and SINEUP prerequisites.

For a contiguous gene deletion, examine each candidate gene separately. Confirm that an intact transcript remains available for SINEUP binding and that increasing its protein level is biologically actionable within the treatment window.

### 2. Compare AAV capsids

Each point in the central design space represents one capsid strategy:

- A point farther to the right has more target-organ exposure relative to weighted off-target exposure.
- A point higher on the chart has a longer modeled effective duration.
- The evidence color indicates the level of literature support.
- The capsid visibility controls can hide candidates or restore all capsids.

Select a point to inspect its route, peak delivery, target-organ exposure share, predicted protein restoration, peak time, episome persistence prior, and primary evidence source. Prefer strategies that combine specificity, persistence, and strong evidence instead of ranking candidates by one metric alone.

### 3. Review expression evidence and regimen combinations

Two supporting sections appear below the design space:

- **Normal-tissue expression evidence** uses GTEx and the Human Protein Atlas to estimate whether the normal target transcript is present. These data do not replace patient-derived or target-cell measurements.
- **Regimen combination screening** shows an exploratory two-agent strategy when one capsid and route cannot adequately cover multiple genes. The combination score is a model-ranking metric, not a probability of clinical success.

### 4. Explore the organ heatmap

Open **Organ Heatmap**, choose the human multiregion or mouse organ model, and then select:

1. administration route;
2. AAV capsid;
3. display metric;
4. time point or early/long-term time window.

Select an organ or region on the anatomy map or ranking list to view its peak concentration, time to peak, exposure share, post-barrier delivery, and full time course. Use the play button and timeline slider to observe how modeled delivery moves through the circulation and tissues.

| Metric | Interpretation |
|---|---|
| Vascular/ISF concentration | Exposure in blood vessels or interstitial fluid at the selected time |
| AUC | Cumulative exposure over the observation period |
| Exposure share | Distribution of total modeled exposure across organs |
| Tmax | Arrival and peak timing across regions |
| Episome | Modeled vector persistence after nuclear delivery |
| Protein | Relative output of the SINEUP–protein pharmacodynamic chain |

## Interpreting the Results

The platform combines published evidence, GTEx/HPA/ClinGen data, and mechanistic modeling. Direct measurements, cross-species priors, and engineering assumptions have different levels of evidence. The outputs are therefore most useful for:

- prioritizing capsids and administration routes for validation;
- comparing target-organ exposure with off-target burden;
- identifying uncertain parameters that strongly affect the result;
- selecting organs and time points for qPCR/ddPCR, protein, and spatial measurements.

Predicted protein restoration, the 65% threshold, human capsid multipliers, and combination scores still require calibration with target-specific, patient-derived cell, and animal data. Always inspect the evidence source and model status, and do not interpret surrogate results as clinical efficacy.

## Reproduce the Model

```bash
python -m pip install -r requirements.txt
python model/export_parameter_registry.py
python model/export_delivery_design_space.py \
  --output sineup-delivery-atlas/public/data/model-results.json
python model/export_safety_screen.py
```

These commands export the parameter-evidence register, the ODE/PBPK trajectories consumed by the frontend, and the research-use safety screen. The complete registry is available at
[`model/data/model_parameter_registry.csv`](model/data/model_parameter_registry.csv); every row records the value, unit, species, evidence class, source, rationale, confidence, and code location.

## Project Structure

| Path | Purpose |
|---|---|
| `model/` | Active ODE/PBPK models, parameter evidence, exporters, and safety screening |
| `model/data/` | Capsid literature data and the machine-readable parameter registry |
| `sineup-delivery-atlas/` | React/vinext disease design space and anatomical heatmap |
| `docs/latex/`, `docs/pdf/` | Main report, technical appendix, and compiled PDFs |
| `docs/wiki/` | Attribution, Engineering, Contribution, Measurement, and Safety drafts |
| `docs/presentation/` | Concise Chinese presentation script |
| `results/` | Curated model figures and tabular outputs |
| `prototypes/` | Earlier standalone visualization prototypes |
| `archive/legacy/` | Superseded scripts retained for provenance, not for current runs |
| `releases/` | Portable demonstration archives |

Parameter evidence is classified as `direct`, `fitted`, `derived`, `scaled`, or `assumed`. Detailed assumptions, limitations, and the proposed wet-lab feedback loop are documented in
[`docs/pdf/aav_sineup_spatial_pk_project_report.pdf`](docs/pdf/aav_sineup_spatial_pk_project_report.pdf) and
[`docs/pdf/ode_report.pdf`](docs/pdf/ode_report.pdf).

## Safety Boundary

The safety score compares modeled scenarios; it is not a clinical safe dose or an adverse-event probability. Any regimen requires histopathology, clinical chemistry, immunogenicity, shedding, and biodistribution validation. See [`docs/wiki/safety.md`](docs/wiki/safety.md).
