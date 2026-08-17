# Contribution: reusable infrastructure for future iGEM teams

Our contribution is not a claim that the current model predicts a clinical AAV
dose. It is a transparent, reusable workflow for turning heterogeneous AAV
evidence into auditable design comparisons.

## 1. Evidence-labelled AAV parameter registry

`model/export_parameter_registry.py` exports every active parameter to CSV,
JSON and a LaTeX table. Each row contains value, unit, species, source,
evidence type, rationale, confidence, code location and calibration priority.

This improves the usual “Value / Reference / Estimated” table in two ways:

- it distinguishes direct, fitted, derived, scaled and assumed values;
- it distinguishes a literature-supported mechanism from an uncalibrated
  numerical value.

Future teams can replace a prior with their own measurement without rewriting
the model narrative.

## 2. A route- and anatomy-aware AAV-to-frontend pipeline

The project provides an executable path from ODE equations to an interactive
atlas:

```text
PBPK + intracellular ODE
        -> capsid/route batch simulation
        -> JSON/CSV with provenance
        -> disease/gene design space and anatomical heat map
```

Unlike a hand-authored score chart, every displayed trajectory can be traced to
a model state and parameter set. Future teams can add a capsid, organ, route or
disease while keeping the same data contract.

## 3. Reusable checks and documentation

The repository includes:

- numerical mass-balance auditing with explicit loss sinks;
- cautious mouse-to-NHP-to-reference-human labels;
- absolute and relative anatomical colour scales;
- a layered-SVG authoring guide;
- Windows and macOS reproduction instructions;
- an exposure-priority safety screen linked to monitoring endpoints;
- archived failed/older models for provenance rather than silent deletion.

## How another team can build on it

1. Add source rows to `model/data/aav_capsid_tropism_literature.csv`.
2. Replace assumed transport or intracellular values in the registry with
   project measurements.
3. Add a target tissue or route in `human_spatial_pbpk.py`.
4. Re-run the exporter and inspect mass balance.
5. Use the generated JSON directly in the React atlas or another wiki frontend.
6. Report the new calibration and uncertainty instead of presenting a single
   deterministic point.

The most valuable future improvement would be a shared, assay-stratified AAV
benchmark dataset containing species, sex, capsid, promoter, dose, route,
sampling time, analyte and uncertainty. The current literature catalog is a
starting schema for that contribution.

## Relationship to earlier iGEM work

We learned from teams that documented parameter calculations and sensitivity
analysis, including GEMS-Taiwan and other 2022 modelling teams. Our extension
is to connect parameter provenance to a multiscale delivery model and a user
interface, while preserving the distinction between estimated and measured
quantities. We also follow the iGEM engineering criterion by documenting
troubleshooting and design changes, not only the final successful plot.

---

# 中文版本：可供未来 iGEM 团队复用的基础设施

本项目的 Contribution 不是声称当前模型已经能够给出临床 AAV 剂量，而是提供一条透明、可审计、可替换参数的工作流，将来源不同的 AAV 证据转化为可以比较的工程方案。

## 1. 带证据标签的 AAV 参数注册表

`model/export_parameter_registry.py` 将所有当前生效的参数导出为 CSV、JSON 和 LaTeX 表。每一行记录参数值、单位、物种、来源、证据类别、选择理由、置信度、代码位置和校准优先级，并提供相应中文字段。

相较于只写“Value / Reference / Estimated”的普通参数表，该注册表有两点改进：

- 区分直接文献值、拟合值、推导值、缩放值与假设值；
- 区分“某一机制有文献支持”和“当前采用的具体数值已经被实验校准”。

未来团队可以只替换一项先验及其证据标签，而不必重写整段模型叙述。

## 2. 识别给药途径与解剖结构的 AAV 到前端管线

仓库提供了从 ODE 方程到交互式图谱的可执行路径：

```text
PBPK + 胞内转运 ODE
        -> 按衣壳/给药途径批量求解
        -> 导出带来源的 JSON/CSV
        -> 疾病/基因设计空间与人体热图
```

与人工填写的递送评分图不同，网页中的每条轨迹都可以追溯到模型状态与参数集合。其他团队可以在保持数据契约不变的情况下增加衣壳、器官、给药途径或疾病。

## 3. 可复用的检查与文档

仓库还包括：

- 带显式损失池的数值质量守恒审计；
- 谨慎区分小鼠、NHP 与参考人体的证据标签；
- 人体热图的绝对颜色与相对颜色模式；
- 分层 SVG 制作指南；
- Windows 与 macOS 的复现说明；
- 将暴露排序映射到监测终点的安全筛查；
- 为追溯保留旧模型和失败版本，而不是静默删除。

## 其他团队如何继续使用

1. 在 `model/data/aav_capsid_tropism_literature.csv` 添加来源记录。
2. 用自己的实验测量替换注册表中的假设性运输或胞内参数。
3. 在 `human_spatial_pbpk.py` 增加目标组织或给药入口。
4. 重新运行导出器并检查质量守恒。
5. 在 React 图谱或自己的 wiki 前端直接读取生成的 JSON。
6. 报告新的校准过程与不确定性，不只展示单个确定性点。

未来最有价值的公共贡献将是按检测体系分层的 AAV 基准数据集，至少包含物种、性别、衣壳、启动子、剂量、途径、取样时间、分析物和不确定性。本项目的文献目录可以作为该数据集的初始 schema。

## 与既往 iGEM 工作的关系

我们参考了 GEMS-Taiwan 等团队对参数计算和敏感性分析的记录方法。我们的扩展是把参数来源、多尺度递送模型与用户界面连接起来，同时始终区分估计量和实测量；Engineering 页面也记录了排错、失败和设计改变，而不只是最后一张成功的图。
