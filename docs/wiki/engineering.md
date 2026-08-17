# Engineering: the model as a Design-Build-Test-Learn system

The engineering problem was not simply to draw AAV concentration curves. We
needed a computational system that could explain where a vector goes, identify
the step limiting expression, expose uncertainty, and return a smaller set of
capsid-route experiments to the wet lab.

## Acceptance criteria

Before comparing biological scenarios, each build had to meet four tests:

1. **Conservation:** administered AAV equals all vector-containing states plus
   explicit cumulative sinks, within numerical tolerance.
2. **Non-negativity and solver success:** no meaningful negative state and no
   failed integration.
3. **Biological separation:** early capsid disappearance, vector genome,
   episome, RNA, protein and therapeutic duration must remain different states.
4. **Reproducibility:** one command regenerates the frontend JSON/CSV, and the
   page displays those exported results without inventing new scores in the
   browser.

Passing these criteria shows internal consistency, not clinical validity.

## Cycle 1: from arbitrary decay to traceable AAV9 kinetics

**Design.** The first implementation used global vascular and ISF half-lives to
make bell-shaped curves. This made the plots readable but did not explain where
capsid went.

**Build.** We added explicit blood loss, organ vascular/RES loss, ISF catabolic
loss and cumulative sink states. Organ half-lives were fitted to Wang et al.
mouse data; the ordinary AAV9 blood half-life was corrected to 5.0 h because
2.4 h referred to a modified AAV9-TC construct.

**Test.** The mass-balance residual fell to numerical precision, and organ
curves retained distinct peaks and declines. Log-fit diagnostics exposed that
muscle was poorly represented by a single exponential.

**Learn.** A measured apparent half-life does not identify endothelial removal
and ISF catabolism separately. We therefore retained the 35:65 rate split as an
explicit structural assumption and ranked it for calibration instead of
presenting it as literature fact.

## Cycle 2: from liver-only uptake to organ- and barrier-specific transduction

**Design.** Organ ISF exposure alone could not answer why two tissues with
similar vector genomes produce different expression.

**Build.** The liver chain became `bound -> EE -> LE -> CY -> Ncap -> Nss ->
Nds -> Epi -> mRNA -> protein`. Kidney received small filtration and
basolateral entry routes; CNS received BBB binding, endothelial trafficking,
transcytosis/recycling and neural-cell uptake.

**Test.** Scenario comparisons changed capsid transport separately from
promoter output. Increasing endosomal escape raised episome more strongly than
protein once the expression module approached saturation.

**Learn.** “More organ exposure” is not equivalent to “more therapeutic
protein”. This supports experiments that measure vector, episome, RNA and
protein rather than a single endpoint. During this audit we also found that
liver transcription was capped at a hard-coded value of 2.0; it now caps at the
scenario's own `k_tx`, so promoter presets behave as declared.

## Cycle 3: from one animal diagram to a route-resolved reference adult

**Design.** Scaling every mouse state by body weight would preserve the wrong
circulation, organ fractions and CSF geometry.

**Build.** We retained the same Q-PS-Kp and intracellular equation family but
introduced explicit right heart, lung, left heart, arterial/venous pools, portal
drainage, 24 regions, 301 states and IV, IT, IM, ICM, ICV and inhaled inputs.

**Test.** Regional blood volumes sum to approximately 5 L; flow fractions close
the circulation; 150 mL CSF with 500 mL/day turnover gives a 4.99 h equivalent
half-life; each route is independently checked for mass balance.

**Learn.** Human anatomy can be mechanistic while capsid parameters remain
uncertain. The interface therefore labels the result “reference-human
projection” and removes the mouse LY6A-dependent PHP.eB CNS gain rather than
calling it a human prediction.

## Cycle 4: from static plots to a disease-facing design tool

**Design.** A wet-lab user should begin with disease and target gene, not with a
list of state-variable names.

**Build.** The exporter solves all capsids, computes organ specificity and
SINEUP-PD duration, then writes JSON/CSV. The React application provides a
disease-to-gene library, CNS-depth profiles, route/capsid selection, a 2D design
space and an anatomical heat map.

**Test.** Frontend tests build the application and verify rendered HTML. The
heat map now supports absolute and within-capsid relative scales; the latter is
the default to reveal spatial differences without claiming cross-capsid dose
equivalence.

**Learn.** Relative colour is useful for pattern recognition but can hide
absolute decline. Both scales remain available and the legend states the
normalisation. We also fixed a hydration mismatch caused by locale-dependent
server/client timestamp formatting.

## Current design decision returned to the wet lab

The model recommends a staged comparison rather than screening every possible
combination:

- choose one broadly distributed reference capsid (AAV9), one candidate with a
  target-organ advantage, and one negative/off-target control;
- compare route only when it changes the biological entry compartment;
- measure early vector/capsid and later episome, RNA and protein;
- prioritise `PS/Kp`, BBB/CSF access, endosomal escape and SINEUP gain because
  they currently drive ranking uncertainty;
- use the safety screen to add liver, platelet/complement, renal, cardiac or
  neurotoxicity endpoints appropriate to the chosen route.

## What failed or remains unfinished

- The 35:65 vascular/ISF loss split is not identifiable from existing data.
- The CNS three-depth model is reduced-order, not an anatomical brain mesh.
- Kidney filtration and receptor rates are hypotheses, not fitted intact-AAV
  measurements.
- Capsid multipliers combine qualitative evidence across incompatible studies;
  they are priors, not meta-analytic effect sizes.
- The 70 kg model has no calibrated human exposure-to-toxicity function.
- A full loop is incomplete until project-specific wet-lab data change at least
  one fitted parameter and the model re-ranks a design.

These limitations define the next Build and Test, rather than being hidden at
the end of the report.

---

# 中文版本：将模型作为 Design-Build-Test-Learn 工程系统

本项目的工程问题并不只是画出 AAV 浓度曲线。我们需要建立一套可以解释载体去向、区分表达瓶颈、暴露不确定性，并把候选衣壳与给药途径缩小到可由湿实验验证范围内的计算系统。

## 验收标准

在比较生物学方案之前，每次模型构建必须通过四项检查：

1. **质量守恒：** 给入的 AAV 必须等于所有含载体状态与显式累计损失池之和，误差应处于数值求解容差内。
2. **非负性与求解成功：** 不出现有生物学意义的负状态，ODE 积分不得失败。
3. **生物过程分离：** 早期衣壳消失、载体基因组、episome、RNA、蛋白和治疗持续时间必须是不同状态，不能共用一个“半衰期”。
4. **可重复性：** 一组命令能够重新生成前端 JSON/CSV；网页只显示导出的模型结果，不在浏览器中临时编造递送分数。

通过这些检查只能说明模型内部一致，并不等于获得临床有效性或安全性验证。

## 第 1 轮：从人为衰减到可追溯的 AAV9 动力学

**Design。** 最初版本使用统一的血管和 ISF 半衰期制造钟形曲线。这样能让图像易读，却没有解释衣壳从系统中去了哪里，也无法区分转运、摄取与降解。

**Build。** 我们加入血液损失、器官血管/RES 损失、ISF 分解和累计损失池。器官半衰期由 Wang 等人的小鼠数据进行对数线性拟合；同时将普通 AAV9 的血液半衰期修正为 5.0 h，因为原先使用的 2.4 h 实际对应经过四半胱氨酸修饰的 AAV9-TC。

**Test。** 质量守恒残差下降到数值精度范围；不同器官保留各自的达峰与下降过程。拟合诊断也显示，肌肉信号先升后降，单指数模型对它的解释能力较弱。

**Learn。** 一个表观半衰期不能单独识别内皮/RES 清除与 ISF 分解。当前将总速率按 35:65 分配只是显式结构假设，因此被列为高优先级校准参数，而不是伪装成文献事实。

## 第 2 轮：从肝脏单器官摄取到器官与屏障特异性转导

**Design。** 仅比较器官 ISF 暴露无法解释为什么两个载体基因组数量相近的组织会产生完全不同的蛋白表达。

**Build。** 肝脏模块扩展为 `bound -> EE -> LE -> CY -> Ncap -> Nss -> Nds -> Epi -> mRNA -> protein`。肾脏增加小比例滤过后顶端摄取与基底外侧摄取两条入口；CNS 增加 BBB 结合、内皮胞内转运、跨胞转运/回收，以及神经细胞摄取。

**Test。** 情景比较能够分别改变衣壳运输和启动子输出。增强内体逃逸首先明显提高 episome；当表达模块接近饱和后，蛋白提升不再与 episome 成正比。代码审计还发现肝转录速率曾被硬编码限制在 2.0，现在改为使用各情景自己的 `k_tx` 上限，保证启动子预设确实按声明生效。

**Learn。** “进入器官更多”不等于“产生治疗蛋白更多”。因此实验不能只测一个终点，而应在同一体系中区分载体/衣壳、episome、RNA 和蛋白。

## 第 3 轮：从小鼠器官图到给药途径分辨的参考成人

**Design。** 把全部小鼠状态按体重线性放大，会保留错误的循环拓扑、器官比例和 CSF 几何，不能代表人体。

**Build。** 人体版本保留同一套 $Q$--$PS$--$K_p$ 与胞内转运方程，但重新建立右心、肺、左心、动静脉池、门静脉回流、24 个解剖分区和 301 个状态，并加入静脉、鞘内、肌内、枕大池、脑室内和吸入六种给药入口。

**Test。** 区域血容量合计约 5 L；血流分数闭合循环；150 mL CSF 与 500 mL/d 周转对应 4.99 h 的等效一阶半衰期；每条给药途径都单独通过质量守恒检查。

**Learn。** 人体解剖结构可以具有机制意义，但衣壳参数仍可能高度不确定。因此界面将输出标记为“参考人体投影”，并在人源结果中取消依赖小鼠 LY6A 的 PHP.eB CNS 增益，而不把它称作人体预测。

## 第 4 轮：从静态图到面向疾病的设计工具

**Design。** 湿实验使用者应从疾病和靶基因出发，而不是先阅读几十个状态变量名。

**Build。** 批量导出器对每种衣壳和给药途径重新求解 ODE，计算器官特异性、CNS 深度信号和 SINEUP-PD 持续时间，再写入 JSON/CSV。React 前端提供疾病到基因的可展开数据库、CNS 深度配置、衣壳/途径选择、二维设计空间和人体解剖热图。

**Test。** 前端测试验证生产构建、SSR HTML、六条给药途径、24 个分区、时间网格、质量守恒和基因表达缓存。热图同时提供绝对映射与衣壳内相对映射；默认相对映射用于看清空间差异，但图例明确说明 100% 是当前衣壳/时点内的归一化最大值。

**Learn。** 相对颜色适合辨认分布图样，却可能掩盖随时间发生的绝对下降；绝对颜色适合比较量级，却会让低暴露分区看起来相同。因此两种映射必须并存。此外，我们修复了由服务器与客户端本地化时间格式不同导致的 hydration mismatch，改用确定性的 UTC 文本。

## 返回湿实验的当前设计建议

- 先比较一个广泛分布的参考衣壳 AAV9、一个在靶器官有优势的候选衣壳，以及一个阴性或脱靶对照，而不是一次筛选所有组合。
- 只有当给药途径改变了真实生物学入口隔室时才将其作为独立变量，例如 IV 进入循环、IT 进入腰段 CSF、IM 先进入局部肌肉仓。
- 早期测量载体基因组/衣壳，后期测量 episome、RNA 和蛋白。
- 优先校准 `PS/Kp`、BBB/CSF 可及性、内体逃逸和 SINEUP 增益，因为它们目前最容易改变方案排序。
- 根据给药途径和预测暴露，为肝脏、血小板/补体、肾脏、心脏或神经毒性增加相应监测终点。

## 已失败或尚未完成的部分

- 35:65 的血管/ISF 损失分配不能由现有总组织数据唯一识别。
- CNS 三级深度模型是降阶隔室模型，不是脑解剖网格或 3D 扩散模型。
- 肾滤过与受体摄取速率是机制假设，尚未由完整 AAV 实验拟合。
- 衣壳倍率来自检测方法、物种、启动子和给药途径不一致的研究，只能作为先验，不能称为荟萃分析效应量。
- 70 kg 模型没有经人体数据校准的暴露-毒性函数。
- 只有当项目自己的湿实验数据真正更新至少一个参数，并使模型重新排序候选设计时，完整的 DBTL 闭环才算完成。

这些限制不是报告末尾的免责声明，而是下一轮 Build 与 Test 的输入。
