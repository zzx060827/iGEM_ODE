# Project structure, web pipeline, and parameter roadmap

This document explains how the current repository works as one system. The
English version is followed by a more detailed Chinese version suitable for
team handover and wiki preparation.

## 1. System boundary

The repository contains three connected products:

1. a mouse PBPK and intracellular AAV model used as the mechanistic reference;
2. a 70 kg, route-resolved reference-human projection using the same equation
   family but different anatomy and inputs;
3. a React/vinext atlas that reads versioned model exports and supports disease,
   gene, capsid, route, time, organ and colour-scale exploration.

The browser is a visualisation and ranking layer. It does not integrate the ODE
or estimate parameters. Model data must be regenerated in Python before a new
parameter set appears on the page.

## 2. Repository map

| Path | Role |
|---|---|
| `README.md` | Bilingual project introduction, reproduction commands and safety boundary |
| `requirements.txt` | Reproducible Python dependency ranges |
| `model/ode1.0.py` | Active 64-state mouse PBPK, liver, kidney, CNS, immune, expression and plotting model |
| `model/human_spatial_pbpk.py` | Active 301-state reference-human multiregion model and six route definitions |
| `model/export_delivery_design_space.py` | Batch-solves capsids/routes, derives specificity and SINEUP-PD metrics, writes frontend data |
| `model/aav_parameter_evidence.py` | AAV9 half-life source data, fit windows, fitted rates and human-projection provenance |
| `model/export_parameter_registry.py` | Generates the bilingual parameter CSV/JSON and LaTeX table from active code |
| `model/export_safety_screen.py` | Converts relative organ AUC into monitoring priorities, not clinical risk probabilities |
| `model/spatial_pk_toy_demo.py` | Separate one-dimensional advection-diffusion-reaction demonstration |
| `model/data/` | Generated parameter registry and curated capsid-tropism literature catalog |
| `sineup-delivery-atlas/` | React/vinext application and static model data consumed by the browser |
| `docs/latex/`, `docs/pdf/` | Scientific report, technical appendix and compiled PDFs |
| `docs/wiki/` | Bilingual iGEM Engineering, Attribution, Contribution, Measurement and Safety copy |
| `docs/presentation/` | Concise Chinese presentation script |
| `results/mouse_pbpk/` | Current mouse-model figures and metrics |
| `results/spatial_pk_1d/` | Current one-dimensional spatial-PK demonstration outputs |
| `results/legacy_mouse_pbpk/` | Older plots retained only for provenance |
| `archive/legacy/python/` | Superseded Python implementations; not current entry points |
| `prototypes/ribosome/` | Independent early HTML prototypes unrelated to the current delivery atlas runtime |
| `releases/` | Portable demonstration archive |

## 3. Frontend data flow

```text
literature and physiology priors
        -> Python parameter dictionaries
        -> solve_ivp for every capsid and route
        -> mass-balance and non-negativity checks
        -> model-results.json / human-spatial-results.json
        -> React import or fetch
        -> disease design space and anatomical heat map
```

`model-results.json` contains organ-level mouse results, design-space points,
CNS profiles and route summaries. `human-spatial-results.json` contains the
large 24-region time-series payload. `gene-expression.json` is a cached GTEx/HPA
evidence snapshot. `safety-screen.json` contains monitoring-priority output.

## 4. Frontend implementation

- `app/page.tsx` mounts the application.
- `app/layout.tsx` defines document metadata, fonts and social preview metadata.
- `app/design-space-app.tsx` controls language, disease/gene selection, capsid
  visibility, 2D point selection, expression evidence and exploratory regimen
  ranking.
- `app/disease-data.ts` is the curated disease-to-gene database and contains
  target organ, CNS profile and evidence links.
- `app/organ-heatmap.tsx` reads mouse data at module load, lazily fetches the
  large human payload, interpolates time-series values, maps model regions to
  SVG paths and supports absolute or within-capsid relative colours.
- `app/globals.css` defines the responsive desktop/mobile layout and scroll
  behaviour.
- `public/reactome-male-body-organs.svg` is the licensed anatomical reference;
  model-coloured paths remain code-controlled overlays.
- `scripts/fetch-gene-expression.mjs` refreshes the GTEx/HPA cache.
- `tests/rendered-html.test.mjs` verifies SSR output, route/region counts,
  mass balance, spatial contrast and expression-cache assumptions.

The `db/`, `drizzle/`, `examples/d1/`, `worker/` and authentication helper are
hosting scaffolding. The current public atlas does not need a live database:
the disease library is TypeScript and the numerical model data are static,
versioned JSON. A database becomes useful only when the team needs accounts,
editable records, saved designs or server-side updates.

## 5. Build and regeneration

From the repository root:

```bash
python -m pip install -r requirements.txt
python model/export_parameter_registry.py
python model/export_delivery_design_space.py \
  --output sineup-delivery-atlas/public/data/model-results.json
python model/export_safety_screen.py
```

If the packages are installed in the existing Conda environment, prefix each
Python command with `conda run -n transformer`, or activate that environment
first. The frontend currently versions both npm and pnpm lockfiles; a team
release should choose one package manager and update only its lockfile.

Then build the frontend:

```bash
cd sineup-delivery-atlas
npm install
npm run lint
npm test
npm run dev
```

Do not manually edit generated JSON, CSV or
`docs/latex/generated_parameter_table.tex`; edit the model/evidence source and
regenerate them.

## 6. Highest-value parameter improvements

| Priority | Parameter group | Best evidence | Why it matters |
|---|---|---|---|
| 1 | Capsid-organ `PS`, `Kp`, and tropism multipliers | Same-study, same-route, same-promoter head-to-head vg/cell or RNA data | These parameters currently determine capsid and organ ranking |
| 1 | BBB/CSF access, CNS depth transfer and cell access | Route-matched NHP/human CSF, region and sorted-cell biodistribution | Dominant uncertainty for CNS disease recommendations |
| 1 | Receptor binding, internalisation, endosomal escape and nuclear entry | Time-resolved cell and subcellular fraction measurements | Separates exposure from functional transduction |
| 1 | SINEUP RNA turnover, EC50, maximum gain and target-protein turnover | Target-specific patient-cell dose-response and time course | Directly determines predicted restoration and duration |
| 2 | Organ vascular/ISF loss split | Dual-label capsid plus genome data with perfusion/fractionation | Current 35:65 split is structurally unidentifiable |
| 2 | Kidney filtration and apical/basolateral uptake | Urine, kidney subregion and proximal-tubule cell data | Current dual-entry rates are mechanistic hypotheses |
| 2 | Route depot release and CSF turnover/access multipliers | Route-specific serial imaging, CSF and tissue sampling | Determines local peak and systemic leakage |
| 3 | Reference physiology volumes and flows | ICRP or subject-specific physiology | Usually better documented and less likely to change capsid ranking |

Sensitivity and identifiability analysis should precede additional fitting.
Correlated parameters such as `PS` and `Kp`, or receptor density and binding
rate, cannot be uniquely estimated from one terminal tissue measurement.

---

# 项目文件、网页构建逻辑与参数优化路线

## 1. 项目整体逻辑

当前仓库不是一张独立网页，也不是单个 ODE 文件，而是三层连接系统：

1. **小鼠机制基线。** 使用 64 个状态描述血液到八类器官、肝脏胞内转运、肾近端小管双入口、BBB/CNS、免疫、episome 和表达。
2. **参考人体投影。** 保留小鼠模型的 $Q$--$PS$--$K_p$--摄取--胞内转运方程族，但重新定义 70 kg 成人的循环拓扑、24 个解剖分区、301 个状态、CSF 和六种给药入口。
3. **前端设计工具。** React/vinext 网页读取已经求解并纳入版本控制的 JSON，让用户按疾病、基因、器官、CNS 深度、衣壳、途径和时间浏览结果。

网页本身不运行 SciPy，也不在用户点击时重新拟合模型。改变模型参数后，必须先运行 Python 导出器，新的数值才会进入网页。这种分层让科学计算与界面展示解耦：ODE 可以独立验证，网页也不会因浏览器性能或随机状态产生另一套结果。

## 2. 根目录文件

| 文件/目录 | 具体作用 |
|---|---|
| `README.md` | 项目入口，包含中英文简介、运行命令、目录图和安全边界 |
| `.gitignore` | 排除 Python/Node/LaTeX 缓存、构建目录和系统文件，避免把本机垃圾文件提交到 Git |
| `requirements.txt` | 限定 NumPy、SciPy、Matplotlib 的兼容版本范围 |
| `model/` | 当前有效的科学计算代码与参数证据，是修改模型时首先进入的目录 |
| `sineup-delivery-atlas/` | 当前有效的 React 前端 |
| `docs/` | 科学报告、wiki 文案、汇报稿和复现指南 |
| `results/` | 经过整理的当前结果与历史对照结果 |
| `archive/` | 旧 Python 模型，仅用于追溯，不应作为新实验入口 |
| `prototypes/` | 与当前运行时分离的早期 HTML 原型 |
| `releases/` | 便于复制到其他电脑的演示压缩包 |

## 3. `model/` 中每个文件的作用

### `ode1.0.py`

这是当前小鼠机制模型主文件。它完成四件事：

- `make_params()` 汇总所有小鼠 PBPK、摄取、胞内、免疫和表达参数；
- `rhs()` 定义 64 个状态的微分方程和所有通量；
- `solve_model()` 使用 `solve_ivp` 分阶段求解给药和给药后过程；
- 绘图与指标函数生成器官浓度、肝/肾/CNS 胞内过程、表达、情景比较和质量守恒结果。

器官通量遵循同一基本结构：血流交换、$PS$ 驱动的血管—ISF 交换、$K_p$ 分配、受体摄取和显式损失。肝、肾与 CNS 在这条全身骨架后加入各自的胞内模块。

### `human_spatial_pbpk.py`

这是参考人体多区域模型。它不是把小鼠结果乘以 70 kg，而是重新定义：

- 右心、肺、左心、动脉和静脉等循环池；
- 门静脉器官回流；
- 24 个组织分区的血管容积、ISF 容积、血流分数、`PS`、`Kp`、摄取与 episome 参数；
- IV、IT、IM、ICM、ICV 和 inhaled 六种给药入口；
- 腰段/颅内 CSF、肌肉注射仓和气道给药仓。

`solve_human_capsid()` 每次只求解一个衣壳和一种途径，`mass_balance_error()` 检查输入剂量是否全部存在于状态或累计损失池中。

### `aav_parameter_evidence.py`

保存 AAV9 半衰期的数据依据，而不只保存最后一个数字：

- 小鼠 125I 时间点和器官均值；
- 各器官的拟合窗口；
- 对数线性拟合得到的速率、半衰期和 $R^2$；
- 5.0 h 普通 AAV9 血液半衰期；
- NHP PET 参数和参考人体投影来源。

这样可以检查“参数从哪里来”和“为什么某些时间点被排除”，而不是只看到人为填写的半衰期。

### `export_delivery_design_space.py`

这是模型与网页之间的核心桥梁：

1. 读取八种衣壳先验；
2. 对小鼠和人体情景逐一重新求解 ODE；
3. 计算峰值、Tmax、AUC、器官暴露份额和特异性；
4. 对 CNS 疾病继续求解浅表层、皮层实质与深部核团的降阶模型；
5. 用 episome 输入 SINEUP-PD，得到蛋白恢复和有效持续时间；
6. 写出 `model-results.json/.csv` 与 `human-spatial-results.json`。

二维图横轴不是直接的衣壳标签，而是靶器官暴露相对于加权脱靶暴露的对数特异性；纵轴来自 episome—SINEUP RNA—蛋白 ODE 超过治疗阈值的持续时间。因此两个轴都来自模型，但仍受参数先验影响。

### `export_parameter_registry.py`

它从实际生效的代码参数生成 467 行中英双语注册表，而不是维护一张容易过期的手工表。英文列在前，中文元数据列在后。变量名、DOI、source 和代码位置保持原文，使参数可以直接回查代码和文献。

### `export_safety_screen.py`

读取人体结果，将不同途径/衣壳的器官 AUC 与 IV AAV9 比较，再产生应优先监测的风险器官和实验终点。它是“测什么”的排序器，不是暴露—毒性模型，也不会输出安全剂量。

### `spatial_pk_toy_demo.py`

这是独立的一维流体/空间 PK 演示，使用对流—扩散—反应思想比较流速、扩散、局部捕获和衣壳增强。它用于展示未来 CFD/3D 空间模型的思路，不参与当前人体热图数值。

### `model/data/`

- `aav_capsid_tropism_literature.csv`：按研究记录物种、途径、衣壳、检测方法和局限性。
- `model_parameter_registry.csv/.json`：由导出器生成的中英双语参数表。

## 4. `docs/` 与 `results/` 的作用

### `docs/`

- `latex/aav_sineup_spatial_pk_project_report.tex`：主科学与 iGEM 报告。
- `latex/ode_report.tex`：完整 ODE、参数和技术附录。
- `latex/generated_parameter_table.tex`：由 Python 自动生成，不应手工修改。
- `pdf/`：上述 TeX 的可展示版本。
- `wiki/engineering.md`：四轮 DBTL、验收标准、失败和下一轮实验。
- `wiki/attributions.md`：仓库作者、外部数据/工具、素材许可和团队待确认项。
- `wiki/contribution.md`：未来团队可复用的参数注册表、导出管线与检查方法。
- `wiki/measurement.md`：每个模型状态对应的实验测量、时间点、对照和拟合/验证划分。
- `wiki/safety.md`：AAV 风险域、当前安全筛查边界和下一步暴露—反应模型。
- `wiki/igem_requirements_traceability.md`：PPT/iGEM 要求与现有证据的对照。
- `presentation/model_engineering_script_zh.md`：简短中文讲稿。
- `aav9_pk_calibration_and_capsid_tropism.md`：AAV9 PK 与衣壳偏好证据综述。
- `aav_spatial_pk_refinement_report.md`：早期模型优化分析，保留用于追溯。
- `windows_download_to_web_demo.md`：从 GitHub 下载到 Windows 网页演示的步骤。

### `results/`

- `mouse_pbpk/` 是当前模型生成的结果。
- `spatial_pk_1d/` 是一维空间 PK 演示输出。
- `legacy_mouse_pbpk/` 是旧图，用于说明模型演化，不能混入当前定量结论。

每个 PNG 对应一个模型观察层：早期器官分布、衣壳衰减、肝/肾/CNS 胞内链、表达、质量守恒或情景比较。CSV 保留绘图背后的数字，便于汇报和后续统计。

## 5. 前端目录和构建方法

### 页面入口与状态管理

- `app/page.tsx` 只负责挂载 `DesignSpaceApp`。
- `app/layout.tsx` 设置 `<html>`、字体、标题、favicon、Open Graph 和 Twitter 图片；所有服务器与客户端都使用确定性元数据，避免 hydration mismatch。
- `app/design-space-app.tsx` 管理中英文、疾病/基因选择、衣壳显示、二维点选、正常组织表达证据和组合方案排序。
- `app/disease-data.ts` 是当前疾病数据库。一个疾病可以展开多个基因；每个基因记录靶器官、细胞类型、CNS profile、位点和证据链接。
- `app/organ-heatmap.tsx` 管理小鼠/人体视图、给药途径、衣壳、指标、时间、播放状态、绝对/相对颜色和 SVG 分区点击。
- `app/globals.css` 管理布局、疾病库独立滚动、响应式移动端和视觉样式。

### 数值数据如何进入页面

`model-results.json` 体积较小，直接随模块导入，用于首屏和小鼠视图；23 MB 的 `human-spatial-results.json` 在进入人体视图时异步 `fetch`，避免首屏 JavaScript 把全部人体轨迹打包。React 只进行线性时间插值、筛选、排序和颜色归一化，不改变原始 ODE 轨迹。

热图中的**绝对映射**在选定时间窗内对所有衣壳使用共同上下界，适合比较量级；**相对映射**将当前衣壳、当前时点最高分区定义为 100%，适合看空间形状。相对值变红不表示绝对浓度上升，只表示该分区在当前归一化集合中更接近最大值。

### SVG 如何被染色

Reactome SVG 是背景解剖参考；可变颜色区域由 `organ-heatmap.tsx` 中与模型 region ID 对应的 path/mask 控制。ODE 输出先按 `region_id` 查找数值，再归一化为 0--1，最后经过五段连续色标转换为 CSS fill。要调整偏移或轮廓，应修改对应的 SVG path 或 `humanRegionShape` 映射，而不是改 ODE 参数。

### 数据库是否必要

当前前端不需要在线数据库：疾病记录保存在 TypeScript，模型轨迹与表达证据保存在版本化静态 JSON，这种方式最适合 iGEM wiki 的只读展示。`db/`、`drizzle/`、`examples/d1/`、`worker/` 和 `chatgpt-auth.ts` 是托管脚手架或未来扩展，目前主页面不依赖它们。

只有在需要以下功能时才值得接入 Supabase、D1 或其他数据库：团队在线编辑疾病记录、用户保存方案、后台重新触发模型、权限管理、版本审批或收集用户反馈。即使接入数据库，23 MB 时间序列仍更适合对象存储/CDN，而不是按页面请求逐行查询关系数据库。

### 前端其余配置文件逐项说明

| 文件 | 当前作用 | 是否参与主页面运行 |
|---|---|---|
| `package.json` | 定义 Node 版本、React/vinext 依赖和 `dev/build/test/lint` 命令 | 是 |
| `pnpm-lock.yaml`、`package-lock.json` | 分别锁定 pnpm 与 npm 依赖；团队应统一选择其中一种 | 构建时二选一 |
| `pnpm-workspace.yaml` | 定义 pnpm 工作区边界 | 使用 pnpm 时参与 |
| `tsconfig.json` | TypeScript 类型检查、路径和 JSX 设置 | 是 |
| `eslint.config.mjs` | ESLint 规则 | 只在 lint/开发阶段 |
| `next.config.ts` | Next/vinext 兼容配置 | 是 |
| `vite.config.ts` | Vite、React Server Components 和托管插件入口 | 是 |
| `postcss.config.mjs` | CSS/PostCSS 处理配置 | 是 |
| `.openai/hosting.json` | Sites/托管平台识别和发布配置 | 托管时参与 |
| `build/sites-vite-plugin.ts` | 将托管环境需要的元数据加入 Vite 构建 | 构建/发布时参与 |
| `worker/index.ts` | Cloudflare/Sites worker 入口 | 托管时参与 |
| `app/chatgpt-auth.ts` | 可选的 ChatGPT 登录头解析和安全返回路径 | 当前公开主页未调用 |
| `db/index.ts`、`db/schema.ts` | 可选 D1/Drizzle 数据库连接和 schema | 当前主页未调用 |
| `drizzle.config.ts`、`drizzle/` | 数据库迁移配置与元数据 | 当前主页未调用 |
| `examples/d1/` | D1 写入/API 示例 | 示例，不参与主页 |
| `scripts/fetch-gene-expression.mjs` | 从 GTEx API 聚合疾病基因的组织表达，并缓存 HPA 链接 | 更新数据时参与 |
| `tests/rendered-html.test.mjs` | SSR、模型数据结构、质量守恒、空间差异和表达缓存测试 | 测试时参与 |
| `start-local.command` | macOS 双击启动本地演示的辅助脚本 | 可选 |
| `public/ASSET_ATTRIBUTION.md` | 解剖 SVG 的来源、许可和修改说明 | 文档/合规 |
| `public/favicon.svg`、`public/og.png` | 浏览器图标与社交分享图 | 是 |
| `public/file.svg`、`globe.svg`、`window.svg` | 模板保留图标，当前核心界面不依赖 | 可删除候选 |
| `public/human-anatomy-organs.svg` | DBCLS 早期人体解剖参考 | 对照/备用 |
| `public/reactome-male-body-organs.svg` | 当前人体背景解剖参考 | 是 |

`dist/`、`.vinext/` 和 `.wrangler/` 是构建或本地托管产生的缓存，不属于科学源文件，也不应手工修改或提交。页面异常时应修改 `app/`、`public/` 或配置源文件后重新构建，而不是修补 `dist/` 中的压缩 JavaScript。

### 本地构建与验证

在仓库根目录重新生成科学数据：

```bash
python -m pip install -r requirements.txt
python model/export_parameter_registry.py
python model/export_delivery_design_space.py \
  --output sineup-delivery-atlas/public/data/model-results.json
python model/export_safety_screen.py
```

如果科学计算依赖安装在本项目现有的 Conda 环境中，可先执行
`conda activate transformer`，或将每条 Python 命令写成
`conda run -n transformer python ...`。前端目前同时保留
`package-lock.json` 和 `pnpm-lock.yaml`；正式团队协作时应统一选择 npm 或
pnpm，只更新对应锁文件，避免两台电脑解析出不同依赖版本。

进入前端并运行：

```bash
cd sineup-delivery-atlas
npm install
npm run lint
npm test
npm run dev
```

`npm test` 会先生产构建，再验证 SSR 文本、六种给药途径、24 个分区、287 个时间点、质量守恒、器官空间差异、Tmax 差异和 GTEx/HPA 缓存。发布前还应在桌面和手机视口进行截图检查。

## 6. 哪些参数最值得继续优化

### 第一优先级：会改变方案排序的参数

1. **衣壳—器官 `PS`、`Kp` 与 tropism multiplier。** 当前倍率来自不同研究的文献约束先验。最好的校准数据是同一研究、同一物种、同一途径、同一启动子和同一 readout 下的多衣壳头对头矩阵。不要直接平均 qPCR、RNA、蛋白和 pooled-barcode 数值。
2. **BBB/CSF 与 CNS 深度参数。** 包括 BBB 结合/跨胞转运、CSF 吸收、腰段到颅内输送、脑表面到深部的层间速率和细胞可及性。优先寻找途径匹配的 NHP 数据、脑区 vg/cell、CSF 时间序列和分选细胞数据。
3. **受体结合、内化、内体逃逸和核转运。** 这些参数决定“暴露”能否转化为 episome。需要同一批样本的表面结合、EE/LE、胞质、核内和 episome 时间序列；只有终点组织 qPCR 无法分别识别它们。
4. **SINEUP-PD 参数。** 当前 0.25 d RNA 半衰期、2 d 蛋白半衰期、EC50、最大翻译增益和 65% 阈值都是透明先验。应对每个靶基因使用患者来源细胞做剂量—时间实验，同时测 SINEUP RNA、靶 mRNA、蛋白和功能终点。

### 第二优先级：会改变曲线形状和安全器官判断的参数

1. **35:65 血管/ISF 损失分配。** 单一总组织放射信号不能区分血管清除与 ISF 分解。需要灌流、双标记、组织分级或成像联合基因组数据。
2. **肾脏滤过与双入口摄取。** 应测尿液、肾皮质/髓质、近端小管细胞和血管残留；最好区分 intact capsid、vector genome 与表达。
3. **给药仓释放。** IM 与 inhaled 的释放半衰期和局部比例决定局部峰值与系统泄漏；IT/ICM/ICV 的区域 access multiplier 决定脑区差异，应由路线特异的连续成像、CSF 和组织数据校准。
4. **episome 持久性。** 与细胞分裂、组织更新、启动子沉默和免疫清除有关，应按器官、细胞类型和年龄分别拟合，而不是全身共用一个值。

### 第三优先级：较容易找到权威数据但通常不决定衣壳排序

成人器官容积、血容量、心输出量、血流分数、CSF 容积和生成速率可以继续用 ICRP、MRI、灌流或生理学数据库细化。这些参数重要且相对容易获得，但在当前模型中通常不如衣壳转运与胞内参数那样改变候选排名。

## 7. 当前结果文件应如何阅读

| 结果文件 | 回答的问题 |
|---|---|
| `01_short_distribution_linear_30min.png` | 给药后最早期各器官血管/ISF 如何分配 |
| `02_short_distribution_log_2h.png` | 用对数纵轴观察低丰度器官的早期到达 |
| `03_bell_shaped_aav_decay_48h.png` | 衣壳输入、分布与显式损失共同形成的峰值和下降 |
| `03_liver_intracellular_56d.png` | 肝 ISF、受体结合、内体、胞质、核与表达链的时间顺序 |
| `04_liver_expression_56d_2.png` | 肝 episome、mRNA 与蛋白的长期输出 |
| `05_antibody_56d.png` | 简化免疫/抗体状态对时间的响应 |
| `06_kidney_multilevel_module.png` | 滤过腔、顶端/基底摄取、胞内转运、episome 与表达 |
| `07_liver_vs_kidney_expression.png` | 相同给药下肝与肾的表达差异 |
| `08_mass_balance_audit.png` | 输入剂量、系统内载体、累计损失与守恒残差 |
| `09_design_scenario_comparison.png/.csv` | 衣壳和启动子情景的定量比较 |
| `10_spatial_pk_1d_demo.png` | 一维流动、扩散、捕获对空间终点分布的影响 |
| `11_cns_bbb_and_transduction.png` | 脑血管、BBB、神经细胞胞内链与表达 |
| `12_cns_design_scenario_comparison.png/.csv` | CNS 衣壳增强和启动子增强是否作用于同一瓶颈 |
| `13_normal_aav9_organ_concentration_comparison_log_axes.png/.csv` | 正常 AAV9 器官浓度、AUC、峰值与 Tmax 的跨量级比较 |
| `14_normal_aav9_organ_concentration_comparison_linear_axes.png` | 在普通线性坐标下检查真实量级差异 |

`results/legacy_mouse_pbpk/` 中同名文件来自旧版本，仅用于比较模型演化；引用结果时应优先使用 `results/mouse_pbpk/`。`sineup-delivery-atlas/public/data/` 则保存网页实际读取的机器数据：图与网页只有在这些 JSON/CSV 来自同一次导出时才可视为同一版本。

## 8. 优化参数时应采用的方法

- 先做局部/全局敏感性分析，识别真正影响 specificity、duration、peak、AUC 和安全排序的参数。
- 检查可辨识性；例如 `PS` 与 `Kp`、受体密度与结合速率常高度相关，单个终点不能同时拟合。
- 将参数分层为 direct、fitted、derived、scaled、assumed，并给出置信区间或分布，而不是只给单点。
- 用一部分剂量/时间/衣壳拟合，保留至少一个条件做外部预测。
- 分物种、给药途径、启动子、检测 readout 和时间窗建立数据矩阵，不混合不兼容研究。
- 每次改变参数后重新导出 JSON，运行质量守恒和前端测试，并记录候选排序是否改变。这才形成真正的下一轮 DBTL。
