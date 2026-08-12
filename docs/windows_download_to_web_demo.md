# iGEM_ODE Windows 下载、建模与网页展示流程

本文档用于在一台新的 Windows 10/11 电脑上，从 GitHub 下载 `zzx060827/iGEM_ODE`，重新运行 AAV PBPK/ODE 模型，生成前端数据，并启动 SINEUP Delivery Atlas 网页。

## 1. 项目数据链

```text
model/ode1.0.py
        +
model/human_spatial_pbpk.py
        |
        v
model/export_delivery_design_space.py
        |
        +--> sineup-delivery-atlas/public/data/model-results.json
        +--> sineup-delivery-atlas/public/data/model-results.csv
        +--> sineup-delivery-atlas/public/data/human-spatial-results.json
                         |
                         v
              React/vinext 前端读取并展示
```

后端和前端不是通过实时服务器 API 连接，而是通过 `public/data` 中的 JSON/CSV 文件连接。因此，修改 ODE 后必须重新执行导出脚本，网页才会显示新结果。

## 2. 推荐软件

安装以下软件：

1. Git for Windows：<https://git-scm.com/download/win>
2. Miniconda：<https://docs.conda.io/projects/miniconda/en/latest/>
3. Node.js 22 LTS：<https://nodejs.org/>
4. Visual Studio Code：<https://code.visualstudio.com/>

安装 Node.js 后确认版本不低于 `22.13.0`。打开 PowerShell：

```powershell
git --version
node --version
npm --version
conda --version
```

如果 PowerShell 找不到 `conda`，从开始菜单打开 **Anaconda Prompt**，或者重新打开终端。

## 3. 从 GitHub 下载

建议把项目放在不含中文和空格的短路径，例如 `D:\igem`：

```powershell
D:
mkdir igem
cd igem
git clone https://github.com/zzx060827/iGEM_ODE.git
cd iGEM_ODE
git status
```

以后获取 GitHub 上的更新：

```powershell
cd D:\igem\iGEM_ODE
git pull origin main
```

不要把网页下载成 ZIP 后长期开发，因为 ZIP 版本不便于同步、比较和提交修改。

## 4. 创建 Python 建模环境

在 Anaconda Prompt 或已经能使用 `conda` 的 PowerShell 中执行：

```powershell
cd D:\igem\iGEM_ODE
conda create -n igem-ode python=3.11 -y
conda activate igem-ode
python -m pip install --upgrade pip
python -m pip install numpy scipy matplotlib
```

检查环境：

```powershell
python -c "import numpy, scipy, matplotlib; print('Python model environment OK')"
```

每次重新打开终端后，都要先执行：

```powershell
conda activate igem-ode
```

## 5. 运行 ODE 并导出网页数据

在仓库根目录执行：

```powershell
cd D:\igem\iGEM_ODE
conda activate igem-ode
python model\export_delivery_design_space.py --model model\ode1.0.py --output sineup-delivery-atlas\public\data\model-results.json
```

成功后应看到类似输出：

```text
Wrote ... modeled design points to ...model-results.json
Wrote ... routes x ...-region human trajectories to ...human-spatial-results.json
```

该命令会生成或更新：

- `model-results.json`：二维设计空间、器官结果、CNS profile 和模型元数据；
- `model-results.csv`：便于 Excel/R/Python 检查的表格；
- `human-spatial-results.json`：人体分区、给药途径和时间轨迹；

可以用以下命令检查文件是否存在：

```powershell
Get-ChildItem sineup-delivery-atlas\public\data
```

## 6. 安装前端依赖

进入前端目录：

```powershell
cd D:\igem\iGEM_ODE\sineup-delivery-atlas
npm ci
```

首次安装可能需要几分钟。`npm ci` 会严格按照 `package-lock.json` 安装，换电脑展示时比 `npm install` 更可重复。

## 7. 在 Windows 上启动网页

当前 `package.json` 的 `npm run dev` 使用了 Unix 风格的环境变量写法。纯 PowerShell 下推荐直接执行：

```powershell
cd D:\igem\iGEM_ODE\sineup-delivery-atlas
npx vinext dev
```

终端会输出本地地址，通常类似：

```text
http://localhost:3000/
```

如果 3000 端口已占用，工具可能自动选择 3001、3012 等端口，以终端实际输出为准。保持终端窗口运行，在浏览器打开该地址。

停止服务器：在运行服务器的终端按 `Ctrl+C`。

### Git Bash 方式

如果使用 Git Bash，可以执行项目原始命令：

```bash
cd /d/igem/iGEM_ODE/sineup-delivery-atlas
npm run dev
```

## 8. 构建正式版本

展示前至少运行一次构建和测试：

```powershell
cd D:\igem\iGEM_ODE\sineup-delivery-atlas
npx vinext build
npm test
```

由于 `npm test` 内部也调用 Unix 风格的 `npm run build`，若在 PowerShell 中失败，可以改在 Git Bash 中运行：

```bash
npm run build
npm test
```

本地演示不要求数据库。当前核心模型数据来自静态 JSON 文件，D1/Drizzle 是可选扩展。

## 9. 一次完整的修改闭环

修改模型后按以下顺序操作：

```powershell
cd D:\igem\iGEM_ODE
conda activate igem-ode

# 1. 修改 model\ode1.0.py、model\human_spatial_pbpk.py 或导出逻辑

# 2. 重新求解并生成前端数据
python model\export_delivery_design_space.py --model model\ode1.0.py --output sineup-delivery-atlas\public\data\model-results.json

# 3. 启动前端
cd sineup-delivery-atlas
npx vinext dev
```

如果服务器已经运行，重新生成 JSON 后通常刷新网页即可；若未更新，停止服务器并重新启动，然后在浏览器执行强制刷新 `Ctrl+F5`。

修改前端时主要关注：

- `app/design-space-app.tsx`：疾病库、二维设计空间和主要交互；
- `app/organ-heatmap.tsx`：人体热图、时间轴、给药途径和 SVG 映射；
- `app/disease-data.ts`：疾病与基因数据；
- `app/globals.css`：布局和视觉样式；
- `public/data/`：后端模型导出的前端输入；

修改模型时以 `model/` 目录为主版本。旧脚本集中保存在
`archive/legacy/python/`，仅用于追溯，不要与现行脚本同步修改。

## 10. 提交到 GitHub

建议每个功能建立分支：

```powershell
cd D:\igem\iGEM_ODE
git switch -c feature\better-human-pbpk
```

Windows/Git 分支名更推荐使用正斜杠：

```powershell
git switch -c feature/better-human-pbpk
```

检查和提交：

```powershell
git status
git diff
git add model sineup-delivery-atlas docs
git commit -m "Improve human PBPK model and atlas visualization"
git push -u origin feature/better-human-pbpk
```

不要提交：

- `node_modules/`；
- Conda 环境目录；
- `.env` 密钥；
- 临时缓存和大型无关输出；

## 11. 换电脑演示的最短流程

电脑已经安装 Git、Node.js 和 Conda 时，只需：

```powershell
git clone https://github.com/zzx060827/iGEM_ODE.git
cd iGEM_ODE

conda create -n igem-ode python=3.11 -y
conda activate igem-ode
python -m pip install numpy scipy matplotlib
python model\export_delivery_design_space.py --model model\ode1.0.py --output sineup-delivery-atlas\public\data\model-results.json

cd sineup-delivery-atlas
npm ci
npx vinext dev
```

如果只是展示已经提交到 GitHub 的计算结果，可以跳过 Python 求解，直接安装前端依赖并启动。为了证明网页数据确实来自当前 ODE，正式汇报前仍建议完整运行一次导出命令。

## 12. 常见问题

### `conda` 不是内部或外部命令

使用 Anaconda Prompt，或执行：

```powershell
conda init powershell
```

然后关闭并重新打开 PowerShell。

### PowerShell 禁止运行脚本

如果出现 `npm.ps1 cannot be loaded`，可以使用：

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

也可以改用 Git Bash。

### `npm run dev` 提示 `WRANGLER_LOG_PATH` 不是命令

这是 Unix/Windows shell 语法差异。PowerShell 使用：

```powershell
npx vinext dev
```

长期方案是让另一个 AI 在前端安装 `cross-env`，并把脚本改为跨平台写法。

### 网页仍显示旧结果

依次检查：

1. 是否运行了 `model\export_delivery_design_space.py`；
2. 输出路径是否为 `sineup-delivery-atlas\public\data\model-results.json`；
3. `human-spatial-results.json` 的修改时间是否更新；
4. 是否在正确的前端目录启动服务器；
5. 浏览器是否执行了 `Ctrl+F5`；

### Hydration mismatch

不要在服务端首屏直接渲染 `Date.now()`、`Math.random()` 或本地时区格式化结果。模型的 `generated_at` 应先按固定 UTC 字符串渲染，或者等组件挂载后再进行本地化显示。

### 页面打开但人体图没有颜色

检查浏览器开发者工具的 Network/Console：

- `/data/human-spatial-results.json` 是否返回 200；
- SVG 区域 ID 是否能映射到模型 region ID；
- JSON 是否因求解中断而只写入了一部分；

## 13. 交给另一个 AI 的接手提示词

可以把下面这段连同 GitHub 地址一起发给另一个 AI：

```text
请接手并完善 GitHub 项目：
https://github.com/zzx060827/iGEM_ODE

这是一个 AAV-SINEUP 全身递送、细胞内转导、CNS 深度分布和长期药效项目。

请先阅读：
1. docs/windows_download_to_web_demo.md
2. docs/aav_spatial_pk_refinement_report.md
3. model/ode1.0.py
4. model/human_spatial_pbpk.py
5. model/export_delivery_design_space.py
6. sineup-delivery-atlas/app/design-space-app.tsx
7. sineup-delivery-atlas/app/organ-heatmap.tsx

工程边界：
- model/ 是当前后端模型主版本；根目录同名 Python 文件可能是历史副本。
- Python 通过 solve_ivp 求解 ODE，并把结果导出到
  sineup-delivery-atlas/public/data/。
- React/vinext 前端读取这些 JSON，不要把模型结果重新写成前端手工常量。
- 修改模型后必须重新运行导出器，并同时验证 JSON、CSV 和网页。
- 人体结果是 reference-human projection，不得表述为临床预测。
- 不同衣壳参数必须保存物种、实验方法、来源和证据等级。
- PHP.eB 的 LY6A 小鼠增益不能直接用于人体。

运行后端：
python model/export_delivery_design_space.py --model model/ode1.0.py --output sineup-delivery-atlas/public/data/model-results.json

运行前端：
cd sineup-delivery-atlas
npm ci
npx vinext dev

开始修改前，请先运行现有模型、构建前端、检查 git status，并说明你验证了哪些输出。每次修改保持后端方程、导出 schema、TypeScript 类型和页面显示同步。不要删除现有科学限制说明或伪造文献参数。
```

## 14. 推荐的下一步工程整理

为了让后续开发更稳定，建议依次完成：

1. 增加根目录 `requirements.txt` 或 Conda `environment.yml`；
2. 用 `cross-env` 修复 Windows 前端脚本；
3. 删除或归档根目录重复 Python 文件，只保留 `model/` 主版本；
4. 为导出 JSON 增加 schema 校验；
5. 为 ODE 增加质量守恒、非负性、参数来源和回归测试；
6. 在前端明确显示模型版本、生成时间和证据等级；
