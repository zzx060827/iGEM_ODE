# AAV 全身递送、入胞与 Spatial PK 建模优化报告

针对现有脚本 `/Volumes/zxzhu/iGEM_ODE/aav_pbpk_liver_plus_multilevel_kidney.py`。

## 1. 你现在已经做对的部分

现有模型已经具备一个很适合 iGEM dry lab 展示的主骨架：

- `blood -> organ vascular -> organ ISF` 的全身 PBPK 分布层。
- 肝脏细胞命运模块：表面结合、早/晚内体、胞质逃逸、入核、单链/双链/episome、mRNA、protein。
- 肾脏近端小管模块：肾血管滤过到 filtrate/PT lumen、apical binding/endocytosis、basolateral ISF uptake、回收内体、溶酶体、表达输出。
- 10 min infusion、半衰期式 apparent clearance、短期分布图和 56 天表达图，适合做故事线。

最大短板不是结构，而是“参数可信度”和“空间层级”。下一阶段应把模型从 demo ODE 变成可校准的 multiscale QSP/spatial-PK 框架。

## 2. 推荐升级路线

### A. 把 PBPK 层改成质量守恒、可校准版本

现在的 `Q_scale` 是为了让曲线好看而缩放心输出量。这在展示中可以用，但答辩或论文式呈现时容易被问到。建议改成两个层次：

1. **真实血流层**：保留文献/解剖血流量、器官体积、血浆/组织分配参数。
2. **有效交换层**：把内皮屏障、窦状隙、肾小球尺寸限制、组织间质扩散写成 `PS` 或 permeability-limited exchange，而不是缩放 `Q`。

建议新增的质量平衡输出：

- cumulative dose input
- extracellular AAV
- intracellular AAV
- urine loss
- RES/nonspecific loss
- lysosomal degradation
- antibody neutralized loss
- numerical mass-balance error

### B. 将“经验半衰期清除”拆成可解释机制

脚本中 `k_clear_blood / k_clear_vascular / k_clear_isf` 目前只是为了产生 bell-shaped curve。建议拆成：

- 肝 Kupffer cell / spleen macrophage uptake
- vascular endothelial nonspecific uptake
- neutralizing antibody/opsonization
- complement/innate immune activation
- renal/urinary loss
- target-cell receptor-mediated uptake

这样同样能得到下降曲线，但每个清除项都能对应一个实验读数。

### C. 用 receptor occupancy 替代单一线性 binding

肝模块现在 `R_tot` 只有一个肝细胞表面受体池。AAV 入胞文献强调 AAVR 是跨 serotype 的关键 receptor，同时不同血清型还依赖 glycan/co-receptor，且 endosomal escape 和 nuclear trafficking 是重要瓶颈。建议扩展成：

```text
surface AAV
  -> glycan/co-receptor binding
  -> AAVR-dependent internalization
  -> endosomal sorting
  -> lysosomal degradation / recycling / escape
  -> nuclear import
  -> uncoating
  -> episome
```

每个 organ/cell type 可以有不同 `R_tot`、`k_int`、`k_escape`、`k_lys`。这会让 capsid engineering、tissue tropism、promoter specificity 都能进入模型。

### D. 肾脏模块建议从“滤过为主”改成“双入口但滤过受限”

AAV 颗粒直径约 25 nm，完整 capsid 通常不应被大量肾小球滤过。你现在已经把 `k_glom_filter` 设得很小，这是合理方向。但后续建议强调：

- intact AAV 的 glomerular filtration 是 size/charge-restricted apparent route。
- 肾表达更可能来自 peritubular capillary/interstitium -> basolateral uptake，或特定 capsid 对肾细胞 tropism 的增强。
- `K_Urine` 可以代表 free genome/capsid fragments 或极少量 intact capsid 的 lumped loss，不要过度解释为完整 AAV 大量尿排。

### E. 加入 serotype/capsid 和 promoter 两层设计变量

建议建立参数表：

| 设计变量 | 模型参数 |
|---|---|
| capsid serotype / engineered capsid | organ `PS`, `Kp`, receptor `R_tot`, `k_on`, `k_int`, `k_escape` |
| promoter | organ/cell-type `k_tx`, `EC50_tx`, max expression |
| miRNA detargeting | liver/spleen expression degradation or translation inhibition |
| dose/infusion route | `dose_vg`, `T_inf`, boundary concentration |
| local fluid control | local velocity, residence time, wall flux, shear-sensitive uptake |

### F. Spatial PK / CFD 的最小可行升级

你的最终目标是“仿真建模 + 流体力学达到好的 spatial PK”。建议三层耦合：

1. **0D PBPK**：给出每个器官入口浓度 `C_in(t)`。
2. **1D/2D organ microtransport**：在血管、肝窦、肾小管或目标组织中解 advection-diffusion-reaction。
3. **cellular fate ODE**：每个空间网格点都有 receptor binding、endosomal trafficking、episome、expression。

核心方程可以先写成：

```text
∂C/∂t + u(x)∂C/∂x = D∂²C/∂x² - k_uptake(C, R)C - k_loss C
dB/dt = k_on C (Bmax - B) - k_off B - k_int B
dE/dt = k_escape I - k_loss_epi E
```

这就是从现在 ODE 走向 spatial PK/CFD 的桥。

## 3. 我建议你下一版代码具体怎么改

- 把参数拆成 `physiology`, `vector`, `cell_fate`, `immune`, `simulation` 五个 dict/dataclass。
- 增加 `run_scenario(params, label)`，一次性比较 serotype、dose、route、flow-control。
- 加 `mass_balance(sol)`，所有答辩图先过质量守恒。
- 加 `fit_or_calibrate()` 接口：先支持手动输入 literature/实验 biodistribution CSV，拟合 `PS`, `k_res`, `k_escape`, `k_tx`。
- 把 plotting 改成输出 metrics CSV + figures，方便交给 wet lab 同学。
- 对肾模块，默认展示两套假设：`filtration_limited` 和 `basolateral_tropism`，避免“完整 AAV 大量滤过”的生物学质疑。

## 4. 建议实验/数据闭环

最小实验设计：

- qPCR/ddPCR：blood、liver、kidney、spleen、urine 的 vector genome time course。
- tissue section 或 reporter imaging：空间分布，尤其肝小叶 periportal/pericentral 或肾皮质/髓质差异。
- cell marker co-localization：proximal tubule、endothelium、immune cells。
- mRNA/protein：验证 episome 到 expression 的滞后。
- neutralizing antibody 或 cytokine markers：给 immune/clearance 模块定量。

## 5. 可引用的代表文献

1. **AAV PBPK/QSP 建模最直接的标杆**  
   “Whole-Body Disposition and Physiologically Based Pharmacokinetic Modeling of Adeno-Associated Viruses and the Transgene Product,” *Journal of Pharmaceutical Sciences*, 2024.  
   https://www.sciencedirect.com/science/article/pii/S0022354923004148

1. **AAVR 作为 AAV 转导核心受体**  
   Pillay et al., *Nature*, 2016, “An essential receptor for adeno-associated virus infection.”  
   https://www.nature.com/articles/nature16465

2. **AAV 入胞和胞内运输综述**  
   Nonnenmacher & Weber, *Gene Therapy*, 2012, “Intracellular transport of recombinant adeno-associated virus vectors.”  
   https://www.nature.com/articles/gt20126.pdf

3. **更新的 AAV intracellular trafficking 综述**  
   Riyad & Weber, *Gene Therapy*, 2021, “Intracellular trafficking of adeno-associated virus (AAV) vectors: challenges and future directions.”  
   https://www.nature.com/articles/s41434-021-00243-z

3. **AAV 免疫与高剂量系统递送风险综述**  
   Verdera et al., *Molecular Therapy*, 2020, “AAV vector immunogenicity in humans: a long journey to successful gene transfer.”  
   https://pmc.ncbi.nlm.nih.gov/articles/PMC7054726/

4. **AAV capsid engineering 和 tropism 优化综述**  
   Li & Samulski, *Nature Reviews Genetics*, 2020, “Engineering adeno-associated virus vectors for gene therapy.”  
   https://www.nature.com/articles/s41576-019-0205-4

5. **高通量 AAV capsid 设计/筛选思路**  
   Bryant et al., *Nature Biotechnology*, 2021, “Deep diversification of an AAV capsid protein by machine learning.”  
   https://www.nature.com/articles/s41587-020-00793-4

6. **肾近端小管 megalin/cubilin endocytosis 基础**  
   Christensen & Birn, *Nature Reviews Molecular Cell Biology*, 2002, “Megalin and cubilin: multifunctional endocytic receptors.”  
   https://www.nature.com/articles/nrm778

7. **药物/颗粒递送中的 transport 与 CFD 思路**  
   Hossain et al., *Advanced Drug Delivery Reviews*, 2019, “Mathematical modelling of tumour drug delivery.”  
   https://pubmed.ncbi.nlm.nih.gov/31325509/

8. **PBPK 建模原则参考**  
   Jones & Rowland-Yeo, *CPT: Pharmacometrics & Systems Pharmacology*, 2013, “Basic concepts in physiologically based pharmacokinetic modeling in drug discovery and development.”  
   https://ascpt.onlinelibrary.wiley.com/doi/full/10.1038/psp.2013.41

9. **AAV whole-body imaging / biodistribution 数据思路**  
   Ballon et al., *Human Gene Therapy*, 2020, “Quantitative Whole-Body Imaging of I-124-Labeled Adeno-Associated Viral Vector Biodistribution in Nonhuman Primates.”  
   https://journals.sagepub.com/doi/10.1089/hum.2020.116

## 6. 现阶段最适合 iGEM 展示的故事线

建议把项目讲成：

> We built a multiscale AAV delivery simulator that connects systemic PBPK, organ exposure, receptor-mediated cellular entry, intracellular trafficking, and spatial fluid transport. The model is not only descriptive: it can test how capsid tropism, infusion route, receptor saturation, renal uptake pathways, and local flow control reshape spatial pharmacokinetics and final expression.

对应图建议：

- 全身 PBPK organ exposure。
- 肝 vs 肾 entry bottleneck。
- receptor saturation 敏感性分析。
- 1D spatial PK demo：同样总剂量下，改变流速/局部输入后 expression spatial uniformity 改变。
- 未来 CFD 框架图：PBPK boundary -> CFD organ flow -> cellular fate。
