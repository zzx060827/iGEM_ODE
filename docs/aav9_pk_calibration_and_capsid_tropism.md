# AAV9 early capsid-PK calibration and capsid-tropism evidence

## 1. What was calibrated

This update calibrates the early decline of circulating and tissue-associated
AAV9 capsid. It does **not** equate capsid disappearance with loss of episome,
mRNA, protein, or therapeutic effect.

The model still separates:

1. extracellular capsid PK: blood, vascular space, ISF, uptake and degradation;
2. intracellular trafficking: bound, EE, LE, CY, Ncap, Nss, Nds and Epi;
3. expression/PD: episome, mRNA, SINEUP RNA and protein.

## 2. Mouse AAV9 fit

The old `2.4 h` blood half-life was replaced by `5.0 h`. Seo et al. reported
`5.0 h` for unmodified AAV9, while `2.4 h` belonged to tetracysteine-modified
AAV9-TC.

For tissues, the means in Wang et al. Table S1 were transcribed at 4, 24, 72
and 168 h. The non-residualizing 125I signal is used as an early
intact/extracellular-capsid proxy. Each selected decay segment is fitted by

```text
ln C_i(t) = alpha_i - k_i t
t_1/2,i = ln(2) / k_i
```

| Tissue | fitted half-life (h) | log-fit R2 | fit points (h) |
|---|---:|---:|---|
| Blood | 5.00 | direct literature estimate | PET circulation |
| Liver | 58.41 | 0.966 | 4, 24, 72, 168 |
| Spleen | 61.66 | 0.948 | 4, 24, 72, 168 |
| Kidney | 21.25 | 0.995 | 4, 24, 72 |
| Heart | 23.17 | 0.988 | 4, 24, 72, 168 |
| Muscle | 32.65 | 0.786 | 4, 24, 72, 168 |
| Lung | 21.42 | 0.942 | 4, 24, 72, 168 |
| Brain | 15.40 | 0.975 | 4, 24, 72 |
| Rest composite | 26.43 | 0.978 | 4, 24, 72, 168 |

The kidney day-7 point was excluded because it forms a near-floor plateau. The
brain day-7 mean was zero and cannot enter a log fit. Muscle has the weakest fit
because its mean first rises and then falls; its single-exponential half-life is
therefore a low-confidence effective parameter.

The raw table values, fit windows, rates and R2 values are stored in
`model/aav_parameter_evidence.py` and are recomputed when the module loads.

## 3. Reference-human projection

No organ-resolved clinical human AAV9 capsid time course suitable for this PBPK
fit was identified. The 70 kg model therefore uses immune-naive NHP I-124-AAV9
PET half-lives from Ballon et al. where reported:

| Compartment | half-life (h) | source status |
|---|---:|---|
| Circulating source | 1.2 | NHP PET |
| Liver | 22.6 | NHP PET |
| Heart | 15.6 | NHP PET |
| Spleen | 22.9 | NHP PET |
| Brain/CSF signal | 24.8 | NHP PET |
| Muscle | 48.7 | NHP body-remainder surrogate |
| Rest | 48.7 | NHP body-remainder surrogate |
| Kidney | 21.25 | provisional mouse fit |
| Lung | 21.42 | provisional mouse fit |

These are a **reference-human projection**, not a clinical prediction. The NHP
study has a very small animal count, the IV brain signal was largely associated
with CSF, and radiolabeled capsid fate is not identical to functional
transduction.

Human physiology now also includes two executable checks:

- central plus regional vascular volume = 5059 mL, within 2% of a 5 L adult
  reference blood volume;
- 150 mL CSF and 500 mL/day production imply an equivalent first-order CSF
  turnover half-life of 4.99 h.

## 4. Remaining structural assumption

The observed organ decline does not separately identify endothelial/RES loss
and ISF catabolism. The current ODE assigns 35% of the apparent organ loss rate
to vascular/endothelial removal and 65% to ISF/catabolic loss:

```text
k_res,i = 0.35 ln(2) / t_1/2,i
k_deg,isf,i = 0.65 ln(2) / t_1/2,i
```

This split is a transparent structural prior, not a literature-fitted result.
It should later be estimated jointly against vascular, extracellular,
internalized and vector-genome observations. A dual-label study is especially
useful for that next step.

## 5. Multi-capsid organ-tropism studies

The machine-readable catalog is
`model/data/aav_capsid_tropism_literature.csv`; the exporter also copies it to
the frontend at `/data/aav_capsid_tropism_literature.csv`.

High-value studies include:

- Zincarelli et al. (2008): AAV1-AAV9, same systemic mouse framework. AAV9 was
  broadest; AAV4 had strong lung/kidney genome delivery; AAV6 favored heart,
  liver and skeletal muscle.
- Walkey et al. (2025): 10 natural serotypes across 22 tissues in male and
  female mice, with qPCR and functional reporter readouts. It exposes important
  sex and readout effects.
- Yang et al. (2025): 21 natural/engineered capsids in C57BL/6, BALB/c and
  cynomolgus macaques. It is currently the most useful cross-species source for
  AAV9, rh10, LK03, PHP.eB and CAP-B10.
- Abele et al. (2025): 34 barcoded capsids after IV versus intraperitoneal
  delivery. It demonstrates that tropism scores must be conditioned on route.
- Liu et al. (2024): AAV8/AAV9 whole-body disposition and PBPK fitting over
  three weeks, useful for the eventual joint fit of vector and transgene data.

The current relative capsid-organ multipliers remain literature-constrained
priors. They were not numerically refitted in this update because values from
different species, promoters, doses, routes, pooled libraries, DNA readouts and
RNA/protein readouts are not directly interchangeable. The next defensible
upgrade is to choose one primary head-to-head matrix per species/route/readout,
normalize within that study, then estimate uncertainty instead of averaging raw
numbers across studies.

## 6. Primary sources

- Wang et al. 2024: https://pmc.ncbi.nlm.nih.gov/articles/PMC11404148/
- Seo et al. 2020: https://pmc.ncbi.nlm.nih.gov/articles/PMC7193641/
- Ballon et al. 2020: https://pmc.ncbi.nlm.nih.gov/articles/PMC7769048/
- Zincarelli et al. 2008: https://pubmed.ncbi.nlm.nih.gov/18414476/
- Walkey et al. 2025: https://pubmed.ncbi.nlm.nih.gov/39863928/
- Yang et al. 2025: https://pmc.ncbi.nlm.nih.gov/articles/PMC11919325/
- Abele et al. 2025: https://pubmed.ncbi.nlm.nih.gov/40337478/
- Liu et al. 2024: https://doi.org/10.1016/j.xphs.2023.10.005
