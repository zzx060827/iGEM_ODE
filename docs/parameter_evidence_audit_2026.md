# Model parameter evidence audit (2026-08-13)

## Scope and decision rule

This audit covers all 467 rows exported by `model/export_parameter_registry.py`.
Evidence was classified by what it can actually constrain:

- **quantitative-direct / quantitative-derived**: a measured value, fitted
  time course, or transparent calculation is comparable to the model quantity;
- **physiology-scaled**: anatomy or flow is supported by physiological data and
  scaled into the model geometry;
- **comparative-prior**: a head-to-head study supports a relative ordering or a
  persistence lower bound, but not an absolute universal constant;
- **mechanistic-only**: the citation supports the modeled step but does not
  identify its numerical value;
- **clinical-context**: a clinical dose or observation supplies context without
  implying product, population, route, or safety equivalence;
- **structural/design**: no external numerical calibration is claimed.

The registry contains 324/467 (69.4%) rows with some external evidence and
197/467 (42.2%) rows with quantitative, physiological, or head-to-head
comparative constraint. The latter is the strict coverage figure; mechanisms
were not counted as numeric support.

For prioritizing calibration experiments, repeated region/capsid rows can also
be grouped into 118 distinct parameter concepts (for example, all 24
`vascular_ml` rows are one physiological concept). Evidence coverage is higher
by row than it first appears because several unsupported mechanism concepts are
repeated across only a few compartments, while the 64 capsid-organ comparison
rows are explicitly supported as comparative priors.

## Parameters changed

| Parameter(s) | Old | New | Evidence and interpretation |
|---|---:|---:|---|
| Mouse `CO` | 25 mL/h | 840 mL/h | Reference unanesthetized mouse cardiac output is about 14 mL/min. The effective exchange multiplier was changed from 0.05 to 1.25/840 so the previously intended 1.25 mL/h pulmonary exchange rate is preserved. |
| Mouse liver vascular / ISF volume | 0.14 / 0.22 mL | 0.095 / 0.190 mL | Replaced with published large-protein PBPK physiological spaces. |
| Mouse kidney vascular / ISF volume | 0.05 / 0.08 mL | 0.030 / 0.101 mL | Same source and definition. |
| Mouse heart vascular / ISF volume | 0.012 / 0.025 mL | 0.007 / 0.019 mL | Same source and definition. |
| Mouse muscle vascular / ISF volume | 0.40 / 1.20 mL | 0.150 / 1.032 mL | Same source and definition. |
| Mouse lung vascular / ISF volume | 0.08 / 0.04 mL | 0.0191 / 0.057 mL | Same source and definition. |
| Human liver episome effective persistence prior | 120 d | 1095 d | Human liver biopsies found transcriptionally competent AAV episomes 2.6-4.1 years post-dose. 1095 d is a conservative prior, not a fitted half-life. |
| Human skeletal-muscle episome effective persistence prior | 365 d | 1460 d | Human muscle biopsies found circular, transcriptionally competent AAV genomes persisting at least four years. 1460 d is a lower-bound prior. |

## Parameters deliberately not changed

- Liu et al. reported steady-state whole-tissue-to-blood ratios of about 50 for
  liver, 10 for heart/muscle, 2 for brain/lung/kidney/adipose/spleen, and at
  most 1 for bone/skin/pancreas. These observations are now recorded as
  calibration targets, but they were not copied into `Kp`: whole-tissue vector
  DNA includes vascular, extracellular, cell-associated, and nuclear material,
  whereas this model's `Kp` acts only in the vascular-to-ISF exchange equation.
- AAV2 internalizes into HeLa cells with a half-time under 10 minutes and
  reaches the nuclear/perinuclear region within about two hours. These data
  reveal large differences from several model trafficking rates, but they are
  not direct AAV9 in-vivo organ estimates; current values remain explicitly
  low-confidence calibration priors.
- Capsid-organ multipliers were not averaged across studies. Zincarelli,
  Walkey, Yang, Abele, Wang, Ballon, and Liu use different species, strains,
  routes, doses, pooled versus single-vector designs, promoters, and DNA versus
  RNA/protein readouts. The registry treats the multipliers as directionally
  constrained relative priors.
- The FDA ZOLGENSMA label supports the clinical context dose of 1.1e14 vg/kg
  and liver-dominant postmortem vector DNA, but does not validate a 70 kg adult
  AAV9 dose or establish equivalence across payloads and populations.

## Core sources

- AAV9 early mouse capsid PK: https://pmc.ncbi.nlm.nih.gov/articles/PMC11404148/
- AAV9/PHP.eB PET and blood half-life: https://pmc.ncbi.nlm.nih.gov/articles/PMC7193641/
- AAV9/rh.10 NHP PET: https://pmc.ncbi.nlm.nih.gov/articles/PMC7769048/
- AAV8/AAV9 PBPK and tissue-to-blood targets: https://doi.org/10.1016/j.xphs.2023.10.005
- Mouse physiological spaces: https://pmc.ncbi.nlm.nih.gov/articles/PMC3727051/
- Mouse cardiac output reference: https://pubmed.ncbi.nlm.nih.gov/8378254/
- Human large-protein PBPK physiology: https://pmc.ncbi.nlm.nih.gov/articles/PMC6890583/
- AAV cell-entry kinetics: https://doi.org/10.1128/JVI.74.6.2777-2785.2000
- Human liver episomes: https://pmc.ncbi.nlm.nih.gov/articles/PMC9018415/
- Human muscle episomes: https://pmc.ncbi.nlm.nih.gov/articles/PMC5374867/
- ZOLGENSMA prescribing information: https://www.fda.gov/media/126109/download
- Human AAV9 postmortem observation: https://www.nejm.org/doi/full/10.1056/NEJMoa2307798
- Multi-capsid sources and limitations are catalogued in
  `model/data/aav_capsid_tropism_literature.csv`.
