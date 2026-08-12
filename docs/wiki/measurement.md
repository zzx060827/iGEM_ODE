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
