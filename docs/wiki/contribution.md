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
