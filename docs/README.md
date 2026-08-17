# Documentation index

The documentation is divided by audience so that the scientific argument,
technical audit and wiki copy do not drift into one oversized file.

| Path | Intended use |
|---|---|
| `latex/aav_sineup_spatial_pk_project_report.tex` | Main scientific and iGEM-facing report |
| `pdf/aav_sineup_spatial_pk_project_report.pdf` | Compiled main report |
| `latex/ode_report.tex` | Full ODE and parameter technical appendix |
| `pdf/ode_report.pdf` | Compiled technical appendix |
| `wiki/attributions.md` | Attribution Form working record and remaining name checks |
| `wiki/engineering.md` | Design-Build-Test-Learn cycles and acceptance criteria |
| `wiki/contribution.md` | Reusable outputs for future teams |
| `wiki/measurement.md` | Assay-to-state mapping and proposed validation |
| `wiki/safety.md` | AAV hazard identification and model boundary |
| `aav_safety_margin_assessment_2026.md` | Route-matched CNS AAV cases and organ exposure-margin calculation |
| `wiki/igem_requirements_traceability.md` | PPT/iGEM requirement-to-evidence review table |
| `presentation/model_engineering_script_zh.md` | Concise Chinese presentation script |
| `windows_download_to_web_demo.md` | Windows reproduction and demonstration workflow |
| `aav9_pk_calibration_and_capsid_tropism.md` | AAV9 PK and capsid evidence review |
| `aav_spatial_pk_refinement_report.md` | Earlier refinement analysis retained for provenance |
| `project_structure_web_pipeline_and_parameter_roadmap.md` | Bilingual file-by-file architecture, frontend data flow and parameter-optimisation roadmap |

The generated parameter table is written to
`latex/generated_parameter_table.tex` by
`python model/export_parameter_registry.py`. Do not edit that table manually;
edit the model, evidence metadata or exporter and regenerate it.

On the current development machine, the numerical dependencies are installed
in the `transformer` Conda environment. Use `conda activate transformer` before
regeneration if the default Python cannot import NumPy or SciPy.

The five main wiki drafts place the English version first and the complete
Chinese version after a horizontal divider. The parameter registry follows the
same principle at column level: English audit fields come first, followed by
Chinese metadata fields. Parameter identifiers, citations and code locations
remain unchanged so that they can be searched directly in the source.
