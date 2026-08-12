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
| `wiki/igem_requirements_traceability.md` | PPT/iGEM requirement-to-evidence review table |
| `presentation/model_engineering_script_zh.md` | Concise Chinese presentation script |
| `windows_download_to_web_demo.md` | Windows reproduction and demonstration workflow |
| `aav9_pk_calibration_and_capsid_tropism.md` | AAV9 PK and capsid evidence review |
| `aav_spatial_pk_refinement_report.md` | Earlier refinement analysis retained for provenance |

The generated parameter table is written to
`latex/generated_parameter_table.tex` by
`python model/export_parameter_registry.py`. Do not edit that table manually;
edit the model, evidence metadata or exporter and regenerate it.
