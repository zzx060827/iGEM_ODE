# Attributions: what we built, what we reused, and what remains to verify

This page is the working record behind the official iGEM Attributions Form. It
should be updated whenever a person, dataset, software package or external
service changes the project. The final form, not this draft alone, is the
competition deliverable.

## Work demonstrably completed in this repository

Git history attributes the current tracked implementation to **Zixuan Zhu**.
The work includes:

- formulation and implementation of the mouse-scale PBPK/ODE system;
- liver, proximal-tubule and BBB/CNS intracellular modules;
- the 70 kg reference-adult regional model and route-specific inputs;
- capsid batch simulation, SINEUP-PD linkage and JSON/CSV export;
- the disease/gene library, 2D design space and anatomical React heat map;
- parameter provenance, mass-balance checks, documentation and local demo;
- iterative visualization changes made in response to internal feedback.

This statement describes repository authorship. It does not imply that the
underlying biological facts, datasets, software libraries or anatomical assets
were created by the team.

## Prior scientific work used

| External work | How it was used | What we changed or added |
|---|---|---|
| Liu et al. (2024) AAV whole-body PBPK | Mechanistic precedent for organ distribution, receptor uptake, intracellular processing and transgene output | Added kidney dual-entry, BBB/CNS, multiple administration routes, SINEUP-PD and a disease-facing frontend |
| Wang et al. (2024) radiolabelled AAV9 mouse data | Fitted early apparent organ capsid half-lives | Kept raw time points, fit windows, log-fit diagnostics and the capsid-versus-episome caveat |
| Ballon et al. (2020) NHP PET | NHP-informed early AAV9 organ priors and CSF/systemic distribution context | Used only as a labelled reference-human projection; kidney/lung gaps remain explicit |
| Zincarelli et al. (2008), Walkey et al. (2025), Yang et al. (2025), Abele et al. (2025) | Head-to-head capsid and route evidence | Built a machine-readable catalog and cautious relative priors rather than pooling incompatible assays |
| GTEx and Human Protein Atlas | Healthy-tissue gene-expression priors | Aggregated tissues into model organs and exposed provenance/limitations in the interface |
| ClinGen | Haploinsufficiency and disease-gene evidence links | Organised disease entries as expandable disease-to-gene records |
| Reactome male-body SVG | Low-opacity anatomical reference | Overlaid independently calculated model regions; source SVG geometry is unchanged |
| DBCLS human anatomy SVG | Earlier anatomical reference retained for comparison | No source-geometry change |

Full citations are in the scientific report and source URLs are stored beside
the parameter/data records.

## Software, services and tools

| Tool | Role | Attribution note |
|---|---|---|
| Python, NumPy, SciPy, Matplotlib | ODE solution, integration, fitting and plots | Open-source scientific software; versions are constrained in `requirements.txt` |
| TypeScript, React, Next/vinext, Vite | Interactive frontend | Open-source software listed in `sineup-delivery-atlas/package.json` |
| Codex / OpenAI tools | Code review, implementation assistance, literature-search assistance, document editing and local testing | AI assistance was supervised by the student author; references and numerical claims require human verification |
| Git and GitHub | Version control and public distribution | Repository history is the audit trail for code authorship |
| LaTeX / Tectonic or TeX Live | Scientific report typesetting | Source and generated PDF are both versioned |

No AlphaFold server or Figma artifact is evidenced in this repository at the
time of writing. Add them to the official form if they were used elsewhere in
the team project.

## Visual asset licences

- Reactome `Male body with organs`, stable identifier `R-ICO-013956`, curated
  by Marija Orlic-Milacic and designed by Cristoffer Sevilla, CC BY 4.0.
- DBCLS `202403 human anatomy organs.svg`, CC BY 4.0.

The full URLs and change descriptions are in
`sineup-delivery-atlas/public/ASSET_ATTRIBUTION.md`.

## Items the full team must verify before submission

The repository cannot determine these contributions. Replace each placeholder
with names, dates and a concrete description in the official Attributions Form.

| Area | Person(s) | Exact contribution | Evidence/status |
|---|---|---|---|
| Project conception and supervision | **TEAM TO VERIFY** | Who framed the AAV-SINEUP therapeutic question and approved scope? | pending |
| Wet-lab design and experiments | **TEAM TO VERIFY** | Constructs, protocols, measurements, analysis and negative results | pending |
| RNA-binder model v2/v3 | **TEAM TO VERIFY** | Architecture, training data, code, compute and interpretation | outside this repository |
| Advisor/expert feedback | **TEAM TO VERIFY** | What decision changed after each consultation? | pending meeting records |
| Wiki integration and visual design | **TEAM TO VERIFY** | Layout, illustrations, copy editing and deployment | pending |
| Institutional facilities and funding | **TEAM TO VERIFY** | Laboratory, computing, reagents, grants and sponsorship | pending |

## Statement of intellectual honesty

Model outputs are the team's calculations, but most parameters are not the
team's measurements. Every output therefore carries a model/evidence label.
Reference-human results are projections rather than clinical predictions, and
the safety screen prioritises measurements rather than declaring a dose safe.
