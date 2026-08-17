# AAV–SINEUP safety-margin assessment (2026-08-13)

## Executive conclusion

The current model does **not** prove that every organ is safe.  It supports a narrower, defensible statement for the default **IV AAV9** scenario: the administered dose of `4.00e+13 vg/kg` is below the marketed ZOLGENSMA dose and below the lowest explicitly reported mouse cardiac-histology signal dose.  Re-solving the same PBPK equations gives organ ISF-AUC margins of `2.75-2.75` relative to the mouse signal-dose context.  This is supportive dose justification, not a human NOAEL.

Of the eight modeled organ groups, 8/8 have some route/capsid-matched contextual support (heart, liver, brain/DRG context, kidney and spleen); 0/8 (lung, muscle and rest-of-body) lack organ-specific toxic exposure thresholds.  Therefore the correct claim is **“lower modeled exposure than selected AAV9 signal/clinical contexts, with residual product-specific risk”**, not “all organs proven safe.”

The current lumbar intrathecal scenario is not supportable at its present dose.  The model applies the same adult total dose `2.80e+15 vg` to CSF delivery, which is `23.3x` the marketed ITVISMA dose (`1.2e14 vg`).  Re-solving at the ITVISMA dose gives reference/current organ-AUC margins of only `0.043-0.043`; equivalently, current modeled exposure is about `23.3-23.3x` the marketed-product context.

## Calculation

For organ `i`, route `r` and capsid `c`:

`Margin_i = AUC_i(reference AAV9, route-matched reference dose) / AUC_i(current route, current capsid)`

- `Margin >= 1` means the current modeled ISF amount AUC is below that reference exposure.
- `Margin < 1` means the current modeled exposure exceeds it.
- A margin is **not** a probability of safety. It inherits all comparator limitations: product, payload, capsid, manufacturing, age, disease, immunity, full/empty ratio and species.

For IV AAV9, the model is nearly dose proportional over this range: the organ margins to ZOLGENSMA's `1.1e14 vg/kg` dose are approximately `2.75`, and margins to the `7.9e13 vg/kg` mouse cardiac-histology signal dose are approximately `1.98`.

## CNS AAV case studies

| Case | Delivery and dose | What it supports | What it does not support |
|---|---|---|---|
| ZOLGENSMA | Single 60-min IV AAV9, `1.1e14 vg/kg` | Systemic AAV9 clinical dose context and monitoring domains | Universal AAV safe dose; label includes fatal/nonfatal liver failure, TMA, thrombocytopenia and troponin signals |
| ITVISMA | Single lumbar IT AAV9, `1.2e14 vg` in 3 mL over about 1-2 min | Best current same-route/same-capsid CNS comparator | No-risk threshold; trial and label include liver enzyme, sensory/DRG-like, platelet, TMA and troponin signals |
| KEBILIDI | AAV2, four stereotactic putaminal infusions, total `1.8e11 vg` in 0.32 mL | Approved local-parenchymal delivery precedent | Comparator for CSF or IV exposure; current model lacks an intraputaminal compartment |
| AAV9-miniSINEUP-GDNF | Mouse striatum, `7.0e9 vg` | Direct SINEUP payload precedent: about two-fold endogenous-protein increase for >=6 months with no reported weight/food-intake signal | Human safety, GLP toxicology or safety of another binding domain/target |

The current model's `2.8e15 vg` local-CNS total dose is also about `15,556x` the KEBILIDI total dose, but this ratio is descriptive only because route, capsid and distribution geometry are different.

## RNA-drug delivery precedents

- SPINRAZA: lumbar intrathecal antisense oligonucleotide; current 2026 label includes `12 mg/5 mL` low-dose injections and a higher-dose regimen. It provides procedure, repeat-dose and monitoring precedent, not an AAV vg conversion.
- QALSODY: lumbar intrathecal antisense oligonucleotide, `100 mg/15 mL`, three loading doses then every 28 days. Myelitis/radiculitis, papilledema/intracranial pressure and aseptic meningitis show that route-specific safety cannot be inferred from AAV biodistribution alone.

## Payload-specific conclusion

The model caps haploinsufficient protein restoration at 100% from a 50% baseline, i.e. a maximum two-fold restoration.  This is mechanism-consistent with reported SINEUP increases (often about 1.5-3-fold) and with the approximately two-fold AAV9-miniSINEUP-GDNF mouse result.  It is not proof that the current plasmid is safe because the repository does not specify one final binding-domain sequence, target transcript/isoform, promoter, cassette sequence, CpG burden, ITR integrity or human off-target screen.

The plasmid claim should therefore be: **“designed for endogenous-transcript-dependent, bounded restoration; existing SINEUP data reduce but do not eliminate overexpression concern.”**  Required construct-specific checks are target-protein upper tolerance, RNA-seq/proteomics off-targets, innate RNA sensing, promoter leakiness, replication-competent AAV, full/empty ratio and long-term integration analysis.

## Recommended safety-model architecture

1. Keep capsid particles, vector genomes, episomes, SINEUP RNA and target protein as separate states.
2. Fit organ vector exposure with PBPK; do not map AUC directly to injury without data.
3. Use a Bayesian interval-censored logistic model for binary histopathology: `logit(P(injury)) = alpha_study + beta*log10(exposure) + route + species + capsid + payload-expression`.
4. Model longitudinal ALT/AST, platelets, creatinine, complement, troponin and neurofilament as continuous submodels; connect them to time-to-event hazards for clinical events.
5. Treat DRG separately from whole brain. The current brain compartment cannot validate DRG safety.
6. Propagate uncertainty and report `P(margin>1)` plus credible intervals, not a single green/red score.
7. Validate on a held-out dose group or second study before human extrapolation.

## Current decision

- **Proceed as a research hypothesis:** IV AAV9 at `4.0e13 vg/kg`, with an explicit ~2-fold margin to the product-specific mouse histology signal context and full liver/TMA/cardiac/DRG monitoring.
- **Do not claim clinical safety:** the evidence is not a product-specific NOAEL and recent marketed-product experience confirms serious AAV liver risk.
- **Do not use the current CSF dose:** reduce and re-solve lumbar IT around `1.2e14 vg total` or lower before interpreting CNS safety; ICM/ICV need route-specific large-animal calibration.
- **Do not claim the plasmid itself safe yet:** sequence- and target-specific evidence is missing.

## Primary sources

- ZOLGENSMA US prescribing information: https://www.fda.gov/media/126109/download
- ITVISMA US prescribing information: https://www.fda.gov/media/193168/download
- KEBILIDI US prescribing information: https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=9d6a6401-c6b5-4f29-af11-67707d249482
- SPINRAZA 2026 US prescribing information: https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/209531s016lbl.pdf
- QALSODY US prescribing information: https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/215887s000lbl.pdf
- AAV9-miniSINEUP-GDNF mouse study: https://pmc.ncbi.nlm.nih.gov/articles/PMC7000958/
- Intrathecal onasemnogene NHP DRG study: https://pmc.ncbi.nlm.nih.gov/articles/PMC9347375/
- AAV8/AAV9 PBPK precedent: https://doi.org/10.1016/j.xphs.2023.10.005
- Antibody/complement-associated AAV TMA study: https://pmc.ncbi.nlm.nih.gov/articles/PMC10760971/
- FDA 2025 ELEVIDYS liver-failure safety action: https://www.fda.gov/news-events/press-announcements/fda-approves-new-safety-warning-and-revised-indication-limits-use-elevidys-following-reports-fatal
