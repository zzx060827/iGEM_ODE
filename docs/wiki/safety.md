# AAV safety research: current evidence and model scope

This is an early hazard-identification framework. It does not demonstrate that
the proposed AAV-SINEUP product is safe, and it must not be used to select a
clinical dose.

## Why safety is product-specific

Risk depends on total capsid and vector-genome dose, full/empty ratio, serotype,
promoter, payload, impurities, route, infusion procedure, age, disease,
pre-existing immunity and concomitant infection or treatment. A dose used by
one licensed product is therefore context, not a universal threshold.

The current 70 kg simulation uses `4.0e13 vg/kg`. For comparison only, the US
ZOLGENSMA label specifies `1.1e14 vg/kg` as a single 60-minute IV infusion in a
defined paediatric SMA population, with corticosteroid prophylaxis and intensive
monitoring. The difference does not prove a safety margin for our construct.

## Clinically important AAV risk domains

| Risk | Evidence context | What should be monitored or modelled |
|---|---|---|
| Hepatotoxicity | A prominent risk after systemic AAV; serious liver injury and failure are included in product warnings | Liver capsid/transgene exposure, ALT, AST, bilirubin, albumin, PT/INR, anti-capsid T cells |
| Thrombocytopenia and TMA | Reported after systemic AAV9; complement activation may accompany microvascular injury | Platelets, hemoglobin, creatinine, urinalysis, LDH/haptoglobin and complement markers |
| Innate/adaptive immunity | Pre-existing NAb can reduce efficacy; capsid and transgene responses can alter safety and persistence | NAb titre, antibody/T-cell response, cytokines, complement, infection status |
| Cardiac signal | Troponin elevation is monitored for onasemnogene abeparvovec | Heart exposure, troponin and cardiac assessment |
| DRG/CNS neurotoxicity | DRG pathology has occurred in nonclinical/clinical CNS-directed AAV programmes; local brain MRI findings are route/procedure dependent | CSF/DRG exposure, neurological/sensory endpoints, MRI, CSF biomarkers and relevant-animal histopathology |
| Unwanted expression | Broad tropism or strong promoter can create on-target/off-tissue effects; excess restoration may also be harmful | Organ-specific promoter activity, target protein dose-response and off-target tissues |
| Integration/tumorigenicity | AAV is mainly episomal but integration risk is not zero | Long-term follow-up and product-specific integration assessment |
| CMC and procedure | Empty particles, residual host-cell material, endotoxin and delivery devices can change risk | Identity, purity, potency, full:empty ratio, endotoxin, sterility and device compatibility |

Primary regulatory context: the FDA's
[AAV clinical-development toxicity examples](https://www.fda.gov/media/167536/download),
[guidance for human gene therapy in neurodegenerative diseases](https://www.fda.gov/media/144886/download),
and the current US
[ZOLGENSMA prescribing information](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=68cd4f06-70e1-40d8-bedb-609ec0afa471).

## What the model now does

`model/export_safety_screen.py` compares each route/capsid organ AUC with IV
AAV9 in the same reference-adult model. It maps relative liver, kidney/spleen,
heart, CNS and systemic exposure to assay priorities. CSF routes automatically
trigger a route-specific CNS/DRG flag.

The screen is useful for questions such as “which organ should receive extra
toxicology measurements?” It does not contain a calibrated relationship from
AUC to ALT elevation, TMA, troponin or neuronal injury, so it cannot answer
“what dose is safe?”.

## Proposed next safety model

The next defensible layer is an exposure-response model fitted to a single
well-characterised vector programme:

1. retain capsid, vector genome, transgene and immune analytes separately;
2. use a saturable or sigmoid dose-exposure model where supported;
3. link liver AUC to longitudinal ALT/AST and immune markers;
4. link platelet/creatinine/complement trajectories to a TMA hazard endpoint;
5. treat DRG/CNS injury as a route- and species-conditioned probability;
6. propagate parameter uncertainty and report probability intervals;
7. validate on held-out dose groups or a second study.

Recent mouse dose-ranging work indicates non-proportional tissue exposure,
sigmoidal transgene expression and greater immune/hepatic signals at high dose.
That supports a nonlinear framework, but mouse AAV8 data cannot directly set a
human AAV9-SINEUP toxicity threshold.

## Risk-reduction decisions available now

- minimise total capsid needed for the target protein window;
- compare local/CSF routes with systemic exposure rather than assuming “local”
  means no peripheral distribution;
- use tissue/cell-selective regulatory elements and evaluate excess expression;
- screen pre-existing immunity and define route-appropriate monitoring;
- characterise full/empty ratio, potency and impurities before interpreting
  dose-response data;
- avoid claims about repeat administration until immune and product-specific
  evidence exists.
