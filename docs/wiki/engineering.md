# Engineering: the model as a Design-Build-Test-Learn system

The engineering problem was not simply to draw AAV concentration curves. We
needed a computational system that could explain where a vector goes, identify
the step limiting expression, expose uncertainty, and return a smaller set of
capsid-route experiments to the wet lab.

## Acceptance criteria

Before comparing biological scenarios, each build had to meet four tests:

1. **Conservation:** administered AAV equals all vector-containing states plus
   explicit cumulative sinks, within numerical tolerance.
2. **Non-negativity and solver success:** no meaningful negative state and no
   failed integration.
3. **Biological separation:** early capsid disappearance, vector genome,
   episome, RNA, protein and therapeutic duration must remain different states.
4. **Reproducibility:** one command regenerates the frontend JSON/CSV, and the
   page displays those exported results without inventing new scores in the
   browser.

Passing these criteria shows internal consistency, not clinical validity.

## Cycle 1: from arbitrary decay to traceable AAV9 kinetics

**Design.** The first implementation used global vascular and ISF half-lives to
make bell-shaped curves. This made the plots readable but did not explain where
capsid went.

**Build.** We added explicit blood loss, organ vascular/RES loss, ISF catabolic
loss and cumulative sink states. Organ half-lives were fitted to Wang et al.
mouse data; the ordinary AAV9 blood half-life was corrected to 5.0 h because
2.4 h referred to a modified AAV9-TC construct.

**Test.** The mass-balance residual fell to numerical precision, and organ
curves retained distinct peaks and declines. Log-fit diagnostics exposed that
muscle was poorly represented by a single exponential.

**Learn.** A measured apparent half-life does not identify endothelial removal
and ISF catabolism separately. We therefore retained the 35:65 rate split as an
explicit structural assumption and ranked it for calibration instead of
presenting it as literature fact.

## Cycle 2: from liver-only uptake to organ- and barrier-specific transduction

**Design.** Organ ISF exposure alone could not answer why two tissues with
similar vector genomes produce different expression.

**Build.** The liver chain became `bound -> EE -> LE -> CY -> Ncap -> Nss ->
Nds -> Epi -> mRNA -> protein`. Kidney received small filtration and
basolateral entry routes; CNS received BBB binding, endothelial trafficking,
transcytosis/recycling and neural-cell uptake.

**Test.** Scenario comparisons changed capsid transport separately from
promoter output. Increasing endosomal escape raised episome more strongly than
protein once the expression module approached saturation.

**Learn.** “More organ exposure” is not equivalent to “more therapeutic
protein”. This supports experiments that measure vector, episome, RNA and
protein rather than a single endpoint. During this audit we also found that
liver transcription was capped at a hard-coded value of 2.0; it now caps at the
scenario's own `k_tx`, so promoter presets behave as declared.

## Cycle 3: from one animal diagram to a route-resolved reference adult

**Design.** Scaling every mouse state by body weight would preserve the wrong
circulation, organ fractions and CSF geometry.

**Build.** We retained the same Q-PS-Kp and intracellular equation family but
introduced explicit right heart, lung, left heart, arterial/venous pools, portal
drainage, 24 regions, 301 states and IV, IT, IM, ICM, ICV and inhaled inputs.

**Test.** Regional blood volumes sum to approximately 5 L; flow fractions close
the circulation; 150 mL CSF with 500 mL/day turnover gives a 4.99 h equivalent
half-life; each route is independently checked for mass balance.

**Learn.** Human anatomy can be mechanistic while capsid parameters remain
uncertain. The interface therefore labels the result “reference-human
projection” and removes the mouse LY6A-dependent PHP.eB CNS gain rather than
calling it a human prediction.

## Cycle 4: from static plots to a disease-facing design tool

**Design.** A wet-lab user should begin with disease and target gene, not with a
list of state-variable names.

**Build.** The exporter solves all capsids, computes organ specificity and
SINEUP-PD duration, then writes JSON/CSV. The React application provides a
disease-to-gene library, CNS-depth profiles, route/capsid selection, a 2D design
space and an anatomical heat map.

**Test.** Frontend tests build the application and verify rendered HTML. The
heat map now supports absolute and within-capsid relative scales; the latter is
the default to reveal spatial differences without claiming cross-capsid dose
equivalence.

**Learn.** Relative colour is useful for pattern recognition but can hide
absolute decline. Both scales remain available and the legend states the
normalisation. We also fixed a hydration mismatch caused by locale-dependent
server/client timestamp formatting.

## Current design decision returned to the wet lab

The model recommends a staged comparison rather than screening every possible
combination:

- choose one broadly distributed reference capsid (AAV9), one candidate with a
  target-organ advantage, and one negative/off-target control;
- compare route only when it changes the biological entry compartment;
- measure early vector/capsid and later episome, RNA and protein;
- prioritise `PS/Kp`, BBB/CSF access, endosomal escape and SINEUP gain because
  they currently drive ranking uncertainty;
- use the safety screen to add liver, platelet/complement, renal, cardiac or
  neurotoxicity endpoints appropriate to the chosen route.

## What failed or remains unfinished

- The 35:65 vascular/ISF loss split is not identifiable from existing data.
- The CNS three-depth model is reduced-order, not an anatomical brain mesh.
- Kidney filtration and receptor rates are hypotheses, not fitted intact-AAV
  measurements.
- Capsid multipliers combine qualitative evidence across incompatible studies;
  they are priors, not meta-analytic effect sizes.
- The 70 kg model has no calibrated human exposure-to-toxicity function.
- A full loop is incomplete until project-specific wet-lab data change at least
  one fitted parameter and the model re-ranks a design.

These limitations define the next Build and Test, rather than being hidden at
the end of the report.
