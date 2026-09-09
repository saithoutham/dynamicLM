# Innovation Lab Notebook

No idea is promoted on plausibility alone. Every verdict below must be backed by
an executable test and a `results/trace.csv` reference.

## IDEA-001: Correlation-aware time-scale selection

Date/commit: PENDING

Hypothesis (falsifiable): A correlation-aware summary delta-AUC test can select
the data-generating time scale with useful power while controlling type-I error.

Why it might matter: It would turn a scientific choice of time origin into a
pre-specified diagnostic using the package's own inferential framework.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: Awaiting the simulation harness and a validated left-truncation-aware
variance calculation.

## IDEA-002: Survey-weighted landmark supermodels

Date/commit: PENDING

Hypothesis (falsifiable): A design-weighted estimator with PSU/stratum-aware
variance changes inference relative to an unweighted landmark analysis in a
public linked-mortality cohort.

Why it might matter: National survey inference requires respecting unequal
selection probabilities and complex sampling.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: Requires a verified public linked-mortality download and a clearly
defined estimand; survey weights must not be bolted onto the existing IID code.

## IDEA-003: Covariate staleness adjustment

Date/commit: PENDING

Hypothesis (falsifiable): Adding time since last measurement improves summary
AUC or Brier score under informative visit timing without harming calibration.

Why it might matter: Last-observation-carried-forward quality depends on how old
the carried value is, and staleness may vary by age.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: Awaiting the PBC and NAFLD application datasets.

## IDEA-004: Inverse-intensity visit weighting

Date/commit: PENDING

Hypothesis (falsifiable): Estimated inverse visit-intensity weights reduce bias
from informative observation times more than a staleness covariate alone.

Why it might matter: Visit timing can be outcome-related and distort covariate
availability.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: This added complexity must beat the simpler staleness adjustment.

## IDEA-005: Trajectory-model evaluation harness

Date/commit: PENDING

Hypothesis (falsifiable): A model-agnostic prediction adapter can evaluate penLM
and trajectory encoders with identical landmark eligibility and summary metrics.

Why it might matter: Fair evaluation should isolate representation quality from
different cohort construction or censoring assumptions.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: Build only after the corrected estimand and data contract are stable.

## IDEA-006: Entry-risk-set shape diagnostic

Date/commit: 2026-09-08 / `49928e4`

Hypothesis (falsifiable): Reporting entries, exits, and eligible counts by
landmark detects age-scale misuse more reliably than checking for monotone risk
sets, which need not be monotone under staggered entry.

Why it might matter: The current shared-origin mental model treats increasing
risk-set counts as suspicious even when they may be valid on attained age.

Test performed: Tabulated unmodified, strict, and delayed eligible counts across
age landmarks in PBC and NAFLD; separately counted complete pre-entry rows.

Result (with numbers + trace refs): At age 50, the unmodified NAFLD stack had
9,472 rows, of which 7,170 were pre-entry and 2,976 pre-entry rows remained
analyzable [TRACE: `03_break_diagnosis::nafld_lm50_stacked_n`,
`::nafld_lm50_preentry_n`, `::nafld_lm50_analyzable_preentry_n`].

Verdict: KEEP

Reasoning: The diagnostic exposed silent misuse even though the naive risk-set
counts decreased monotonically. It belongs in every age-scale run.

## IDEA-007: Entry-conditional censoring ratios

Date/commit: 2026-09-08 / PENDING correction commit

Hypothesis (falsifiable): Dividing the left-truncated censoring survival at the
evaluation time by its value at subject entry yields weights that reduce to the
existing right-censored implementation when everyone starts at the landmark.

Why it might matter: A raw inverse `1/G(t)` grants late entrants censor-free
time before they could be observed.

Test performed: Evaluated a hand-checkable staggered-entry risk set and compared
the strict-entry special case directly with development `riskRegression::Score()`.

Result (with numbers + trace refs): Early and late entrants had conditional
inverse horizon weights 2 and 1; at the strict PBC validation landmark, AUC and
Brier differences from `riskRegression` were 1.1102230e-16 and 2.7755576e-17
[TRACE: `04_correction::toy_inverse_weight_horizon_early_entry`,
`::toy_inverse_weight_horizon_late_entry`,
`::strict_validation_max_auc_difference`,
`::strict_validation_max_brier_difference`].

Verdict: KEEP

Reasoning: The finite-sample implementation and limiting special case are
verified. Consistency and its analytic influence function remain unproved.

## IDEA-008: Within-window delayed-entry augmentation

Date/commit: 2026-09-08 / PENDING correction commit

Hypothesis (falsifiable): Admitting subjects at `max(landmark, entry)` improves
precision without materially increasing coefficient bias relative to strict
landmarking.

Why it might matter: Strict eligibility discards people entering shortly after
a landmark, which can be severe with coarse age grids.

Test performed: Fitted corrected strict and delayed-entry supermodels to PBC and
NAFLD.

Result (with numbers + trace refs): Delayed-entry robust SEs were smaller for
the reported PBC bilirubin and NAFLD BMI coefficients [TRACE:
`04_correction::pbcseq_death_strict_log_bili_z_robust_se`,
`::pbcseq_death_delayed_log_bili_z_robust_se`,
`::nafld_strict_bmi_z_robust_se`, `::nafld_delayed_bmi_z_robust_se`].

Verdict: PARK

Reasoning: The populations differ, so public-data SE reductions do not establish
better efficiency at a common estimand. Known-truth simulation is required.

## IDEA-009: Raw attained-age polynomial in the linear predictor

Date/commit: 2026-09-08 / `49928e4`

Hypothesis (falsifiable): Existing raw `LM1` terms can be reused unchanged when
landmarks are attained ages.

Why it might matter: Reuse would minimize code and tutorial changes.

Test performed: Fitted unmodified PBC and NAFLD age-scale models with raw `LM1`.

Result (with numbers + trace refs): Every fitted cause/model returned one
non-finite coefficient, identified as `LM1` in
`results/age_model_diagnostics.csv` [TRACE:
`03_break_diagnosis::pbcseq_cause1_nonfinite_coefficients`,
`::pbcseq_cause2_nonfinite_coefficients`,
`::nafld_cause1_nonfinite_coefficients`].

Verdict: KILL

Reasoning: The raw term was degenerate in both stress tests. Centered/scaled
interaction terms are retained; the unidentified landmark main term is not.

## IDEA-010: Missing covariates as an implicit entry filter

Date/commit: 2026-09-08 / `49928e4`

Hypothesis (falsifiable): Long-format pre-entry rows will always receive missing
time-varying values and therefore cannot contaminate a fitted model.

Why it might matter: If true, explicit entry filtering would be less urgent.

Test performed: Counted complete pre-entry rows after the unmodified stacker in
PBC and NAFLD.

Result (with numbers + trace refs): PBC had zero analyzable pre-entry rows at
the tested landmarks, but NAFLD had 2,240 at age 40 and 2,976 at age 50 [TRACE:
`03_break_diagnosis::pbcseq_lm50_analyzable_preentry_n`,
`::nafld_lm40_analyzable_preentry_n`,
`::nafld_lm50_analyzable_preentry_n`].

Verdict: KILL

Reasoning: NA masking is dataset-dependent and historical measurements defeat
it. Eligibility must be explicit.
