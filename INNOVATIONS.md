# Innovation Lab Notebook

No idea is promoted on plausibility alone. Every verdict below must be backed by
an executable test and a `results/trace.csv` reference.

## IDEA-001: Correlation-aware time-scale selection

Date/commit: 2026-09-08 / PENDING simulation commit

Hypothesis (falsifiable): A correlation-aware summary delta-AUC test can select
the data-generating time scale with useful power while controlling type-I error.

Why it might matter: It would turn a scientific choice of time origin into a
pre-specified diagnostic using the package's own inferential framework.

Test performed: Compared naive and strict summary AUC with subject-paired
bootstrap resamples across landmarks in 1,000 independently seeded datasets per
entry-spread condition.

Result (with numbers + trace refs): At entry-age SD 4 and 10, rejection
probabilities were 0.045 (MCSE 0.00656) and 0.046 (MCSE 0.00662), respectively;
mean delta AUCs were 0.000273 and 0.000425 [TRACE:
`05_simulation::diagnostic_entry_sd_4_rejection_probability`,
`::diagnostic_entry_sd_4_rejection_mcse`,
`::diagnostic_entry_sd_10_rejection_probability`,
`::diagnostic_entry_sd_10_rejection_mcse`,
`::diagnostic_entry_sd_4_mean_delta_auc`,
`::diagnostic_entry_sd_10_mean_delta_auc`].

Verdict: KILL

Reasoning: In this prespecified independent-entry DGP, the wrong eligibility
rule barely changed discrimination. The test behaved like a level test, not a
diagnostic with power to identify the wrong time origin. Delta AUC must not be
used as a general time-scale selector.

## IDEA-002: Survey-weighted landmark supermodels

Date/commit: 2026-09-08 / PENDING application commit

Hypothesis (falsifiable): A design-weighted estimator with PSU/stratum-aware
variance changes inference relative to an unweighted landmark analysis in a
public linked-mortality cohort.

Why it might matter: National survey inference requires respecting unequal
selection probabilities and complex sampling.

Test performed: Linked the public-use NHANES 1999--2000 examination file to its
public-use mortality file, constructed strict and delayed age-scale stacks, and
compared subject-clustered `coxph()` with a nested PSU/participant
`survey::svycoxph()` using MEC examination weights and survey strata.

Result (with numbers + trace refs): In the strict stack the unweighted male
coefficient was 0.525765 (SE 0.086943), versus 0.481217 (SE 0.075297) under the
survey design [TRACE:
`06_application::nhanes_1999_2000_strict_coxph_cluster_1_male_estimate`,
`::nhanes_1999_2000_strict_coxph_cluster_1_male_standard_error`,
`::nhanes_1999_2000_strict_svycoxph_1_male_estimate`,
`::nhanes_1999_2000_strict_svycoxph_1_male_standard_error`].

Verdict: PARK

Reasoning: The machinery runs and weighting changes the empirical estimate, but
this is not a proof that ordinary survey Cox variance remains valid for an
overlapping landmark stack. A design-based summary-AUC/Brier estimand and its
influence function are still PENDING, so the idea is not promoted.

## IDEA-003: Covariate staleness adjustment

Date/commit: 2026-09-08 / PENDING application commit

Hypothesis (falsifiable): Adding time since last measurement improves summary
AUC or Brier score under informative visit timing without harming calibration.

Why it might matter: Last-observation-carried-forward quality depends on how old
the carried value is, and staleness may vary by age.

Test performed: Added time since last measurement to five-fold
subject-disjoint landmark Cox fits in NAFLD and PBC, then compared out-of-fold
summary AUC with paired subject bootstraps across landmarks.

Result (with numbers + trace refs): Delta AUC was 0.003890 (bootstrap SE
0.012467; p=0.7550) in NAFLD and -0.002384 (SE 0.013043; p=0.8550) in PBC
[TRACE: `06_application::nafld_staleness_delta_auc`,
`::nafld_staleness_bootstrap_se`, `::nafld_staleness_p_value`,
`::pbcseq_staleness_delta_auc`, `::pbcseq_staleness_bootstrap_se`,
`::pbcseq_staleness_p_value`].

Verdict: KILL

Reasoning: The simple additive staleness term showed no useful discrimination
gain in either public-data test. This kills the proposed off-the-shelf fix, not
the broader possibility that visit timing contains outcome information.

## IDEA-004: Inverse-intensity visit weighting

Date/commit: 2026-09-08 / PENDING application commit

Hypothesis (falsifiable): Estimated inverse visit-intensity weights reduce bias
from informative observation times more than a staleness covariate alone.

Why it might matter: Visit timing can be outcome-related and distort covariate
availability.

Test performed: Estimated a Poisson visit-rate approximation in each training
fold from sex, entry age, and landmark, capped inverse predicted-rate weights at
the training-fold 99th percentile, normalized them to mean one, and evaluated
weighted landmark Cox fits out of fold.

Result (with numbers + trace refs): Delta AUC was -0.022056 (bootstrap SE
0.012824; p=0.08545) in NAFLD and 0.001217 (SE 0.002319; p=0.5997) in PBC
[TRACE: `06_application::nafld_iiw_delta_auc`,
`::nafld_iiw_bootstrap_se`, `::nafld_iiw_p_value`,
`::pbcseq_iiw_delta_auc`, `::pbcseq_iiw_bootstrap_se`,
`::pbcseq_iiw_p_value`].

Verdict: KILL

Reasoning: The approximation did not improve PBC and moved NAFLD AUC in the
wrong direction. A properly specified recurrent-event intensity model could be
revisited, but this landmark-row approximation should not ship.

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

## IDEA-011: Paired subject bootstrap as an immediate IF substitute

Date/commit: 2026-09-08 / PENDING simulation commit

Hypothesis (falsifiable): A subject bootstrap shared across landmarks provides
near-nominal summary-AUC coverage over the entry-heterogeneity stress test.

Why it might matter: It would permit honest inference while the analytic
left-truncation influence function remains open.

Test performed: Used 100 paired bootstrap resamples in each of 1,000 outer
datasets at each of three entry-age spreads, for naive and strict eligibility.

Result (with numbers + trace refs): Strict coverage was 0.932, 0.940, and 0.954
as entry-age SD increased; corresponding MCSEs were 0.00796, 0.00751, and
0.00662 [TRACE: `05_simulation::auc_entry_sd_0p25_strict_coverage`,
`::auc_entry_sd_4_strict_coverage`, `::auc_entry_sd_10_strict_coverage`,
`::auc_entry_sd_0p25_strict_coverage_mcse`,
`::auc_entry_sd_4_strict_coverage_mcse`,
`::auc_entry_sd_10_strict_coverage_mcse`].

Verdict: PARK

Reasoning: Coverage was not uniformly nominal and was low in the near-shared
condition. The bootstrap is exposed as empirical inference, not endorsed as a
drop-in theoretical correction.

## IDEA-012: Validated scalar Cox engine for large method grids

Date/commit: 2026-09-08 / PENDING simulation commit

Hypothesis (falsifiable): A direct Breslow score/Newton solver with cluster
sandwich variance can reproduce `coxph()` closely enough to run the full grid.

Why it might matter: Repeated general-purpose formula parsing would dominate
the required 324,000 fits.

Test performed: Compared coefficient and robust SE outputs with `coxph()` for
naive, strict, and delayed stacks, then required every Newton score residual to
be near zero.

Result (with numbers + trace refs): Maximum coefficient and robust-SE
differences were 1.1709017e-09 and 6.8008849e-11; all 324,000 grid fits
converged [TRACE: `05_simulation::engine_validation_max_estimate_difference`,
`::engine_validation_max_robust_se_difference`,
`::total_model_fits`, `::fit_failures`].

Verdict: KEEP

Reasoning: It matched the reference to a tolerance far below the Monte Carlo
precision and made the mandated grid feasible. It remains deliberately limited
to the one-covariate validation DGP.
