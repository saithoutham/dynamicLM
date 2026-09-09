# Innovation Lab Notebook

No idea is promoted on plausibility alone. Every verdict below must be backed by
an executable test and a `results/trace.csv` reference.

## IDEA-001: Correlation-aware time-scale selection

Date/commit: 2026-09-09 / `700fcc2`

Hypothesis (falsifiable): A correlation-aware summary delta-AUC test can select
the data-generating time scale with useful power while controlling type-I error.

Why it might matter: It would turn a scientific choice of time origin into a
pre-specified diagnostic using the package's own inferential framework.

Test performed: Compared naive and strict summary AUC with subject-paired
bootstrap resamples across landmarks in 1,000 independently seeded datasets per
entry-spread condition. Phase 6 then re-audited strict-target AUC bias under
independent, measured covariate-dependent, and shared-frailty entry using the
oracle generating score.

Result (with numbers + trace refs): At entry-age SD 4 and 10, rejection
probabilities were 0.045 (MCSE 0.00656) and 0.046 (MCSE 0.00662), respectively;
mean delta AUCs were 0.000273 and 0.000425 [TRACE:
`05_simulation::diagnostic_entry_sd_4_rejection_probability`,
`::diagnostic_entry_sd_4_rejection_mcse`,
`::diagnostic_entry_sd_10_rejection_probability`,
`::diagnostic_entry_sd_10_rejection_mcse`,
`::diagnostic_entry_sd_4_mean_delta_auc`,
`::diagnostic_entry_sd_10_mean_delta_auc`]. Under wide measured entry at gamma
`-4`, naive and strict AUC biases were 0.020885 and -0.000816, whereas the
original independent-entry test remained null [TRACE:
`07_informative_entry::auc_cell_003_naive_bias`,
`::auc_cell_003_strict_bias`].

Verdict: KILL

Reasoning: The selector can separate methods in a strong measured-entry
mechanism, but that does not make it general: the same eligibility defect was
invisible under independent entry, frailty can preserve oracle ranking while
breaking coefficient estimation, and delayed entry changes the target
population. Delta AUC must not choose a scientific time origin. Direct
eligibility auditing is retained instead.

## IDEA-002: Survey-weighted landmark supermodels

Date/commit: 2026-09-08 / `67eb695`

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

Date/commit: 2026-09-08 / `67eb695`

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

Date/commit: 2026-09-08 / `67eb695`

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

Date/commit: 2026-09-09 / `700fcc2`

Hypothesis (falsifiable): Reporting entries, exits, and eligible counts by
landmark detects age-scale misuse more reliably than checking for monotone risk
sets, which need not be monotone under staggered entry.

Why it might matter: The current shared-origin mental model treats increasing
risk-set counts as suspicious even when they may be valid on attained age.

Test performed: Tabulated unmodified, strict, and delayed eligible counts across
age landmarks in PBC and NAFLD; separately counted complete pre-entry rows. In
Phase 6, repeated the pre-entry fraction, entry--covariate correlation, and
naive-versus-pre-entry covariate-distribution tests in 1,000 independently
seeded datasets for each independent, covariate-dependent, and frailty-dependent
entry mechanism.

Result (with numbers + trace refs): At age 50, the unmodified NAFLD stack had
9,472 rows, of which 7,170 were pre-entry and 2,976 pre-entry rows remained
analyzable [TRACE: `03_break_diagnosis::nafld_lm50_stacked_n`,
`::nafld_lm50_preentry_n`, `::nafld_lm50_analyzable_preentry_n`]. In Phase 6,
structural detection was 1.000 for the wide-entry gamma-zero control as well as
the strongest measured informative-entry condition. At wide entry spread, the
observed-covariate correlation test detected gamma `-4` with probability 1.000
but detected shared frailty entry with probability only 0.065 [TRACE:
`07_informative_entry::diagnostic_cell_012_lm1_structural_detection_probability`,
`::diagnostic_cell_003_lm1_structural_detection_probability`,
`::diagnostic_cell_003_lm1_entry_x_detection_probability`,
`::diagnostic_cell_021_lm1_entry_x_detection_probability`].

Verdict: KEEP

Reasoning: The diagnostic exposed silent misuse even though the naive risk-set
counts decreased monotonically. It belongs in every age-scale run as a
**structural eligibility audit**. It does not identify whether entry is
informative: pre-entry rows also occur under independent entry, observed
covariates miss unmeasured frailty, and conditioning on survival to entry can
induce entry--covariate association.

## IDEA-007: Entry-conditional censoring ratios

Date/commit: 2026-09-08 / `eb7ab56`

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

Date/commit: 2026-09-08 / `eb7ab56`

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

Date/commit: 2026-09-08 / `67eb695`

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

Date/commit: 2026-09-08 / `67eb695`

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

## IDEA-013: Observed entry--covariate correlation as a general screen

Date/commit: 2026-09-09 / `700fcc2`

Hypothesis (falsifiable): Testing entry age against modeled covariates detects
informative entry while controlling false positives under independent entry.

Why it might matter: Unlike prediction performance, it interrogates the entry
mechanism directly and is inexpensive to run before model fitting.

Test performed: In each of 1,000 seeded datasets per mechanism and entry spread,
tested the observed entry--`X` correlation and compared the disjoint pre-entry
and strict risk-set `X` distributions at every landmark.

Result (with numbers + trace refs): At entry SD 10, the correlation screen's
detection probability was 0.156 even for gamma `0`, 1.000 for gamma `-4`, and
0.065 for shared-frailty entry (`theta = 0.5`, `delta = -2`). The corresponding
frailty risk-set-shift detection probability was 0.049 [TRACE:
`07_informative_entry::diagnostic_cell_012_lm1_entry_x_detection_probability`,
`::diagnostic_cell_003_lm1_entry_x_detection_probability`,
`::diagnostic_cell_021_lm1_entry_x_detection_probability`,
`::diagnostic_cell_021_lm1_any_riskset_shift_detection_probability`].

Verdict: KILL

Reasoning: It is useful evidence when it fires strongly for measured entry, but
it is neither a calibrated general test after survival-to-entry selection nor
sensitive to unmeasured causes. It must not be sold as an informative-entry
diagnostic.

## IDEA-014: Eligibility correction under measured informative entry

Date/commit: 2026-09-09 / `700fcc2`

Hypothesis (falsifiable): When entry depends on modeled `X` but is conditionally
independent of failure, strict and delayed eligibility remove the naive
coefficient bias.

Why it might matter: This is the mechanism in which an EHR analyst can measure
and adjust the variable driving differential presentation.

Test performed: Swept gamma over `0`, `-1`, `-2`, and `-4`, then confirmed the
largest naive-versus-strict separation over the complete 108-cell Phase 4 grid,
with 1,000 replicates and three methods in every cell.

Result (with numbers + trace refs): At gamma `-4`, pooled full-grid bias was
0.043250 for naive eligibility, 0.001762 for strict eligibility, and 0.002134
for delayed eligibility. Corresponding coverage was 0.92494, 0.94370, and
0.94464 [TRACE: `07_informative_entry::pooled_confirmation_naive_bias`,
`::pooled_confirmation_naive_coverage`,
`::pooled_confirmation_strict_bias`,
`::pooled_confirmation_strict_coverage`,
`::pooled_confirmation_delayed_bias`, and
`::pooled_confirmation_delayed_coverage`].

Verdict: KEEP

Reasoning: The correction removed most measured-entry bias throughout a large
known-truth grid. This is empirical evidence, not a general consistency proof,
and it does not license claims under unmeasured dependent truncation.

## IDEA-015: Eligibility correction is sufficient under shared frailty

Date/commit: 2026-09-09 / `700fcc2`

Hypothesis (falsifiable): Correct risk-set eligibility alone restores unbiased
coefficient estimation when an unmeasured frailty drives both entry and failure.

Why it might matter: If true, the code patch would be enough for EHR cohorts
even when presentation depends on latent health.

Test performed: Crossed `theta` in `{0, 0.5}` and `delta` in `{0, -2}` over all
three entry spreads, with 1,000 replicates and naive, strict, and delayed fits.

Result (with numbers + trace refs): With `theta = 0.5` and `delta = -2`, pooled
bias was -0.054841, -0.052436, and -0.053453 for naive, strict, and delayed
fits; coverage was 0.9187, 0.9197, and 0.9207 [TRACE:
`07_informative_entry::pooled_frailty_theta_0p5_delta_m2_naive_bias`,
`::pooled_frailty_theta_0p5_delta_m2_naive_coverage`,
`::pooled_frailty_theta_0p5_delta_m2_strict_bias`,
`::pooled_frailty_theta_0p5_delta_m2_strict_coverage`,
`::pooled_frailty_theta_0p5_delta_m2_delayed_bias`, and
`::pooled_frailty_theta_0p5_delta_m2_delayed_coverage`].

> **Superseded verdict and reasoning (2026-09-09):** `KILL`. “All methods
> failed similarly. Because `theta > 0` also makes the one-covariate model a
> marginal frailty mixture, this experiment cannot isolate dependent truncation
> from omitted-variable model misspecification. It does establish that
> eligibility correction alone is not sufficient for this stress test.”

Re-analysis (2026-09-09): The factorial design does isolate the `delta`
increment within the simulation. For all three methods, both `delta`-contrast
Monte Carlo intervals and the `theta`-by-`delta` interaction interval included
zero; both `theta` main-effect intervals excluded zero [TRACE:
`08_frailty_decomposition::{method}_{contrast}_{contrast_estimate,ci_lower,ci_upper}`].

Verdict: KILL AS A TEST OF CORRECTION SUFFICIENCY

Reasoning: The hypothesis asked eligibility correction to recover a
one-covariate conditional coefficient after an omitted frailty was added to the
hazard. The re-analysis attributes the measured bias to that `theta` change,
not to the added entry dependence through `delta`. The original hypothesis was
therefore a misspecified test of the patch rather than evidence that the patch
failed under dependent truncation. This does not prove adequacy under other
dependent-entry mechanisms.

## IDEA-016: Manifest-regenerated raw result artifacts

Date/commit: 2026-09-09 / `PENDING`

Hypothesis (falsifiable): The two tracked result files above 25 MB can be
removed without losing a reproducible inferential artifact if a clean
manifest-driven replay recreates each file byte-identically.

Why it might matter: Keeping the files accepts a GitHub size warning and makes
every clone carry reproducible rows. Gzip would reduce transfer size but retain
a large generated binary and require path changes for existing consumers.
Verified regeneration would keep seed-level provenance while avoiding both.

Test performed: Evaluated all three options. Option 1 retained the existing
warning. Option 2 was feasible but would leave generated artifacts under
version control. For option 3, rebuilt both raw tables from the committed seed
manifests with the production generator and fitting engine, wrote temporary
CSVs, and compared the complete byte vectors with `identical()` before any
tracked file was removed.

Result (with numbers + trace refs): The regenerated Phase 4 table contained
324,000 rows and 35,068,735 bytes; the informative-entry table contained
396,000 rows and 77,294,318 bytes. Both byte comparisons returned one [TRACE:
`08_hygiene::simulation_raw_rows`, `::simulation_raw_bytes`,
`::simulation_raw_byte_identical`, `::informative_entry_raw_rows`,
`::informative_entry_raw_bytes`,
`::informative_entry_raw_byte_identical`]. Full MD5 hashes and the measured
807.106-second replay are in
`results/check/raw_regeneration_verification.csv` [TRACE:
`08_hygiene::regeneration_runtime_seconds`].

Verdict: KEEP OPTION 3

Reasoning: The prerequisite exact test passed for both oversized files. They
are now ignored working artifacts regenerated by
`analysis/regenerate_raw.R --write`; the committed verification record and
seed manifests remain. The removed bytes remain recoverable from prior Git
history as well as by deterministic replay. This decision does not rewrite
repository history, so existing clones may still contain the old blob.
