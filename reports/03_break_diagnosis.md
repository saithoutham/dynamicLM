# Phase 2 — naive age port and break diagnosis

## Bottom line

The left-truncation hypothesis is verified against both the package source and two
public-data stress tests. The unmodified stacker admits rows satisfying
`landmark < exit_age` without checking `entry_age <= landmark`. This has two
different observable consequences:

1. In `pbcseq`, the long-format covariate lookup masks the error accidentally:
   every pre-entry row has a missing time-varying covariate and is dropped during
   fitting. The fit still returns an object with a non-finite `LM1` coefficient
   for each cause, and prediction then stops with the exact error
   `all rows of newdata have missing values`.
2. In `nafld`, laboratory records dated before the cohort index date make many
   pre-entry rows complete. Those rows survive model fitting, so this is the more
   dangerous silent-wrong-answer case.

The public-data coefficient comparisons below are **discrepancies, not bias
estimates**: these data have no known truth, and strict landmarking and within-window
delayed entry target different populations. Bias is reserved for the known-truth
simulation in Phase 4.

## Construction and asserted eligibility

For observation interval `[E_i, X_i)`, the eligible strict landmark set is
`E_i <= a < X_i`. The unmodified package checks only `a < X_i`; see the quoted
source in `01_recon.md`. Both datasets were passed through `stack_data()`,
`add_interactions()`, and `dynamic_lm()` without age-entry modifications.

The prediction window was ten years. The PBC landmarks were ages 40, 50, 60,
and 70; the NAFLD landmarks were ages 40, 50, 60, 70, and 80. These parameters,
and every resulting count, are recorded by
`analysis/03_break_diagnosis.R` under the `03_break_diagnosis` trace namespace.
[TRACE: `03_break_diagnosis::prediction_window_years`,
`::pbc_landmark_1`–`::pbc_landmark_4`, and
`::nafld_landmark_1`–`::nafld_landmark_5`.]

The complete diagnostic table is
[`age_naive_stack_diagnostics.csv`](../results/age_naive_stack_diagnostics.csv).
Selected rows show the failure clearly:

| data / age | unmodified stack | pre-entry | strict eligible | analyzable pre-entry |
|---|---:|---:|---:|---:|
| PBC / 50 | 221 | 154 | 67 | 0 |
| NAFLD / 50 | 9,472 | 7,170 | 2,302 | 2,976 |
| NAFLD / 70 | 3,262 | 1,693 | 1,569 | 1,209 |

Trace refs: `03_break_diagnosis::pbcseq_lm50_stacked_n`,
`::pbcseq_lm50_preentry_n`, `::pbcseq_lm50_strict_n`,
`::pbcseq_lm50_analyzable_preentry_n`,
`::nafld_lm50_stacked_n`, `::nafld_lm50_preentry_n`,
`::nafld_lm50_strict_n`, `::nafld_lm50_analyzable_preentry_n`,
`::nafld_lm70_stacked_n`, `::nafld_lm70_preentry_n`,
`::nafld_lm70_strict_n`, and `::nafld_lm70_analyzable_preentry_n`.

This also refutes a useful but unsafe diagnostic: the unmodified stack counts
decline with age, even though the risk-set membership is wrong. On an age scale,
a correct strict risk set can first grow as people enter and later shrink as they
exit. Monotonicity of the naive counts is therefore falsely reassuring.

## Minimal failure taxonomy

### Fundamental assumption violation: no observation-entry filter

The package has no argument representing subject-specific observation entry.
The incorrect rows are not a malformed input under its current API; they expose
the shared-origin assumption. Classification: **fundamental assumption
violation**, with a silent-wrong-answer manifestation when historical covariates
are present.

### Downstream hard failure: PBC prediction

The PBC cause-specific model used 189 stacked observations. The two cause fits
used 12 and 52 events, respectively, and each returned one non-finite
coefficient, `LM1` [TRACE:
`03_break_diagnosis::pbcseq_cause1_observations_used`,
`::pbcseq_cause1_events_used`, `::pbcseq_cause2_events_used`,
`::pbcseq_cause1_nonfinite_coefficients`,
`::pbcseq_cause2_nonfinite_coefficients`; term names are in
`age_model_diagnostics.csv`]. Model construction itself did not error. Prediction
did, so `score()` was not attempted [TRACE:
`03_break_diagnosis::pbc_model_failed`,
`::pbc_prediction_failed`, `::pbc_score_attempted`]. Classification:
**hard downstream failure caused by a silently degenerate fit**.

### Silent degenerate fit: NAFLD

The NAFLD model used 16,890 stacked observations and 762 events, yet returned one
non-finite coefficient, again `LM1` [TRACE:
`03_break_diagnosis::nafld_cause1_observations_used`,
`::nafld_cause1_events_used`,
`::nafld_cause1_nonfinite_coefficients`; term name in
`age_model_diagnostics.csv`]. `dynamic_lm()` did not signal an error [TRACE:
`03_break_diagnosis::nafld_model_failed`]. Classification: **silent degenerate
fit**. The absolute-age landmark polynomial is numerically and structurally a
poor parameterization here; centering landmarks is a correction candidate, not
evidence that the entry problem is solved.

### Independent package bug: formula symbols

Passing an inline `as.formula(...)` works, while passing the same formula through
a local symbol fails with `object 'formula_object' not found`. Classification:
**package call-evaluation bug**, unrelated to left truncation. The exact condition
and reproducer are in
[`age_break_conditions.csv`](../results/age_break_conditions.csv) and
`analysis/03_break_diagnosis.R`.

The run captured two errors and no warnings [TRACE:
`03_break_diagnosis::error_count`, `::warning_count`]. No condition was silently
coerced.

## Single-landmark reference fits

At age 50 with a ten-year window, three plain `coxph()` analyses were fitted with
counting-process `Surv(entry, exit, event)`:

- **naive:** include everyone whose exit exceeds the landmark and set analysis
  entry to the landmark;
- **strict:** additionally require observation entry no later than the landmark;
- **delayed:** admit people who enter during the window, with analysis entry
  `max(landmark, entry_age)`.

All retained intervals were asserted to have `entry < exit`, all event indicators
were binary, and every coefficient and standard error was asserted finite.

For NAFLD, the naive BMI log-hazard coefficient was 0.2414276. The strict estimate
was 0.0614864, a naive-minus-strict discrepancy of 0.1799411; the delayed-entry
estimate was 0.2295852, a discrepancy of 0.0118424 [TRACE:
`03_break_diagnosis::nafld_naive_minus_strict_bmi_z_naive`,
`::nafld_naive_minus_strict_bmi_z_reference`,
`::nafld_naive_minus_strict_bmi_z_difference`,
`::nafld_naive_minus_delayed_bmi_z_reference`,
`::nafld_naive_minus_delayed_bmi_z_difference`].

For PBC death, the naive standardized log-bilirubin coefficient was 0.7868795.
The strict estimate was 1.8219446, a discrepancy of -1.0350651; the delayed-entry
estimate was 1.1414514, a discrepancy of -0.3545718 [TRACE:
`03_break_diagnosis::pbcseq_death_naive_minus_strict_log_bili_z_naive`,
`::pbcseq_death_naive_minus_strict_log_bili_z_reference`,
`::pbcseq_death_naive_minus_strict_log_bili_z_difference`,
`::pbcseq_death_naive_minus_delayed_log_bili_z_reference`,
`::pbcseq_death_naive_minus_delayed_log_bili_z_difference`].

The full estimates, standard errors, sample sizes, and event counts are in
[`age_coefficient_reference_fits.csv`](../results/age_coefficient_reference_fits.csv)
and [`age_coefficient_differences.csv`](../results/age_coefficient_differences.csv).

## Verdict before correction

- **PROVEN (code/set logic):** current risk-set construction cannot enforce
  `E_i <= a`; it implements only the exit-side condition.
- **EMPIRICAL-ONLY:** both public datasets contain many rows that the unmodified
  stacker assigns before cohort entry; NAFLD retains complete covariates on a
  substantial subset of them.
- **EMPIRICAL-ONLY:** naive, strict, and within-window delayed-entry coefficients
  differ in these applications. The data do not identify which discrepancy is
  bias.
- **PENDING:** known-truth bias, interval-estimator coverage, and the validity of
  any left-truncation influence-function correction. These require Phase 4.
