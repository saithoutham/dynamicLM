# Phase 3 — age-scale correction

## Implemented result

The `age-scale` branch now has an explicit observation-entry contract.
`stack_data()` accepts `entry` and `entry_mode`:

- `shared` is the original, backward-compatible behavior;
- `strict` retains a subject at landmark `a` only when `E_i <= a < X_i`;
- `delayed` also admits subjects first observed before `a + w`, sets the
  analysis start to `max(a, E_i)`, and retrieves a long-format covariate at that
  subject-specific start.

Delayed stacks carry `.LM_entry` separately from `LM`. Penalized responses now
use `lmdata$entry_col`, and unpenalized fitting stops if the formula uses the
wrong counting-process start. The original wide-data path no longer asks for an
irrelevant `rtime` argument. A pre-existing formula-call bug was also fixed:
formula objects stored in symbols now reach the fitter.

**PROVEN (implementation invariant):** every corrected output interval is
asserted to satisfy `analysis_entry < analysis_exit`; strict rows are asserted
to satisfy `entry_age <= LM`; delayed starts are asserted equal to
`max(LM, entry_age)`.

## Risk-set correction

The full computed table is
[`corrected_risk_sets.csv`](../results/corrected_risk_sets.csv). At the age-50
landmark, strict PBC retained 67 rows and delayed entry retained 164 [TRACE:
`04_correction::pbcseq_death_strict_lm_50_rows`,
`::pbcseq_death_delayed_lm_50_rows`]. In NAFLD, the corresponding counts were
2,324 and 5,311 [TRACE: `04_correction::nafld_strict_lm_50_rows`,
`::nafld_delayed_lm_50_rows`].

These counts are not interchangeable sample-size gains. Strict landmarking
answers risk conditional on already being observed at age 50. Delayed mode
augments the window with people first observed later and therefore targets a
different population and covariate-assessment rule.

## Partial-likelihood propagation

Both modes were fitted through the patched `dynamic_lm()` path. Landmark values
were centered and scaled before constructing interactions; the unidentified raw
`LM1` main effect exposed in Phase 2 was not included.

| data | mode | stacked observations | events | term | estimate | robust SE |
|---|---|---:|---:|---|---:|---:|
| PBC death | strict | 167 | 38 | standardized log bilirubin | 1.436461 | 0.204402 |
| PBC death | delayed | 410 | 111 | standardized log bilirubin | 1.030305 | 0.093049 |
| NAFLD | strict | 7,993 | 429 | standardized BMI | 0.128007 | 0.075931 |
| NAFLD | delayed | 17,217 | 890 | standardized BMI | 0.184504 | 0.044780 |

Trace refs follow the pattern
`04_correction::<dataset>_<mode>_<term>_{observations,events,estimate,robust_se}`;
the exact values and subject counts are in
[`corrected_model_comparison.csv`](../results/corrected_model_comparison.csv).

**EMPIRICAL-ONLY:** delayed mode produced smaller robust standard errors in
these fits. This is not yet an efficiency theorem or a bias result because the
modes target different observed populations. Known-truth bias and repeated-run
efficiency are evaluated in Phase 4.

## Left-truncation-aware IPCW

For one landmark, `lt_censoring_weights()` fits the reverse product-limit curve

`survfit(Surv(analysis_entry, analysis_exit, censor_event) ~ 1)`.

If `G` denotes that fitted censoring survival, an observed event at `T_i` uses
inverse weight `G(A_i) / G(T_i-)`; a subject known to be event-free at horizon
`h` uses `G(A_i) / G(h)`. Thus censoring survival is conditional on the
subject-specific analysis entry `A_i`, rather than granting a late entrant
uncensored observation time before entry. Administrative exits at the horizon
are not counted as censoring events before that horizon.

The scorer implements cumulative/dynamic AUC with half credit for tied risks
and cause-specific Brier score. It asserts finite non-negative weights, positive
case/control totals, valid probabilities, and one row per subject-landmark.

In a hand-checkable delayed-entry risk set, the fitted horizon censoring
survival was 0.5 for both an early and a late entrant, but their conditional
inverse horizon weights were 2 and 1 [TRACE:
`04_correction::toy_censor_survival_horizon_early_entry`,
`::toy_censor_survival_horizon_late_entry`,
`::toy_inverse_weight_horizon_early_entry`,
`::toy_inverse_weight_horizon_late_entry`]. The resulting AUC was 1 and Brier
score was 0.045 [TRACE: `04_correction::toy_auc`, `::toy_brier`].

As a compatibility check, at the strict PBC age-50 landmark the new scorer and
the development `riskRegression::Score()` returned the same AUC to an absolute
difference of 1.1102230e-16 and the same Brier score to
2.7755576e-17 [TRACE:
`04_correction::strict_validation_max_auc_difference`,
`::strict_validation_max_brier_difference`; full values in
[`corrected_score_validation.csv`](../results/corrected_score_validation.csv)].
This is **EMPIRICAL-ONLY**, not proof for all inputs.

The apparent, in-sample PBC summary scores were AUC 0.901556 and Brier 0.164284
under strict entry, versus AUC 0.816041 and Brier 0.151394 under delayed entry
[TRACE: `04_correction::pbcseq_death_strict_summary_auc`,
`::pbcseq_death_strict_summary_brier`,
`::pbcseq_death_delayed_summary_auc`,
`::pbcseq_death_delayed_summary_brier`]. They are implementation checks, not
claims of out-of-sample performance.

## Inference boundary

`score_left_truncated()` deliberately does not reuse the original analytic IID
decomposition. Its object records:

`PENDING: not derived for staggered entry`.

When requested, it supplies a subject-level nonparametric bootstrap using the
same resampled subjects across landmarks, so the empirical resampling preserves
cross-landmark dependence. This is an implementation of a bootstrap procedure,
not an asserted asymptotic result.

- **PROVEN:** the code computes the stated counting-process risk sets and the
  stated finite-sample IPCW formulas.
- **CONJECTURED:** under suitable independent censoring/truncation conditions,
  the conditional weight ratio is the appropriate ingredient for a consistent
  age-scale landmark metric.
- **EMPIRICAL-ONLY:** point estimates agree with the established right-censored
  implementation in the strict-entry special case tested above.
- **PENDING:** a formal left-truncation Hájek/influence-function derivation and
  analytic summary-metric confidence intervals.

## Regression protection

`tests/testthat/test-age-scale.R` locks the strict and delayed risk sets,
subject-specific long-format covariate lookup, penalized-response entry column,
wrong-formula rejection, conditional IPCW weights, exact toy scores, wide-data
compatibility, and formula-symbol fix. `tests/testthat/test-baseline.R` reruns
the original PBC pipeline and compares it with the frozen enrollment-time
artifact. The complete suite passed after these changes.
