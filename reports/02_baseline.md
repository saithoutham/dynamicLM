# Baseline Reproduction

## Outcome

The current PBC competing-risks tutorial ran end to end through the
unpenalized supermodel, `pen_lm()`, subject-level cross-validation with
`cv.pen_lm()`, penalized model reconstruction, prediction, and `score()`.
The fixed seed was 20260908 (all Phase 1 traces use that seed).

The tutorial data contained 1,589 long-format rows from 280 subjects (TRACE
`02_baseline::pbc_long_rows` and `02_baseline::pbc_unique_ids`). Stacking five
landmarks produced 1,175 rows (TRACE `02_baseline::landmark_count` and
`02_baseline::stacked_rows`). The landmark risk-set counts were 280, 261, 247,
215, and 172 (TRACE `02_baseline::risk_set_lm{0,1,2,3,4}`); the asserted
non-increasing shared-origin property held.

An initial strict-interval assertion correctly stopped the run because 24
terminal visit rows have `time == tstart` after `get_pbc_long()` rounds days to
one decimal year (TRACE `02_baseline::terminal_visit_ties`). Inspection showed
no `time < tstart` rows. The assertion was therefore corrected to the actual
data invariant, `time >= tstart`; no row was silently changed or removed.

## Penalized fits and warnings

Each cause-specific ridge path contained 100 candidate lambdas (TRACE
`02_baseline::path_lambda_count_cause1` and
`02_baseline::path_lambda_count_cause2`). Subject-level cross-validation chose
minimum-error lambdas 0.025579 and 0.140643, and one-standard-error lambdas
52.603026 and 1.439528 (TRACE `02_baseline::cv_lambda_min_cause{1,2}` and
`02_baseline::cv_lambda_1se_cause{1,2}`). The tutorial's one-standard-error
choice was used for prediction.

The run emitted 488 instances of the exact warning `cox.fit: algorithm did not
converge`: 45 during `pen_lm()` and 443 during `cv.pen_lm()` (TRACE
`02_baseline::warning_count`, `02_baseline::pen_lm_warning_count`, and
`02_baseline::cv_pen_lm_warning_count`). There were two unique stage/message
patterns (TRACE `02_baseline::unique_warning_patterns`). The complete condition
log is `results/baseline_conditions.csv`. These warnings are a reproducibility
finding, not dismissed because the README calls them common.

## Landmark-specific performance

These are apparent, same-data estimates for the tutorial's transplant outcome;
they must not be interpreted as externally validated performance.

| Landmark | Model | AUC | SE | Confidence interval |
|---:|---|---:|---:|---:|
| 0 | LM | 0.729712 | 0.058066 | 0.615906–0.843519 |
| 0 | penLM | 0.434736 | 0.045243 | 0.346062–0.523411 |
| 2 | LM | 0.761378 | 0.053854 | 0.655827–0.866930 |
| 2 | penLM | 0.363519 | 0.072250 | 0.221912–0.505126 |
| 4 | LM | 0.860045 | 0.074571 | 0.713888–1.000000 |
| 4 | penLM | 0.205293 | 0.089878 | 0.029135–0.381452 |

Trace refs follow
`02_baseline::AUC_<model>_lm<landmark>_{AUC,se,lower,upper}`.

| Landmark | Model | Brier | SE | Confidence interval |
|---:|---|---:|---:|---:|
| 0 | Null model | 0.043060 | 0.011625 | 0.020276–0.065844 |
| 0 | LM | 0.042218 | 0.011257 | 0.020155–0.064281 |
| 0 | penLM | 0.043031 | 0.011631 | 0.020235–0.065826 |
| 2 | Null model | 0.068254 | 0.016761 | 0.035403–0.101104 |
| 2 | LM | 0.064972 | 0.015541 | 0.034512–0.095432 |
| 2 | penLM | 0.068879 | 0.016962 | 0.035634–0.102124 |
| 4 | Null model | 0.062902 | 0.020782 | 0.022170–0.103635 |
| 4 | LM | 0.051229 | 0.016482 | 0.018925–0.083534 |
| 4 | penLM | 0.064541 | 0.021535 | 0.022334–0.106749 |

Trace refs follow
`02_baseline::Brier_<model>_lm<landmark>_{Brier,se,lower,upper}`.

## Summary metrics

| Metric | Model | Estimate | SE | Confidence interval |
|---|---|---:|---:|---:|
| Summary AUC | LM | 0.783712 | 0.042085 | 0.701227–0.866197 |
| Summary AUC | penLM | 0.334516 | 0.041553 | 0.253073–0.415959 |
| Summary Brier | Null model | 0.058072 | 0.009963 | 0.038544–0.077600 |
| Summary Brier | LM | 0.052806 | 0.008795 | 0.035569–0.070044 |
| Summary Brier | penLM | 0.058817 | 0.010178 | 0.038868–0.078766 |

Trace refs follow
`02_baseline::{AUC_summary,Brier_summary}_<model>_{AUC/Brier,se,lower,upper}`.

The penalized tutorial model is a negative result in this run: its summary AUC
was lower than the unpenalized model's and below the no-discrimination value.
That result is consistent with its narrow prediction range, 0.016923–0.081403,
versus 0.002004–0.422060 for the unpenalized model (TRACE
`02_baseline::prediction_{min,max}_{penLM,LM}`). This does not diagnose why the
ridge path behaved poorly; the convergence warnings, apparent evaluation, and
one-standard-error penalty choice are all plausible contributors. No causal or
clinical conclusion is drawn.

## What the summary CI computes

**PROVEN (implementation algebra, not asymptotic validity).** At each landmark,
`riskRegression::Score()` supplies a subject-level IID contribution for the
IPCW AUC or Brier estimate. `summary_metric()` reshapes those contributions into
a subject-by-landmark matrix, inserts zero for subjects absent at a later
landmark, and rescales each column by its landmark sample fraction. It then
computes the cross-landmark sample covariance matrix. For an equal-weighted
average over (L) landmarks, the code uses gradient
(g=(1/L,\ldots,1/L)^T) and reports

\[
\widehat{SE}(\bar\theta) =
\sqrt{g^T\widehat\Sigma g/N}.
\]

[SOURCE: local `R/summary_metric.R`, lines 79–170]

This retains off-diagonal covariance between landmark estimates. Averaging the
per-landmark standard errors or confidence limits would discard those terms and
is not what the package does. The paper describes this as a multivariate
extension of the landmark-specific asymptotic IID decomposition.[^1]

**CONJECTURED outside the baseline setting.** This implementation does not
establish validity with random left truncation. That extension remains open.

**EMPIRICAL-ONLY for this run.** The displayed intervals are computed outputs,
but coverage was not assessed in Phase 1. Moreover, the tutorial evaluates the
same subjects used for fitting. The penLM paper states that its inferential
properties require independent testing data and independent censoring; those
conditions are not satisfied by this apparent-performance reproduction.[^1]

## Regression lock

`tests/testthat/test-baseline.R` reruns the full competing-risks pipeline with
the recorded seed and compares cohort dimensions, risk sets, event counts,
penalty paths, selected lambdas, coefficients, prediction ranges, score tables,
and captured warnings against `results/baseline_expected.rds`. The numerical
tolerance is declared in the test. This fixture is the enrollment-scale
backward-compatibility lock for later refactoring.

## Sources

[^1]: Fries AH, Choi E, Han SS. “[Penalized landmark supermodels (penLM) for dynamic prediction for time-to-event outcomes in high-dimensional data](https://link.springer.com/article/10.1186/s12874-024-02418-9).” *BMC Medical Research Methodology*. See the summary-metric methods and stated independent-censoring/independent-test-data assumptions.
[^2]: The Han Lab. “[dynamicLM PBC tutorial](https://github.com/thehanlab/dynamicLM).” The run follows the competing-risks PBC path through `pen_lm()` and `cv.pen_lm()`.
