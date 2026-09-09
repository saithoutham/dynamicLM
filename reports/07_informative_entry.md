# Phase 6 — informative entry

## Result in one paragraph

**EMPIRICAL-ONLY:** measured covariate-dependent entry made the documented
eligibility defect consequential. At the anchor condition with the strongest
tested dependence (`gamma = -4`) and the widest entry spread, naive coefficient
bias was 0.070566, versus -0.000180 for strict eligibility and -0.000697 for
within-window delayed entry. In the full-factorial confirmation, pooled naive
bias was 0.043250 and coverage was 0.92494; strict and delayed biases were
0.001762 and 0.002134, with coverage 0.94370 and 0.94464. With shared unmeasured
frailty, all methods were biased, but the Phase 7 factorial decomposition
attributes the measured increment to `theta`, not `delta`: intervals for all
three methods' two `delta` contrasts and `theta`-by-`delta` interactions include
zero, whereas both `theta` contrasts exclude zero. The patch therefore corrects
observable eligibility; this experiment does not demonstrate an additional
dependent-truncation failure [TRACE:
`07_informative_entry::focused_cell_003_naive_bias`,
`::focused_cell_003_strict_bias`, `::focused_cell_003_delayed_bias`, and
`::pooled_confirmation_naive_bias`, `::pooled_confirmation_naive_coverage`,
`::pooled_confirmation_strict_bias`, `::pooled_confirmation_strict_coverage`,
`::pooled_confirmation_delayed_bias`,
`::pooled_confirmation_delayed_coverage`, and
`08_frailty_decomposition::{method}_{contrast}_{contrast_estimate,ci_lower,ci_upper}`].

## Why the independent-entry null needs a qualification

The Phase 4 generator samples source-population entry age `E` independently of
the measured covariate `X`, samples event age `T` conditional on `X`, and then
retains only subjects for whom `T > E`. The brief proposed that future entrants
therefore form an unbiased source-population sample. That statement is not
generally correct after retention.

**PROVEN for this data-generating mechanism:** let `A = {T > E}`. Before
retention, `f(X,E) = f_X(X) f_E(E)`. Bayes' rule gives

```text
f(x,e | A)
  = f_X(x) f_E(e) Pr(T > e | X=x) / Pr(A)
  = f_X(x) f_E(e) S_T(e | x) / Pr(A).
```

When the event hazard depends on `x`, `S_T(e | x)` is generally not separable
into a function of `e` times a function of `x`; `X` and `E` can therefore be
associated in the retained cohort even though they were independent in the
source population. More specifically, the measured-covariate density among
retained future entrants at landmark `a` is proportional to

```text
f_X(x) integral[e > a] f_E(e) S_T(e | x) de,
```

whereas the strict risk set has density proportional to

```text
f_X(x) Pr(E <= a) S_T(a | x).
```

They need not match. This refutes the proposed *unbiased-sample explanation*;
it does not refute the Phase 4 computation. The earlier null remains an
**EMPIRICAL-ONLY** robustness/cancellation result under that generator, not a
proof that admitting pre-entry subjects is harmless.

For Regime B, entry and failure are independent conditional on measured `X` in
the source mechanism. Hence, for `t >= a`,

```text
Pr(T in [t,t+dt) | T >= t, E <= a, X=x) / dt
  -> lambda0 exp(beta x).
```

This conditional-hazard equality is **PROVEN for the stated generator**.
Consistency and asymptotic normality of the fitted strict landmark
supermodel are **NOT PROVEN**; the simulation evidence below supports, but
does not establish, consistency.

## Implementation and exact backward compatibility

`generate_left_truncated_cohort()` now accepts `gamma`, `theta`, and `delta`.
The default `gamma = theta = delta = 0` branch preserves the old calls and old
random-number order. A hard gate loaded the engine at commit `f6bc7ec` and used
every Phase 4 seed: all 108 design cells, 1,000 replicates per cell, and all
three methods produced `identical()` cohorts, censored cohorts, stacks, and fit
objects. No numerical tolerance was used [TRACE:
`07_informative_entry::regression_cells`,
`::regression_replicates_per_cell`, `::methods`,
`::regression_exact_cohort_cells`, `::regression_exact_fit_cells`]. Cell-level
results are in
[`regime_a_regression_check.csv`](../results/regime_a_regression_check.csv).

The informative-entry sampler uses inverse-CDF draws from the requested
truncated normal. Every generated cohort asserts `event_age > entry`; every
stack asserts valid counting-process intervals, binary events, and unique
subject-landmark rows. Censoring was recalibrated independently for every
mechanism using a 50,000-subject pilot. Achieved and target censoring are both
retained in
[`informative_entry_calibration.csv`](../results/informative_entry_calibration.csv)
[TRACE: `07_informative_entry::calibration_pilot_n`].

## Simulation design

The focused sweep crossed four gamma values (`0`, `-1`, `-2`, `-4`) with entry
SDs `0.25`, `4`, and `10` at landmark spacing `4`, censoring target `0.30`, and
sample size `750`. Each cell used 1,000 replicates and all three methods. The
largest mean absolute naive-versus-strict bias separation selected `gamma =
-4` for the complete 108-cell Phase 4 grid. The mean separation used for that
selection was 0.041039 [TRACE: `07_informative_entry::focused_cells`,
`::replicates_per_cell`, `::selected_gamma`,
`::selected_gamma_mean_absolute_separation`, `::confirmation_cells`].

The frailty sweep crossed `theta` in `{0, 0.5}`, `delta` in `{0, -2}`, and all
three entry spreads at the anchor condition. Across focused, confirmation, and
frailty analyses, 396,000 fits were attempted and none failed [TRACE:
`07_informative_entry::frailty_cells`, `::total_fits`, `::fit_failures`]. Every
seed is in
[`informative_entry_seed_manifest.csv`](../results/informative_entry_seed_manifest.csv).

## Measured covariate-dependent entry

The mechanism fired: at `gamma = -4`, the mean observed entry--`X`
correlations were -0.997973, -0.705448, and -0.358284 as entry SD increased.
The corresponding pre-entry fractions at the first landmark were 0.14909,
0.24253, and 0.37183 [TRACE:
`07_informative_entry::focused_cell_001_naive_mean_entry_x_correlation`
through `::focused_cell_003_naive_mean_entry_x_correlation`, plus
`::focused_cell_001_naive_mean_preentry_fraction_lm1` through
`::focused_cell_003_naive_mean_preentry_fraction_lm1`].

Primary anchor results are below. Coverage Monte Carlo standard errors use the
binomial formula with the cell's successful replicate count.

| entry SD | method | bias | empirical SE | mean SE | SE ratio | coverage (MCSE) | power |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0.25 | naive | 0.016876 | 0.124506 | 0.119808 | 0.9623 | 0.941 (0.00745) | 0.936 |
| 0.25 | strict | 0.002205 | 0.130351 | 0.126914 | 0.9736 | 0.942 (0.00739) | 0.896 |
| 0.25 | delayed | 0.002117 | 0.127661 | 0.122977 | 0.9633 | 0.944 (0.00727) | 0.911 |
| 4 | naive | 0.025247 | 0.114827 | 0.120380 | 1.0484 | 0.961 (0.00612) | 0.953 |
| 4 | strict | -0.012453 | 0.122738 | 0.129524 | 1.0553 | 0.964 (0.00589) | 0.861 |
| 4 | delayed | -0.012105 | 0.119174 | 0.125002 | 1.0489 | 0.964 (0.00589) | 0.887 |
| 10 | naive | 0.070566 | 0.128658 | 0.127963 | 0.9946 | 0.927 (0.00823) | 0.965 |
| 10 | strict | -0.000180 | 0.139836 | 0.137891 | 0.9861 | 0.937 (0.00768) | 0.831 |
| 10 | delayed | -0.000697 | 0.131766 | 0.130822 | 0.9928 | 0.938 (0.00763) | 0.871 |

All table entries are mapped in `trace.csv` by the label template
`focused_cell_{001..003}_{method}_{field}`. The remaining gamma cells, including the
independent-entry control, and every requested row/event/subject mean and
landmark-specific pre-entry fraction are in
[`informative_entry_summary.csv`](../results/informative_entry_summary.csv).

The full-factorial confirmation pooled equally replicated cell-level estimates:

| method | estimates | bias | coverage (MCSE) | power (MCSE) |
|---|---:|---:|---:|---:|
| naive | 108,000 | 0.043250 | 0.92494 (0.00080) | 0.92991 (0.00078) |
| strict | 108,000 | 0.001762 | 0.94370 (0.00070) | 0.85323 (0.00108) |
| delayed | 108,000 | 0.002134 | 0.94464 (0.00070) | 0.87613 (0.00100) |

These values are mapped in `trace.csv` by the label template
`pooled_confirmation_{method}_{field}` and are stored in
[`informative_entry_pooled_summary.csv`](../results/informative_entry_pooled_summary.csv).
The full cell results—not only the pooled summary—remain the inferential record.

## Shared-frailty limitation

> **Superseded interpretation (2026-09-09):** “The eligibility patch is not
> sufficient under this shared frailty stress test. This experiment does not
> isolate dependent truncation: with `theta > 0`, omitting `U` also makes the
> fitted one-covariate hazard model misspecified through frailty mixing. The
> similar results for `delta = 0` and `delta = -2` show that the observed bias
> cannot honestly be attributed solely to informative entry.”
>
> **Revision (2026-09-09):** the replicate-level `delta` contrasts do isolate
> the increment associated with making entry depend on `U` within this
> factorial design. The original paragraph correctly warned against attributing
> all bias to entry, but it was too conservative about what the `delta`
> comparison can identify.

When `theta = 0.5`, all three one-covariate fits were biased regardless of
whether entry also depended on frailty through `delta`. Pooling the three entry
spreads:

| theta | delta | method | bias | coverage (MCSE) |
|---:|---:|---|---:|---:|
| 0 | -2 | naive | 0.000489 | 0.9490 (0.00402) |
| 0 | -2 | strict | -0.000179 | 0.9503 (0.00397) |
| 0 | -2 | delayed | -0.000119 | 0.9507 (0.00395) |
| 0 | 0 | naive | 0.002654 | 0.9440 (0.00420) |
| 0 | 0 | strict | 0.001865 | 0.9460 (0.00413) |
| 0 | 0 | delayed | 0.001759 | 0.9437 (0.00421) |
| 0.5 | -2 | naive | -0.054841 | 0.9187 (0.00499) |
| 0.5 | -2 | strict | -0.052436 | 0.9197 (0.00496) |
| 0.5 | -2 | delayed | -0.053453 | 0.9207 (0.00493) |
| 0.5 | 0 | naive | -0.050872 | 0.9277 (0.00473) |
| 0.5 | 0 | strict | -0.051890 | 0.9290 (0.00469) |
| 0.5 | 0 | delayed | -0.051292 | 0.9267 (0.00476) |

All table entries are mapped in `trace.csv` by the label template
`pooled_frailty_theta_{theta}_delta_{delta}_{method}_{field}`. The formal
factorial re-analysis is in
[`08_frailty_decomposition.md`](08_frailty_decomposition.md).

**EMPIRICAL-ONLY:** for naive, strict, and delayed fits, the Monte Carlo
intervals for the `delta` contrast included zero at `theta = 0` and `theta =
0.5`, and the `theta`-by-`delta` interaction intervals also included zero. The
two `theta` main-effect intervals excluded zero for every method [TRACE:
`08_frailty_decomposition::naive_delta_at_theta_0_estimate`,
`::naive_delta_at_theta_0p5_estimate`, `::naive_theta_x_delta_estimate`,
`::strict_delta_at_theta_0_estimate`,
`::strict_delta_at_theta_0p5_estimate`, `::strict_theta_x_delta_estimate`,
`::delayed_delta_at_theta_0_estimate`,
`::delayed_delta_at_theta_0p5_estimate`, and
`::delayed_theta_x_delta_estimate`, with matching `_ci_low` and `_ci_high`
fields; the `theta_at_delta_*` fields trace the main effects]. Within this
design, the observed coefficient bias is driven by including an omitted hazard
frailty (`theta`), not by the added entry dependence (`delta`). This does not
prove that dependent truncation is harmless under other mechanisms.

## Entry-risk-set diagnostic

The promoted diagnostic reports (i) the pre-entry fraction, (ii) observed
entry--covariate correlation, and (iii) a Welch test comparing `X` among the
strict risk set and the disjoint pre-entry rows at each landmark. Landmark tests
use a predeclared Bonferroni threshold; all detection probabilities below come
from 1,000 datasets [TRACE: `07_informative_entry::diagnostic_alpha`,
`::diagnostic_bonferroni_alpha`, `::diagnostic_replicates_per_cell`].

| regime/mechanism | entry SD | corr test | structural pre-entry | any `X` shift | combined |
|---|---:|---:|---:|---:|---:|
| A: gamma 0 | 0.25 | 0.054 | 0.000 | 0.000 | 0.054 |
| A: gamma 0 | 10 | 0.156 | 1.000 | 0.041 | 0.176 |
| B: gamma -4 | 0.25 | 1.000 | 1.000 | 1.000 | 1.000 |
| B: gamma -4 | 10 | 1.000 | 1.000 | 1.000 | 1.000 |
| C: theta 0.5, delta -2 | 0.25 | 0.514 | 1.000 | 0.021 | 0.525 |
| C: theta 0.5, delta -2 | 10 | 0.065 | 1.000 | 0.049 | 0.104 |

These first-landmark values are mapped in `trace.csv` under diagnostic cells
`010`, `012`, `001`, `003`, `019`, and `021`. All landmarks and mechanisms are in
[`entry_diagnostic_summary.csv`](../results/entry_diagnostic_summary.csv).

Verdict: retain IDEA-006 as a **structural eligibility audit**, not as a test of
informative entry. The pre-entry fraction correctly detects when naive and
strict risk sets differ, including under independent entry. Observed-covariate
screens detected Regime B but had increasing false positives in the retained
gamma-zero cohort and low sensitivity to unmeasured frailty. A general
correlation-based informativeness test is killed.

## Prediction versus estimation

The AUC experiment reused the focused and frailty mechanisms at the anchor
condition. It used the true generating score—`beta*X` in Regimes A/B and
`beta*X + theta*U` in Regime C—so fitted-model error is absent. Each of 24 cells
used 1,000 outer datasets, three methods, and 100 subject bootstrap resamples
shared across landmarks: 72,000 method datasets and 7,200,000 resamples. No
scoring task failed [TRACE:
`07_informative_entry::auc_design_cells`,
`::auc_outer_replicates_per_cell`, `::auc_methods`,
`::auc_bootstrap_replicates`, `::auc_total_method_datasets`,
`::auc_total_bootstrap_resamples`, `::auc_failures`].

### Strict-target numerical truth

Let `q_a(x,u) = Pr(E <= a | x,u)`, `eta = beta*x + theta*u`, and
`S(t|eta) = exp{-lambda0*t*exp(eta)}`. Among the source normal density
`f(x,u)`, the unnormalized eligible case and control masses at landmark `a` and
horizon `a+w` are

```text
C_a(x,u) = f(x,u) q_a(x,u) [S(a|eta) - S(a+w|eta)]
D_a(x,u) = f(x,u) q_a(x,u) S(a+w|eta).
```

Normalizing these masses and integrating concordance of the generating score
gives the cumulative/dynamic AUC among people enrolled and event-free at the
landmark. **PROVEN for this finite DGP:** these masses follow directly from the
joint event/entry probabilities. The reported integrals are numerical, not
closed-form. Fine versus coarse quadrature differed by at most 9.502184e-06,
below the predeclared `1e-5` criterion [TRACE:
`07_informative_entry::auc_max_quadrature_difference`,
`::auc_quadrature_tolerance`].

Naive and delayed results below are differences from this **strict scientific
target**. In particular, delayed entry admits a different population and uses a
shorter interval for post-landmark entrants; its difference must not be called
bias for the delayed method's own estimand.

### Coefficient and AUC bias in the same cells

| entry SD | method | coefficient bias | coefficient coverage | strict-target AUC bias | AUC coverage |
|---:|---|---:|---:|---:|---:|
| 0.25 | naive | 0.016876 | 0.941 | 0.010075 | 0.924 |
| 0.25 | strict | 0.002205 | 0.942 | -0.000365 | 0.946 |
| 0.25 | delayed | 0.002117 | 0.944 | 0.009611 | 0.929 |
| 4 | naive | 0.025247 | 0.961 | 0.012072 | 0.931 |
| 4 | strict | -0.012453 | 0.964 | -0.003855 | 0.962 |
| 4 | delayed | -0.012105 | 0.964 | 0.011644 | 0.931 |
| 10 | naive | 0.070566 | 0.927 | 0.020885 | 0.904 |
| 10 | strict | -0.000180 | 0.937 | -0.000816 | 0.939 |
| 10 | delayed | -0.000697 | 0.938 | 0.029105 | 0.845 |

The coefficient fields map to focused cells `001`--`003`; the AUC fields map to
the `trace.csv` label template `auc_cell_{001..003}_{method}_{field}`. The exact joined table is
[`informative_prediction_estimation_comparison.csv`](../results/informative_prediction_estimation_comparison.csv).

**EMPIRICAL-ONLY:** at the widest spread, naive eligibility distorted both the
hazard coefficient and discrimination relative to the strict target. Strict
eligibility removed most of both differences. Delayed coefficient behavior was
close to strict, but its AUC difference grew because its evaluation population
and risk interval differ; that is an estimand warning, not evidence that its
oracle score ranks its own target population poorly.

The wide-entry shared-frailty cell separates estimation and prediction even
more sharply:

| method | coefficient bias | coefficient coverage | strict-target AUC bias | AUC coverage |
|---|---:|---:|---:|---:|
| naive | -0.058213 | 0.911 | 0.007122 | 0.934 |
| strict | -0.053566 | 0.909 | 0.000031 | 0.943 |
| delayed | -0.055533 | 0.916 | 0.018159 | 0.888 |

These values trace to coefficient frailty cell `009` and AUC cell `021` fields.
The strict coefficient is biased because the fitted model omits frailty, while
the strict oracle AUC uses that frailty and remains close to its generating
target. This is not attainable analyst performance when `U` is unmeasured.

The strict bootstrap did not have uniformly nominal coverage. In the three
`gamma = -4` cells its coverage was 0.946, 0.962, and 0.939, with Monte Carlo
standard errors 0.00715, 0.00605, and 0.00757 [TRACE:
`07_informative_entry::auc_cell_001_strict_coverage`,
`::auc_cell_002_strict_coverage`, `::auc_cell_003_strict_coverage`, and matching
`_coverage_mcse` fields]. The bootstrap remains **EMPIRICAL-ONLY**, not a
replacement for the missing analytic influence function. Full results are in
[`informative_auc_coverage_summary.csv`](../results/informative_auc_coverage_summary.csv).

## Tolerances

One Phase 6 tolerance was loosened and recorded [TRACE:
`07_informative_entry::phase6_tolerance_changes`]. Recomputing the deterministic
AUC truth table in a fresh process differed from the saved table by at most
1.113554e-13, so an attempted exact-equality check was replaced by a `1e-12`
threshold [TRACE: `07_informative_entry::auc_truth_reproduction_max_abs_difference`,
`::tolerance_4_original`, `::tolerance_4_new`]. This remains far below the
predeclared `1e-5` fine-versus-coarse quadrature tolerance; it is nevertheless
a real post-run relaxation and is not described as exact [TRACE:
`07_informative_entry::auc_quadrature_tolerance`].

The three Phase 4 changes are documented retroactively in
[`tolerance_changes.csv`](../results/tolerance_changes.csv): censoring validation
`1e-5` to `0.01`, AUC bounds `0` to `sqrt(.Machine$double.eps)`, and Cox-engine
validation `1e-10` to `1e-8`. The exact original failure magnitudes for the
first two were not retained and remain `PENDING`; the ledger does not invent
them [TRACE: `07_informative_entry::tolerance_1_original` through
`::tolerance_4_new`].

## CRAN incoming-feasibility check

The correct tarball workflow completed:

```text
R CMD build .
_R_CHECK_CRAN_INCOMING_=TRUE _R_CHECK_FORCE_SUGGESTS_=false \
  R CMD check --as-cran --no-manual dynamicLM_1.0.0.tar.gz
```

The check process exited successfully, but its package status was not clean: 1
WARNING and 5 NOTEs. The WARNING was a pre-existing data-documentation mismatch
for `splc` and `splc_test`; NOTES covered incoming metadata/URL issues, inability
to verify current time, non-standard top-level files, visible global bindings,
and long examples. Tests passed inside the check. This completes the previously
stalled feasibility item but is **not** a claim of CRAN readiness [TRACE:
`07_informative_entry::cran_check_exit_code`, `::cran_check_warnings`,
`::cran_check_notes`; details in
[`cran_check_summary.csv`](../results/cran_check_summary.csv)].

A direct `Rscript tests/testthat.R` invocation from the repository root failed
with `No test files found`; that entry point expects the installed-package
layout created by `R CMD check`. The README now uses `devtools::test()`, which
ran 36 tests with zero failures, warnings, or skips [TRACE:
`07_informative_entry::project_test_passes`, `::project_test_failures`,
`::project_test_warnings`, `::project_test_skips`; exact commands in
[`project-tests.csv`](../results/check/project-tests.csv)].

## Theory status and remaining boundary

- **PROVEN:** the retention-conditioned density calculation above, the
  conditional hazard identity within the stated Regime B DGP, and finite code
  assertions/exact backward compatibility.
- **EMPIRICAL-ONLY:** coefficient/AUC bias, SE, coverage, power, diagnostic
  sensitivity, and numerical quadrature behavior.
- **CONJECTURED:** strict/delayed landmark estimators remain consistent whenever
  entry is independent of failure conditional on all modeled predictors.
- **NOT PROVEN / PENDING:** the left-truncation Hájek influence function and an
  analytic standard error for staggered entry. `score_left_truncated()` still
  returns `"PENDING: not derived for staggered entry"`. The subject bootstrap
  is evaluated empirically and is not described as a replacement.
