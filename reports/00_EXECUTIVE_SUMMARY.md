# Executive summary

## What I proved

The source audit proved a concrete eligibility defect in the original
shared-origin implementation: the age-scale stacker checked that exit followed
the landmark but did not require observation entry to precede it. The corrected
contract is `entry <= landmark < exit` for strict landmarking, or counting-
process entry at `max(entry, landmark)` for delayed entry. These are finite
risk-set definitions, not asymptotic claims; the source quotations and minimal
reproducers are in [01_recon.md](01_recon.md) and
[03_break_diagnosis.md](03_break_diagnosis.md).

The correction is backward compatible in the tested enrollment-scale and
gamma-zero simulation paths. For Phase 6, the new generator was compared with
the engine at commit `f6bc7ec` using all 108 Phase 4 cells and all 1,000 saved
seeds per cell. Cohorts, censored cohorts, three method-specific stacks, and fit
objects were exactly `identical()`—no tolerance [TRACE:
`07_informative_entry::regression_cells`,
`::regression_replicates_per_cell`, `::regression_exact_cohort_cells`,
`::regression_exact_fit_cells`]. Package-level enrollment-scale regression
checks are described in [04_correction.md](04_correction.md).

For the stated simulation DGP, I proved the retention-conditioning identity

```text
f(x,e | T > E) = f_X(x) f_E(e) S_T(e | x) / Pr(T > E).
```

It shows that source-population independence of `X` and entry `E` need not
survive the study's `T > E` retention rule. I also derived the post-landmark
conditional hazard under Regime B: when entry and failure are independent
conditional on measured `X`, conditioning on `E <= a` does not change the
specified hazard after `a`. These derivations are shown step by step in
[07_informative_entry.md](07_informative_entry.md). They do **not** prove Cox
estimator consistency or a variance formula.

## What I showed empirically but did not prove

Under independent entry, the complete Phase 4 grid falsified the proposed
coverage-degradation headline: naive, strict, and delayed fits had similar
coefficient behavior, and the targeted summary-AUC study likewise showed
little separation. This was a negative result, not a failed analysis; see
[05_simulation.md](05_simulation.md).

Under measured covariate-dependent entry, the defect mattered. At the
predeclared strongest-separation mechanism (`gamma = -4`), the full-factorial
confirmation pooled 108,000 estimates per method. Naive bias was 0.043250 with
coverage 0.92494, versus strict bias 0.001762 with coverage 0.94370 and delayed
bias 0.002134 with coverage 0.94464 [TRACE:
`07_informative_entry::pooled_confirmation_naive_bias`,
`::pooled_confirmation_naive_coverage`,
`::pooled_confirmation_strict_bias`, `::pooled_confirmation_strict_coverage`,
`::pooled_confirmation_delayed_bias`, and
`::pooled_confirmation_delayed_coverage`]. This is strong simulation evidence that explicit
eligibility corrects most measured-entry coefficient bias in the tested DGP;
it is not a general theorem.

Under shared frailty, eligibility correction was not sufficient. With `theta =
0.5` and `delta = -2`, pooled coefficient biases were -0.054841, -0.052436, and
-0.053453 for naive, strict, and delayed fits [TRACE:
`07_informative_entry::pooled_frailty_theta_0p5_delta_m2_naive_bias`,
`::pooled_frailty_theta_0p5_delta_m2_strict_bias`, and
`::pooled_frailty_theta_0p5_delta_m2_delayed_bias`]. Because
frailty also makes the one-covariate fitted model a marginal mixture, the design
does not isolate dependent truncation from omitted-variable misspecification.

The entry-risk-set audit reliably exposed structural pre-entry inclusion, and
observed-covariate tests detected strong Regime B dependence. Those tests were
not general informativeness diagnostics: they had elevated false positives
after survival-to-entry selection and low sensitivity to unmeasured frailty.
IDEA-006 remains `KEEP` only as a structural audit; the general correlation
screen is `KILL`.

Prediction and coefficient conclusions were not interchangeable. In the wide
`gamma = -4` cell, naive coefficient bias was 0.070566 and strict-target AUC
bias was 0.020885; strict values were -0.000180 and -0.000816. Delayed entry had
coefficient bias -0.000697 but differed from the strict-target AUC by 0.029105
because it evaluates a different entry population and interval [TRACE:
`07_informative_entry::focused_cell_003_naive_bias`,
`::focused_cell_003_strict_bias`, `::focused_cell_003_delayed_bias`,
`::auc_cell_003_naive_bias`, `::auc_cell_003_strict_bias`, and
`::auc_cell_003_delayed_bias`]. Under wide shared frailty, strict coefficient
bias was -0.053566 while strict oracle-score AUC bias was 0.000031 [TRACE:
`07_informative_entry::frailty_cell_009_strict_bias`,
`::auc_cell_021_strict_bias`]. Thus a risk-factor failure can coexist with
intact oracle discrimination.

The Phase 6 AUC study completed 7,200,000 subject-bootstrap resamples with no
captured scoring failure [TRACE:
`07_informative_entry::auc_total_bootstrap_resamples`, `::auc_failures`]. Strict
coverage was not uniformly nominal, so these intervals remain empirical
procedures rather than a replacement for the missing analytic variance.

## What I could not do, and exactly what is needed to do it

I did not derive the left-truncation Hájek projection or influence function for
the summary AUC/Brier estimator. `score_left_truncated()` therefore still
returns `"PENDING: not derived for staggered entry"` for analytic influence-
function inference. Completing this requires a senior-author-reviewed
counting-process derivation that includes estimation of the left-truncated
censoring law, followed by a proof of the cross-landmark covariance and a
simulation verification independent of the derivation.

I did not prove consistency, asymptotic normality, or nominal coverage for the
strict or delayed supermodel under informative entry. That requires a formal
estimating-equation argument under explicit conditional-independence,
positivity, censoring, and overlapping-landmark conditions. The simulations
support only the mechanisms that were run.

I could not identify dependent truncation from observed data when its common
cause was unmeasured. Addressing that requires measured proxies, a defensible
joint entry/failure model, an instrumental or sensitivity-analysis strategy,
and validation against a setting where the entry mechanism is substantively
known. The current structural audit cannot supply this information.

I did not establish design-based influence-function inference for the NHANES
landmark metrics or validate a trajectory encoder; those innovation tracks
remain `PARK`. They require, respectively, a survey-design estimand/variance
derivation and a frozen prediction-adapter/data contract before model work.

The CRAN incoming-feasibility check was completed, but the package is not
CRAN-clean: the tarball check ended with 1 WARNING and 5 NOTEs [TRACE:
`07_informative_entry::cran_check_warnings`, `::cran_check_notes`]. Resolving
them requires correcting the `splc` documentation mismatch, metadata/URL and
namespace diagnostics, top-level packaging exclusions, and slow examples, then
rerunning the tarball check in an environment able to verify current time.

## Decision

Use explicit entry on attained-age analyses. Strict eligibility is the simpler
primary analysis; delayed entry can be a prespecified efficiency analysis when
its different at-risk population is scientifically acceptable. Always publish
the pre-entry fraction by landmark. Do not claim that this patch solves
unmeasured dependent truncation, and do not report an analytic summary-metric
CI until its influence function is derived and verified.
