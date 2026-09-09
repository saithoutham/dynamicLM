# Phase 7 — frailty decomposition

## Revised conclusion

**EMPIRICAL-ONLY:** the saved Regime C factorial separates the entry effect
(`delta`) from the omitted-hazard-frailty effect (`theta`). Across naive,
strict, and delayed fits, both `delta` contrasts and the `theta × delta`
interaction were compatible with zero at Monte Carlo resolution. Both `theta`
contrasts were negative and their Monte Carlo intervals excluded zero. The
observed coefficient bias is therefore driven almost entirely by including
`U` in the hazard while omitting it from the fitted model, not by making entry
depend on `U`.

This sharpens the conservative Phase 6 statement. It does **not** prove that
dependent truncation is harmless generally. **PROVEN for the stated
generator:** in the executed `theta = 0`, `delta = -2` arm, entry depends on
unmeasured `U`, but event time does not; entry and failure are therefore
independent conditional on measured `X`. When `theta = 0.5`, `delta = 0`, `U`
affects failure but not entry, yet the one-covariate fit is already biased.
**EMPIRICAL-ONLY:** adding `delta = -2` contributed little resolved additional
bias.

## Data and estimands

No cohorts were generated and no model was refit. The analysis read the 36,000
successful Regime C method rows in
[`informative_entry_raw.csv`](../results/informative_entry_raw.csv) [TRACE:
`08_frailty_decomposition::raw_frailty_rows`]. For each method and each
`theta`/`delta` arm, the three entry-spread cells contributed 3,000 replicate
estimates [TRACE: `08_frailty_decomposition::arm_replicates`].

The four arm means for method `m` are denoted
`b_m(theta, delta)`. Because the generating `beta` cancels in differences,
contrasts of mean estimates equal contrasts of bias:

```text
delta at theta:  b_m(theta, -2) - b_m(theta, 0)
theta at delta:  b_m(0.5, delta) - b_m(0, delta)
interaction:     [b_m(0.5, -2) - b_m(0.5, 0)]
               - [b_m(0, -2)   - b_m(0, 0)].
```

The theta/delta arms have disjoint cell-specific seed ranges. The observed
seed-overlap count was zero [TRACE:
`08_frailty_decomposition::seed_overlap_count`]. They were therefore treated as
independent. Methods within an arm share cohort seeds, but no cross-method
contrast is reported here.

For a contrast `C = sum_k w_k b_k`, the estimated Monte Carlo variance was

```text
Var_MC(C) = sum_k w_k^2 s_k^2 / R_k,
```

where `s_k^2` is the replicate-level estimate variance and `R_k` the arm
replicate count. Intervals are `C ± qnorm(0.975) SE_MC(C)`. These are 95%
Monte Carlo precision intervals for the simulation contrasts, not confidence
intervals for patient data and not evidence of estimator asymptotics [TRACE:
`08_frailty_decomposition::confidence_level`].

## Results

| method | contrast | estimate | MCSE | 95% Monte Carlo interval | includes zero |
|---|---|---:|---:|---:|---:|
| naive | delta at theta 0 | -0.002166 | 0.003187 | [-0.008412, 0.004081] | yes |
| naive | delta at theta 0.5 | -0.003969 | 0.003168 | [-0.010178, 0.002240] | yes |
| naive | theta at delta 0 | -0.053526 | 0.003167 | [-0.059733, -0.047318] | no |
| naive | theta at delta -2 | -0.055330 | 0.003188 | [-0.061578, -0.049081] | no |
| naive | theta × delta | -0.001804 | 0.004494 | [-0.010612, 0.007004] | yes |
| strict | delta at theta 0 | -0.002044 | 0.003272 | [-0.008457, 0.004369] | yes |
| strict | delta at theta 0.5 | -0.000546 | 0.003246 | [-0.006909, 0.005816] | yes |
| strict | theta at delta 0 | -0.053755 | 0.003240 | [-0.060106, -0.047404] | no |
| strict | theta at delta -2 | -0.052258 | 0.003278 | [-0.058682, -0.045833] | no |
| strict | theta × delta | 0.001497 | 0.004609 | [-0.007536, 0.010531] | yes |
| delayed | delta at theta 0 | -0.001879 | 0.003191 | [-0.008134, 0.004376] | yes |
| delayed | delta at theta 0.5 | -0.002161 | 0.003172 | [-0.008378, 0.004056] | yes |
| delayed | theta at delta 0 | -0.053051 | 0.003171 | [-0.059266, -0.046837] | no |
| delayed | theta at delta -2 | -0.053334 | 0.003192 | [-0.059591, -0.047077] | no |
| delayed | theta × delta | -0.000282 | 0.004500 | [-0.009101, 0.008537] | yes |

Every table value is stored in
[`frailty_decomposition.csv`](../results/frailty_decomposition.csv) and mapped
to `results/trace.csv` by
`08_frailty_decomposition::{method}_{contrast}_{field}`. The full-precision
fields are `contrast_estimate`, `monte_carlo_se`, `ci_lower`, `ci_upper`, and
`interval_includes_zero`.

The `delta` estimates ranged from -0.003969 to -0.000546, and all six intervals
included zero. The two method-specific `theta` effects ranged from -0.055330 to
-0.052258, and all six intervals excluded zero. All three interaction intervals
included zero. These are **EMPIRICAL-ONLY** comparisons; absence of a detectable
`delta` effect is not proof of exact equality.

## Interpretation correction

The Phase 6 report originally said that the design could not identify a
dependent-truncation increment because `theta > 0` also induced frailty mixing.
That statement was too conservative: the factorial `delta` contrast directly
estimates the increment within each `theta` level. The increment was small
relative to its Monte Carlo error for every method, while the `theta` effect was
large and stable across `delta`.

The corrected conclusion is narrower and more useful. The code patch addresses
eligibility. Regime C does not show that the patch fails under informative
truncation; it shows attenuation when a hazard frailty is omitted from the
fitted model. A mechanism in which an unmeasured variable jointly affects entry
and failure **without** independently misspecifying the hazard model was not
isolated here and remains untested.

## Claim labels

- **PROVEN:** contrast algebra, disjoint seed sets, and the finite calculations
  asserted in `analysis/08_frailty_decomposition.R`.
- **EMPIRICAL-ONLY:** the magnitudes, Monte Carlo intervals, and decomposition
  interpretation for the executed simulation arms.
- **NOT PROVEN:** consistency, asymptotic normality, general robustness to
  dependent truncation, or any analytic summary-metric variance result.
