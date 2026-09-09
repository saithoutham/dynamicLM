# Phase 4 — simulation with known truth

## Headline: the proposed degradation pattern did not occur

The prespecified headline was falsified under independent staggered entry.
Naive coefficient coverage did **not** deteriorate as entry-age heterogeneity
increased, and strict correction did **not** restore a pattern that was absent.
The same was true for the targeted summary-AUC experiment: the naive and strict
procedures had similar bias and coverage.

This negative result is prominent because it narrows the claim supported by the
public-data break diagnosis. The unmodified stacker demonstrably includes
unobserved person-time, but entry heterogeneity alone is not sufficient to
guarantee detectable coefficient or AUC bias under every data-generating
mechanism.

## Data-generating mechanism

Each source subject received a standard-normal covariate `X`, an entry age from
a truncated normal distribution, and an event age from an exponential
proportional-hazards model. Only subjects with event age greater than entry age
were retained. The true log hazard ratio was 0.4054651, corresponding to a
hazard ratio of 1.5, and the baseline hazard was 0.01 per year [TRACE:
`05_simulation::truth_beta`, `::hazard_ratio`, `::lambda0`]. Censor time after
entry was generated independently; rates were calibrated with a separate
50,000-subject pilot for each design condition [TRACE:
`05_simulation::calibration_pilot_n`].

Every generated cohort was asserted to satisfy `event_age > entry`, every
counting-process row satisfied `start < stop`, event indicators were binary,
and subject-landmark rows were unique. Failed fits were captured with their
exact messages rather than dropped.

The grid crossed the following traced design values:

- landmark spacings 2, 4, and 6 years;
- censoring targets 0, 0.15, 0.30, and 0.50;
- sample sizes 500, 750, and 1,500;
- entry-age standard deviations 0.25, 4, and 10 years.

Trace refs: `05_simulation::design_landmark_interval_1` through `_3`,
`::design_censoring_target_1` through `_4`, `::design_sample_size_1` through
`_3`, and `::design_entry_sd_1` through `_3`.

This produced 108 cells, each with 1,000 replicates and three methods, for
324,000 fitted models. No fit failed [TRACE: `05_simulation::grid_cells`,
`::replicates_per_cell`, `::methods`, `::total_model_fits`,
`::fit_failures`]. Eight processes were requested; measured simulation wall
time was 409.345 seconds [TRACE: `05_simulation::parallel_cores`,
`::wall_runtime_seconds`]. The exact seed for every cell-replicate pair is in
[`simulation_seed_manifest.csv`](../results/simulation_seed_manifest.csv).

## Fitting-engine verification

The grid used a scalar Newton solver for the Breslow partial likelihood and a
subject-clustered sandwich variance. Before use, it was compared with
`survival::coxph()` on naive, strict, and delayed stacks. The maximum absolute
coefficient difference was 1.1709017e-09 and the maximum robust-SE difference
was 6.8008849e-11 [TRACE:
`05_simulation::engine_validation_max_estimate_difference`,
`::engine_validation_max_robust_se_difference`]. The exact comparison is in
[`simulation_engine_validation.csv`](../results/simulation_engine_validation.csv).

This establishes numerical agreement for the tested scalar models; it is not a
general reimplementation of `coxph()`.

## Coefficient results over the complete grid

Pooled over the 36 cells within each entry-spread condition [TRACE:
`05_simulation::design_cells_per_entry_sd`]:

| entry-age SD | method | bias | mean estimated SE | coverage (MCSE) | power (MCSE) |
|---:|---|---:|---:|---:|---:|
| 0.25 | naive | 0.001819 | 0.114842 | 0.94572 (0.00119) | 0.90633 (0.00154) |
| 0.25 | strict | 0.001819 | 0.114842 | 0.94572 (0.00119) | 0.90633 (0.00154) |
| 0.25 | delayed | 0.001819 | 0.114842 | 0.94572 (0.00119) | 0.90633 (0.00154) |
| 4 | naive | 0.000338 | 0.115656 | 0.94319 (0.00122) | 0.90117 (0.00157) |
| 4 | strict | 0.000414 | 0.117152 | 0.94308 (0.00122) | 0.89506 (0.00162) |
| 4 | delayed | 0.000270 | 0.115654 | 0.94342 (0.00122) | 0.90058 (0.00158) |
| 10 | naive | 0.003436 | 0.123246 | 0.94397 (0.00121) | 0.87639 (0.00173) |
| 10 | strict | 0.001805 | 0.129711 | 0.94297 (0.00122) | 0.85186 (0.00187) |
| 10 | delayed | 0.001622 | 0.123188 | 0.94297 (0.00122) | 0.87600 (0.00174) |

Every cell in this pooled table has 1,000 successful replicates and zero
failures. Trace refs follow
`05_simulation::entry_sd_<spread>_<method>_<field>` and the exact pooled values
are in
[`simulation_heterogeneity_summary.csv`](../results/simulation_heterogeneity_summary.csv).

Bias, empirical SE, mean estimated SE, coverage, coverage Monte Carlo error, and
power are reported for every individual grid cell in
[`simulation_summary.csv`](../results/simulation_summary.csv). For example, in
the spacing-4, 30%-censoring, `n=750` conditions, cell-level results are traceable
under `cell_020_*`, `cell_056_*`, and `cell_092_*` for increasing entry spread.

At the widest spread in that slice, naive empirical SE and mean estimated SE
were 0.126309 and 0.127165, with coverage 0.948 (MCSE 0.00715) [TRACE:
`05_simulation::cell_092_naive_empirical_se`,
`::cell_092_naive_mean_estimated_se`, `::cell_092_naive_coverage`,
`::cell_092_naive_coverage_mcse`]. Strict values were 0.132774 and 0.132789,
with coverage 0.935 (MCSE 0.00780) [TRACE:
`05_simulation::cell_092_strict_empirical_se`,
`::cell_092_strict_mean_estimated_se`, `::cell_092_strict_coverage`,
`::cell_092_strict_coverage_mcse`].

**EMPIRICAL-ONLY:** within-window delayed entry recovered much of the precision
lost by strict eligibility at wide spread, with similarly small pooled bias.
Because delayed and strict samples are not identical, this does not establish a
universal efficiency advantage.

## Numerical AUC truth

The true cumulative/dynamic AUC was computed by deterministic normal-grid
quadrature under the generating hazard. Repeating the calculation at half the
grid resolution changed each value by less than 2.4e-08 [TRACE fields ending in
`quadrature_difference`].

For spacing 4, landmark AUCs were 0.613348 at age 50, 0.613002 at age 54, and
0.612664 at age 58; the true summary AUC was 0.613005 [TRACE:
`05_simulation::auc_truth_interval_4_landmark_50_auc`,
`::auc_truth_interval_4_landmark_54_auc`,
`::auc_truth_interval_4_landmark_58_auc`,
`::auc_truth_interval_4_landmark_50_summary_auc`]. These are **numerical truth
values**, not closed-form proofs. All spacings are in
[`simulation_auc_truth.csv`](../results/simulation_auc_truth.csv).

## Summary-AUC interval coverage

A separate experiment isolated metric inference by using the true generating
risk score instead of an estimated model. It fixed spacing 4, target censoring
at 0.30, and sample size 750 [TRACE:
`05_simulation::auc_design_landmark_interval`,
`::auc_design_censoring_target`, `::auc_design_sample_size`], then ran 1,000
outer datasets per entry spread.
Each interval used 100 subject bootstrap resamples shared across landmarks.
That is 600,000 bootstrap resamples across naive and strict analyses; none of
the outer scoring tasks failed [TRACE:
`05_simulation::auc_outer_replicates_per_cell`,
`::auc_bootstrap_replicates`, `::auc_total_bootstrap_resamples`,
`::auc_coverage_failures`].

| entry-age SD | method | AUC bias | empirical SE | mean bootstrap SE | coverage (MCSE) |
|---:|---|---:|---:|---:|---:|
| 0.25 | naive | 0.000323 | 0.034452 | 0.033465 | 0.930 (0.00807) |
| 0.25 | strict | 0.000323 | 0.034452 | 0.033446 | 0.932 (0.00796) |
| 4 | naive | -0.000009 | 0.034123 | 0.033445 | 0.937 (0.00768) |
| 4 | strict | -0.000282 | 0.034658 | 0.033780 | 0.940 (0.00751) |
| 10 | naive | -0.000804 | 0.034312 | 0.035476 | 0.951 (0.00683) |
| 10 | strict | -0.001229 | 0.036211 | 0.037617 | 0.954 (0.00662) |

Trace refs follow `05_simulation::auc_entry_sd_<spread>_<method>_<field>`;
full results are in
[`simulation_auc_coverage_summary.csv`](../results/simulation_auc_coverage_summary.csv).

The strict subject bootstrap did not achieve nominal coverage uniformly: its
near-shared coverage was below 0.95 by multiple reported Monte Carlo
standard errors, while wide-entry coverage was compatible with 0.95. Therefore
the bootstrap remains **EMPIRICAL-ONLY** and `PARK`, not a validated replacement
for the missing analytic influence function.

## Time-scale selection diagnostic: killed

The seeded proposal was tested with subject-paired bootstrap delta AUC. At entry
SD 4, rejection probability was 0.045 (MCSE 0.00656); at SD 10 it was 0.046
(MCSE 0.00662). Mean naive-minus-strict delta AUCs were only 0.000273 and
0.000425 [TRACE:
`05_simulation::diagnostic_entry_sd_4_rejection_probability`,
`::diagnostic_entry_sd_4_rejection_mcse`,
`::diagnostic_entry_sd_10_rejection_probability`,
`::diagnostic_entry_sd_10_rejection_mcse`,
`::diagnostic_entry_sd_4_mean_delta_auc`,
`::diagnostic_entry_sd_10_mean_delta_auc`].

The test controlled rejection near its nominal level but had no useful signal
for the wrong eligibility rule in this DGP. Verdict: **KILL as a general
time-scale selector**. Discrimination can remain unchanged even when risk-set
construction is scientifically invalid.

## Claims and boundaries

- **PROVEN:** eligibility formulas, the scalar score/Hessian calculations, and
  the executed assertions are finite-sample code properties.
- **EMPIRICAL-ONLY:** all bias, SE, coverage, power, AUC, and runtime results.
- **CONJECTURED:** the subject bootstrap may be a serviceable interim inference
  method in some regimes; the grid refutes a claim of uniformly nominal
  coverage.
- **NOT PROVEN:** consistency, asymptotic normality, or a left-truncation
  Hájek/influence-function result.
