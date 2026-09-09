# Phase 5 — public-data application

## Outcome

The corrected stacker ran on all four installed public datasets and on a
downloaded public-use NHANES linked-mortality cycle. The application confirms
the structural point that motivated the patch: age-scale risk sets can grow
when new subjects enter. It does **not** turn differences between strict and
delayed estimates into bias estimates, because these datasets do not provide a
known truth.

Two seeded visit-timing ideas were tested and killed in their implemented
forms. An additive staleness covariate did not improve out-of-fold summary AUC
in NAFLD or PBC. A Poisson inverse-visit-rate approximation also failed to
improve PBC and moved NAFLD discrimination in the wrong direction.

## Data and eligibility contract

The installed `survival` datasets were read directly from the recorded R
environment:

- `nafld1` supplied entry age, follow-up, death, sex, and BMI; post-index
  `nafld2` HDL visits supplied the longitudinal covariate.
- `pbcseq` supplied irregular bilirubin and albumin measurements. Death was the
  event and transplant was treated as non-death censoring for these
  application fits.
- `mgus2` supplied time to plasma-cell progression and time to death in months;
  the first event defined the competing-risks outcome.
- `flchain` supplied baseline free-light-chain measurements and mortality.

The variable definitions were checked against the installed help files rather
than inferred from the mission brief [SOURCE: the
[`survival` package manual](https://spout.ussg.indiana.edu/CRAN/web/packages/survival/survival.pdf)].

For strict landmarking, every row was asserted to satisfy
`entry_age <= landmark < exit_age`. For delayed landmarking, every row was
asserted to start at `max(landmark, entry_age)` and to have positive interval
length. Subject-landmark keys were unique. Model inputs with unavailable
longitudinal values were removed explicitly rather than relying on `coxph()` to
drop them silently.

Selected risk-set counts show why monotonicity is not a valid age-scale check:

| dataset | mode | age landmarks | eligible rows at successive landmarks |
|---|---|---|---|
| NAFLD | strict | 50, 60, 70 | 2,145; 2,080; 1,518 |
| NAFLD | delayed | 50, 60, 70 | 4,705; 4,017; 2,529 |
| PBC | strict | 40, 50, 60 | 38; 67; 62 |
| PBC | delayed | 40, 50, 60 | 86; 116; 92 |
| MGUS | strict | 60, 70, 80, 90 | 160; 319; 393; 148 |
| FLCHAIN | strict | 60, 70, 80, 90 | 2,403; 2,212; 1,464; 423 |

Trace refs follow
`06_application::<dataset>_<mode>_<landmark>_<cause>_rows`; the complete table,
including event counts and delayed-entry results, is in
[`application_risk_sets.csv`](../results/application_risk_sets.csv).
The PBC enrollment-scale risk sets, by contrast, decreased from 312 to 278 to
225 across its shared-origin landmarks [TRACE:
`06_application::pbcseq_enrollment_shared_0_1_rows`,
`::pbcseq_enrollment_shared_2_1_rows`,
`::pbcseq_enrollment_shared_4_1_rows`].

## Corrected Cox fits

All estimates below use Breslow ties and subject-clustered sandwich standard
errors. Landmark interactions used centered, ten-year-scaled attained age
[TRACE: `06_application::landmark_interaction_scale_years`].
They are descriptive application results, not causal effects or bias
estimates.

| dataset/outcome | mode | selected standardized term | estimate | robust SE |
|---|---|---|---:|---:|
| NAFLD death | strict | HDL | -0.298886 | 0.136155 |
| NAFLD death | delayed | HDL | -0.299282 | 0.135639 |
| PBC death | strict | log bilirubin | 1.526521 | 0.239153 |
| PBC death | delayed | log bilirubin | 1.249496 | 0.158983 |
| MGUS progression | strict | M-spike | 0.552080 | 0.157939 |
| MGUS progression | delayed | M-spike | 0.449481 | 0.140803 |
| FLCHAIN death | strict | log lambda | 0.269715 | 0.059964 |
| FLCHAIN death | delayed | log lambda | 0.294898 | 0.052595 |

Trace refs follow
`06_application::<dataset>_<mode>_coxph_cluster_<cause>_<term>_estimate` and
`..._standard_error`; all terms, events, observations, and unique-subject
counts are in [`application_models.csv`](../results/application_models.csv).

The PBC enrollment-scale log-bilirubin estimate was 1.035086 (SE 0.095648)
[TRACE:
`06_application::pbcseq_enrollment_shared_coxph_cluster_1_log_bili_z_estimate`,
`::pbcseq_enrollment_shared_coxph_cluster_1_log_bili_z_standard_error`]. Its
attained-age strict estimate is shown above. This contrast is **not** a
time-scale selection test: landmark populations and risk sets differ.

Delayed entry added little analyzable NAFLD longitudinal information: the
strict and delayed fits used 4,467 and 4,492 rows, respectively [TRACE:
`06_application::nafld_strict_coxph_cluster_1_hdl_z_observations`,
`::nafld_delayed_coxph_cluster_1_hdl_z_observations`]. Many otherwise eligible
within-window entrants had no HDL value available at their analysis entry. This
is a practical constraint on the theoretically more inclusive design.

## Informative visit timing experiments

Subjects, rather than rows, were split across five folds [TRACE:
`06_application::cv_folds`]. Each test subject was absent from its training
folds. Summary AUC used the corrected
entry-conditional censoring weights and a paired subject bootstrap shared over
landmarks. These are **EMPIRICAL-ONLY** uncertainty estimates; they are not the
unproved analytic influence function.

The staleness term was years since the carried measurement. Median and
90th-percentile staleness were 0.713 and 2.984 years in NAFLD, and 0.593 and
3.309 years in PBC [TRACE:
`06_application::nafld_staleness_median_staleness`,
`::nafld_staleness_p90_staleness`,
`::pbcseq_staleness_median_staleness`,
`::pbcseq_staleness_p90_staleness`].

| dataset | variant | base AUC | variant AUC | delta AUC | bootstrap SE | p-value |
|---|---|---:|---:|---:|---:|---:|
| NAFLD | staleness | 0.599327 | 0.603217 | 0.003890 | 0.012467 | 0.7550 |
| NAFLD | inverse visit rate | 0.599327 | 0.577272 | -0.022056 | 0.012824 | 0.08545 |
| PBC | staleness | 0.876120 | 0.873736 | -0.002384 | 0.013043 | 0.8550 |
| PBC | inverse visit rate | 0.876120 | 0.877337 | 0.001217 | 0.002319 | 0.5997 |

Trace refs follow `06_application::<dataset>_<variant>_<field>` and exact
values are in
[`application_irregular_visits.csv`](../results/application_irregular_visits.csv).
The inverse-rate implementation was deliberately modest: a training-fold
Poisson visit model used sex, entry age, and landmark, after which inverse
predicted rates were upper-tail capped and normalized to mean one. It is an
approximation, not a derived recurrent-event weighting estimator.

The PBC staleness and weighting comparisons had 993 and 994 evaluable
bootstrap resamples; seven and six lacked cases or controls at at least one
landmark [TRACE: `06_application::pbcseq_staleness_bootstrap_successful`,
`::pbcseq_iiw_bootstrap_successful`,
`::pbcseq_staleness_bootstrap_failed`,
`::pbcseq_iiw_bootstrap_failed`]. No failed resample was silently replaced.

**Verdicts:** the additive staleness proposal is `KILL` for this specification,
and the Poisson inverse-rate approximation is `KILL`. Neither result rules out
a prespecified model for nonlinear staleness or a rigorously estimated
recurrent-visit intensity process.

## A deliberately retained evaluation failure

The initial PBC design used a ten-year prediction window. At the age-40
landmark, 38 analytic rows contained seven events but only one horizon control
[TRACE: `06_application::40_rows`, `::40_events`,
`::40_horizon_controls`]. With a fixed 100-resample minimal reproducer, only 49
resamples were evaluable, below the predeclared minimum of 90 [TRACE:
`06_application::pbc_w10_bootstrap_attempted`,
`::pbc_w10_bootstrap_successful`,
`::pbc_w10_bootstrap_minimum_required`]. The captured error was:

```text
Fewer than 90% of bootstrap resamples were evaluable: 49/100; Fast censoring product-limit estimate reached zero. | AUC requires cases and controls.
```

The production PBC comparison therefore uses a five-year window, which also
matches its enrollment-scale comparison [TRACE:
`06_application::pbc_application_window_years`]. The reproducer and exact
condition are retained in
[`application_pbc_w10_design_failure.csv`](../results/application_pbc_w10_design_failure.csv)
and [`application_conditions.csv`](../results/application_conditions.csv).

## NHANES and NHIS public-use attempt

Network access was available. The official public-use NHANES 1999--2000
mortality and examination files were downloaded, checksummed, and linked by
`SEQN`. The analytic merge contained 4,973 eligible adults and 1,454 deaths
[TRACE: `06_application::nhanes_analysis_rows`,
`::nhanes_analysis_events`]. The mortality file, examination file, and NHIS
mortality file all downloaded successfully; their exact byte counts and MD5
hashes are in
[`application_downloads.csv`](../results/application_downloads.csv).

The current CDC page distinguishes restricted linked mortality updated through
the end of 2022 from earlier public-use files with mortality follow-up through
2019. It also warns that selected public-use follow-up and cause-of-death values
are perturbed or synthetic, while mortality status is not altered [SOURCE:
[CDC Linked Mortality Files](https://www.cdc.gov/nchs/linked-data/mortality-files/index.html);
[CDC public-use file description](https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf)].
Those limitations apply to this application.

Age-scale NHANES strict risk sets contained 1,473, 1,227, 1,189, and 907
subjects at ages 50, 60, 70, and 80; delayed entry contained 1,969, 1,882,
1,636, and 1,178 [TRACE refs beginning
`06_application::nhanes_1999_2000_strict_` and
`::nhanes_1999_2000_delayed_`].

For the strict stack, the unweighted male coefficient was 0.525765 (SE
0.086943); the nested survey-design fit using MEC weights, strata, survey PSU,
and participant clusters gave 0.481217 (SE 0.075297) [TRACE:
`06_application::nhanes_1999_2000_strict_coxph_cluster_1_male_estimate`,
`::nhanes_1999_2000_strict_coxph_cluster_1_male_standard_error`,
`::nhanes_1999_2000_strict_svycoxph_1_male_estimate`,
`::nhanes_1999_2000_strict_svycoxph_1_male_standard_error`].

This establishes that the design-weighted fitting path executes and changes
the empirical estimate. It does **not** prove valid design-based inference for
overlapping landmark rows, and no survey-weighted summary AUC/Brier CI is
claimed. That influence-function problem remains `PENDING`; IDEA-002 stays
`PARK`.

The NHIS 2018 mortality file was also downloaded successfully [TRACE:
`06_application::nhis_2018_mort_2019_public_dat_success`], but no NHIS model is
reported. `PENDING`: obtain and verify the matching public survey person file,
merge its age/covariate/weight/stratum/PSU fields, and then repeat the asserted
NHANES workflow. A mortality link file by itself is not a valid application
dataset.

NHANES and NHIS are single-exam sources. They test attained-age stacking and
survey machinery, not within-person covariate updating; only the NAFLD and PBC
analyses above exercise that feature.

## Claim status

- **PROVEN:** the executed finite-sample data assertions, eligibility
  identities, unique keys, positive counting-process intervals, and seed/file
  manifests.
- **EMPIRICAL-ONLY:** all coefficients, standard errors, AUC comparisons,
  bootstrap results, and survey-weighted application results.
- **CONJECTURED:** a fully specified recurrent-visit-process weighting method
  could behave differently from the killed approximation.
- **NOT PROVEN / PENDING:** the left-truncation Hájek influence function and the
  design-based influence function for survey-weighted summary metrics.
