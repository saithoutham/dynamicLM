# Landmark supermodels on the attained-age scale with staggered cohort entry

[AUTHOR: ...]

## Abstract

### Background

Landmark supermodels support repeated dynamic risk prediction by stacking risk
sets defined at several prediction times. Existing implementations naturally
encode a shared cohort origin. In cohorts without a meaningful common origin,
attained age can be the scientifically relevant time scale, but observation
then begins at subject-specific ages.

### Methods

We audited and extended an existing landmark-supermodel implementation to
enforce the age-scale eligibility condition `E <= a < X`, where `E` is entry
age, `a` is landmark age, and `X` is observed exit age. We implemented strict
eligibility at the landmark and within-window delayed entry, together with an
entry-conditional inverse-probability-of-censoring estimator. Known-truth
experiments compared unmodified, strict, and delayed construction under
independent entry, entry dependent on a measured predictor, and a factorial
shared-frailty mechanism. Every empirical result used saved seeds and a
number-level trace ledger.

### Results

**EMPIRICAL-ONLY:** independent staggered entry did not produce the proposed
degradation pattern in the executed grid. Under the strongest measured-entry
mechanism, pooled coefficient bias was 0.043250 for unmodified eligibility,
0.001762 for strict eligibility, and 0.002134 for delayed entry [TRACE:
`07_informative_entry::pooled_confirmation_naive_bias`;
`07_informative_entry::pooled_confirmation_strict_bias`;
`07_informative_entry::pooled_confirmation_delayed_bias`]. In the frailty
factorial, every method-specific Monte Carlo interval for the entry-dependence
contrast and its interaction with hazard frailty included zero, whereas every
hazard-frailty contrast excluded zero [TRACE:
`08_frailty_decomposition::{method}_{contrast}_{contrast_estimate,ci_lower,ci_upper}`].
Coefficient and discrimination results diverged: in one wide-entry
measured-entry cell, unmodified coefficient bias and strict-target summary-AUC
bias were 0.070566 and 0.020885, while the strict values were -0.000180 and
-0.000816 [TRACE: `07_informative_entry::focused_cell_003_naive_bias`;
`07_informative_entry::auc_cell_003_naive_bias`;
`07_informative_entry::focused_cell_003_strict_bias`;
`07_informative_entry::auc_cell_003_strict_bias`].

### Conclusions

Attained-age landmarking requires explicit observation-entry eligibility.
**EMPIRICAL-ONLY:** that correction removed most measured-entry coefficient
bias in the simulated mechanisms where the defect mattered, but it did not
resolve omitted-predictor misspecification or identify informativeness from
observed data. The analytic influence function for left-truncated summary
metrics remains open.

## Background

Dynamic risk prediction updates an individual's prognosis as follow-up and
covariate history accumulate. Landmark analysis performs this update at a
collection of landmark times by retaining people at risk, carrying their
available history to the landmark, and predicting over a fixed horizon.
Landmark supermodels pool these landmark-specific datasets and allow covariate
effects to vary smoothly with landmark time, avoiding unstable separate fits
late in follow-up (van Houwelingen and Putter, 2011; Nicolaie et al., 2013).
Penalized landmark supermodels extend this construction to larger predictor
sets and overlapping prediction windows (Fries et al., 2025).

The usual presentation assumes a meaningful shared time origin, such as a
diagnosis or enrollment date. That origin is arbitrary for a healthy cohort
assembled from longitudinal health records: a person's first recorded contact
need not mark biological onset or comparable observation density. Attained age
instead poses a clinically interpretable question, such as risk conditional on
history available by a specified age. It also changes eligibility. A person
whose record begins after landmark age `a` was not observable at `a`, even if a
later event age is present in the assembled cohort.

This work is framed as a methodological extension, not a software defect
report. Applying landmark supermodels to attained age introduces an eligibility
condition that shared-origin implementations do not need to encode. We ask
when violating that condition affects coefficient estimation or
discrimination, when it does not in known-truth experiments, and what a
corrected implementation does and does not solve.

## Methods

### General framework for landmark supermodels

Let `a` denote a landmark age, `w` the prediction horizon, `Z(a)` the covariate
history summarized at `a`, `E` the age at which observation begins, and `X` the
observed event or censoring age. A landmark supermodel stacks subject-landmark
records and fits a cause-specific Cox model whose coefficients may interact
with landmark age. Repeated contributions from the same person are handled by
a subject-clustered sandwich calculation. Overlapping windows allow one event
to contribute to more than one landmark record, as in penLM (Fries et al.,
2025).

### Attained-age eligibility and two entry modes

The strict risk set is

```text
R_strict(a) = {i : E_i <= a < X_i}.
```

Each retained interval begins at `a` and is administratively stopped at
`min(X_i, a + w)`. The alternative delayed-entry construction admits a person
who first becomes observable within the prediction window. Its counting-process
interval starts at `max(E_i, a)`, provided this is earlier than
`min(X_i, a + w)`. Strict and delayed modes therefore answer different
population questions; the latter is not merely a more efficient version of the
former.

The audited shared-origin stacker required only that exit follow the landmark.
On attained age, that admits rows with `a < E_i < X_i`. The corrected interface
adds explicit `entry` and `entry_mode` arguments and asserts the relevant
interval and subject-landmark invariants. Enrollment-scale calls retain their
original behavior; exact saved-seed regression tests found identical cohorts,
stacks, and fit objects in all 108 tested cells and 1,000 replicates per cell
[TRACE: `07_informative_entry::regression_cells`;
`07_informative_entry::regression_replicates_per_cell`;
`07_informative_entry::regression_exact_cohort_cells`;
`07_informative_entry::regression_exact_fit_cells`].

### Entry-conditional censoring weights and the open inference problem

For analysis entry `A_i` and analysis exit `Y_i`, we estimate the reverse
product-limit curve with counting-process records
`Surv(A_i, Y_i, censor_event_i)`. If `G` is the resulting censoring survival,
an observed case at time `T_i` receives the inverse weight
`G(A_i)/G(T_i-)`; a subject known to be event-free at horizon `h` receives
`G(A_i)/G(h)`. This is the finite algorithm implemented for cumulative/dynamic
AUC and cause-specific Brier score.

The left-truncation Hájek projection and cross-landmark influence function are
**NOT PROVEN / PENDING**. The implementation deliberately returns
`"PENDING: not derived for staggered entry"` where an analytic standard error
would be reported. Subject-level bootstrap intervals were evaluated as an
empirical procedure, not as a substitute for the missing derivation. Formal
work must extend time-dependent IPCW metric theory (for example, Blanche et
al., 2013) using left-truncation counting-process arguments and must account for
estimation of `G` and correlation across landmarks.

### Retention conditioning

The independent-entry experiment generates `X_cov` and `E` independently in a
source population and retains only people with event time `T > E`. Let
`A = {T > E}`. Bayes' rule gives, step by step,

```text
f(x,e | A)
  = f(x,e) Pr(A | x,e) / Pr(A)
  = f_X(x) f_E(e) Pr(T > e | X_cov=x) / Pr(A)
  = f_X(x) f_E(e) S_T(e | x) / Pr(A).
```

This identity is **PROVEN for the stated generator**. Unless the survival term
separates into functions of `x` and `e`, source-population independence does
not survive retention. Thus the independent-entry null result cannot be
explained by treating pre-entry rows as an unbiased sample; it remains an
empirical robustness or cancellation result.

Under the measured-entry generator, entry and failure are independent
conditional on the measured predictor in the source mechanism. Conditioning
on `E <= a` therefore leaves the specified post-landmark conditional hazard
unchanged for that mechanism. This finite conditional-hazard identity is
**PROVEN for the stated generator**; no Cox-estimator consistency or
asymptotic-normality result is claimed.

### Simulation studies

All simulations generated a measured standard-normal predictor, an entry age
from a bounded truncated-normal law, and an exponential proportional-hazards
event age, retaining only draws with event age after entry. Censoring was
generated after entry and calibrated separately by mechanism. Explicit seeds
were saved for every cell-replicate pair.

Regime A made entry independent in the source population. Its full grid crossed
three landmark spacings, four censoring targets, three sample sizes, and three
entry-age spreads; 108 cells were run with 1,000 replicates and three methods
[TRACE: `05_simulation::design_landmark_interval_1`--`_3`;
`05_simulation::design_censoring_target_1`--`_4`;
`05_simulation::design_sample_size_1`--`_3`;
`05_simulation::design_entry_sd_1`--`_3`;
`05_simulation::grid_cells`; `05_simulation::replicates_per_cell`;
`05_simulation::methods`].

Regime B shifted mean entry age by `gamma X_cov`, with larger risk values
entering earlier. A focused sweep crossed four dependence strengths and three
entry spreads at the anchor design; the strongest separation was then tested
through the full Regime A grid [TRACE:
`07_informative_entry::focused_cells`;
`07_informative_entry::selected_gamma`;
`07_informative_entry::confirmation_cells`]. Regime C crossed the presence or
absence of hazard frailty (`theta`) with the presence or absence of an entry
shift by the same latent variable (`delta`), over three entry spreads [TRACE:
`07_informative_entry::frailty_cells`].

For each cell and method we summarized coefficient bias, empirical dispersion,
mean estimated standard error, interval coverage, rejection of a zero
coefficient, stack rows, events, subjects, entry-predictor correlation, and
pre-entry fractions. Separate oracle-risk-score experiments summarized
landmark AUC with entry-conditional weights and a paired subject bootstrap.
Oracle scoring intentionally removes fitted-model error; in Regime C it uses a
latent variable unavailable to an analyst.

### Factorial frailty decomposition

No new cohort was generated and no model was refit for the decomposition. We
pooled the three entry-spread cells within each `theta`/`delta` arm, giving
3,000 replicate estimates per method-arm [TRACE:
`08_frailty_decomposition::arm_replicates`]. The two `delta` contrasts, two
`theta` contrasts, and their interaction were differences of arm-specific mean
biases. Arms used disjoint seed sets and were treated as independent [TRACE:
`08_frailty_decomposition::seed_overlap_count`]. Monte Carlo variances were the
sum of arm estimate variances divided by their replicate counts, with squared
contrast weights.

### Structural and exploratory diagnostics

The structural audit reports the fraction of unmodified landmark rows whose
entry follows the landmark. It also reports entry-predictor correlations and
Welch comparisons between the strict risk set and the disjoint pre-entry rows.
Only the pre-entry fraction directly audits eligibility; the other statistics
cannot detect dependence through unmeasured risk and can be distorted by
conditioning on survival to entry.

Installed public cohorts from the `survival` package were used to exercise the
attained-age machinery, including irregular longitudinal records in NAFLD and
PBC data and competing outcomes in MGUS data. These are implementation and
structural stress tests, not known-truth validation and not proxies for the
intended healthy pre-cancer electronic-record cohort.

## Results

### Source audit and public-data behavior

The original stack construction admitted pre-entry rows on attained age. In
the NAFLD age-50 risk set, 7,170 of 9,472 unmodified rows were pre-entry, and
2,976 of those rows retained complete analysis covariates [TRACE:
`03_break_diagnosis::nafld_lm50_stacked_n`;
`03_break_diagnosis::nafld_lm50_preentry_n`;
`03_break_diagnosis::nafld_lm50_analyzable_preentry_n`]. Thus missing
longitudinal history does not reliably remove invalid rows. Corrected strict
and delayed fits executed on all public cohorts studied, but their differences
cannot be interpreted as bias because no population truth is known.

### Independent and measured predictor-dependent entry

**EMPIRICAL-ONLY:** Regime A falsified the prespecified degradation hypothesis.
Unmodified, strict, and delayed estimates had similar coefficient behavior
across increasing entry spread; the targeted AUC experiment was likewise
similar. This negative result limits the claim: incorrect eligibility is not
sufficient by itself to guarantee detectable inferential distortion under
every entry mechanism.

**EMPIRICAL-ONLY:** Regime B made the defect consequential. At the strongest
tested dependence and widest entry spread, coefficient bias was 0.070566 for
unmodified construction, -0.000180 for strict eligibility, and -0.000697 for
delayed entry [TRACE:
`07_informative_entry::focused_cell_003_naive_bias`;
`07_informative_entry::focused_cell_003_strict_bias`;
`07_informative_entry::focused_cell_003_delayed_bias`]. Across the full
confirmation grid, the corresponding pooled biases were 0.043250, 0.001762,
and 0.002134 [TRACE:
`07_informative_entry::pooled_confirmation_naive_bias`;
`07_informative_entry::pooled_confirmation_strict_bias`;
`07_informative_entry::pooled_confirmation_delayed_bias`]. These are simulation
results for the executed data-generating process, not general estimator
properties.

### Factorial shared-frailty decomposition

**EMPIRICAL-ONLY:** the entry-dependence increment was small at Monte Carlo
resolution. For strict eligibility, the `delta` contrast was -0.002044 when
`theta` was absent and -0.000546 when it was present; the interaction was
0.001497. Their Monte Carlo intervals all included zero [TRACE:
`08_frailty_decomposition::strict_delta_at_theta_0_contrast_estimate`;
`08_frailty_decomposition::strict_delta_at_theta_0p5_contrast_estimate`;
`08_frailty_decomposition::strict_theta_by_delta_interaction_contrast_estimate`;
matching `_ci_lower` and `_ci_upper` trace fields]. The strict `theta` effects
were -0.053755 and -0.052258 at the two `delta` levels, and both intervals
excluded zero [TRACE:
`08_frailty_decomposition::strict_theta_at_delta_0_contrast_estimate`;
`08_frailty_decomposition::strict_theta_at_delta_m2_contrast_estimate`;
matching `_ci_lower` and `_ci_upper` trace fields]. Naive and delayed methods
showed the same interval pattern [TRACE:
`08_frailty_decomposition::{naive,delayed}_{contrast}_{field}`].

The Phase 6 description that the design could not identify a truncation
increment was too conservative. The factorial `delta` comparison identifies
that increment within the simulation. The measured bias is driven by adding an
omitted hazard frailty through `theta`, rather than by a detected additional
entry-dependence effect. This does not establish absence of such an effect in
other dependent-entry mechanisms.

### Coefficient estimation and discrimination can diverge

**EMPIRICAL-ONLY:** in the wide measured-entry cell, unmodified coefficient
bias and strict-target summary-AUC bias were 0.070566 and 0.020885; strict
values were -0.000180 and -0.000816 [TRACE:
`07_informative_entry::focused_cell_003_naive_bias`;
`07_informative_entry::auc_cell_003_naive_bias`;
`07_informative_entry::focused_cell_003_strict_bias`;
`07_informative_entry::auc_cell_003_strict_bias`]. Delayed coefficient bias was
-0.000697, while its difference from the strict-population AUC target was
0.029105 [TRACE:
`07_informative_entry::focused_cell_003_delayed_bias`;
`07_informative_entry::auc_cell_003_delayed_bias`]. The delayed AUC comparison
contains a deliberate estimand mismatch because delayed entry changes both the
eligible population and interval.

Under wide shared frailty, strict one-covariate coefficient bias was -0.053566
while strict oracle-score AUC bias was 0.000031 [TRACE:
`07_informative_entry::frailty_cell_009_strict_bias`;
`07_informative_entry::auc_cell_021_strict_bias`]. The oracle score contains
the latent frailty, so this is a decomposition of model fitting from metric
construction, not attainable prediction performance.

### Eligibility diagnostic and negative innovation results

The pre-entry fraction detected structural disagreement whenever unmodified
and strict risk sets differed. It does not test whether entry is informative.
Observed entry-predictor tests detected strong measured dependence but had
elevated rejections after retention under a source-independent control and low
sensitivity to latent entry dependence. The structural audit was retained; the
general correlation screen was rejected.

The proposed change in summary AUC was also rejected as a time-scale selector:
discrimination could remain stable despite invalid eligibility and could
change because strict and delayed modes target different populations. In
subject-disjoint public-data experiments, neither an additive measurement-
staleness term nor a simple inverse-visit-rate weighting approximation improved
discrimination in both datasets. The staleness delta AUCs were 0.003890 for
NAFLD and -0.002384 for PBC [TRACE:
`06_application::nafld_staleness_delta_auc`;
`06_application::pbcseq_staleness_delta_auc`]. The inverse-rate delta AUCs were
-0.022056 and 0.001217 [TRACE:
`06_application::nafld_iiw_delta_auc`;
`06_application::pbcseq_iiw_delta_auc`]. These negative results reject the
implemented shortcuts, not all models of observation timing.

## Discussion

The central requirement is simple: attained-age risk sets must respect when a
person became observable. A shared-origin implementation can omit entry because
all participants are observed from the common landmarking origin. Reusing that
logic on age silently assigns event and covariate information to prediction
times before observation. Explicit strict or counting-process delayed entry
makes the scientific population visible in the model contract.

The simulations show why code correctness and inferential consequence must be
separated. Independent staggered entry did not yield the proposed worsening
pattern, even though invalid rows were present. Measured risk-dependent entry
did produce coefficient and AUC distortion, and strict eligibility removed
most of it in that generator. Neither result should be universalized. The
retention identity explains why source independence alone is not a proof of
harmlessness, while the Regime B conditional-hazard calculation explains why
strict eligibility is coherent for that particular mechanism.

The frailty decomposition further narrows the claim. Bias appeared when the
latent variable entered the hazard even if it did not enter the entry model.
Adding latent entry dependence produced no resolved increment in these runs.
Eligibility construction was therefore not the source of the measured frailty
bias. This is an omitted-predictor and marginal-mixture warning, not evidence
that the correction rescues all dependent truncation.

Several attractive diagnostics or quick fixes failed. A delta-AUC selector
confused predictive ranking with scientific time origin. A general entry-
covariate screen could neither control retention-induced association nor see an
unmeasured common cause. A linear staleness term and a Poisson visit-rate
approximation did not provide portable improvements. Reporting these failures
is part of the contribution: each limits a tempting but unsupported workflow.

Strict entry is the clearer default when the target is risk among people
already observed at the landmark. Delayed entry can use more within-window
observation but changes the target population and prediction interval. It
should be prespecified rather than selected because it produces smaller
standard errors or a favorable apparent score.

## Limitations

The analytic left-truncation Hájek influence function is **NOT PROVEN /
PENDING**. Consequently, the proposed implementation does not supply analytic
summary-AUC or Brier intervals. The subject bootstrap is not a replacement: its
empirical coverage was 0.930 in one independent-entry condition and 0.845 for
delayed entry in one wide measured-entry condition [TRACE:
`05_simulation::auc_entry_sd_0p25_naive_coverage`;
`07_informative_entry::auc_cell_003_delayed_coverage`]. A formal derivation must
be reviewed independently and tested separately from the derivation.

No proof of consistency, asymptotic normality, or target coverage is supplied
for strict or delayed supermodels. The known-truth generator is scalar and
parametric; it omits nonlinear predictors, time-varying effects, competing
events, measurement error, and realistic observation processes. The Regime C
result covers one normal frailty mechanism and does not establish general
robustness to dependent truncation.

The delayed AUC experiments use a strict-population target and therefore do not
estimate a delayed-population AUC bias. Regime C oracle scores include an
unmeasured frailty and cannot be attained by an analyst. Public datasets verify
execution and reveal structural risk-set behavior, but they do not represent a
healthy pre-cancer electronic-record cohort and have no known inferential
truth. NHANES supplies only a single examination and therefore cannot validate
within-person covariate updating.

The structural pre-entry fraction detects eligibility mismatch, not
informativeness. Observed entry-predictor associations can be induced by the
retention rule, and absence of such an association cannot rule out latent
dependence. Survey-design variance for overlapping landmark stacks and summary
metrics also remains undeveloped.

## Conclusions

Attained age is a useful landmarking scale when no meaningful common origin
exists, but it requires the explicit eligibility condition `E <= a < X` or a
prespecified within-window delayed-entry analogue. **EMPIRICAL-ONLY:** enforcing
that condition removed most bias in the tested measured-entry mechanism and
did not create a claimed benefit where the independent-entry experiment was
null. It does not solve omitted-predictor bias or supply inference for
left-truncated summary metrics. A structural pre-entry audit should accompany
every age-scale analysis, while scientific time-scale selection remains a
design decision.

## Availability of data and materials

All data used here are public. The implementation and analysis are available
on the `age-scale` branch of `github.com/saithoutham/dynamicLM`. The Phase 7
analysis commit is `[COMMIT: PENDING UNTIL FINAL PHASE 7 COMMIT]`. Cell-level
seed manifests and `results/trace.csv` provide deterministic inputs and
number-level provenance. Large regenerateable raw result tables are governed by
the verified policy documented in `INNOVATIONS.md` and `README.md`.

## Abbreviations

AUC: area under the receiver operating characteristic curve; IPCW: inverse
probability of censoring weighting; NAFLD: nonalcoholic fatty liver disease;
PBC: primary biliary cholangitis.

## References

Blanche P, Dartigues J-F, Jacqmin-Gadda H. Estimating and comparing
time-dependent areas under receiver operating characteristic curves for
censored event times with competing risks. *Statistics in Medicine*. 2013.

Fries AH, Choi E, Wu JT, et al. Software Application Profile: dynamicLM.
*International Journal of Epidemiology*. 2023; dyad122.

Fries AH, Choi E, Han SS. Penalized landmark supermodels. *BMC Medical Research
Methodology*. 2025;25:22. doi:10.1186/s12874-024-02418-9.

Nicolaie MA, van Houwelingen HC, de Witte TM, Putter H. Dynamic prediction in
competing risks models. *Statistics in Medicine*. 2013.

van Houwelingen HC, Putter H. *Dynamic Prediction in Clinical Survival
Analysis*. Chapman & Hall/CRC; 2011.
