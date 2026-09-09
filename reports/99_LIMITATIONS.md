# Limitations

## The central variance theory is missing

The project does **not** provide a proved influence function, Hájek projection,
or analytic standard error for summary AUC/Brier under staggered entry. The
implemented entry-conditional IPCW weights passed finite checks and special-
case comparisons, but that is not an asymptotic derivation. The exposed return
value remains `"PENDING: not derived for staggered entry"`.

The paired subject bootstrap is not a validated substitute. In Phase 4, one
condition had coverage 0.930 with Monte Carlo standard error 0.00807 [TRACE:
`05_simulation::auc_entry_sd_0p25_naive_coverage`,
`::auc_entry_sd_0p25_naive_coverage_mcse`]. Any summary-metric interval in this
repository must be read as **EMPIRICAL-ONLY**.

## The simulations do not prove estimator validity

All bias, empirical SE, estimated SE, coverage, power, and AUC findings are
Monte Carlo results under specified generators. None proves consistency or
asymptotic normality. The focused and full-factorial grids are broad in sample
size, censoring, landmark spacing, entry spread, and entry dependence, but they
still use an exponential proportional-hazards event mechanism with a single
normally distributed measured covariate. They omit nonlinear effects,
time-varying coefficients, multiple measured predictors, competing events,
model selection, and covariate measurement error.

The exact gamma-zero regression check proves software equivalence to the Phase
4 engine for saved seeds. It does not prove that either engine is scientifically
correct outside the separately validated scalar setting.

## The frailty experiment is deliberately confounded

Regime C makes `U` affect both event risk and, when `delta != 0`, entry. The
fitted coefficient model omits `U`. Consequently, its bias can reflect marginal
frailty mixing, selection through entry, or both. Similar bias when `delta = 0`
and `delta = -2` means the experiment cannot identify a dependent-truncation
increment. It supports only the narrower conclusion that eligibility correction
alone is insufficient under this shared-frailty stress test.

The summary-AUC analysis uses the oracle generating score, including frailty in
Regime C. That isolates metric construction from fitted-model error but is not
available to an analyst with unmeasured `U`. Oracle discrimination results must
not be presented as attainable prediction performance.

The Phase 6 AUC study used the oracle generating score. In the wide shared-
frailty cell, strict oracle AUC bias was 0.000031 even though strict fitted-
coefficient bias was -0.053566 [TRACE:
`07_informative_entry::auc_cell_021_strict_bias`,
`::frailty_cell_009_strict_bias`]. This is a useful separation of discrimination
from coefficient estimation, but oracle access to `U` makes the AUC result
unattainable in the stated analyst setting. It must not be used to claim that an
observed-predictor model is unbiased or well calibrated.

The delayed oracle AUC was compared with a strict-population target. Its large
wide-entry differences—0.029105 under `gamma = -4` and 0.018159 in the shared-
frailty cell—partly encode a deliberate estimand mismatch [TRACE:
`07_informative_entry::auc_cell_003_delayed_bias`,
`::auc_cell_021_delayed_bias`]. They are not estimates of delayed-method bias
for a delayed-entry target.

## Strict and delayed entry do not target identical populations

Strict landmarking conditions on observation by the landmark; delayed entry
admits people who enter later within the prediction window. Delayed entry can
use more rows, but its risk-set composition and prediction interval differ.
Lower SE or similar bias does not by itself establish superior efficiency at a
common estimand. Any comparison must state the target population and interval.
Phase 6 AUC differences for naive or delayed rows are reported against the
strict scientific target and are therefore target differences, not necessarily
bias for each method's own estimand.

## Entry diagnostics do not establish informativeness

A nonzero pre-entry fraction proves that naive and strict eligibility differ;
it does not prove informative entry or inferential bias. Conversely, zero
observed entry--covariate correlation does not exclude entry dependence through
unmeasured risk. Conditioning the source population on `T > E` can itself
associate entry and measured risk, so ordinary correlation p-values are not
calibrated tests of the source entry mechanism.

In the Phase 6 wide-entry gamma-zero control, the entry--`X` screen rejected in
0.156 of datasets; in the wide-entry shared-frailty condition it rejected in
only 0.065 [TRACE:
`07_informative_entry::diagnostic_cell_012_lm1_entry_x_detection_probability`,
`::diagnostic_cell_021_lm1_entry_x_detection_probability`]. The structural
audit is retained; a general informative-entry diagnostic is not.

## Public data are machinery tests, not the intended EHR validation

The public cohorts do not reproduce the lab's healthy pre-cancer EHR target.
NAFLD and PBC provide irregular longitudinal measurements but different
disease, referral, and follow-up processes. `mgus2` and `flchain` are useful
age-scale and competing-risk checks with limited longitudinal updating. NHANES
is a single-exam survey: it tests age-scale and survey-weighted mechanics, not
within-person covariate trajectories. No restricted EHR, SEER, Medicare, MHOS,
PLCO, NLST, UK Biobank, All of Us, or MIMIC data were acquired or simulated as
if observed.

The NHANES mortality linkage covers the executed public cycle and variables;
it does not establish transportability to later cycles or a national dynamic-
prediction estimand. Survey-weighted overlapping-landmark variance remains
theoretical work, not a supported feature.

## Censoring, visits, and missingness remain assumptions

Simulation censoring is independent after entry. Real EHR loss to follow-up can
depend on latent health, insurance, migration, and care intensity. Correct
left-truncated Kaplan--Meier mechanics do not repair dependent censoring.

Last-observation-carried-forward remains vulnerable to informative visit timing
and covariate staleness. The tested additive staleness term and approximate
inverse-intensity weighting did not improve the public-data AUC comparisons;
those negative tests kill the specific fixes, not the problem. Missing
pre-entry covariates sometimes mask invalid rows and sometimes do not, so
missingness cannot serve as an eligibility rule.

## Multiplicity and selection affect interpretation

The full confirmation gamma was selected using the largest focused
naive-versus-strict bias separation. The confirmation grid then reused that
mechanism across new factor combinations, but it is not an independent
mechanism-discovery sample. Cell-level Monte Carlo errors quantify simulation
noise, not uncertainty due to selecting gamma. The many diagnostic tests and
simulation summaries are descriptive; no family-wise claim is made beyond the
within-dataset Bonferroni adjustment used for landmark risk-set shifts.

## Numerical and tolerance history is not pristine

Phase 4 loosened three tolerances during development. Phase 6 loosened an exact
cross-process AUC truth-table comparison to `1e-12` after observing a maximum
difference of 1.113554e-13 [TRACE:
`07_informative_entry::tolerance_4_new`,
`::auc_truth_reproduction_max_abs_difference`]. They are documented in
[`tolerance_changes.csv`](../results/tolerance_changes.csv), but the exact
original failure magnitudes for censoring validation and AUC bounds were not
retained and remain `PENDING` [TRACE:
`07_informative_entry::phase6_tolerance_changes`]. This limits retrospective
auditability even though the final values and justifications are explicit.

Deterministic AUC truth is numerical quadrature, not closed-form truth. Fine--
coarse grid agreement checks discretization sensitivity only for the grids
used; it is not a proof of exact integration.

## Package and reproducibility limitations

The tarball `R CMD check --as-cran` completed but was not clean: 1 WARNING and
5 NOTEs [TRACE: `07_informative_entry::cran_check_warnings`,
`::cran_check_notes`]. Tests passed inside the check, but a data-documentation
WARNING, packaging/metadata and namespace NOTES, slow examples, and an
unverifiable current time remain. Passing tests must not be summarized as
passing CRAN checks.

The large grids are reproducible from explicit seed manifests, but wall time and
floating-point results can vary by R, BLAS, compiler, and platform. Exact
Regime A identity is guaranteed only for the recorded environment and compared
commit. `results/session.txt` records that environment.
