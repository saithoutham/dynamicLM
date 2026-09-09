# Figure plan

No figures are generated in Phase 7. Every panel below has a named input file
and should be produced by a later scripted analysis with assertions and trace
entries before manuscript submission.

## Figure 1 — Why attained age changes eligibility

- **Form:** three-panel schematic, not a data plot.
- **Panels:** shared-origin risk set; attained-age strict eligibility
  (`E <= a < X`); within-window delayed entry beginning at `max(E, a)`.
- **Purpose:** make the population difference between strict and delayed entry
  visible before presenting results.
- **Input:** definitions implemented in `R/stack_data.R` and asserted in
  `tests/testthat/test-age-scale.R`; no numerical result file.

## Figure 2 — Coefficient bias under measured predictor-dependent entry

- **Form:** coefficient-bias point/interval plot faceted by entry-age spread,
  with unmodified, strict, and delayed methods aligned within each mechanism.
- **Primary input:** `results/informative_entry_summary.csv` for the focused
  sweep.
- **Confirmation input:** `results/informative_entry_pooled_summary.csv` for a
  compact full-grid inset.
- **Purpose:** show that method separation grows with the executed dependence
  mechanism and that the independent-entry control remains visible.
- **Trace family:** `07_informative_entry::focused_cell_*` and
  `07_informative_entry::pooled_confirmation_*`.

## Figure 3 — Frailty factorial decomposition

- **Form:** forest plot of `delta`, `theta`, and interaction contrasts, faceted
  by method, with a zero reference line.
- **Input:** `results/frailty_decomposition.csv`.
- **Purpose:** distinguish the small entry-dependence contrasts from the
  omitted-hazard-frailty contrasts and display Monte Carlo precision directly.
- **Trace family:** `08_frailty_decomposition::{method}_{contrast}_{field}`.

## Figure 4 — Coefficient bias versus strict-target summary-AUC bias

- **Form:** paired outcome plot with coefficient bias and AUC bias on separate
  aligned panels; do not combine them on one numerical axis.
- **Inputs:** `results/informative_entry_summary.csv`,
  `results/informative_auc_coverage_summary.csv`, and
  `results/prediction_estimation_comparison.csv`.
- **Purpose:** show that an eligibility mechanism may affect coefficient
  estimation and ranking differently, and flag the delayed-entry target
  mismatch in the caption.
- **Trace family:** `07_informative_entry::focused_cell_*`,
  `::frailty_cell_*`, and `::auc_cell_*`.

## Supplementary Figure S1 — Structural eligibility diagnostic

- **Form:** pre-entry fraction by landmark and mechanism, with separate panels
  for entry-predictor correlation and risk-set-shift rejection frequency.
- **Input:** `results/entry_diagnostic_summary.csv`.
- **Purpose:** demonstrate why the pre-entry fraction is retained as a
  structural audit while the correlation screen is not an informativeness
  test.
- **Trace family:** `07_informative_entry::diagnostic_cell_*`.

## Supplementary Figure S2 — Public-data risk-set shapes

- **Form:** connected counts across attained-age landmarks, faceted by public
  dataset and strict/delayed mode.
- **Input:** `results/application_risk_sets.csv`.
- **Purpose:** show that valid attained-age risk-set counts can increase as new
  people enter and that monotonicity is not a valid age-scale assertion.
- **Trace family:** `06_application::<dataset>_<mode>_<landmark>_<cause>_rows`.

## Supplementary Figure S3 — Negative visit-timing experiments

- **Form:** delta-AUC forest plot for additive staleness and inverse visit-rate
  variants in NAFLD and PBC.
- **Input:** `results/application_irregular_visits.csv`.
- **Purpose:** report killed ideas with the same visual prominence as positive
  simulation findings.
- **Trace family:** `06_application::<dataset>_<variant>_<field>`.
