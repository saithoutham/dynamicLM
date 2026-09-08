# Innovation Lab Notebook

No idea is promoted on plausibility alone. Every verdict below must be backed by
an executable test and a `results/trace.csv` reference.

## IDEA-001: Correlation-aware time-scale selection

Date/commit: PENDING

Hypothesis (falsifiable): A correlation-aware summary delta-AUC test can select
the data-generating time scale with useful power while controlling type-I error.

Why it might matter: It would turn a scientific choice of time origin into a
pre-specified diagnostic using the package's own inferential framework.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: Awaiting the simulation harness and a validated left-truncation-aware
variance calculation.

## IDEA-002: Survey-weighted landmark supermodels

Date/commit: PENDING

Hypothesis (falsifiable): A design-weighted estimator with PSU/stratum-aware
variance changes inference relative to an unweighted landmark analysis in a
public linked-mortality cohort.

Why it might matter: National survey inference requires respecting unequal
selection probabilities and complex sampling.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: Requires a verified public linked-mortality download and a clearly
defined estimand; survey weights must not be bolted onto the existing IID code.

## IDEA-003: Covariate staleness adjustment

Date/commit: PENDING

Hypothesis (falsifiable): Adding time since last measurement improves summary
AUC or Brier score under informative visit timing without harming calibration.

Why it might matter: Last-observation-carried-forward quality depends on how old
the carried value is, and staleness may vary by age.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: Awaiting the PBC and NAFLD application datasets.

## IDEA-004: Inverse-intensity visit weighting

Date/commit: PENDING

Hypothesis (falsifiable): Estimated inverse visit-intensity weights reduce bias
from informative observation times more than a staleness covariate alone.

Why it might matter: Visit timing can be outcome-related and distort covariate
availability.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: This added complexity must beat the simpler staleness adjustment.

## IDEA-005: Trajectory-model evaluation harness

Date/commit: PENDING

Hypothesis (falsifiable): A model-agnostic prediction adapter can evaluate penLM
and trajectory encoders with identical landmark eligibility and summary metrics.

Why it might matter: Fair evaluation should isolate representation quality from
different cohort construction or censoring assumptions.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: Build only after the corrected estimand and data contract are stable.

## IDEA-006: Entry-risk-set shape diagnostic

Date/commit: PENDING

Hypothesis (falsifiable): Reporting entries, exits, and eligible counts by
landmark detects age-scale misuse more reliably than checking for monotone risk
sets, which need not be monotone under staggered entry.

Why it might matter: The current shared-origin mental model treats increasing
risk-set counts as suspicious even when they may be valid on attained age.

Test performed: PENDING

Result (with numbers + trace refs): PENDING

Verdict: PARK

Reasoning: Awaiting empirical risk-set inventories and minimal reproducers.
