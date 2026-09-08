# Reconnaissance

## Bottom line

The left-truncation hypothesis is **confirmed for cohort construction**, with an
important qualification: the current source already uses each landmark as a
counting-process entry time during model fitting. What it cannot represent is a
subject-specific observation-entry age distinct from the landmark. The current
stacker selects `exit > landmark` and never checks `observation entry <=
landmark`. A deterministic two-subject reproducer therefore returns two rows
when only one subject is observed at the landmark; the extra row is a pre-entry
subject (TRACE `01_recon::toy_current_rows`, `01_recon::toy_correct_rows`, and
`01_recon::toy_preentry_rows`).

This is not one isolated defect. The installed development source of
`riskRegression::Score()` also removes entry from the censoring-model response,
and `dynamicLM::summary_metric()` assumes the first landmark contains every ID.
The latter fails with the exact error `Subsample size mismatch.` in a minimal
staggered-entry example. Exact captured conditions are in
`results/recon_errors.csv`.

## Source versions inspected

- `dynamicLM` upstream commit `e9d6feaeabcf35fdc32850de47ef35be2bc2d4cc`
  on the local `age-scale` branch. [SOURCE: local Git commit]
- `riskRegression` development commit
  `b3eea672eab976daf60585cf690d32cb50cdade6`, installed package version
  `2026.05.21`. [SOURCE: installed package `DESCRIPTION` fields `RemoteSha` and
  `Version`]
- `survival` package version `3.8-3`. [SOURCE: `results/session.txt`]
- `dynpred` is PENDING: CRAN did not provide a build for the installed R
  release. It is not an imported dependency of the inspected `dynamicLM`
  source; `get_lm_data()` only credits its earlier `cutLM()` implementation.

The GitHub README requires the development `riskRegression` build, consistent
with the installation route used here.[^1]

## Function-by-function audit

### `stack_data()` and `get_lm_data()`

`stack_data()` loops over landmarks and delegates every slice to
`get_lm_data()` (`R/stack_data.R`, lines 145–168). The decisive eligibility
line is:

```r
lmdata <- lmdata[lmdata[[outcome$time]] > lm, ]
```

[SOURCE: local `R/get_lm_data.R`, line 101]

There is no entry-time argument or filter. In long format, the function locates
the most recent row at or before the landmark. If none exists, it copies the
first row and sets only the time-varying fields and running time to `NA`
(`R/get_lm_data.R`, lines 85–96). Consequently:

- with any requested varying covariate, downstream complete-case behavior can
  hide the bad eligibility row by dropping it;
- with fixed covariates only, the pre-entry row remains fully analyzable and is
  silently included;
- administrative censoring is then applied at `landmark + window`
  (`R/get_lm_data.R`, lines 107–110), so a future entrant can contribute
  manufactured event-free time before observation.

There is a separate hard bug in wide format. `stack_data()` evaluates `rtime`
and uses `(id, rtime)` before branching on format (`R/stack_data.R`, lines
90–109). Omitting `rtime` as documented for wide data produces the exact error
`argument "rtime" is missing, with no default` (captured in
`results/recon_errors.csv`).

### `add_interactions()`

This function is agnostic to entry. It transforms the landmark column and
multiplies each selected covariate by the requested landmark basis functions
(`R/add_interactions.R`, lines 82–123). It will work on an age scale, but raw
ages make polynomial terms poorly scaled. Centered landmark bases should be
considered as a numerical-stability option, not as a left-truncation fix.

### `dynamic_lm()`

The unpenalized dispatcher sends the supplied formula directly to
`survival::coxph()` or `riskRegression::CSC()` (`R/dynamic_lm_helper.R`, lines
80–96). The formula can therefore represent counting-process data.

The penalized reconstruction already hard-codes the landmark as entry:

```r
LHS_surv <- paste0("Surv(", entry, ",", exit, ",", status, ")")
```

[SOURCE: local `R/dynamic_lm.R`, lines 448–460]

For competing risks it analogously creates `Hist(exit, status, entry)`. This is
correct for the ordinary stacked landmark interval, whose risk period begins at
the landmark, but it cannot use `max(landmark, observation_entry)` for a
within-window delayed-entry design because `entry` is fixed to `lmdata$lm_col`.

### `pen_lm()` and `cv.pen_lm()`

Both call `check_penlm_inputs()`, which constructs
`Hist(exit, status, landmark)` and converts it to one or more start/stop `Surv`
responses (`R/checks.R`, lines 450–505). `pen_lm()` then calls
`glmnet(..., family = "cox")`; `cv.pen_lm()` does the same through
`cv.glmnet()` and assigns folds at the subject level (`R/pen_lm.R`, lines
70–97; `R/cv.pen_lm.R`, lines 88–135). Neither pathway has a distinct
observation-entry field.

### `predict.dynamicLM()`

Prediction computes the cause-specific linear predictor, obtains baseline
survival curves from the fitted Cox models, and integrates hazard increments
between the requested landmark and `landmark + window` (`R/predict.R`, lines
163–235). Entry eligibility is therefore inherited from the supplied prediction
rows; prediction does not validate that the subject was observed by the
landmark.

### `score()` and IPCW

`dynamicLM::score()` subsets prediction rows separately at each landmark and
calls `riskRegression::Score()` at `landmark + window` while requesting the IID
decomposition (`R/score.R`, lines 204–275). It then sends the landmark-specific
scores and IID contributions to `summary_metric()` (`R/score.R`, lines
359–369).

Inside the inspected `riskRegression` development source, marginal censoring
weights are fitted with:

```r
sFormula <- update(formula,
                   "Surv(riskRegression_time,riskRegression_status)~1")
fit <- prodlim::prodlim(sFormula, data = data, reverse = TRUE)
```

[SOURCE: `riskRegression` commit
`b3eea672eab976daf60585cf690d32cb50cdade6`, `R/getCensoringWeights.R`,
lines 17–25]

The Cox censoring branch likewise uses `Surv(time, status)` without an entry
term (lines 35–38). Thus a delayed-entry formula can be parsed for the outcome,
but entry is discarded when the censoring distribution is estimated.

This omission has two different implications:

- under **strict landmark eligibility**, every included subject is already
  observed at the landmark, so estimating censoring conditionally in that
  landmark subset does not require heterogeneous post-landmark entry;
- under **entry within the prediction window**, subjects begin risk at different
  ages and the current censoring weights are wrong because pre-entry time is
  included in the censoring risk set.

### Summary IID calculation

`summary_metric()` averages the landmark-specific point estimates. For its CI,
it reshapes subject-level IID contributions across landmarks, inserts zero for a
subject absent at a landmark, rescales by each landmark subsample fraction,
estimates the cross-landmark covariance matrix, and applies the delta method to
the equal-weighted mean (`R/summary_metric.R`, lines 79–170). This is why the CI
is not an average of landmark-specific confidence intervals: it uses covariance
between landmarks induced by repeated subjects.[^2]

The current implementation defines total sample size as all unique IDs, then
asserts that the first landmark's contribution count equals that total
(`R/summary_metric.R`, lines 61–63 and 101–114). Staggered entry violates this
implementation assumption even when the eligibility sets themselves are
correct. The minimal reproducer stops with `Subsample size mismatch.`

## Eligibility claim

**PROVEN (set-theoretic, conditional on the observation-interval definition).**
Let subject (i) be observable on the half-open age interval
([E_i, X_i)). By definition, the subject is observable at landmark age (a)
if and only if (a \in [E_i, X_i)). Expanding membership in a half-open interval
gives (E_i \le a) and (a < X_i). Therefore the eligible landmark set is

\[
R(a) = \{i : E_i \le a < X_i\}.
\]

The existing filter implements only (a < X_i). Whenever (a < E_i < X_i),
that filter includes an individual who is not observable at (a). The local
reproducer demonstrates exactly this case at landmark age 55 with a ten-year
window (TRACE `01_recon::toy_landmark` and `01_recon::toy_window`).

**CONJECTURED.** The existing summary-IID variance cannot be made generally
valid under random left truncation merely by replacing the first-landmark size
assertion. A formal influence-function derivation for the desired conditional
or source-population estimand is not provided here. The paper states that the
published CI relies on independent censoring and independent testing data; it
does not establish the left-truncated extension.[^2]

**EMPIRICAL-ONLY.** Coverage of any candidate correction remains PENDING until
the simulation grid is run.

## Public-data inventory

All values below were computed from the installed `survival` datasets, not
copied from their help pages.

| Dataset | Rows | Columns | Unique subjects | Variables |
|---|---:|---:|---:|---|
| `nafld1` | 17,549 | 9 | 17,549 | `id`, `age`, `male`, `weight`, `height`, `bmi`, `case.id`, `futime`, `status` |
| `nafld2` | 400,123 | 4 | 15,666 | `id`, `days`, `test`, `value` |
| `nafld3` | 34,327 | 3 | 12,453 | `id`, `days`, `event` |
| `pbcseq` | 1,945 | 19 | 312 | `id`, `futime`, `status`, `trt`, `age`, `sex`, `day`, `ascites`, `hepato`, `spiders`, `edema`, `bili`, `chol`, `albumin`, `alk.phos`, `ast`, `platelet`, `protime`, `stage` |
| `mgus2` | 1,384 | 11 | 1,384 | `id`, `age`, `sex`, `dxyr`, `hgb`, `creat`, `mspike`, `ptime`, `pstat`, `futime`, `death` |
| `flchain` | 7,874 | 11 | 7,874 | `age`, `sex`, `sample.yr`, `kappa`, `lambda`, `flc.grp`, `creatinine`, `mgus`, `futime`, `death`, `chapter` |
| `rotterdam` | 2,982 | 15 | 2,982 | `pid`, `year`, `age`, `meno`, `size`, `grade`, `nodes`, `pgr`, `er`, `hormon`, `chemo`, `rtime`, `recur`, `dtime`, `death` |
| `colon` | 1,858 | 16 | 929 | `id`, `study`, `rx`, `sex`, `age`, `obstruct`, `perfor`, `adhere`, `nodes`, `status`, `differ`, `extent`, `surg`, `node4`, `time`, `etype` |

Trace refs, in table order: `01_recon::nafld1_{rows,columns,unique_ids}`;
`01_recon::nafld2_{rows,columns,unique_ids}`;
`01_recon::nafld3_{rows,columns,unique_ids}`;
`01_recon::pbcseq_{rows,columns,unique_ids}`;
`01_recon::mgus2_{rows,columns,unique_ids}`;
`01_recon::flchain_{rows,columns,unique_ids}`;
`01_recon::rotterdam_{rows,columns,unique_ids}`; and
`01_recon::colon_{rows,columns,unique_ids}`.

The installed help text is inconsistent with the installed NAFLD objects: it
describes ten variables for `nafld1` but the object has nine columns, and it
describes 34,340 `nafld3` observations while the object has 34,327 rows.
[SOURCE: installed `survival` help topic `nafld`; computed object dimensions are
traced above.] The data object, not the prose documentation, is used downstream.

### Entry-age distributions

Values are minimum / first quartile / median / third quartile / maximum, followed
by standard deviation, in years.

| Dataset | Distribution | SD |
|---|---:|---:|
| `nafld1` | 18 / 42 / 53 / 63 / 98 | 14.723 |
| `pbcseq` | 26.278 / 42.239 / 49.795 / 56.715 / 78.439 | 10.581 |
| `mgus2` | 24 / 63 / 72 / 79 / 96 | 12.172 |
| `flchain` | 50 / 55 / 63 / 72 / 101 | 10.463 |
| `rotterdam` | 24 / 45 / 54 / 65 / 90 | 12.953 |
| `colon` | 18 / 53 / 61 / 69 / 85 | 11.949 |

Trace refs, by row: `01_recon::<dataset>_entry_age_{minimum,q1,median,q3,maximum,sd}`.
The detailed table, including means at full precision, is
`results/recon_entry_age_summary.csv`.

### Event structure

- `nafld1`: 1,364 deaths and 16,185 censored records (TRACE
  `01_recon::nafld1_event_1` and `01_recon::nafld1_event_0`). `nafld2` contains
  longitudinal laboratory measurements; `nafld3` contains named dated events.
- `pbcseq`: 29 transplants, 140 deaths, and 143 censored subjects (TRACE
  `01_recon::pbcseq_event_1`, `01_recon::pbcseq_event_2`, and
  `01_recon::pbcseq_event_0`). Its help page explicitly notes irregular extra
  visits and potentially informative missing laboratory values.[^3]
- `mgus2`: using time to first progression or death, 115 progressions, 860 deaths
  without prior progression, and 409 censored subjects (TRACE
  `01_recon::mgus2_first_event_event_{1,2,0}`). The event construction matches
  the `survival` manual example.[^3]
- `flchain`: 2,169 deaths and 5,705 censored subjects (TRACE
  `01_recon::flchain_event_1` and `01_recon::flchain_event_0`). It is a single
  measurement per subject, not longitudinal covariate follow-up.[^3]
- `rotterdam`: the constructed first-event endpoint contains 1,518 recurrences,
  195 deaths without prior recurrence, and 1,269 censored subjects (TRACE
  `01_recon::rotterdam_first_event_event_{1,2,0}`).
- `colon`: exactly two endpoint rows per subject, one for recurrence and one for
  death. It is not directly a single time-to-first-competing-event table.

## Dataset suitability decisions

- **KEEP:** `nafld1/2/3` for the strongest public age-entry and irregular-visit
  stress test; `pbcseq` for the tutorial-aligned two-time-scale comparison;
  `mgus2` for competing risks; `flchain` for scale and static-covariate testing.
- **PARK:** `rotterdam` as a secondary competing-risk application because it is
  not longitudinal and recurrence/death censoring has a documented ambiguity.
- **KILL for primary application:** `colon` because its paired endpoint rows do
  not directly encode a single first-event process and it adds no longitudinal
  age-scale stress beyond stronger candidates.

## Classification of findings

| Finding | Classification | Reason |
|---|---|---|
| No `E_i <= landmark` filter | Fundamental assumption violation exposed as a package limitation | The current API has no observation-entry argument; shared-origin data do not reveal the gap. |
| Pre-entry fixed-covariate row retained | Silent wrong answer / bug | The row is analyzable and contributes manufactured risk time. |
| Pre-entry varying-covariate row converted to `NA` | Silent masking / documented LOCF behavior with unintended age-scale consequence | Complete-case deletion can conceal the eligibility failure. |
| Wide format requires `rtime` | Bug | It contradicts the documented “long format only” argument. |
| IPCW response drops entry | Fundamental assumption violation for within-window delayed entry | Both marginal and Cox censoring branches reconstruct a right-censored response. |
| First landmark must contain every ID | Fundamental assumption embedded as a hard assertion | Valid staggered-entry landmark sets need not contain all eventual IDs at the first age. |

## Sources

[^1]: The Han Lab. “[dynamicLM GitHub repository and README](https://github.com/thehanlab/dynamicLM).” Accessed during this audit.
[^2]: Fries AH, Choi E, Han SS. “[Penalized landmark supermodels (penLM) for dynamic prediction for time-to-event outcomes in high-dimensional data](https://link.springer.com/article/10.1186/s12874-024-02418-9).” *BMC Medical Research Methodology*. The methods describe the multivariate IID extension and state its independent-censoring/test-data assumptions.
[^3]: Therneau T. “[survival package reference manual](https://spout.ussg.indiana.edu/CRAN/web/packages/survival/survival.pdf).” Dataset documentation and examples for `pbcseq`, `mgus2`, `flchain`, `nafld`, `rotterdam`, and `colon`.
