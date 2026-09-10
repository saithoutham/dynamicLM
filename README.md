# Attained-age landmarking with `dynamicLM`

I use this branch to study what changes when a landmark supermodel is indexed
by attained age rather than time since a shared baseline. The practical issue
is staggered observation: at a landmark age, a person should not enter the risk
set before their record begins.

This is a research branch of
[`thehanlab/dynamicLM`](https://github.com/thehanlab/dynamicLM), not a separate
package release. I kept the original package guide in `README.Rmd`; this file
describes the attained-age work and how I reproduced it. All analyses use
public data.

## What I changed

I added two explicit entry rules to the landmark-data construction:

- `strict` requires observation by the landmark (`entry <= landmark < exit`);
- `delayed` allows entry during the prediction window and uses a
  counting-process start time.

I also added entry-conditional censoring weights for AUC and Brier-score point
estimates. The analytic influence function for staggered entry is still open,
so `score_left_truncated()` deliberately reports
`"PENDING: not derived for staggered entry"` instead of presenting an
unproved standard error.

The simulations separate three questions that are easy to conflate: whether a
risk set is constructed correctly, whether a fitted coefficient is biased, and
whether discrimination changes. The independent-entry experiment was largely
null. Measured risk-dependent entry exposed the eligibility problem, while the
frailty decomposition pointed to omitted-predictor bias rather than an added
truncation effect in that particular design. I report the negative results
alongside the positive ones.

## Where I would start

- [`reports/00_EXECUTIVE_SUMMARY.md`](reports/00_EXECUTIVE_SUMMARY.md) gives the
  short version of the project.
- [`reports/99_LIMITATIONS.md`](reports/99_LIMITATIONS.md) lists the unresolved
  statistical and computational limitations.
- [`manuscript/age_scale_landmark.md`](manuscript/age_scale_landmark.md) is the
  current manuscript draft.
- [`INNOVATIONS.md`](INNOVATIONS.md) records ideas I tested, including the ones
  I rejected.
- [`results/README.md`](results/README.md) explains the result files and the
  intentionally header-only failure logs.

## Local setup

I run the analysis scripts from the repository root. Paths are resolved with
`here::here()`.

```sh
Rscript analysis/00_setup.R
```

The setup script audits the local R environment; it does not install packages.
The analyses require the GitHub development version of `riskRegression`. On
the recorded R installation, `dynpred` remains unavailable but is not imported
by the source used here. Apple silicon systems with a missing Fortran library
path can use `tools/Makevars.macos-arm64`; that file is not intended for other
platforms.

## Reproducing the analyses

The numbered scripts are meant to be run in order:

```sh
Rscript analysis/01_recon.R
Rscript analysis/02_baseline.R
Rscript analysis/03_break_diagnosis.R
Rscript analysis/04_correction.R
Rscript analysis/05_engine_validation.R
Rscript analysis/05_simulation.R
Rscript analysis/05_auc_coverage.R
Rscript analysis/05_timescale_diagnostic.R
Rscript analysis/06_application.R
Rscript analysis/07_informative_entry.R
Rscript analysis/07_entry_diagnostic.R
Rscript analysis/07_auc_informative.R
Rscript analysis/08_frailty_decomposition.R
```

The simulation scripts accept `SIM_CORES`. Their committed summaries already
come from the full runs; most readers will only need the reports, result tables,
and seed manifests.

Two large replicate-level tables are reproducible but not tracked:
`results/simulation_raw.csv` and `results/informative_entry_raw.csv`. I removed
them only after rebuilding them from the committed seed manifests and comparing
the files byte for byte. The retained hashes and file sizes are in
[`results/check/raw_regeneration_verification.csv`](results/check/raw_regeneration_verification.csv).
To recreate the raw tables:

```sh
REGENERATE_CORES=8 Rscript analysis/regenerate_raw.R --write
```

The regeneration script refuses to overwrite an existing raw table. Its
`--verify` mode compares existing tables with a fresh manifest-driven rebuild.

## Provenance and tests

Numeric results cited in the reports are indexed in `results/trace.csv`, with
the script, function, seed, runtime, and commit that produced it. Seed manifests
are stored beside the corresponding result tables. The trace-reference audit
is:

```sh
Rscript tools/audit_trace.R
```

I run the package tests from the repository root with:

```sh
Rscript -e 'devtools::test()'
```

`tests/testthat.R` is the installed-package entry point used by `R CMD check`;
running that file directly from the repository root does not select the test
directory correctly.

For a CRAN-style check, I build the source package first:

```sh
R CMD build .
_R_CHECK_CRAN_INCOMING_=TRUE _R_CHECK_FORCE_SUGGESTS_=false \
  R CMD check --as-cran --no-manual dynamicLM_1.0.0.tar.gz
```

The latest recorded check is not clean. Its remaining warning and notes are
also present upstream; the side-by-side attribution and preserved logs are in
[`results/check/cran_attribution.csv`](results/check/cran_attribution.csv).
