# Age-Scale Landmark Supermodels

This branch is a source-audited extension of `dynamicLM` for attained-age time
scales with staggered cohort entry. It uses public data only.

Run scripts from the repository root. All analysis paths use `here::here()`;
computed report values are cross-referenced in `results/trace.csv`, and every
stochastic analysis writes its seed manifest under `results/`.

## Current prerequisite status

Run:

```sh
Rscript analysis/00_setup.R
```

The required GitHub development build of `riskRegression` must be installed.
On Apple silicon systems whose R configuration points to absent CRAN Fortran
libraries, the project includes `tools/Makevars.macos-arm64` as an explicit,
platform-specific workaround. Do not use it on other platforms.

The original upstream README is retained in `README.Rmd` while the research
workflow is under construction.

## Reproducible raw simulation artifacts

The two replicate-level tables above 25 MB are intentionally not tracked:
`results/simulation_raw.csv` and `results/informative_entry_raw.csv`. Before
removal, both were rebuilt from their committed seed manifests and verified
byte-for-byte against the originals. The row counts, byte counts, MD5 hashes,
and equality results are in
[`raw_regeneration_verification.csv`](results/check/raw_regeneration_verification.csv).

To recreate both tables from a clean checkout without rerunning the analysis
reports, run:

```sh
REGENERATE_CORES=8 Rscript analysis/regenerate_raw.R --write
```

The script refuses to overwrite either raw table. To audit existing tables
against fresh manifest-driven rebuilds, use `--verify` instead. Phase 7's
frailty decomposition requires `results/informative_entry_raw.csv`; regenerate
it first if the table is absent.

## Reproduce the analyses

After setup, run:

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

The simulation scripts accept `SIM_CORES` to control local parallelism. Phase 6
also accepts explicit small-run environment switches for development, but the
committed result artifacts were produced with the required full settings; a
small run must not be reported as the study result. Raw regeneration uses
`REGENERATE_CORES` and is itself a full manifest replay, not a smoke test.

Run package and project tests with:

```sh
Rscript -e 'devtools::test()'
```

`tests/testthat.R` is the installed-package entry point used by `R CMD check`;
invoking it directly from the repository root does not use the correct test
directory.

For the CRAN-style check, build the tarball before checking it:

```sh
R CMD build .
_R_CHECK_CRAN_INCOMING_=TRUE _R_CHECK_FORCE_SUGGESTS_=false \
  R CMD check --as-cran --no-manual dynamicLM_1.0.0.tar.gz
```

The final Phase 7 check completed with inherited non-clean diagnostics. See
[`cran_attribution.csv`](results/check/cran_attribution.csv) and the preserved
post-fix log rather than treating test success as CRAN readiness.
