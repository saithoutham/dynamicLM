# Age-Scale Landmark Supermodels

This branch is a source-audited extension of `dynamicLM` for attained-age time
scales with staggered cohort entry. It uses public data only.

Reproduction is being built phase by phase. Until the final audit is complete,
run scripts from the repository root in numeric order. All analysis paths use
`here::here()`; computed report values are cross-referenced in
`results/trace.csv`.

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
