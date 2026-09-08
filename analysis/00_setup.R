#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
})
source(here::here("analysis", "_helpers.R"))

required <- c(
  "devtools", "survival", "prodlim", "dynpred", "glmnet", "data.table",
  "pec", "msm", "testthat", "here", "riskRegression", "dynamicLM"
)
status <- data.frame(
  package = required,
  installed = vapply(required, requireNamespace, logical(1L), quietly = TRUE),
  version = vapply(required, function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return("PENDING")
    as.character(utils::packageVersion(package))
  }, character(1L)),
  stringsAsFactors = FALSE
)
stopifnot(nrow(status) == length(required), !anyDuplicated(status$package))

session_path <- here::here("results", "session.txt")
dir.create(dirname(session_path), recursive = TRUE, showWarnings = FALSE)
sink(session_path)
cat("Package status\n")
print(status, row.names = FALSE)
cat("\nSession\n")
print(sessionInfo())
sink()
stopifnot(file.exists(session_path), file.info(session_path)$size > 0L)

runtime <- timed({
  installed_count <- sum(status$installed)
  missing_count <- sum(!status$installed)
  c(installed_count = installed_count, missing_count = missing_count)
})
trace_values(
  report = "setup",
  labels = names(runtime$value),
  values = unname(runtime$value),
  units = "packages",
  script = "analysis/00_setup.R",
  function_name = "base::sum",
  runtime_seconds = runtime$runtime_seconds,
  notes = "Package availability at setup audit; versions are recorded in results/session.txt."
)

if (!all(status$installed)) {
  message("PENDING packages: ", paste(status$package[!status$installed], collapse = ", "))
}
