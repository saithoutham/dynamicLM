#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(parallel)
})
source(here::here("analysis", "_helpers.R"))
source(here::here("analysis", "_simulation_engine.R"))

analysis_seed <- 20260908L
assert_single_seed(analysis_seed)
outer_replicates <- as.integer(Sys.getenv("AUC_REPLICATES", "1000"))
bootstrap_replicates <- as.integer(Sys.getenv("AUC_BOOTSTRAPS", "100"))
allow_small <- identical(Sys.getenv("SIM_ALLOW_SMALL", "0"), "1")
if (!allow_small) stopifnot(outer_replicates >= 1000L)
stopifnot(outer_replicates >= 1L, bootstrap_replicates >= 2L)
cores <- as.integer(Sys.getenv("SIM_CORES", as.character(min(8L, parallel::detectCores()))))
stopifnot(cores >= 1L)

truth_beta <- log(1.5)
lambda0 <- 0.01
entry_mean <- 45
entry_sds <- c(0.25, 4, 10)
landmark_interval <- 4
landmarks <- 50 + c(0, landmark_interval, 2 * landmark_interval)
prediction_window <- 10
study_end <- max(landmarks) + prediction_window
censoring_target <- 0.30
sample_size <- 750L
nominal_coverage <- 0.95
methods <- c("naive", "strict")

calibrations <- data.table::fread(project_path("results", "simulation_calibration.csv"))
target_interval <- landmark_interval
target_censoring <- censoring_target
calibrations <- calibrations[
  landmark_interval == target_interval & censoring_target == target_censoring &
    entry_sd %in% entry_sds
]
stopifnot(nrow(calibrations) == length(entry_sds),
          !anyDuplicated(calibrations$entry_sd))
truth_table <- data.table::fread(project_path("results", "simulation_auc_truth.csv"))
truth_summary_auc <- unique(truth_table[
  landmark_interval == target_interval, summary_AUC
])
stopifnot(length(truth_summary_auc) == 1L, is.finite(truth_summary_auc))

tasks <- CJ(entry_sd_index = seq_along(entry_sds),
            replicate = seq_len(outer_replicates))
tasks[, entry_sd := entry_sds[entry_sd_index]]
tasks[, outer_seed := 710000000L + entry_sd_index * 2000L + replicate]
tasks[, bootstrap_seed_naive := 810000000L + entry_sd_index * 4000L + replicate]
tasks[, bootstrap_seed_strict := 910000000L + entry_sd_index * 4000L + replicate]
stopifnot(!anyDuplicated(tasks$outer_seed),
          !anyDuplicated(tasks$bootstrap_seed_naive),
          !anyDuplicated(tasks$bootstrap_seed_strict))

run_task <- function(i) {
  task <- tasks[i]
  set.seed(task$outer_seed)
  cohort <- generate_left_truncated_cohort(
    sample_size, truth_beta, lambda0, entry_mean, task$entry_sd, study_end
  )
  rate <- calibrations[entry_sd == task$entry_sd]$censor_rate
  stopifnot(length(rate) == 1L, is.finite(rate), rate > 0)
  cohort <- apply_independent_censoring(cohort, rate, study_end)
  results <- lapply(methods, function(method) {
    stack <- make_simulation_stack(cohort, landmarks, prediction_window, method)
    true_risk <- 1 - exp(-lambda0 * prediction_window * exp(truth_beta * stack$x))
    bootstrap_seed <- if (method == "naive") {
      task$bootstrap_seed_naive
    } else {
      task$bootstrap_seed_strict
    }
    fit <- tryCatch(
      bootstrap_summary_auc_fast(
        stack, true_risk, prediction_window, bootstrap_replicates, bootstrap_seed
      ),
      error = function(e) structure(conditionMessage(e), class = "auc_error")
    )
    if (inherits(fit, "auc_error")) {
      return(data.table(
        entry_sd = task$entry_sd, replicate = task$replicate,
        outer_seed = task$outer_seed, bootstrap_seed = bootstrap_seed,
        method = method, estimate = NA_real_, estimated_se = NA_real_,
        covered = NA_integer_, reject_half = NA_integer_, converged = 0L,
        error = as.character(fit)
      ))
    }
    z <- stats::qnorm(1 - (1 - nominal_coverage) / 2)
    data.table(
      entry_sd = task$entry_sd, replicate = task$replicate,
      outer_seed = task$outer_seed, bootstrap_seed = bootstrap_seed,
      method = method, estimate = fit$estimate, estimated_se = fit$se,
      covered = as.integer(abs(fit$estimate - truth_summary_auc) <= z * fit$se),
      reject_half = as.integer(abs(fit$estimate - 0.5) > z * fit$se),
      converged = 1L, error = ""
    )
  })
  rbindlist(results)
}

wall_start <- proc.time()[["elapsed"]]
runs <- parallel::mclapply(
  seq_len(nrow(tasks)), run_task, mc.cores = cores,
  mc.preschedule = TRUE, mc.set.seed = FALSE
)
wall_runtime <- proc.time()[["elapsed"]] - wall_start
if (any(vapply(runs, inherits, logical(1L), "try-error")))
  stop("At least one AUC coverage task failed outside the captured scorer.")
raw <- rbindlist(runs)
stopifnot(nrow(raw) == nrow(tasks) * length(methods),
          raw[, all(.N == outer_replicates), by = .(entry_sd, method)]$V1)
failures <- raw[converged == 0L, .(failure_count = .N),
                by = .(entry_sd, method, error)]
successful <- raw[converged == 1L]
summary <- successful[, .(
  successful_replicates = .N,
  bias = mean(estimate - truth_summary_auc),
  empirical_se = stats::sd(estimate),
  mean_estimated_se = mean(estimated_se),
  se_ratio = mean(estimated_se) / stats::sd(estimate),
  coverage = mean(covered),
  coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / .N),
  power_against_half = mean(reject_half),
  power_mcse = sqrt(mean(reject_half) * (1 - mean(reject_half)) / .N)
), by = .(entry_sd, method)]
failure_counts <- raw[, .(failed_replicates = sum(converged == 0L)),
                      by = .(entry_sd, method)]
summary <- merge(summary, failure_counts, by = c("entry_sd", "method"))
setorder(summary, entry_sd, method)
stopifnot(nrow(summary) == length(entry_sds) * length(methods),
          all(summary$successful_replicates + summary$failed_replicates ==
                outer_replicates),
          all(is.finite(as.matrix(summary[, setdiff(names(summary),
                                                    c("entry_sd", "method")),
                                          with = FALSE]))))

fwrite(raw, project_path("results", "simulation_auc_coverage_raw.csv"))
fwrite(summary, project_path("results", "simulation_auc_coverage_summary.csv"))
fwrite(failures, project_path("results", "simulation_auc_coverage_failures.csv"))
fwrite(tasks, project_path("results", "simulation_auc_seed_manifest.csv"))

trace_values(
  report = "05_simulation",
  labels = c("auc_outer_replicates_per_cell", "auc_bootstrap_replicates",
             "auc_coverage_cells", "auc_total_outer_datasets",
             "auc_total_bootstrap_resamples", "auc_truth_summary_interval_4",
             "auc_coverage_wall_runtime_seconds", "auc_coverage_failures"),
  values = c(outer_replicates, bootstrap_replicates, nrow(summary), nrow(tasks),
             nrow(tasks) * length(methods) * bootstrap_replicates,
             truth_summary_auc, wall_runtime, nrow(raw) - nrow(successful)),
  units = c("replicates", "bootstrap replicates", "cells", "datasets",
            "resamples", "AUC", "seconds", "datasets"),
  script = "analysis/05_auc_coverage.R",
  function_name = "bootstrap_summary_auc_fast",
  seed = "results/simulation_auc_seed_manifest.csv",
  runtime_seconds = wall_runtime
)

fields <- c("successful_replicates", "failed_replicates", "bias",
            "empirical_se", "mean_estimated_se", "se_ratio", "coverage",
            "coverage_mcse", "power_against_half", "power_mcse")
for (i in seq_len(nrow(summary))) {
  row <- summary[i]
  prefix <- paste0("auc_entry_sd_", gsub("\\.", "p", row$entry_sd), "_", row$method)
  trace_values(
    report = "05_simulation", labels = paste0(prefix, "_", fields),
    values = as.numeric(row[, ..fields]),
    units = c("replicates", "replicates", "AUC", "standard error",
              "standard error", "ratio", "proportion",
              "Monte Carlo standard error", "proportion",
              "Monte Carlo standard error"),
    script = "analysis/05_auc_coverage.R",
    function_name = "bootstrap_summary_auc_fast/aggregation",
    seed = "results/simulation_auc_seed_manifest.csv",
    runtime_seconds = wall_runtime
  )
}

stopifnot(file.exists(project_path("results", "simulation_auc_coverage_summary.csv")),
          file.exists(project_path("results", "simulation_auc_seed_manifest.csv")))
