#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(parallel)
  library(survival)
})
source(here::here("analysis", "_helpers.R"))
source(here::here("analysis", "_simulation_engine.R"))

analysis_seed <- 20260908L
assert_single_seed(analysis_seed)

replicates <- as.integer(Sys.getenv("SIM_REPLICATES", "1000"))
allow_small <- identical(Sys.getenv("SIM_ALLOW_SMALL", "0"), "1")
if (!allow_small) stopifnot(replicates >= 1000L)
stopifnot(replicates >= 1L)
cores <- as.integer(Sys.getenv("SIM_CORES", as.character(min(8L, parallel::detectCores()))))
stopifnot(cores >= 1L)

truth_beta <- log(1.5)
lambda0 <- 0.01
entry_mean <- 45
prediction_window <- 10
nominal_coverage <- 0.95
landmark_intervals <- c(2, 4, 6)
censoring_targets <- c(0, 0.15, 0.30, 0.50)
sample_sizes <- c(500L, 750L, 1500L)
entry_sds <- c(0.25, 4, 10)
methods <- c("naive", "strict", "delayed")

grid <- as.data.table(expand.grid(
  landmark_interval = landmark_intervals,
  censoring_target = censoring_targets,
  n = sample_sizes,
  entry_sd = entry_sds,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
))
grid[, cell_id := .I]
setcolorder(grid, c("cell_id", "landmark_interval", "censoring_target", "n", "entry_sd"))
stopifnot(nrow(grid) == length(landmark_intervals) * length(censoring_targets) *
            length(sample_sizes) * length(entry_sds), !anyDuplicated(grid$cell_id))

calibration_grid <- unique(grid[, .(landmark_interval, censoring_target, entry_sd)])
setorder(calibration_grid, landmark_interval, entry_sd, censoring_target)
calibration_grid[, calibration_id := .I]

calibrations <- lapply(seq_len(nrow(calibration_grid)), function(i) {
  row <- calibration_grid[i]
  seed <- 510000000L + i
  set.seed(seed)
  landmarks <- 50 + c(0, row$landmark_interval, 2 * row$landmark_interval)
  study_end <- max(landmarks) + prediction_window
  pilot <- generate_left_truncated_cohort(
    50000L, truth_beta, lambda0, entry_mean, row$entry_sd, study_end
  )
  if (row$censoring_target == 0) {
    rate <- 0
    achieved <- 0
  } else {
    rate <- calibrate_censor_rate(row$censoring_target, pilot, study_end)
    # Validate on fresh censoring uniforms from the same fixed pilot cohort.
    u <- stats::runif(nrow(pilot))
    censor_age <- pilot$entry - log(u) / rate
    achieved <- mean(censor_age < pilot$event_age & censor_age < study_end)
  }
  stopifnot(is.finite(rate), rate >= 0, is.finite(achieved))
  if (abs(achieved - row$censoring_target) >= 0.01)
    stop("Censoring calibration missed target: id=", i,
         ", target=", row$censoring_target, ", achieved=", achieved)
  data.table(
    calibration_id = row$calibration_id,
    landmark_interval = row$landmark_interval,
    censoring_target = row$censoring_target,
    entry_sd = row$entry_sd,
    study_end = study_end,
    censor_rate = rate,
    pilot_achieved_censoring = achieved,
    pilot_n = nrow(pilot),
    seed = seed
  )
})
calibrations <- rbindlist(calibrations)
stopifnot(nrow(calibrations) == nrow(calibration_grid),
          !anyDuplicated(calibrations[, .(landmark_interval, censoring_target, entry_sd)]))

grid <- merge(
  grid, calibrations[, .(landmark_interval, censoring_target, entry_sd,
                         study_end, censor_rate)],
  by = c("landmark_interval", "censoring_target", "entry_sd"),
  all.x = TRUE, sort = FALSE
)
setorder(grid, cell_id)
stopifnot(nrow(grid) == max(grid$cell_id), !anyNA(grid))

truth_auc <- rbindlist(lapply(landmark_intervals, function(interval) {
  landmarks <- 50 + c(0, interval, 2 * interval)
  auc_fine <- vapply(
    landmarks, numerical_true_auc, numeric(1L),
    w = prediction_window, beta = truth_beta, lambda0 = lambda0,
    grid_step = 0.00125
  )
  auc_coarse <- vapply(
    landmarks, numerical_true_auc, numeric(1L),
    w = prediction_window, beta = truth_beta, lambda0 = lambda0,
    grid_step = 0.0025
  )
  data.table(
    landmark_interval = interval,
    landmark = landmarks,
    AUC = auc_fine,
    coarse_AUC = auc_coarse,
    quadrature_difference = auc_fine - auc_coarse,
    summary_AUC = mean(auc_fine)
  )
}))
stopifnot(all(abs(truth_auc$quadrature_difference) < 1e-5),
          all(truth_auc$AUC >= 0 & truth_auc$AUC <= 1))

seed_manifest <- grid[, .(replicate = seq_len(replicates)), by = cell_id]
seed_manifest[, seed := 410000000L + cell_id * 2000L + replicate]
stopifnot(!anyDuplicated(seed_manifest$seed), all(seed_manifest$seed > 0))

simulate_cell <- function(i) {
  row <- grid[i]
  cell_seeds <- seed_manifest[cell_id == row$cell_id]
  landmarks <- 50 + c(0, row$landmark_interval, 2 * row$landmark_interval)
  started <- proc.time()[["elapsed"]]
  cell_results <- vector("list", replicates * length(methods))
  cursor <- 1L
  for (replicate_id in seq_len(replicates)) {
    replicate_seed <- cell_seeds$seed[replicate_id]
    set.seed(replicate_seed)
    cohort <- generate_left_truncated_cohort(
      row$n, truth_beta, lambda0, entry_mean, row$entry_sd, row$study_end
    )
    cohort <- apply_independent_censoring(cohort, row$censor_rate, row$study_end)
    observed_entry_sd <- stats::sd(cohort$entry)
    observed_censoring <- mean(cohort$censored_before_end)
    for (method in methods) {
      stack <- make_simulation_stack(cohort, landmarks, prediction_window, method)
      fitted <- tryCatch(
        fit_fast_cluster_cox(stack),
        error = function(e) structure(conditionMessage(e), class = "simulation_fit_error")
      )
      if (inherits(fitted, "simulation_fit_error")) {
        cell_results[[cursor]] <- data.table(
          cell_id = row$cell_id, replicate = replicate_id, seed = replicate_seed,
          method = method, estimate = NA_real_, estimated_se = NA_real_,
          covered = NA_integer_, reject_null = NA_integer_, rows = nrow(stack),
          events = sum(stack$event), subjects = uniqueN(stack$id),
          observed_entry_sd = observed_entry_sd,
          observed_censoring = observed_censoring,
          converged = 0L, error = as.character(fitted)
        )
      } else {
        covered <- abs(fitted$estimate - truth_beta) <=
          stats::qnorm(1 - (1 - nominal_coverage) / 2) * fitted$robust_se
        reject_null <- abs(fitted$estimate / fitted$robust_se) >
          stats::qnorm(1 - (1 - nominal_coverage) / 2)
        cell_results[[cursor]] <- data.table(
          cell_id = row$cell_id, replicate = replicate_id, seed = replicate_seed,
          method = method, estimate = fitted$estimate,
          estimated_se = fitted$robust_se,
          covered = as.integer(covered), reject_null = as.integer(reject_null),
          rows = fitted$rows, events = fitted$events, subjects = fitted$subjects,
          observed_entry_sd = observed_entry_sd,
          observed_censoring = observed_censoring,
          converged = 1L, error = ""
        )
      }
      cursor <- cursor + 1L
    }
  }
  results <- rbindlist(cell_results)
  stopifnot(nrow(results) == replicates * length(methods),
            all(results$cell_id == row$cell_id),
            all(table(results$method) == replicates))
  list(results = results,
       runtime_seconds = proc.time()[["elapsed"]] - started)
}

wall_start <- proc.time()[["elapsed"]]
cell_runs <- parallel::mclapply(
  seq_len(nrow(grid)), simulate_cell,
  mc.cores = cores, mc.preschedule = TRUE, mc.set.seed = FALSE
)
wall_runtime <- proc.time()[["elapsed"]] - wall_start
if (any(vapply(cell_runs, inherits, logical(1L), "try-error")))
  stop("At least one parallel simulation cell failed.")
raw <- rbindlist(lapply(cell_runs, `[[`, "results"))
cell_runtime <- data.table(
  cell_id = grid$cell_id,
  runtime_seconds = vapply(cell_runs, `[[`, numeric(1L), "runtime_seconds")
)
stopifnot(
  nrow(raw) == nrow(grid) * replicates * length(methods),
  raw[, all(.N == replicates), by = .(cell_id, method)]$V1,
  all(raw$seed == seed_manifest$seed[
    match(paste(raw$cell_id, raw$replicate),
          paste(seed_manifest$cell_id, seed_manifest$replicate))
  ])
)

failures <- raw[converged == 0L,
                .(failure_count = .N), by = .(cell_id, method, error)]
successful <- raw[converged == 1L]
simulation_summary <- successful[, .(
  successful_replicates = .N,
  bias = mean(estimate - truth_beta),
  empirical_se = stats::sd(estimate),
  mean_estimated_se = mean(estimated_se),
  se_ratio = mean(estimated_se) / stats::sd(estimate),
  coverage = mean(covered),
  coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / .N),
  power = mean(reject_null),
  power_mcse = sqrt(mean(reject_null) * (1 - mean(reject_null)) / .N),
  mean_rows = mean(rows),
  mean_events = mean(events),
  mean_subjects = mean(subjects),
  mean_observed_entry_sd = mean(observed_entry_sd),
  mean_observed_censoring = mean(observed_censoring)
), by = .(cell_id, method)]
simulation_summary <- merge(simulation_summary, grid, by = "cell_id", all.x = TRUE)
failure_counts <- raw[, .(failed_replicates = sum(converged == 0L)),
                      by = .(cell_id, method)]
simulation_summary <- merge(simulation_summary, failure_counts,
                            by = c("cell_id", "method"), all.x = TRUE)
setorder(simulation_summary, cell_id, method)
stopifnot(nrow(simulation_summary) == nrow(grid) * length(methods),
          all(simulation_summary$successful_replicates +
                simulation_summary$failed_replicates == replicates),
          all(is.finite(as.matrix(simulation_summary[, .(
            bias, empirical_se, mean_estimated_se, se_ratio, coverage,
            coverage_mcse, power, power_mcse, mean_rows, mean_events,
            mean_subjects, mean_observed_entry_sd, mean_observed_censoring
          )]))))

successful_with_grid <- merge(
  successful, grid[, .(cell_id, entry_sd)], by = "cell_id", all.x = TRUE
)
heterogeneity_summary <- successful_with_grid[, .(
  successful_replicates = .N,
  bias = mean(estimate - truth_beta),
  mean_estimated_se = mean(estimated_se),
  coverage = mean(covered),
  coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / .N),
  power = mean(reject_null),
  power_mcse = sqrt(mean(reject_null) * (1 - mean(reject_null)) / .N),
  mean_observed_entry_sd = mean(observed_entry_sd),
  mean_observed_censoring = mean(observed_censoring)
), by = .(entry_sd, method)]
# Count failures without a cartesian join.
heterogeneity_failures <- merge(
  raw[, .(failed_replicates = sum(converged == 0L)), by = .(cell_id, method)],
  grid[, .(cell_id, entry_sd)], by = "cell_id"
)[, .(failed_replicates = sum(failed_replicates)), by = .(entry_sd, method)]
heterogeneity_summary <- merge(heterogeneity_summary, heterogeneity_failures,
                               by = c("entry_sd", "method"), all.x = TRUE,
                               suffixes = c("", ".counted"))
if ("failed_replicates.counted" %in% names(heterogeneity_summary)) {
  heterogeneity_summary[, failed_replicates := failed_replicates.counted]
  heterogeneity_summary[, failed_replicates.counted := NULL]
}
setorder(heterogeneity_summary, entry_sd, method)

data.table::fwrite(raw, project_path("results", "simulation_raw.csv"))
data.table::fwrite(seed_manifest, project_path("results", "simulation_seed_manifest.csv"))
data.table::fwrite(grid, project_path("results", "simulation_grid.csv"))
data.table::fwrite(calibrations, project_path("results", "simulation_calibration.csv"))
data.table::fwrite(truth_auc, project_path("results", "simulation_auc_truth.csv"))
data.table::fwrite(simulation_summary, project_path("results", "simulation_summary.csv"))
data.table::fwrite(heterogeneity_summary,
                   project_path("results", "simulation_heterogeneity_summary.csv"))
data.table::fwrite(failures, project_path("results", "simulation_failures.csv"))
data.table::fwrite(cell_runtime, project_path("results", "simulation_cell_runtime.csv"))

trace_values(
  report = "05_simulation",
  labels = c("truth_beta", "hazard_ratio", "lambda0", "entry_mean",
             "prediction_window", "nominal_coverage", "replicates_per_cell",
             "grid_cells", "methods", "total_model_fits", "parallel_cores",
             "wall_runtime_seconds", "fit_failures"),
  values = c(truth_beta, exp(truth_beta), lambda0, entry_mean,
             prediction_window, nominal_coverage, replicates, nrow(grid),
             length(methods), nrow(raw), cores, wall_runtime, nrow(raw) - nrow(successful)),
  units = c("log hazard ratio", "hazard ratio", "hazard per year", "years",
            "years", "proportion", "replicates", "cells", "methods", "fits",
            "cores", "seconds", "fits"),
  script = "analysis/05_simulation.R", function_name = "simulation grid",
  seed = "results/simulation_seed_manifest.csv",
  runtime_seconds = wall_runtime
)

for (i in seq_len(nrow(heterogeneity_summary))) {
  row <- heterogeneity_summary[i]
  prefix <- paste0("entry_sd_", gsub("\\.", "p", row$entry_sd), "_", row$method)
  fields <- c("successful_replicates", "failed_replicates", "bias",
              "mean_estimated_se", "coverage", "coverage_mcse", "power",
              "power_mcse", "mean_observed_entry_sd", "mean_observed_censoring")
  trace_values(
    report = "05_simulation", labels = paste0(prefix, "_", fields),
    values = as.numeric(row[, ..fields]),
    units = c("replicates", "replicates", "log hazard ratio", "standard error",
              "proportion", "Monte Carlo standard error", "proportion",
              "Monte Carlo standard error", "years", "proportion"),
    script = "analysis/05_simulation.R",
    function_name = "fit_fast_cluster_cox/heterogeneity aggregation",
    seed = "results/simulation_seed_manifest.csv", runtime_seconds = wall_runtime
  )
}

for (i in seq_len(nrow(truth_auc))) {
  row <- truth_auc[i]
  prefix <- paste0("auc_truth_interval_", row$landmark_interval,
                   "_landmark_", row$landmark)
  trace_values(
    report = "05_simulation",
    labels = paste0(prefix, c("_auc", "_summary_auc", "_quadrature_difference")),
    values = c(row$AUC, row$summary_AUC, row$quadrature_difference),
    units = c("AUC", "AUC", "AUC"), script = "analysis/05_simulation.R",
    function_name = "numerical_true_auc", seed = NA_integer_,
    runtime_seconds = wall_runtime
  )
}

# Every cell-level quantity that may be quoted is traceable to its exact seed range.
cell_fields <- c("successful_replicates", "failed_replicates", "bias",
                 "empirical_se", "mean_estimated_se", "se_ratio", "coverage",
                 "coverage_mcse", "power", "power_mcse", "mean_rows",
                 "mean_events", "mean_subjects", "mean_observed_entry_sd",
                 "mean_observed_censoring")
for (i in seq_len(nrow(simulation_summary))) {
  row <- simulation_summary[i]
  prefix <- paste0("cell_", sprintf("%03d", row$cell_id), "_", row$method, "_")
  seeds <- seed_manifest[cell_id == row$cell_id]$seed
  trace_values(
    report = "05_simulation", labels = paste0(prefix, cell_fields),
    values = as.numeric(row[, ..cell_fields]),
    units = c("replicates", "replicates", "log hazard ratio", "standard error",
              "standard error", "ratio", "proportion", "Monte Carlo standard error",
              "proportion", "Monte Carlo standard error", "rows", "events",
              "subjects", "years", "proportion"),
    script = "analysis/05_simulation.R",
    function_name = "fit_fast_cluster_cox/cell aggregation",
    seed = paste0(min(seeds), "-", max(seeds)), runtime_seconds = cell_runtime[
      cell_id == row$cell_id
    ]$runtime_seconds
  )
}

stopifnot(
  file.exists(project_path("results", "simulation_raw.csv")),
  file.exists(project_path("results", "simulation_summary.csv")),
  file.exists(project_path("results", "simulation_seed_manifest.csv")),
  all(cell_runtime$runtime_seconds > 0)
)
