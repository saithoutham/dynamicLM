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
outer_replicates <- as.integer(Sys.getenv("DIAGNOSTIC_REPLICATES", "1000"))
bootstrap_replicates <- as.integer(Sys.getenv("DIAGNOSTIC_BOOTSTRAPS", "100"))
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
alpha <- 0.05

calibrations <- fread(project_path("results", "simulation_calibration.csv"))
calibrations <- calibrations[
  landmark_interval == 4 & censoring_target == 0.30 & entry_sd %in% entry_sds
]
stopifnot(nrow(calibrations) == length(entry_sds))

tasks <- CJ(entry_sd_index = seq_along(entry_sds),
            replicate = seq_len(outer_replicates))
tasks[, entry_sd := entry_sds[entry_sd_index]]
tasks[, outer_seed := 710000000L + entry_sd_index * 2000L + replicate]
tasks[, bootstrap_seed := 1010000000L + entry_sd_index * 2000L + replicate]
stopifnot(!anyDuplicated(tasks$outer_seed), !anyDuplicated(tasks$bootstrap_seed))

run_task <- function(i) {
  task <- tasks[i]
  set.seed(task$outer_seed)
  cohort <- generate_left_truncated_cohort(
    sample_size, truth_beta, lambda0, entry_mean, task$entry_sd, study_end
  )
  rate <- calibrations[entry_sd == task$entry_sd]$censor_rate
  stopifnot(length(rate) == 1L, is.finite(rate), rate > 0)
  cohort <- apply_independent_censoring(cohort, rate, study_end)
  naive <- make_simulation_stack(cohort, landmarks, prediction_window, "naive")
  strict <- make_simulation_stack(cohort, landmarks, prediction_window, "strict")
  naive_risk <- 1 - exp(-lambda0 * prediction_window * exp(truth_beta * naive$x))
  strict_risk <- 1 - exp(-lambda0 * prediction_window * exp(truth_beta * strict$x))
  fit <- tryCatch(
    paired_bootstrap_auc_delta_fast(
      naive, naive_risk, strict, strict_risk, cohort$id,
      prediction_window, bootstrap_replicates, task$bootstrap_seed
    ),
    error = function(e) structure(conditionMessage(e), class = "diagnostic_error")
  )
  if (inherits(fit, "diagnostic_error")) {
    return(data.table(
      entry_sd = task$entry_sd, replicate = task$replicate,
      outer_seed = task$outer_seed, bootstrap_seed = task$bootstrap_seed,
      naive_auc = NA_real_, strict_auc = NA_real_, delta_auc = NA_real_,
      estimated_se = NA_real_, reject_equal_auc = NA_integer_,
      degenerate_equal_stacks = NA_integer_, converged = 0L,
      error = as.character(fit)
    ))
  }
  reject <- if (fit$se == 0) {
    fit$delta != 0
  } else {
    abs(fit$delta / fit$se) > stats::qnorm(1 - alpha / 2)
  }
  data.table(
    entry_sd = task$entry_sd, replicate = task$replicate,
    outer_seed = task$outer_seed, bootstrap_seed = task$bootstrap_seed,
    naive_auc = fit$estimate_a, strict_auc = fit$estimate_b,
    delta_auc = fit$delta, estimated_se = fit$se,
    reject_equal_auc = as.integer(reject),
    degenerate_equal_stacks = as.integer(fit$se == 0 && fit$delta == 0),
    converged = 1L, error = ""
  )
}

wall_start <- proc.time()[["elapsed"]]
runs <- parallel::mclapply(
  seq_len(nrow(tasks)), run_task, mc.cores = cores,
  mc.preschedule = TRUE, mc.set.seed = FALSE
)
wall_runtime <- proc.time()[["elapsed"]] - wall_start
if (any(vapply(runs, inherits, logical(1L), "try-error")))
  stop("At least one time-scale diagnostic task failed outside capture.")
raw <- rbindlist(runs)
stopifnot(nrow(raw) == nrow(tasks))
failures <- raw[converged == 0L, .(failure_count = .N), by = .(entry_sd, error)]
successful <- raw[converged == 1L]
summary <- successful[, .(
  successful_replicates = .N,
  mean_delta_auc = mean(delta_auc),
  empirical_se = stats::sd(delta_auc),
  mean_estimated_se = mean(estimated_se),
  rejection_probability = mean(reject_equal_auc),
  rejection_mcse = sqrt(mean(reject_equal_auc) *
                          (1 - mean(reject_equal_auc)) / .N),
  equal_stack_probability = mean(degenerate_equal_stacks)
), by = entry_sd]
failure_counts <- raw[, .(failed_replicates = sum(converged == 0L)), by = entry_sd]
summary <- merge(summary, failure_counts, by = "entry_sd")
setorder(summary, entry_sd)
stopifnot(nrow(summary) == length(entry_sds),
          all(summary$successful_replicates + summary$failed_replicates ==
                outer_replicates),
          all(is.finite(as.matrix(summary))))

fwrite(raw, project_path("results", "timescale_diagnostic_raw.csv"))
fwrite(summary, project_path("results", "timescale_diagnostic_summary.csv"))
fwrite(failures, project_path("results", "timescale_diagnostic_failures.csv"))
fwrite(tasks, project_path("results", "timescale_diagnostic_seed_manifest.csv"))

trace_values(
  report = "05_simulation",
  labels = c("diagnostic_replicates_per_entry_sd", "diagnostic_bootstrap_replicates",
             "diagnostic_wall_runtime_seconds", "diagnostic_failures"),
  values = c(outer_replicates, bootstrap_replicates, wall_runtime,
             nrow(raw) - nrow(successful)),
  units = c("replicates", "bootstrap replicates", "seconds", "datasets"),
  script = "analysis/05_timescale_diagnostic.R",
  function_name = "paired_bootstrap_auc_delta_fast",
  seed = "results/timescale_diagnostic_seed_manifest.csv",
  runtime_seconds = wall_runtime
)
fields <- c("successful_replicates", "failed_replicates", "mean_delta_auc",
            "empirical_se", "mean_estimated_se", "rejection_probability",
            "rejection_mcse", "equal_stack_probability")
for (i in seq_len(nrow(summary))) {
  row <- summary[i]
  prefix <- paste0("diagnostic_entry_sd_", gsub("\\.", "p", row$entry_sd), "_")
  trace_values(
    report = "05_simulation", labels = paste0(prefix, fields),
    values = as.numeric(row[, ..fields]),
    units = c("replicates", "replicates", "AUC", "standard error",
              "standard error", "proportion", "Monte Carlo standard error",
              "proportion"),
    script = "analysis/05_timescale_diagnostic.R",
    function_name = "paired_bootstrap_auc_delta_fast/aggregation",
    seed = "results/timescale_diagnostic_seed_manifest.csv",
    runtime_seconds = wall_runtime
  )
}
