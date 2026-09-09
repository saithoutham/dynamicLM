#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(parallel)
})
source(here::here("analysis", "_helpers.R"))
source(here::here("analysis", "_simulation_engine.R"))

analysis_seed <- 20260909L
assert_single_seed(analysis_seed)
outer_replicates <- as.integer(Sys.getenv("AUC_INFORMATIVE_REPLICATES", "1000"))
bootstrap_replicates <- as.integer(Sys.getenv("AUC_INFORMATIVE_BOOTSTRAPS", "100"))
allow_small <- identical(Sys.getenv("SIM_ALLOW_SMALL", "0"), "1")
if (!allow_small) stopifnot(outer_replicates >= 1000L,
                            bootstrap_replicates == 100L)
stopifnot(outer_replicates >= 1L, bootstrap_replicates >= 2L)
cores <- as.integer(Sys.getenv(
  "SIM_CORES", as.character(min(8L, parallel::detectCores()))
))
stopifnot(cores >= 1L)

truth_beta <- log(1.5)
lambda0 <- 0.01
entry_mean <- 45
prediction_window <- 10
landmarks <- c(50, 54, 58)
study_end <- 68
sample_size <- 750L
nominal_coverage <- 0.95
methods <- c("naive", "strict", "delayed")
quadrature_tolerance <- 1e-5

logspace_subtract <- function(a, b) {
  stopifnot(length(a) == length(b), all(a >= b))
  a + log1p(-exp(b - a))
}

normal_interval_log_probability <- function(lower_z, upper_z) {
  stopifnot(length(lower_z) == length(upper_z),
            all(is.finite(lower_z)), all(is.finite(upper_z)),
            all(lower_z < upper_z))
  out <- numeric(length(lower_z))
  left <- upper_z <= 0
  right <- lower_z >= 0
  middle <- !(left | right)
  if (any(left)) {
    upper_log <- stats::pnorm(upper_z[left], log.p = TRUE)
    lower_log <- stats::pnorm(lower_z[left], log.p = TRUE)
    out[left] <- logspace_subtract(upper_log, lower_log)
  }
  if (any(right)) {
    lower_log <- stats::pnorm(
      lower_z[right], lower.tail = FALSE, log.p = TRUE
    )
    upper_log <- stats::pnorm(
      upper_z[right], lower.tail = FALSE, log.p = TRUE
    )
    out[right] <- logspace_subtract(lower_log, upper_log)
  }
  if (any(middle)) {
    probability <- stats::pnorm(upper_z[middle]) -
      stats::pnorm(lower_z[middle])
    stopifnot(all(is.finite(probability)), all(probability > 0))
    out[middle] <- log(probability)
  }
  stopifnot(all(is.finite(out)))
  out
}

truncated_entry_probability <- function(landmark, conditional_mean, entry_sd,
                                        lower = 25,
                                        upper = study_end - 1e-6) {
  stopifnot(landmark > lower, landmark < upper, entry_sd > 0,
            all(is.finite(conditional_mean)))
  lower_z <- (lower - conditional_mean) / entry_sd
  landmark_z <- (landmark - conditional_mean) / entry_sd
  upper_z <- (upper - conditional_mean) / entry_sd
  log_numerator <- normal_interval_log_probability(lower_z, landmark_z)
  log_denominator <- normal_interval_log_probability(lower_z, upper_z)
  probability <- exp(log_numerator - log_denominator)
  numerical_tolerance <- sqrt(.Machine$double.eps)
  stopifnot(all(is.finite(probability)),
            all(probability >= -numerical_tolerance),
            all(probability <= 1 + numerical_tolerance))
  pmax(0, pmin(1, probability))
}

numerical_strict_target_auc <- function(landmark, entry_sd, gamma, theta,
                                        delta, grid_step) {
  stopifnot(landmark %in% landmarks, entry_sd > 0, grid_step > 0,
            all(is.finite(c(gamma, theta, delta))))
  if (gamma == 0 && theta == 0 && delta == 0) {
    # With no entry predictor, the landmark entry CDF is a common factor in
    # case and control masses and cancels exactly.  Regime-C delta values stay
    # on the executed two-dimensional quadrature path so reruns reproduce the
    # committed numerical truth table.
    return(numerical_true_auc(
      landmark, prediction_window, truth_beta, lambda0,
      grid_step = grid_step, limit = 6
    ))
  }
  x <- seq(-6, 6, by = grid_step)
  if (theta == 0 && delta == 0) {
    eta <- truth_beta * x
    entry_probability <- truncated_entry_probability(
      landmark, entry_mean + gamma * x, entry_sd
    )
    source_mass <- stats::dnorm(x) * grid_step
  } else {
    # Frailty cases require two-dimensional deterministic normal-grid
    # quadrature.  expand.grid ordering is immaterial because the AUC helper
    # sorts on the risk score.
    grid <- CJ(x = x, u = x)
    eta <- truth_beta * grid$x + theta * grid$u
    entry_probability <- truncated_entry_probability(
      landmark, entry_mean + gamma * grid$x + delta * grid$u, entry_sd
    )
    source_mass <- stats::dnorm(grid$x) * stats::dnorm(grid$u) * grid_step^2
  }
  hazard <- lambda0 * exp(eta)
  survival_landmark <- exp(-hazard * landmark)
  survival_horizon <- exp(-hazard * (landmark + prediction_window))
  case_mass <- source_mass * entry_probability *
    (survival_landmark - survival_horizon)
  control_mass <- source_mass * entry_probability * survival_horizon
  stopifnot(all(is.finite(case_mass)), all(is.finite(control_mass)),
            all(case_mass >= 0), all(control_mass >= 0),
            sum(case_mass) > 0, sum(control_mass) > 0)
  weighted_auc_fast(eta, case_mass, control_mass)
}

calibration <- fread(project_path("results", "informative_entry_calibration.csv"))
calibration <- calibration[
  (phase == "focused" | phase == "frailty") & landmark_interval == 4 &
    censoring_target == 0.30
]
design <- unique(calibration[, .(
  phase, gamma, theta, delta, entry_sd, censor_rate, study_end
)])
design[, regime := fifelse(phase == "frailty", "C",
                           fifelse(gamma == 0, "A", "B"))]
setorder(design, phase, gamma, theta, delta, entry_sd)
design[, auc_cell_id := .I]
stopifnot(nrow(design) == 24L, all(design$study_end == study_end),
          !anyDuplicated(design[, .(phase, gamma, theta, delta, entry_sd)]))

truth_rows <- lapply(seq_len(nrow(design)), function(i) {
  row <- design[i]
  rbindlist(lapply(landmarks, function(landmark) {
    fine <- numerical_strict_target_auc(
      landmark, row$entry_sd, row$gamma, row$theta, row$delta,
      grid_step = if (row$theta == 0 && row$delta == 0) 0.00125 else 0.025
    )
    coarse <- numerical_strict_target_auc(
      landmark, row$entry_sd, row$gamma, row$theta, row$delta,
      grid_step = if (row$theta == 0 && row$delta == 0) 0.0025 else 0.05
    )
    data.table(
      auc_cell_id = row$auc_cell_id, phase = row$phase, regime = row$regime,
      gamma = row$gamma, theta = row$theta, delta = row$delta,
      entry_sd = row$entry_sd, landmark = landmark,
      strict_target_auc = fine, coarse_auc = coarse,
      quadrature_difference = fine - coarse
    )
  }))
})
truth <- rbindlist(truth_rows)
truth[, strict_target_summary_auc := mean(strict_target_auc), by = auc_cell_id]
if (any(abs(truth$quadrature_difference) >= quadrature_tolerance))
  stop("Informative-entry AUC quadrature failed its predeclared tolerance.")
stopifnot(nrow(truth) == nrow(design) * length(landmarks),
          all(truth$strict_target_auc >= 0 & truth$strict_target_auc <= 1),
          all(is.finite(truth$strict_target_summary_auc)))

main_seeds <- fread(project_path("results", "informative_entry_seed_manifest.csv"))
phase_cell_map <- unique(fread(
  project_path("results", "informative_entry_summary.csv")
)[phase %in% c("focused", "frailty"), .(
  phase, cell_id, gamma, theta, delta, entry_sd
)])
tasks <- merge(
  phase_cell_map, main_seeds[replicate <= outer_replicates],
  by = c("phase", "cell_id"), all.x = TRUE
)
tasks <- merge(
  tasks, design,
  by = c("phase", "gamma", "theta", "delta", "entry_sd"), all.x = TRUE,
  suffixes = c("", ".design")
)
tasks <- merge(
  tasks, unique(truth[, .(auc_cell_id, strict_target_summary_auc)]),
  by = "auc_cell_id", all.x = TRUE
)
setorder(tasks, auc_cell_id, replicate)
tasks[, bootstrap_seed_naive := as.integer(
  1800000000 + auc_cell_id * 5000L + 1000L + replicate
)]
tasks[, bootstrap_seed_strict := as.integer(
  1800000000 + auc_cell_id * 5000L + 2000L + replicate
)]
tasks[, bootstrap_seed_delayed := as.integer(
  1800000000 + auc_cell_id * 5000L + 3000L + replicate
)]
stopifnot(nrow(tasks) == nrow(design) * outer_replicates,
          !anyNA(tasks),
          !anyDuplicated(tasks[, .(auc_cell_id, replicate)]),
          !anyDuplicated(tasks$bootstrap_seed_naive),
          !anyDuplicated(tasks$bootstrap_seed_strict),
          !anyDuplicated(tasks$bootstrap_seed_delayed))

run_task <- function(i) {
  task <- tasks[i]
  set.seed(task$seed)
  cohort <- generate_left_truncated_cohort(
    sample_size, truth_beta, lambda0, entry_mean, task$entry_sd, study_end,
    gamma = task$gamma, theta = task$theta, delta = task$delta
  )
  cohort <- apply_independent_censoring(cohort, task$censor_rate, study_end)
  results <- lapply(methods, function(method) {
    stack <- make_simulation_stack(cohort, landmarks, prediction_window, method)
    frailty <- if ("u" %in% names(stack)) stack$u else 0
    eta <- truth_beta * stack$x + task$theta * frailty
    interval_length <- stack$landmark + prediction_window - stack$start
    stopifnot(all(interval_length > 0), all(interval_length <= prediction_window))
    true_risk <- 1 - exp(-lambda0 * interval_length * exp(eta))
    bootstrap_seed <- task[[paste0("bootstrap_seed_", method)]]
    fitted <- tryCatch(
      bootstrap_summary_auc_fast(
        stack, true_risk, prediction_window, bootstrap_replicates,
        bootstrap_seed
      ),
      error = function(error)
        structure(conditionMessage(error), class = "informative_auc_error")
    )
    common <- data.table(
      auc_cell_id = task$auc_cell_id, phase = task$phase,
      regime = task$regime, gamma = task$gamma, theta = task$theta,
      delta = task$delta, entry_sd = task$entry_sd,
      replicate = task$replicate, outer_seed = task$seed,
      bootstrap_seed = bootstrap_seed, method = method,
      strict_target_summary_auc = task$strict_target_summary_auc,
      rows = nrow(stack), events = sum(stack$event),
      subjects = uniqueN(stack$id)
    )
    if (inherits(fitted, "informative_auc_error")) {
      return(cbind(common, data.table(
        estimate = NA_real_, estimated_se = NA_real_, covered = NA_integer_,
        reject_half = NA_integer_, converged = 0L, error = as.character(fitted)
      )))
    }
    z <- stats::qnorm(1 - (1 - nominal_coverage) / 2)
    cbind(common, data.table(
      estimate = fitted$estimate, estimated_se = fitted$se,
      covered = as.integer(abs(fitted$estimate -
        task$strict_target_summary_auc) <= z * fitted$se),
      reject_half = as.integer(abs(fitted$estimate - 0.5) > z * fitted$se),
      converged = 1L, error = ""
    ))
  })
  rbindlist(results)
}

started <- proc.time()[["elapsed"]]
runs <- mclapply(
  seq_len(nrow(tasks)), run_task, mc.cores = cores,
  mc.preschedule = TRUE, mc.set.seed = FALSE
)
wall_runtime <- proc.time()[["elapsed"]] - started
if (any(vapply(runs, inherits, logical(1L), "try-error")))
  stop("At least one informative AUC worker failed outside condition capture.")
raw <- rbindlist(runs)
stopifnot(nrow(raw) == nrow(tasks) * length(methods),
          raw[, all(.N == outer_replicates), by = .(auc_cell_id, method)]$V1)

grouping <- c("auc_cell_id", "phase", "regime", "gamma", "theta", "delta",
              "entry_sd", "method", "strict_target_summary_auc")
successful <- raw[converged == 1L]
summary <- successful[, .(
  successful_replicates = .N,
  bias = mean(estimate - strict_target_summary_auc),
  empirical_se = stats::sd(estimate),
  mean_estimated_se = mean(estimated_se),
  se_ratio = mean(estimated_se) / stats::sd(estimate),
  coverage = mean(covered),
  coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / .N),
  power_against_half = mean(reject_half),
  power_mcse = sqrt(mean(reject_half) * (1 - mean(reject_half)) / .N),
  mean_rows = mean(rows), mean_events = mean(events),
  mean_subjects = mean(subjects)
), by = grouping]
failure_counts <- raw[, .(failed_replicates = sum(converged == 0L)),
                      by = grouping]
summary <- merge(summary, failure_counts, by = grouping, all.x = TRUE)
setorder(summary, auc_cell_id, method)
numeric_fields <- setdiff(
  names(summary), c(grouping, "successful_replicates", "failed_replicates")
)
stopifnot(nrow(summary) == nrow(design) * length(methods),
          all(summary$successful_replicates + summary$failed_replicates ==
                outer_replicates),
          all(is.finite(as.matrix(summary[, ..numeric_fields]))))
failures <- raw[converged == 0L,
                .(failure_count = .N), by = .(auc_cell_id, method, error)]

coefficient <- fread(
  project_path("results", "informative_entry_summary.csv")
)[phase %chin% c("focused", "frailty")]
comparison <- merge(
  coefficient[, .(
    coefficient_cell_id = cell_id,
    phase, regime, gamma, theta, delta, entry_sd, method,
    coefficient_bias = bias,
    coefficient_empirical_se = empirical_se,
    coefficient_mean_estimated_se = mean_estimated_se,
    coefficient_coverage = coverage,
    coefficient_coverage_mcse = coverage_mcse
  )],
  summary[, .(
    auc_cell_id,
    phase, regime, gamma, theta, delta, entry_sd, method,
    auc_target = strict_target_summary_auc,
    auc_bias = bias,
    auc_empirical_se = empirical_se,
    auc_mean_bootstrap_se = mean_estimated_se,
    auc_coverage = coverage,
    auc_coverage_mcse = coverage_mcse
  )],
  by = c("phase", "regime", "gamma", "theta", "delta", "entry_sd", "method"),
  all = FALSE
)
stopifnot(nrow(comparison) == nrow(summary), !anyNA(comparison))
setorder(comparison, auc_cell_id, method)

fwrite(raw, project_path("results", "informative_auc_coverage_raw.csv"))
fwrite(summary,
       project_path("results", "informative_auc_coverage_summary.csv"))
fwrite(comparison,
       project_path("results", "informative_prediction_estimation_comparison.csv"))
fwrite(truth, project_path("results", "informative_auc_truth.csv"))
fwrite(failures,
       project_path("results", "informative_auc_coverage_failures.csv"))
fwrite(tasks[, .(
  auc_cell_id, phase, regime, gamma, theta, delta, entry_sd, replicate,
  outer_seed = seed, bootstrap_seed_naive, bootstrap_seed_strict,
  bootstrap_seed_delayed
)], project_path("results", "informative_auc_seed_manifest.csv"))

trace_values(
  report = "07_informative_entry",
  labels = c(
    "auc_outer_replicates_per_cell", "auc_bootstrap_replicates",
    "auc_design_cells", "auc_methods", "auc_total_method_datasets",
    "auc_total_bootstrap_resamples", "auc_failures",
    "auc_quadrature_tolerance", "auc_max_quadrature_difference",
    "auc_runtime_seconds"
  ),
  values = c(
    outer_replicates, bootstrap_replicates, nrow(design), length(methods),
    nrow(raw), nrow(raw) * bootstrap_replicates,
    nrow(raw) - nrow(successful), quadrature_tolerance,
    max(abs(truth$quadrature_difference)), wall_runtime
  ),
  units = c(
    "replicates", "resamples", "cells", "methods", "datasets",
    "resamples", "datasets", "AUC", "AUC", "seconds"
  ),
  script = "analysis/07_auc_informative.R",
  function_name = "bootstrap_summary_auc_fast",
  seed = "results/informative_auc_seed_manifest.csv",
  runtime_seconds = wall_runtime
)

fields <- c(
  "strict_target_summary_auc", "successful_replicates", "failed_replicates",
  "bias", "empirical_se", "mean_estimated_se", "se_ratio", "coverage",
  "coverage_mcse", "power_against_half", "power_mcse", "mean_rows",
  "mean_events", "mean_subjects"
)
for (i in seq_len(nrow(summary))) {
  row <- summary[i]
  prefix <- paste0("auc_cell_", sprintf("%03d", row$auc_cell_id), "_",
                   row$method, "_")
  trace_values(
    report = "07_informative_entry", labels = paste0(prefix, fields),
    values = as.numeric(row[, ..fields]), units = fields,
    script = "analysis/07_auc_informative.R",
    function_name = "bootstrap_summary_auc_fast/aggregation",
    seed = "results/informative_auc_seed_manifest.csv",
    runtime_seconds = wall_runtime
  )
}

for (i in seq_len(nrow(truth))) {
  row <- truth[i]
  prefix <- paste0("auc_truth_cell_", sprintf("%03d", row$auc_cell_id),
                   "_landmark_", row$landmark, "_")
  trace_values(
    report = "07_informative_entry",
    labels = paste0(prefix, c("auc", "summary_auc", "quadrature_difference")),
    values = c(row$strict_target_auc, row$strict_target_summary_auc,
               row$quadrature_difference),
    units = "AUC", script = "analysis/07_auc_informative.R",
    function_name = "numerical_strict_target_auc", seed = NA_integer_,
    runtime_seconds = wall_runtime
  )
}

stopifnot(
  file.exists(project_path("results", "informative_auc_coverage_summary.csv")),
  file.exists(project_path("results", "informative_auc_seed_manifest.csv")),
  file.exists(project_path("results", "informative_auc_truth.csv"))
)
