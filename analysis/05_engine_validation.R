#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(devtools)
  library(here)
  library(survival)
})
source(here::here("analysis", "_helpers.R"))
source(here::here("analysis", "_simulation_engine.R"))
devtools::load_all(here::here(), quiet = TRUE)

analysis_seed <- 20260909L
assert_single_seed(analysis_seed)
run <- timed({
  set.seed(analysis_seed)
  cohort <- generate_left_truncated_cohort(
    n = 500L, beta = log(1.5), lambda0 = 0.01,
    entry_mean = 45, entry_sd = 8, study_end = 72
  )
  cohort <- apply_independent_censoring(cohort, censor_rate = 0.03,
                                        study_end = 72)
  comparisons <- rbindlist(lapply(c("naive", "strict", "delayed"), function(method) {
    stack <- make_simulation_stack(cohort, c(50, 56, 62), 10, method)
    fast <- fit_fast_cluster_cox(stack)
    reference <- survival::coxph(
      survival::Surv(start, stop, event) ~ x + cluster(id),
      data = stack, ties = "breslow", robust = TRUE
    )
    reference_estimate <- as.numeric(stats::coef(reference))
    reference_se <- as.numeric(sqrt(stats::vcov(reference)))
    stopifnot(is.finite(reference_estimate), is.finite(reference_se))
    data.table(
      method = method,
      fast_estimate = fast$estimate,
      coxph_estimate = reference_estimate,
      estimate_difference = fast$estimate - reference_estimate,
      fast_robust_se = fast$robust_se,
      coxph_robust_se = reference_se,
      robust_se_difference = fast$robust_se - reference_se,
      score_at_solution = fast$score,
      rows = fast$rows,
      events = fast$events
    )
  }))
  stopifnot(max(abs(comparisons$estimate_difference)) < 1e-8,
            max(abs(comparisons$robust_se_difference)) < 1e-8,
            max(abs(comparisons$score_at_solution)) < 1e-5)

  # Validate the fast product-limit/AUC path against the package implementation.
  strict_stack <- make_simulation_stack(cohort, c(50, 56, 62), 10, "strict")
  risk <- 1 - exp(-0.01 * 10 * exp(log(1.5) * strict_stack$x))
  fast_auc <- summary_auc_fast(strict_stack, risk, 10)
  package_auc <- dynamicLM::score_left_truncated(
    risk, strict_stack, time_col = "stop", status_col = "event",
    entry_col = "start", landmark_col = "landmark", id_col = "id",
    cause = 1, w = 10
  )$summary$AUC
  stopifnot(abs(fast_auc - package_auc) < 1e-12)
  list(comparisons = comparisons, fast_auc = fast_auc,
       package_auc = package_auc)
})

utils::write.csv(run$value$comparisons,
                 project_path("results", "simulation_engine_validation.csv"),
                 row.names = FALSE)
trace_values(
  report = "05_simulation",
  labels = c("engine_validation_max_estimate_difference",
             "engine_validation_max_robust_se_difference",
             "engine_validation_max_score_residual",
             "engine_validation_auc_difference"),
  values = c(max(abs(run$value$comparisons$estimate_difference)),
             max(abs(run$value$comparisons$robust_se_difference)),
             max(abs(run$value$comparisons$score_at_solution)),
             abs(run$value$fast_auc - run$value$package_auc)),
  units = c("log hazard ratio", "standard error", "score", "AUC"),
  script = "analysis/05_engine_validation.R",
  function_name = "fit_fast_cluster_cox/score_left_truncated validation",
  seed = analysis_seed, runtime_seconds = run$runtime_seconds
)
simulation_grid <- data.table::fread(project_path("results", "simulation_grid.csv"))
calibration <- data.table::fread(project_path("results", "simulation_calibration.csv"))
trace_values(
  report = "05_simulation",
  labels = c(paste0("design_landmark_interval_", seq_along(sort(unique(simulation_grid$landmark_interval)))),
             paste0("design_censoring_target_", seq_along(sort(unique(simulation_grid$censoring_target)))),
             paste0("design_sample_size_", seq_along(sort(unique(simulation_grid$n)))),
             paste0("design_entry_sd_", seq_along(sort(unique(simulation_grid$entry_sd)))),
             "design_cells_per_entry_sd", "calibration_pilot_n",
             "auc_design_landmark_interval", "auc_design_censoring_target",
             "auc_design_sample_size"),
  values = c(sort(unique(simulation_grid$landmark_interval)),
             sort(unique(simulation_grid$censoring_target)),
             sort(unique(simulation_grid$n)), sort(unique(simulation_grid$entry_sd)),
             nrow(simulation_grid) / length(unique(simulation_grid$entry_sd)),
             unique(calibration$pilot_n), 4, 0.30, 750),
  units = c(rep("years", 3), rep("proportion", 4), rep("subjects", 3),
            rep("years", 3), "cells", "subjects", "years", "proportion",
            "subjects"),
  script = "analysis/05_engine_validation.R",
  function_name = "simulation design audit", seed = NA_integer_,
  runtime_seconds = run$runtime_seconds
)
stopifnot(file.exists(project_path("results", "simulation_engine_validation.csv")))
