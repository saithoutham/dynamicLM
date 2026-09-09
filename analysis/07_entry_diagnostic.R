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
replicates <- as.integer(Sys.getenv("PHASE6_REPLICATES", "1000"))
allow_small <- identical(Sys.getenv("SIM_ALLOW_SMALL", "0"), "1")
if (!allow_small) stopifnot(replicates >= 1000L)
stopifnot(replicates >= 1L)
cores <- as.integer(Sys.getenv(
  "SIM_CORES", as.character(min(8L, parallel::detectCores()))
))
stopifnot(cores >= 1L)

truth_beta <- log(1.5)
lambda0 <- 0.01
entry_mean <- 45
prediction_window <- 10
landmarks <- c(50, 54, 58)
alpha <- 0.05
bonferroni_alpha <- alpha / length(landmarks)

calibration <- fread(project_path("results", "informative_entry_calibration.csv"))
calibration <- calibration[
  (phase == "focused" | phase == "frailty") & landmark_interval == 4 &
    censoring_target == 0.30
]
design_columns <- c("phase", "gamma", "theta", "delta", "entry_sd",
                    "study_end", "censor_rate")
design <- unique(calibration[, ..design_columns])
design[, regime := fifelse(phase == "frailty", "C",
                           fifelse(gamma == 0, "A", "B"))]
setorder(design, phase, gamma, theta, delta, entry_sd)
design[, diagnostic_cell_id := .I]
stopifnot(nrow(design) == 24L, !anyDuplicated(design[, ..design_columns]))

main_seeds <- fread(project_path("results", "informative_entry_seed_manifest.csv"))
phase_cell_map <- unique(fread(
  project_path("results", "informative_entry_summary.csv")
)[phase %in% c("focused", "frailty"), .(
  phase, cell_id, regime, gamma, theta, delta, entry_sd
)])
tasks <- merge(
  phase_cell_map, main_seeds[replicate <= replicates],
  by = c("phase", "cell_id"), all.x = TRUE
)
tasks <- merge(
  tasks, design,
  by = c("phase", "regime", "gamma", "theta", "delta", "entry_sd"),
  all.x = TRUE
)
setorder(tasks, diagnostic_cell_id, replicate)
stopifnot(nrow(tasks) == nrow(design) * replicates,
          !anyNA(tasks),
          !anyDuplicated(tasks[, .(diagnostic_cell_id, replicate)]))

safe_correlation_test <- function(x, y) {
  stopifnot(length(x) == length(y), length(x) >= 3L,
            all(is.finite(x)), all(is.finite(y)))
  if (stats::sd(x) == 0 || stats::sd(y) == 0)
    return(list(estimate = 0, p.value = 1))
  result <- suppressWarnings(stats::cor.test(x, y, method = "pearson"))
  stopifnot(is.finite(result$estimate), is.finite(result$p.value),
            result$p.value >= 0, result$p.value <= 1)
  list(estimate = unname(result$estimate), p.value = result$p.value)
}

run_task <- function(i) {
  task <- tasks[i]
  set.seed(task$seed)
  cohort <- generate_left_truncated_cohort(
    750L, truth_beta, lambda0, entry_mean, task$entry_sd, task$study_end,
    gamma = task$gamma, theta = task$theta, delta = task$delta
  )
  cohort <- apply_independent_censoring(
    cohort, task$censor_rate, task$study_end
  )
  entry_x <- safe_correlation_test(cohort$entry, cohort$x)
  entry_u_correlation <- if ("u" %in% names(cohort))
    stats::cor(cohort$entry, cohort$u) else NA_real_
  stopifnot(is.na(entry_u_correlation) || is.finite(entry_u_correlation))

  landmark_rows <- lapply(seq_along(landmarks), function(j) {
    landmark <- landmarks[j]
    naive <- cohort$exit > landmark
    strict <- cohort$entry <= landmark & cohort$exit > landmark
    preentry <- cohort$entry > landmark & cohort$exit > landmark
    stopifnot(all((strict | preentry) == naive),
              sum(naive) == sum(strict) + sum(preentry), sum(naive) > 0L)
    fraction <- sum(preentry) / sum(naive)
    naive_strict_difference <- mean(cohort$x[naive]) - mean(cohort$x[strict])
    if (sum(preentry) >= 2L && sum(strict) >= 2L &&
        stats::sd(cohort$x[preentry]) > 0 && stats::sd(cohort$x[strict]) > 0) {
      # The two tested groups are disjoint.  Under equal component means,
      # mean(naive)-mean(strict) is zero because naive=strict union pre-entry.
      test <- stats::t.test(cohort$x[preentry], cohort$x[strict],
                            var.equal = FALSE)
      p_value <- test$p.value
    } else {
      p_value <- 1
    }
    stopifnot(is.finite(fraction), fraction >= 0, fraction <= 1,
              is.finite(naive_strict_difference), is.finite(p_value),
              p_value >= 0, p_value <= 1)
    data.table(
      landmark_index = j, landmark = landmark, naive_rows = sum(naive),
      strict_rows = sum(strict), preentry_rows = sum(preentry),
      preentry_fraction = fraction,
      naive_minus_strict_mean_x = naive_strict_difference,
      preentry_vs_strict_p = p_value,
      riskset_shift_detected = as.integer(p_value < bonferroni_alpha)
    )
  })
  rows <- rbindlist(landmark_rows)
  rows[, `:=`(
    diagnostic_cell_id = task$diagnostic_cell_id,
    phase = task$phase, regime = task$regime, replicate = task$replicate,
    seed = task$seed, gamma = task$gamma, theta = task$theta,
    delta = task$delta, entry_sd = task$entry_sd,
    entry_x_correlation = entry_x$estimate,
    entry_x_p = entry_x$p.value,
    entry_x_detected = as.integer(entry_x$p.value < alpha),
    entry_u_correlation_oracle = entry_u_correlation
  )]
  rows[, structural_difference_detected :=
         as.integer(any(preentry_fraction > 0))]
  rows[, any_riskset_shift_detected :=
         as.integer(any(riskset_shift_detected == 1L))]
  rows[, combined_observed_diagnostic := as.integer(
    entry_x_detected == 1L | any_riskset_shift_detected == 1L
  )]
  rows
}

started <- proc.time()[["elapsed"]]
runs <- mclapply(
  seq_len(nrow(tasks)), run_task, mc.cores = cores,
  mc.preschedule = TRUE, mc.set.seed = FALSE
)
wall_runtime <- proc.time()[["elapsed"]] - started
if (any(vapply(runs, inherits, logical(1L), "try-error")))
  stop("At least one entry-diagnostic worker failed.")
raw <- rbindlist(runs)
stopifnot(nrow(raw) == nrow(tasks) * length(landmarks),
          raw[, all(.N == length(landmarks)),
              by = .(diagnostic_cell_id, replicate)]$V1)

cell_group <- c("diagnostic_cell_id", "phase", "regime", "gamma", "theta",
                "delta", "entry_sd")
subject_level <- unique(raw[, c(
  cell_group, "replicate", "entry_x_correlation", "entry_x_p",
  "entry_x_detected", "entry_u_correlation_oracle",
  "structural_difference_detected", "any_riskset_shift_detected",
  "combined_observed_diagnostic"
), with = FALSE])
stopifnot(nrow(subject_level) == nrow(tasks))

overall <- subject_level[, .(
  replicates = .N,
  mean_entry_x_correlation = mean(entry_x_correlation),
  entry_x_detection_probability = mean(entry_x_detected),
  entry_x_detection_mcse = sqrt(mean(entry_x_detected) *
    (1 - mean(entry_x_detected)) / .N),
  mean_entry_u_correlation_oracle = if (all(is.na(entry_u_correlation_oracle)))
    NA_real_ else mean(entry_u_correlation_oracle, na.rm = TRUE),
  structural_detection_probability = mean(structural_difference_detected),
  structural_detection_mcse = sqrt(mean(structural_difference_detected) *
    (1 - mean(structural_difference_detected)) / .N),
  any_riskset_shift_detection_probability = mean(any_riskset_shift_detected),
  any_riskset_shift_detection_mcse = sqrt(mean(any_riskset_shift_detected) *
    (1 - mean(any_riskset_shift_detected)) / .N),
  combined_detection_probability = mean(combined_observed_diagnostic),
  combined_detection_mcse = sqrt(mean(combined_observed_diagnostic) *
    (1 - mean(combined_observed_diagnostic)) / .N)
), by = cell_group]

landmark_summary <- raw[, .(
  mean_naive_rows = mean(naive_rows),
  mean_strict_rows = mean(strict_rows),
  mean_preentry_rows = mean(preentry_rows),
  mean_preentry_fraction = mean(preentry_fraction),
  mean_naive_minus_strict_mean_x = mean(naive_minus_strict_mean_x),
  riskset_shift_detection_probability = mean(riskset_shift_detected),
  riskset_shift_detection_mcse = sqrt(mean(riskset_shift_detected) *
    (1 - mean(riskset_shift_detected)) / .N)
), by = c(cell_group, "landmark_index", "landmark")]

summary <- merge(overall, landmark_summary, by = cell_group, all = TRUE)
setorder(summary, diagnostic_cell_id, landmark_index)
stopifnot(nrow(summary) == nrow(design) * length(landmarks),
          all(summary$replicates == replicates),
          all(is.finite(as.matrix(summary[, .(
            mean_entry_x_correlation, entry_x_detection_probability,
            entry_x_detection_mcse, structural_detection_probability,
            structural_detection_mcse,
            any_riskset_shift_detection_probability,
            any_riskset_shift_detection_mcse,
            riskset_shift_detection_probability,
            riskset_shift_detection_mcse, combined_detection_probability,
            combined_detection_mcse, mean_naive_rows, mean_strict_rows,
            mean_preentry_rows, mean_preentry_fraction,
            mean_naive_minus_strict_mean_x
          )]))))

fwrite(raw, project_path("results", "entry_diagnostic_raw.csv"))
fwrite(summary, project_path("results", "entry_diagnostic_summary.csv"))
fwrite(tasks[, .(diagnostic_cell_id, phase, regime, gamma, theta, delta,
                 entry_sd, replicate, seed)],
       project_path("results", "entry_diagnostic_seed_manifest.csv"))

trace_values(
  report = "07_informative_entry",
  labels = c("diagnostic_replicates_per_cell", "diagnostic_cells",
             "diagnostic_landmarks", "diagnostic_alpha",
             "diagnostic_bonferroni_alpha", "diagnostic_runtime_seconds"),
  values = c(replicates, nrow(design), length(landmarks), alpha,
             bonferroni_alpha, wall_runtime),
  units = c("replicates", "cells", "landmarks", "probability",
            "probability", "seconds"),
  script = "analysis/07_entry_diagnostic.R",
  function_name = "entry-risk-set diagnostic",
  seed = "results/entry_diagnostic_seed_manifest.csv",
  runtime_seconds = wall_runtime
)

fields <- c(
  "replicates", "mean_entry_x_correlation", "entry_x_detection_probability",
  "entry_x_detection_mcse", "structural_detection_probability",
  "structural_detection_mcse", "any_riskset_shift_detection_probability",
  "any_riskset_shift_detection_mcse", "riskset_shift_detection_probability",
  "riskset_shift_detection_mcse", "combined_detection_probability",
  "combined_detection_mcse", "mean_naive_rows", "mean_strict_rows",
  "mean_preentry_rows", "mean_preentry_fraction",
  "mean_naive_minus_strict_mean_x"
)
for (i in seq_len(nrow(summary))) {
  row <- summary[i]
  prefix <- paste0("diagnostic_cell_", sprintf("%03d", row$diagnostic_cell_id),
                   "_lm", row$landmark_index, "_")
  cell_seeds <- tasks[diagnostic_cell_id == row$diagnostic_cell_id]$seed
  trace_values(
    report = "07_informative_entry", labels = paste0(prefix, fields),
    values = as.numeric(row[, ..fields]), units = fields,
    script = "analysis/07_entry_diagnostic.R",
    function_name = "entry-risk-set diagnostic aggregation",
    seed = paste0(min(cell_seeds), "-", max(cell_seeds)),
    runtime_seconds = wall_runtime
  )
  if (is.finite(row$mean_entry_u_correlation_oracle)) {
    trace_values(
      report = "07_informative_entry",
      labels = paste0(prefix, "mean_entry_u_correlation_oracle"),
      values = row$mean_entry_u_correlation_oracle, units = "correlation",
      script = "analysis/07_entry_diagnostic.R",
      function_name = "entry-risk-set diagnostic aggregation",
      seed = paste0(min(cell_seeds), "-", max(cell_seeds)),
      runtime_seconds = wall_runtime
    )
  }
}

stopifnot(file.exists(project_path("results", "entry_diagnostic_summary.csv")),
          file.exists(project_path("results", "entry_diagnostic_seed_manifest.csv")))
