#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
  library(here)
})
source(here::here("analysis", "_helpers.R"))
devtools::load_all(here::here(), quiet = TRUE)
source(here::here("analysis", "_baseline_pipeline.R"))

baseline_seed <- 20260908L
run <- timed(run_pbc_baseline(seed = baseline_seed))
baseline <- run$value
stopifnot(identical(baseline$seed, baseline_seed))

score_output <- data.table::rbindlist(lapply(names(baseline$score_tables), function(section) {
  table <- as.data.frame(baseline$score_tables[[section]])
  table$section <- section
  table
}), fill = TRUE)
score_output <- as.data.frame(score_output)
row.names(score_output) <- NULL
utils::write.csv(score_output, project_path("results", "baseline_metrics.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(baseline$conditions,
                 project_path("results", "baseline_conditions.csv"),
                 row.names = FALSE, na = "")
saveRDS(baseline, project_path("results", "baseline_expected.rds"), version = 3)

trace_values(
  report = "02_baseline",
  labels = c("pbc_long_rows", "pbc_unique_ids", "terminal_visit_ties", "stacked_rows",
             "landmark_count", "prediction_window", "warning_count",
             "unique_warning_patterns", "pen_lm_warning_count",
             "cv_pen_lm_warning_count"),
  values = c(baseline$pbc_rows, baseline$pbc_ids, baseline$terminal_visit_ties,
             baseline$lm_rows,
             length(baseline$landmarks), baseline$window,
             nrow(baseline$conditions), nrow(unique(baseline$conditions)),
             sum(baseline$conditions$stage == "pen_lm"),
             sum(baseline$conditions$stage == "cv.pen_lm")),
  units = c("rows", "subjects", "rows", "rows", "landmarks", "years",
            "warnings", "patterns", "warnings", "warnings"),
  script = "analysis/02_baseline.R",
  function_name = "run_pbc_baseline",
  seed = baseline_seed,
  runtime_seconds = run$runtime_seconds
)
trace_values(
  report = "02_baseline",
  labels = paste0("risk_set_lm", baseline$landmarks),
  values = baseline$risk_counts,
  units = "rows",
  script = "analysis/02_baseline.R",
  function_name = "dynamicLM::stack_data",
  seed = baseline_seed,
  runtime_seconds = run$runtime_seconds
)
trace_values(
  report = "02_baseline",
  labels = c("path_lambda_count_cause1", "path_lambda_count_cause2",
             "cv_lambda_min_cause1", "cv_lambda_min_cause2",
             "cv_lambda_1se_cause1", "cv_lambda_1se_cause2",
             "prediction_min_LM", "prediction_max_LM",
             "prediction_min_penLM", "prediction_max_penLM"),
  values = c(baseline$path_lambda_counts, baseline$cv_lambda_min,
             baseline$cv_lambda_1se, as.vector(t(baseline$prediction_ranges))),
  units = c("lambdas", "lambdas", rep("lambda", 4L), rep("risk", 4L)),
  script = "analysis/02_baseline.R",
  function_name = "pen_lm/cv.pen_lm/predict.dynamicLM",
  seed = baseline_seed,
  runtime_seconds = run$runtime_seconds
)

metric_fields <- intersect(c("AUC", "Brier", "se", "lower", "upper"),
                           names(score_output))
for (i in seq_len(nrow(score_output))) {
  fields <- metric_fields[!is.na(unlist(score_output[i, metric_fields, drop = FALSE]))]
  if (length(fields) == 0L) next
  model <- gsub("[^[:alnum:]]+", "_", as.character(score_output$model[i]))
  landmark <- if ("tLM" %in% names(score_output) && !is.na(score_output$tLM[i])) {
    paste0("_lm", score_output$tLM[i])
  } else {
    ""
  }
  trace_values(
    report = "02_baseline",
    labels = paste0(score_output$section[i], "_", model, landmark, "_", fields),
    values = as.numeric(score_output[i, fields, drop = TRUE]),
    units = ifelse(fields == "se", "standard error", "probability scale"),
    script = "analysis/02_baseline.R",
    function_name = "dynamicLM::score",
    seed = baseline_seed,
    runtime_seconds = run$runtime_seconds
  )
}

stopifnot(
  file.exists(project_path("results", "baseline_metrics.csv")),
  file.exists(project_path("results", "baseline_expected.rds")),
  file.exists(project_path("results", "baseline_conditions.csv"))
)
