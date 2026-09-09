#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(devtools)
  library(here)
  library(riskRegression)
  library(survival)
})
source(here::here("analysis", "_helpers.R"))
source(here::here("analysis", "_age_data.R"))
devtools::load_all(here::here(), quiet = TRUE)

analysis_seed <- 20260908L
assert_single_seed(analysis_seed)
set.seed(analysis_seed)

add_centered_interactions <- function(lmdata, variables, center, scale = 10) {
  stopifnot(inherits(lmdata, "LMdataframe"), length(center) == 1L,
            length(scale) == 1L, is.finite(center), is.finite(scale), scale > 0)
  lmdata$data$LM_centered <- (lmdata$data$LM - center) / scale
  for (variable in variables) {
    stopifnot(variable %in% names(lmdata$data))
    lmdata$data[[paste0(variable, "_LM_centered")]] <-
      lmdata$data[[variable]] * lmdata$data$LM_centered
  }
  lmdata
}

fit_corrected_model <- function(lmdata, dataset, interaction_variables,
                                base_terms, center) {
  lmdata <- add_centered_interactions(lmdata, interaction_variables, center)
  interaction_terms <- paste0(interaction_variables, "_LM_centered")
  entry <- lmdata$entry_col
  formula <- stats::as.formula(paste0(
    "Surv(", entry, ", exit_age, event) ~ ",
    paste(c(base_terms, interaction_terms, "cluster(id)"), collapse = " + ")
  ))
  fit <- dynamicLM::dynamic_lm(lmdata, formula, "coxph", x = TRUE)
  coefficients <- stats::coef(fit$model)
  robust_se <- sqrt(diag(stats::vcov(fit$model)))
  stopifnot(inherits(fit, "LMcoxph"), all(is.finite(coefficients)),
            all(is.finite(robust_se)), identical(names(coefficients), names(robust_se)),
            identical(fit$entry_col, lmdata$entry_col),
            identical(fit$entry_mode, lmdata$entry_mode))
  list(model = fit, lmdata = lmdata, coefficients = coefficients,
       robust_se = robust_se, dataset = dataset)
}

summarize_model <- function(fitted, terms) {
  fit <- fitted$model$model
  stopifnot(all(terms %in% names(fitted$coefficients)))
  data.frame(
    dataset = fitted$dataset,
    entry_mode = fitted$model$entry_mode,
    observations_used = fit$n,
    events_used = fit$nevent,
    unique_subjects = length(unique(fitted$lmdata$data$id)),
    term = terms,
    estimate = unname(fitted$coefficients[terms]),
    robust_se = unname(fitted$robust_se[terms]),
    stringsAsFactors = FALSE
  )
}

compare_strict_to_riskregression <- function(predictions, data, w,
                                             max_rows = 200L) {
  rows <- lapply(sort(unique(data$LM)), function(landmark) {
    index <- data$LM == landmark
    di <- data[index, , drop = FALSE]
    risk <- predictions[index]
    if (nrow(di) > max_rows) {
      event_rows <- which(di$event > 0L)
      nonevent_rows <- which(di$event == 0L)
      nonevent_rows <- nonevent_rows[order(di$exit_age[nonevent_rows],
                                           decreasing = TRUE)]
      selected <- c(event_rows,
                    head(nonevent_rows, max_rows - length(event_rows)))
      di <- di[selected, , drop = FALSE]
      risk <- risk[selected]
    }
    ours <- dynamicLM::score_left_truncated(
      risk = risk, data = di, time_col = "exit_age", status_col = "event",
      entry_col = "LM", landmark_col = "LM", id_col = "id", cause = 1, w = w
    )$score
    rr <- riskRegression::Score(
      list(model = risk), formula = Surv(exit_age, event) ~ 1,
      data = di, metrics = c("auc", "brier"),
      times = landmark + w - 1e-5, cause = 1, se.fit = FALSE
    )
    stopifnot(nrow(ours) == 1L)
    data.frame(
      landmark = landmark,
      our_auc = ours$AUC,
      riskregression_auc = rr$AUC$score$AUC[1L],
      our_brier = ours$Brier,
      riskregression_brier = rr$Brier$score$Brier[
        rr$Brier$score$model == "model"
      ][1L],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$auc_difference <- out$our_auc - out$riskregression_auc
  out$brier_difference <- out$our_brier - out$riskregression_brier
  stopifnot(all(is.finite(as.matrix(out[, -1L]))))
  out
}

run <- timed({
  pbc <- make_pbc_age()
  pbc$long$log_bili_z <- as.numeric(scale(pbc$long$log_bili))
  pbc$long$albumin_z <- as.numeric(scale(pbc$long$albumin_locf))
  stopifnot(!anyNA(pbc$long[c("log_bili_z", "albumin_z")]))
  nafld <- make_nafld_age(include_preindex_labs = FALSE)
  pbc_landmarks <- c(40, 50, 60)
  nafld_landmarks <- c(40, 50, 60, 70, 80)
  prediction_window <- 10

  pbc_stacks <- lapply(c("strict", "delayed"), function(mode) {
    dynamicLM::stack_data(
      pbc$long, outcome = list(time = "exit_age", status = "event"),
      lms = pbc_landmarks, w = prediction_window,
      covs = list(fixed = c("entry_age", "male", "stage0", "trt0"),
                  varying = c("log_bili_z", "albumin_z")),
      format = "long", id = "id", rtime = "visit_age",
      entry = "entry_age", entry_mode = mode
    )
  })
  names(pbc_stacks) <- c("strict", "delayed")

  nafld_stacks <- lapply(c("strict", "delayed"), function(mode) {
    dynamicLM::stack_data(
      nafld$baseline, outcome = list(time = "exit_age", status = "event"),
      lms = nafld_landmarks, w = prediction_window,
      covs = list(fixed = c("entry_age", "male", "bmi_z"), varying = NULL),
      format = "wide", id = "id", entry = "entry_age", entry_mode = mode
    )
  })
  names(nafld_stacks) <- c("strict", "delayed")

  for (stack in c(pbc_stacks, nafld_stacks)) {
    stopifnot(
      all(stack$data[[stack$entry_col]] >= stack$data$LM),
      all(stack$data[[stack$entry_col]] < stack$data$exit_age),
      !anyDuplicated(stack$data[, c("id", "LM")])
    )
    if (stack$entry_mode == "strict")
      stopifnot(all(stack$data$entry_age <= stack$data$LM),
                identical(stack$entry_col, "LM"))
    if (stack$entry_mode == "delayed")
      stopifnot(all(stack$data$.LM_entry ==
                      pmax(stack$data$LM, stack$data$entry_age)),
                identical(stack$entry_col, ".LM_entry"))
  }

  pbc_models <- lapply(pbc_stacks, fit_corrected_model,
                       dataset = "pbcseq_death",
                       interaction_variables = c("log_bili_z", "albumin_z"),
                       base_terms = c("male", "stage0", "log_bili_z", "albumin_z"),
                       center = mean(pbc_landmarks))
  nafld_models <- lapply(nafld_stacks, fit_corrected_model,
                         dataset = "nafld",
                         interaction_variables = "bmi_z",
                         base_terms = c("male", "bmi_z"),
                         center = mean(nafld_landmarks))

  model_summary <- rbind(
    do.call(rbind, lapply(pbc_models, summarize_model,
                          terms = c("male", "log_bili_z"))),
    do.call(rbind, lapply(nafld_models, summarize_model,
                          terms = c("male", "bmi_z")))
  )
  rownames(model_summary) <- NULL

  scores <- list()
  predictions <- list()
  # The package predictor performs a row-by-row baseline-hazard integration.
  # PBC is sufficient to validate the corrected scoring path here; the much
  # larger NAFLD prediction benchmark is handled separately in Phase 5.
  all_models <- pbc_models
  model_names <- paste0("pbcseq_death_", names(pbc_models))
  names(all_models) <- model_names
  for (name in model_names) {
    fitted <- all_models[[name]]
    pred <- stats::predict(fitted$model)
    stopifnot(inherits(pred, "LMpred"), nrow(pred$preds) == nrow(pred$data),
              all(is.finite(pred$preds$risk)),
              all(pred$preds$risk >= 0 & pred$preds$risk <= 1))
    predictions[[name]] <- list(risk = pred$preds$risk, data = pred$data)
    scores[[name]] <- dynamicLM::score_left_truncated(
      risk = setNames(data.frame(pred$preds$risk), name),
      data = pred$data,
      time_col = "exit_age", status_col = "event",
      entry_col = fitted$lmdata$entry_col, landmark_col = "LM",
      id_col = "id", cause = 1, w = prediction_window
    )
  }
  score_summary <- do.call(rbind, lapply(names(scores), function(name) {
    transform(scores[[name]]$summary,
              dataset = sub("_(strict|delayed)$", "", name),
              entry_mode = sub("^.*_(strict|delayed)$", "\\1", name))
  }))
  rownames(score_summary) <- NULL

  strict_validation <- rbind(
    transform(compare_strict_to_riskregression(
      predictions$pbcseq_death_strict$risk[
        predictions$pbcseq_death_strict$data$LM == 50
      ],
      predictions$pbcseq_death_strict$data[
        predictions$pbcseq_death_strict$data$LM == 50, , drop = FALSE
      ], prediction_window
    ), dataset = "pbcseq_death")
  )
  stopifnot(max(abs(strict_validation$auc_difference)) < 1e-10,
            max(abs(strict_validation$brier_difference)) < 1e-5)

  toy_weights <- dynamicLM::lt_censoring_weights(
    entry = c(0, 0, 3, 3), exit = c(2, 5, 4, 5),
    status = c(0, 0, 1, 0), horizon = 5
  )
  stopifnot(all.equal(toy_weights$survival_at_horizon, rep(0.5, 4)),
            all.equal(toy_weights$inverse_at_horizon, c(2, 2, 1, 1)))
  toy_score <- dynamicLM::score_left_truncated(
    risk = c(0.1, 0.2, 0.9, 0.3),
    data = data.frame(id = 1:4, LM = 0, analysis_entry = c(0, 0, 3, 3),
                      observed_exit = c(2, 5, 4, 5), event = c(0, 0, 1, 0)),
    time_col = "observed_exit", status_col = "event",
    entry_col = "analysis_entry", id_col = "id", cause = 1, w = 5
  )
  stopifnot(all.equal(toy_score$score$AUC, 1),
            all.equal(toy_score$score$Brier, 0.045),
            grepl("PENDING", toy_score$analytic_influence_function))

  risk_sets <- rbind(
    do.call(rbind, lapply(names(pbc_stacks), function(mode) {
      aggregate(id ~ LM, pbc_stacks[[mode]]$data, length) |>
        transform(dataset = "pbcseq_death", entry_mode = mode)
    })),
    do.call(rbind, lapply(names(nafld_stacks), function(mode) {
      aggregate(id ~ LM, nafld_stacks[[mode]]$data, length) |>
        transform(dataset = "nafld", entry_mode = mode)
    }))
  )
  names(risk_sets)[names(risk_sets) == "id"] <- "rows"
  risk_sets <- risk_sets[, c("dataset", "entry_mode", "LM", "rows")]

  list(
    model_summary = model_summary,
    score_summary = score_summary,
    strict_validation = strict_validation,
    risk_sets = risk_sets,
    toy_weights = toy_weights,
    toy_score = toy_score
  )
})

out <- run$value
drop_trace_values("04_correction", "test_assertions_passed")
utils::write.csv(out$model_summary,
                 project_path("results", "corrected_model_comparison.csv"),
                 row.names = FALSE)
utils::write.csv(out$score_summary,
                 project_path("results", "corrected_score_comparison.csv"),
                 row.names = FALSE)
utils::write.csv(out$strict_validation,
                 project_path("results", "corrected_score_validation.csv"),
                 row.names = FALSE)
utils::write.csv(out$risk_sets,
                 project_path("results", "corrected_risk_sets.csv"),
                 row.names = FALSE)

for (i in seq_len(nrow(out$risk_sets))) {
  row <- out$risk_sets[i, ]
  trace_values(
    report = "04_correction",
    labels = paste(row$dataset, row$entry_mode, "lm", row$LM, "rows", sep = "_"),
    values = row$rows, units = "rows",
    script = "analysis/04_correction.R",
    function_name = "dynamicLM::stack_data", seed = analysis_seed,
    runtime_seconds = run$runtime_seconds
  )
}

for (i in seq_len(nrow(out$model_summary))) {
  row <- out$model_summary[i, ]
  prefix <- paste(row$dataset, row$entry_mode, row$term, sep = "_")
  trace_values(
    report = "04_correction",
    labels = paste0(prefix, c("_observations", "_events", "_subjects",
                              "_estimate", "_robust_se")),
    values = unlist(row[c("observations_used", "events_used", "unique_subjects",
                          "estimate", "robust_se")], use.names = FALSE),
    units = c("rows", "events", "subjects", "log hazard ratio", "standard error"),
    script = "analysis/04_correction.R", function_name = "dynamic_lm/survival::coxph",
    seed = analysis_seed, runtime_seconds = run$runtime_seconds
  )
}

for (i in seq_len(nrow(out$score_summary))) {
  row <- out$score_summary[i, ]
  prefix <- paste(row$dataset, row$entry_mode, sep = "_")
  trace_values(
    report = "04_correction",
    labels = paste0(prefix, c("_summary_auc", "_summary_brier")),
    values = c(row$AUC, row$Brier), units = c("AUC", "Brier score"),
    script = "analysis/04_correction.R",
    function_name = "dynamicLM::score_left_truncated", seed = analysis_seed,
    runtime_seconds = run$runtime_seconds
  )
}

trace_values(
  report = "04_correction",
  labels = c("strict_validation_max_auc_difference",
             "strict_validation_max_brier_difference",
             "toy_censor_survival_horizon_early_entry",
             "toy_censor_survival_horizon_late_entry",
             "toy_inverse_weight_horizon_early_entry",
             "toy_inverse_weight_horizon_late_entry",
             "toy_auc", "toy_brier"),
  values = c(max(abs(out$strict_validation$auc_difference)),
             max(abs(out$strict_validation$brier_difference)),
             out$toy_weights$survival_at_horizon[c(1, 3)],
             out$toy_weights$inverse_at_horizon[c(1, 3)],
             out$toy_score$score$AUC, out$toy_score$score$Brier),
  units = c("AUC", "Brier score", rep("probability", 2), rep("inverse probability", 2),
            "AUC", "Brier score"),
  script = "analysis/04_correction.R",
  function_name = "score validation/devtools::test", seed = analysis_seed,
  runtime_seconds = run$runtime_seconds
)

stopifnot(
  file.exists(project_path("results", "corrected_model_comparison.csv")),
  file.exists(project_path("results", "corrected_score_comparison.csv")),
  file.exists(project_path("results", "corrected_score_validation.csv")),
  file.exists(project_path("results", "corrected_risk_sets.csv"))
)
