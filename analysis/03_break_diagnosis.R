#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(devtools)
  library(here)
  library(survival)
})
source(here::here("analysis", "_helpers.R"))
devtools::load_all(here::here(), quiet = TRUE)

analysis_seed <- 20260908L
assert_single_seed(analysis_seed)
set.seed(analysis_seed)

condition_log <- data.frame(
  dataset = character(), stage = character(), type = character(),
  message = character(), stringsAsFactors = FALSE
)

run_stage <- function(dataset, stage, expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(force(expr), error = function(e) structure(conditionMessage(e), class = "captured_error")),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warnings) > 0L) {
    condition_log <<- rbind(
      condition_log,
      data.frame(dataset = dataset, stage = stage, type = "warning",
                 message = warnings, stringsAsFactors = FALSE)
    )
  }
  if (inherits(value, "captured_error")) {
    condition_log <<- rbind(
      condition_log,
      data.frame(dataset = dataset, stage = stage, type = "error",
                 message = as.character(value), stringsAsFactors = FALSE)
    )
  }
  value
}

stack_diagnostics <- function(lmdata, observation_entry, covariates) {
  data <- lmdata$data
  assert_required_columns(data, c("id", "LM", observation_entry, covariates))
  preentry <- data[[observation_entry]] > data$LM
  complete <- stats::complete.cases(data[, covariates, drop = FALSE])
  out <- data.table(
    landmark = sort(unique(data$LM))
  )
  out[, stacked_n := vapply(landmark, function(a) sum(data$LM == a), integer(1L))]
  out[, preentry_n := vapply(landmark, function(a) sum(data$LM == a & preentry), integer(1L))]
  out[, strict_n := stacked_n - preentry_n]
  out[, delayed_n := vapply(
    landmark,
    function(a) sum(data$LM == a & data[[observation_entry]] < a + lmdata$w),
    integer(1L)
  )]
  out[, analyzable_n := vapply(landmark, function(a) sum(data$LM == a & complete), integer(1L))]
  out[, analyzable_preentry_n := vapply(
    landmark, function(a) sum(data$LM == a & preentry & complete), integer(1L)
  )]
  stopifnot(
    all(out$preentry_n <= out$stacked_n),
    all(out$strict_n <= out$delayed_n),
    all(out$delayed_n <= out$stacked_n),
    all(out$analyzable_n <= out$stacked_n),
    all(out$analyzable_preentry_n <= out$preentry_n)
  )
  out
}

model_diagnostics <- function(object, dataset) {
  if (inherits(object, "captured_error")) return(data.frame())
  models <- if (object$type == "coxph") list(object$model) else object$model$models
  out <- do.call(rbind, lapply(seq_along(models), function(cause) {
    coefficients <- stats::coef(models[[cause]])
    variances <- diag(stats::vcov(models[[cause]]))
    standard_errors <- sqrt(variances)
    data.frame(
      dataset = dataset,
      cause = cause,
      observations_used = models[[cause]]$n,
      events_used = models[[cause]]$nevent,
      coefficient_count = length(coefficients),
      nonfinite_coefficients = sum(!is.finite(coefficients)),
      nonfinite_coefficient_terms = paste(names(coefficients)[!is.finite(coefficients)],
                                          collapse = ";"),
      nonfinite_standard_errors = sum(!is.finite(standard_errors)),
      nonfinite_standard_error_terms = paste(names(standard_errors)[!is.finite(standard_errors)],
                                             collapse = ";"),
      maximum_absolute_finite_coefficient = if (any(is.finite(coefficients))) {
        max(abs(coefficients[is.finite(coefficients)]))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  stopifnot(nrow(out) == length(models), all(out$coefficient_count > 0L))
  out
}

make_reference_rows <- function(data, landmark, window, method) {
  stopifnot(method %in% c("naive", "strict", "delayed"))
  horizon <- landmark + window
  if (method == "naive") {
    keep <- data$exit_age > landmark
    entry <- rep(landmark, sum(keep))
  } else if (method == "strict") {
    keep <- data$entry_age <= landmark & data$exit_age > landmark
    entry <- rep(landmark, sum(keep))
  } else {
    keep <- data$entry_age < horizon & data$exit_age > pmax(landmark, data$entry_age)
    entry <- pmax(landmark, data$entry_age[keep])
  }
  out <- data[keep, , drop = FALSE]
  out$analysis_entry <- entry
  out$analysis_exit <- pmin(out$exit_age, horizon)
  out$analysis_event <- as.integer(out$event == 1L & out$exit_age <= horizon)
  stopifnot(
    nrow(out) == length(entry), nrow(out) > 0L,
    all(out$analysis_entry < out$analysis_exit),
    all(out$analysis_event %in% 0:1),
    sum(out$analysis_event) <= nrow(out)
  )
  out
}

fit_reference_set <- function(dataset, data, landmark, window, formula_rhs) {
  methods <- c("naive", "strict", "delayed")
  fits <- lapply(methods, function(method) {
    rows <- make_reference_rows(data, landmark, window, method)
    formula <- stats::as.formula(
      paste0("survival::Surv(analysis_entry, analysis_exit, analysis_event) ~ ",
             formula_rhs)
    )
    fit <- survival::coxph(formula, data = rows, ties = "breslow", x = TRUE)
    stopifnot(inherits(fit, "coxph"), all(is.finite(stats::coef(fit))))
    se <- sqrt(diag(stats::vcov(fit)))
    stopifnot(identical(names(se), names(stats::coef(fit))), all(is.finite(se)))
    data.frame(
      dataset = dataset,
      landmark = landmark,
      window = window,
      method = method,
      n = nrow(rows),
      events = sum(rows$analysis_event),
      term = names(stats::coef(fit)),
      estimate = unname(stats::coef(fit)),
      se = unname(se),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, fits)
  stopifnot(nrow(out) > 0L)
  references <- out[out$method != "naive", c("dataset", "term", "method", "estimate")]
  names(references)[names(references) == "method"] <- "reference"
  names(references)[names(references) == "estimate"] <- "reference_estimate"
  naive <- out[out$method == "naive", c("dataset", "term", "estimate")]
  names(naive)[names(naive) == "estimate"] <- "naive_estimate"
  differences <- merge(references, naive, by = c("dataset", "term"))
  differences$naive_minus_reference <- differences$naive_estimate - differences$reference_estimate
  stopifnot(nrow(differences) == 2L * nrow(naive),
            all(is.finite(differences$naive_minus_reference)))
  list(fits = out, differences = differences)
}

make_pbc_age <- function() {
  data <- as.data.table(survival::pbcseq)
  setorder(data, id, day)
  baseline <- data[, .SD[1L], by = id]
  baseline <- baseline[, .(
    id, entry_age = age, exit_age = age + futime / 365.25,
    event = as.integer(status == 2L), cr_status = status,
    male = as.integer(sex == "m"),
    stage0 = stage, trt0 = trt, bili0 = bili
  )]
  baseline <- baseline[complete.cases(baseline)]
  baseline[, log_bili_z := as.numeric(scale(log(bili0)))]
  stopifnot(!anyNA(baseline), all(baseline$entry_age < baseline$exit_age))

  data[, bili_locf := nafill(bili, type = "locf"), by = id]
  data[, albumin_locf := nafill(albumin, type = "locf"), by = id]
  data <- merge(
    data[, .(id, day, bili_locf, albumin_locf)],
    baseline[, .(id, entry_age, exit_age, event, cr_status, male, stage0, trt0)],
    by = "id", all = FALSE
  )
  data <- data[complete.cases(data)]
  data[, visit_age := entry_age + day / 365.25]
  data <- data[visit_age < exit_age]
  data[, log_bili := log(bili_locf)]
  stopifnot(!anyDuplicated(data[, .(id, visit_age)]), !anyNA(data),
            all(data$visit_age >= data$entry_age))
  list(long = as.data.frame(data), baseline = as.data.frame(baseline))
}

make_nafld_age <- function() {
  baseline <- as.data.table(survival::nafld1)
  baseline <- baseline[complete.cases(baseline[, .(id, age, futime, status, male, bmi)])]
  baseline <- baseline[, .(
    id, entry_age = age, exit_age = age + futime / 365.25,
    event = as.integer(status == 1L), male, bmi
  )]
  baseline[, bmi_z := as.numeric(scale(bmi))]
  stopifnot(!anyDuplicated(baseline$id), !anyNA(baseline),
            all(baseline$entry_age < baseline$exit_age))

  labs <- as.data.table(survival::nafld2)
  labs <- labs[test %in% c("hdl", "chol")]
  labs <- labs[, .(value = mean(value)), by = .(id, days, test)]
  labs <- dcast(labs, id + days ~ test, value.var = "value")
  labs <- labs[complete.cases(labs)]
  data <- merge(labs, baseline, by = "id", all = FALSE)
  data[, visit_age := entry_age + days / 365.25]
  data <- data[visit_age < exit_age]
  data[, hdl_z := as.numeric(scale(hdl))]
  stopifnot(!anyDuplicated(data[, .(id, visit_age)]), !anyNA(data))
  list(long = as.data.frame(data), baseline = as.data.frame(baseline))
}

run <- timed({
  pbc <- make_pbc_age()
  nafld <- make_nafld_age()
  pbc_landmarks <- c(40, 50, 60, 70)
  nafld_landmarks <- c(40, 50, 60, 70, 80)
  prediction_window <- 10

  pbc_stack <- dynamicLM::stack_data(
    pbc$long,
    outcome = list(time = "exit_age", status = "cr_status"),
    lms = pbc_landmarks, w = prediction_window,
    covs = list(fixed = c("entry_age", "male", "stage0", "trt0"),
                varying = c("log_bili", "albumin_locf")),
    format = "long", id = "id", rtime = "visit_age"
  )
  pbc_diag <- stack_diagnostics(
    pbc_stack, "entry_age", c("male", "stage0", "trt0", "log_bili", "albumin_locf")
  )
  pbc_diag$dataset <- "pbcseq"
  pbc_stack <- dynamicLM::add_interactions(
    pbc_stack, lm_covs = c("log_bili", "albumin_locf"),
    func_covars = "linear", func_lms = "linear"
  )
  pbc_model <- run_stage(
    "pbcseq", "dynamic_lm",
    dynamicLM::dynamic_lm(
      pbc_stack,
      stats::as.formula("Hist(exit_age, cr_status, LM) ~ male + stage0 + log_bili + log_bili_LM1 + albumin_locf + albumin_locf_LM1 + LM1 + cluster(id)"),
      "CSC", x = TRUE
    )
  )
  pbc_score <- NULL
  pbc_prediction <- NULL
  if (!inherits(pbc_model, "captured_error")) {
    pbc_prediction <- run_stage("pbcseq", "predict", stats::predict(pbc_model, cause = 2))
    if (!inherits(pbc_prediction, "captured_error")) {
      pbc_score <- run_stage(
        "pbcseq", "score_summary",
        dynamicLM::score(list(naive_age = pbc_prediction),
                         times = pbc_landmarks, cause = 2, seed = analysis_seed,
                         metrics = c("auc", "brier"), summary = TRUE)
      )
    }
  }
  pbc_model_diagnostics <- model_diagnostics(pbc_model, "pbcseq")

  nafld_stack <- dynamicLM::stack_data(
    nafld$long,
    outcome = list(time = "exit_age", status = "event"),
    lms = nafld_landmarks, w = prediction_window,
    covs = list(fixed = c("entry_age", "male", "bmi_z"), varying = "hdl_z"),
    format = "long", id = "id", rtime = "visit_age"
  )
  nafld_diag <- stack_diagnostics(
    nafld_stack, "entry_age", c("male", "bmi_z", "hdl_z")
  )
  nafld_diag$dataset <- "nafld"
  nafld_stack <- dynamicLM::add_interactions(
    nafld_stack, lm_covs = "hdl_z", func_covars = "linear",
    func_lms = "linear"
  )
  nafld_model <- run_stage(
    "nafld", "dynamic_lm",
    dynamicLM::dynamic_lm(
      nafld_stack,
      stats::as.formula("Surv(LM, exit_age, event) ~ male + bmi_z + hdl_z + hdl_z_LM1 + LM1 + cluster(id)"),
      "coxph", x = TRUE
    )
  )
  nafld_model_diagnostics <- model_diagnostics(nafld_model, "nafld")

  # Non-age-specific call-evaluation bug discovered while reproducing the tutorial.
  formula_symbol_error <- local({
    formula_object <- stats::as.formula(
      "Surv(LM, exit_age, event) ~ male + bmi_z + hdl_z + cluster(id)"
    )
    run_stage(
      "nafld", "dynamic_lm_formula_symbol",
      dynamicLM::dynamic_lm(nafld_stack, formula_object, "coxph")
    )
  })
  stopifnot(inherits(formula_symbol_error, "captured_error"))

  nafld_reference <- fit_reference_set(
    "nafld", nafld$baseline, landmark = 50, window = prediction_window,
    formula_rhs = "male + bmi_z"
  )
  pbc_reference <- fit_reference_set(
    "pbcseq_death", pbc$baseline, landmark = 50, window = prediction_window,
    formula_rhs = "male + log_bili_z"
  )

  diagnostics <- rbindlist(list(pbc_diag, nafld_diag), fill = TRUE)
  setcolorder(diagnostics, c("dataset", "landmark", "stacked_n", "preentry_n",
                            "strict_n", "delayed_n", "analyzable_n",
                            "analyzable_preentry_n"))
  stopifnot(
    any(diagnostics$preentry_n > 0L),
    any(diagnostics$analyzable_preentry_n > 0L),
    all(diagnostics$stacked_n == diagnostics$preentry_n +
          (diagnostics$stacked_n - diagnostics$preentry_n))
  )

  list(
    diagnostics = diagnostics,
    model_diagnostics = rbind(pbc_model_diagnostics, nafld_model_diagnostics),
    coefficient_fits = rbind(nafld_reference$fits, pbc_reference$fits),
    coefficient_differences = rbind(nafld_reference$differences,
                                    pbc_reference$differences),
    condition_log = condition_log,
    pbc_score_failed = inherits(pbc_score, "captured_error"),
    pbc_score_attempted = !is.null(pbc_score),
    pbc_prediction_failed = inherits(pbc_prediction, "captured_error"),
    pbc_model_failed = inherits(pbc_model, "captured_error"),
    nafld_model_failed = inherits(nafld_model, "captured_error")
  )
})

out <- run$value
utils::write.csv(out$diagnostics, project_path("results", "age_naive_stack_diagnostics.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(out$coefficient_fits,
                 project_path("results", "age_coefficient_reference_fits.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(out$coefficient_differences,
                 project_path("results", "age_coefficient_differences.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(out$condition_log,
                 project_path("results", "age_break_conditions.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(out$model_diagnostics,
                 project_path("results", "age_model_diagnostics.csv"),
                 row.names = FALSE, na = "")

for (i in seq_len(nrow(out$diagnostics))) {
  row <- as.data.frame(out$diagnostics[i, ])
  trace_values(
    report = "03_break_diagnosis",
    labels = paste0(row$dataset, "_lm", row$landmark, "_",
                    c("stacked_n", "preentry_n", "strict_n", "delayed_n",
                      "analyzable_n",
                      "analyzable_preentry_n")),
    values = as.numeric(unlist(
      row[1L, c("stacked_n", "preentry_n", "strict_n", "delayed_n", "analyzable_n",
                "analyzable_preentry_n")], use.names = FALSE
    )),
    units = "rows",
    script = "analysis/03_break_diagnosis.R",
    function_name = "dynamicLM::stack_data/stack_diagnostics",
    seed = analysis_seed,
    runtime_seconds = run$runtime_seconds
  )
}

for (i in seq_len(nrow(out$model_diagnostics))) {
  row <- out$model_diagnostics[i, ]
  fields <- c("observations_used", "events_used", "coefficient_count",
              "nonfinite_coefficients", "nonfinite_standard_errors",
              "maximum_absolute_finite_coefficient")
  trace_values(
    report = "03_break_diagnosis",
    labels = paste0(row$dataset, "_cause", row$cause, "_", fields),
    values = as.numeric(unlist(row[, fields], use.names = FALSE)),
    units = c("rows", "events", "coefficients", "coefficients",
              "standard errors", "log hazard ratio"),
    script = "analysis/03_break_diagnosis.R",
    function_name = "dynamicLM::dynamic_lm/model_diagnostics",
    seed = analysis_seed,
    runtime_seconds = run$runtime_seconds
  )
}

for (i in seq_len(nrow(out$coefficient_fits))) {
  row <- out$coefficient_fits[i, ]
  prefix <- paste(row$dataset, row$method, gsub("[^[:alnum:]]+", "_", row$term), sep = "_")
  trace_values(
    report = "03_break_diagnosis",
    labels = paste0(prefix, c("_n", "_events", "_estimate", "_se")),
    values = as.numeric(unlist(row[, c("n", "events", "estimate", "se")],
                               use.names = FALSE)),
    units = c("subjects", "events", "log hazard ratio", "standard error"),
    script = "analysis/03_break_diagnosis.R",
    function_name = "survival::coxph",
    seed = analysis_seed,
    runtime_seconds = run$runtime_seconds
  )
}

for (i in seq_len(nrow(out$coefficient_differences))) {
  row <- out$coefficient_differences[i, ]
  prefix <- paste(row$dataset, "naive_minus", row$reference,
                  gsub("[^[:alnum:]]+", "_", row$term), sep = "_")
  trace_values(
    report = "03_break_diagnosis",
    labels = paste0(prefix, c("_naive", "_reference", "_difference")),
    values = as.numeric(unlist(
      row[, c("naive_estimate", "reference_estimate", "naive_minus_reference")],
      use.names = FALSE
    )),
    units = "log hazard ratio",
    script = "analysis/03_break_diagnosis.R",
    function_name = "fit_reference_set",
    seed = analysis_seed,
    runtime_seconds = run$runtime_seconds
  )
}

trace_values(
  report = "03_break_diagnosis",
  labels = c("prediction_window_years",
             paste0("pbc_landmark_", seq_along(c(40, 50, 60, 70))),
             paste0("nafld_landmark_", seq_along(c(40, 50, 60, 70, 80)))),
  values = c(10, c(40, 50, 60, 70), c(40, 50, 60, 70, 80)),
  units = "years of age",
  script = "analysis/03_break_diagnosis.R",
  function_name = "run",
  seed = analysis_seed,
  runtime_seconds = run$runtime_seconds
)

trace_values(
  report = "03_break_diagnosis",
  labels = c("condition_count", "error_count", "warning_count",
             "pbc_prediction_failed", "pbc_score_attempted",
             "pbc_score_failed", "pbc_model_failed", "nafld_model_failed"),
  values = c(nrow(out$condition_log), sum(out$condition_log$type == "error"),
             sum(out$condition_log$type == "warning"),
             as.integer(out$pbc_prediction_failed),
             as.integer(out$pbc_score_attempted),
             as.integer(out$pbc_score_failed), as.integer(out$pbc_model_failed),
             as.integer(out$nafld_model_failed)),
  units = c("conditions", "errors", "warnings", rep("indicator", 5L)),
  script = "analysis/03_break_diagnosis.R",
  function_name = "run_stage",
  seed = analysis_seed,
  runtime_seconds = run$runtime_seconds
)

stopifnot(
  file.exists(project_path("results", "age_naive_stack_diagnostics.csv")),
  file.exists(project_path("results", "age_coefficient_reference_fits.csv")),
  file.exists(project_path("results", "age_coefficient_differences.csv")),
  file.exists(project_path("results", "age_break_conditions.csv")),
  file.exists(project_path("results", "age_model_diagnostics.csv"))
)
