run_pbc_baseline <- function(seed = 20260908L) {
  assert_single_seed(seed)
  set.seed(seed)
  conditions <- data.frame(stage = character(), message = character(),
                           call = character(), stringsAsFactors = FALSE)
  run_stage <- function(stage, expr) {
    withCallingHandlers(
      force(expr),
      warning = function(w) {
        warning_call <- conditionCall(w)
        conditions <<- rbind(
          conditions,
          data.frame(
            stage = stage,
            message = conditionMessage(w),
            call = if (is.null(warning_call)) "" else paste(deparse(warning_call), collapse = " "),
            stringsAsFactors = FALSE
          )
        )
        invokeRestart("muffleWarning")
      }
    )
  }

  pbc_df <- dynamicLM::get_pbc_long()
  assert_required_columns(
    pbc_df,
    c("id", "time", "status", "tstart", "male", "stage", "trt", "age",
      "albumin", "alk.phos", "ascites", "ast", "bili", "edema", "hepato",
      "platelet", "protime", "spiders"),
    "pbc_df"
  )
  stopifnot(
    nrow(pbc_df) > 0L,
    length(unique(pbc_df$id)) > 0L,
    !anyNA(pbc_df),
    all(pbc_df$tstart >= 0),
    all(pbc_df$time >= pbc_df$tstart),
    all(pbc_df$status %in% 0:2)
  )

  outcome <- list(time = "time", status = "status")
  fixed_variables <- c("male", "stage", "trt", "age")
  varying_variables <- c(
    "albumin", "alk.phos", "ascites", "ast", "bili", "edema", "hepato",
    "platelet", "protime", "spiders"
  )
  covars <- list(fixed = fixed_variables, varying = varying_variables)
  window <- 5
  landmarks <- seq(0, 4, by = 1)

  lmdata <- dynamicLM::stack_data(
    pbc_df, outcome, landmarks, window, covars, format = "long",
    id = "id", rtime = "tstart"
  )
  stopifnot(inherits(lmdata, "LMdataframe"), nrow(lmdata$data) > 0L)
  risk_counts <- as.integer(table(factor(lmdata$data$LM, levels = landmarks)))
  stopifnot(length(risk_counts) == length(landmarks), all(risk_counts > 0L),
            all(diff(risk_counts) <= 0L), sum(risk_counts) == nrow(lmdata$data))

  lmdata$data$age <- lmdata$data$age + lmdata$data$LM
  lmdata <- dynamicLM::add_interactions(
    lmdata,
    func_covars = c("linear", "quadratic"),
    func_lms = c("linear", "quadratic")
  )
  stopifnot(
    all(c("bili_LM1", "bili_LM2", "albumin_LM1", "albumin_LM2",
          "LM1", "LM2") %in% names(lmdata$data)),
    !anyNA(lmdata$data)
  )

  supermodel <- run_stage("dynamic_lm_unpenalized", dynamicLM::dynamic_lm(
    lmdata,
    stats::as.formula("Hist(time, status, LM) ~ stage + bili + bili_LM1 + bili_LM2 + albumin + albumin_LM1 + albumin_LM2 + LM1 + LM2 + cluster(id)"),
    "CSC", x = TRUE
  ))
  stopifnot(inherits(supermodel, "LMCSC"), length(supermodel$model$models) == 2L)

  path <- run_stage("pen_lm", dynamicLM::pen_lm(lmdata, alpha = 0))
  stopifnot(inherits(path, "pen_lm"), length(path) == 2L)

  cv_model <- run_stage(
    "cv.pen_lm",
    dynamicLM::cv.pen_lm(lmdata, alpha = 0, seed = seed)
  )
  stopifnot(inherits(cv_model, "cv.pen_lm"), length(cv_model) == 2L)
  supermodel_pen <- run_stage(
    "dynamic_lm_penalized",
    dynamicLM::dynamic_lm(cv_model, lambda = "lambda.1se", x = TRUE)
  )
  stopifnot(inherits(supermodel_pen, "penLMCSC"),
            length(supermodel_pen$model$models) == 2L)

  predictions <- run_stage("predict", list(
    LM = stats::predict(supermodel),
    penLM = stats::predict(supermodel_pen)
  ))
  for (prediction in predictions) {
    stopifnot(inherits(prediction, "LMpred"), nrow(prediction$preds) > 0L)
    risk <- prediction$preds$risk
    stopifnot(!anyNA(risk), all(is.finite(risk)), all(risk >= 0), all(risk <= 1))
  }
  stopifnot(
    identical(predictions$LM$data$id, predictions$penLM$data$id),
    identical(predictions$LM$preds$LM, predictions$penLM$preds$LM)
  )

  scores <- run_stage("score", dynamicLM::score(
    predictions, times = c(0, 2, 4), seed = seed,
    metrics = c("auc", "brier"), summary = TRUE, se.fit = TRUE
  ))
  stopifnot(inherits(scores, "LMScore"))
  score_tables <- list(
    AUC = scores$AUC$score,
    Brier = scores$Brier$score,
    AUC_summary = scores$AUC_summary$score,
    Brier_summary = scores$Brier_summary$score
  )
  stopifnot(all(vapply(score_tables, nrow, integer(1L)) > 0L))
  for (table in score_tables) {
    numeric_columns <- vapply(table, is.numeric, logical(1L))
    table_df <- as.data.frame(table)
    stopifnot(all(is.finite(as.matrix(table_df[, numeric_columns, drop = FALSE]))))
  }

  event_counts <- as.data.frame(
    table(LM = factor(lmdata$data$LM, levels = landmarks),
          status = factor(lmdata$data$status, levels = 0:2)),
    stringsAsFactors = FALSE
  )
  stopifnot(sum(event_counts$Freq) == nrow(lmdata$data))

  list(
    seed = seed,
    conditions = conditions,
    pbc_rows = nrow(pbc_df),
    pbc_ids = length(unique(pbc_df$id)),
    terminal_visit_ties = sum(pbc_df$time == pbc_df$tstart),
    lm_rows = nrow(lmdata$data),
    landmarks = landmarks,
    window = window,
    risk_counts = risk_counts,
    event_counts = event_counts,
    unpenalized_coefficients = unlist(lapply(supermodel$model$models, stats::coef)),
    path_lambda_counts = vapply(path, function(model) length(model$lambda), integer(1L)),
    cv_lambda_min = vapply(cv_model, function(model) model$lambda.min, numeric(1L)),
    cv_lambda_1se = vapply(cv_model, function(model) model$lambda.1se, numeric(1L)),
    prediction_ranges = rbind(
      LM = range(predictions$LM$preds$risk),
      penLM = range(predictions$penLM$preds$risk)
    ),
    score_tables = score_tables
  )
}
