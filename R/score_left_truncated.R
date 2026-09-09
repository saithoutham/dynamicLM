#' Left-truncation-aware censoring weights
#'
#' Estimate the marginal censoring survival curve on the analysis time scale
#' using counting-process risk sets, then return inverse probabilities
#' conditional on each subject's analysis entry. Administrative censoring at
#' `horizon` is not counted as a censoring event before the prediction horizon.
#'
#' @param entry Numeric vector of counting-process entry times.
#' @param exit Numeric vector of observed or administratively censored exit
#'   times.
#' @param status Integer event code, with zero denoting censoring.
#' @param horizon Scalar prediction horizon.
#'
#' @return A list containing the fitted censoring curve and inverse censoring
#'   weights at each subject's observed exit (left limit) and at `horizon`.
#' @export
lt_censoring_weights <- function(entry, exit, status, horizon) {
  n <- length(entry)
  if (length(exit) != n || length(status) != n || n == 0L)
    stop("entry, exit, and status must have the same positive length.")
  if (length(horizon) != 1L || !is.finite(horizon))
    stop("horizon must be one finite number.")
  if (!is.numeric(entry) || !is.numeric(exit) ||
      any(!is.finite(entry)) || any(!is.finite(exit)))
    stop("entry and exit must contain finite numeric values.")
  if (any(entry >= exit))
    stop("Every counting-process interval must satisfy entry < exit.")
  if (any(exit > horizon))
    stop("exit cannot exceed the prediction horizon.")
  if (anyNA(status) || any(status < 0) || any(status != as.integer(status)))
    stop("status must contain non-negative integer event codes.")

  censor_event <- as.integer(status == 0L & exit < horizon)
  fit <- survival::survfit(
    survival::Surv(entry, exit, censor_event) ~ 1,
    type = "kaplan-meier"
  )

  step_survival <- function(at, before = FALSE) {
    vapply(at, function(value) {
      eligible <- if (before) fit$time < value else fit$time <= value
      if (any(eligible)) utils::tail(fit$surv[eligible], 1L) else 1
    }, numeric(1L))
  }

  survival_at_entry <- step_survival(entry, before = FALSE)
  survival_before_exit <- step_survival(exit, before = TRUE)
  survival_at_horizon <- step_survival(rep(horizon, n), before = FALSE)
  if (any(!is.finite(c(survival_at_entry, survival_before_exit,
                       survival_at_horizon))) ||
      any(c(survival_at_entry, survival_before_exit,
            survival_at_horizon) <= 0))
    stop("Estimated censoring survival is zero or non-finite; IPCW is undefined.")

  inverse_at_exit <- survival_at_entry / survival_before_exit
  inverse_at_horizon <- survival_at_entry / survival_at_horizon
  if (any(!is.finite(c(inverse_at_exit, inverse_at_horizon))) ||
      any(c(inverse_at_exit, inverse_at_horizon) < 1))
    stop("Internal error: invalid conditional inverse censoring weights.")

  list(
    fit = fit,
    censor_event = censor_event,
    survival_at_entry = survival_at_entry,
    survival_before_exit = survival_before_exit,
    survival_at_horizon = survival_at_horizon,
    inverse_at_exit = inverse_at_exit,
    inverse_at_horizon = inverse_at_horizon
  )
}


weighted_dynamic_auc <- function(risk, case_weight, control_weight) {
  if (length(risk) != length(case_weight) ||
      length(risk) != length(control_weight))
    stop("risk and AUC weights must have equal lengths.")
  if (any(!is.finite(c(risk, case_weight, control_weight))) ||
      any(case_weight < 0) || any(control_weight < 0))
    stop("risk and AUC weights must be finite and weights non-negative.")
  total_cases <- sum(case_weight)
  total_controls <- sum(control_weight)
  if (total_cases <= 0 || total_controls <= 0)
    stop("AUC requires positive weighted case and control totals.")

  grouped <- stats::aggregate(
    cbind(case_weight, control_weight),
    by = list(risk = risk), FUN = sum
  )
  grouped <- grouped[order(grouped$risk), , drop = FALSE]
  control_below <- c(0, utils::head(cumsum(grouped$control_weight), -1L))
  numerator <- sum(grouped$case_weight *
                     (control_below + 0.5 * grouped$control_weight))
  auc <- numerator / (total_cases * total_controls)
  tolerance <- sqrt(.Machine$double.eps)
  if (!is.finite(auc) || auc < -tolerance || auc > 1 + tolerance)
    stop("Internal error: AUC outside [0, 1].")
  pmin(1, pmax(0, auc))
}


score_lt_once <- function(risk, data, time_col, status_col, entry_col,
                          landmark_col, cause, w) {
  landmarks <- sort(unique(data[[landmark_col]]))
  models <- colnames(risk)
  if (is.null(models) || any(!nzchar(models)))
    stop("risk must have non-empty model column names.")

  rows <- lapply(landmarks, function(landmark) {
    index <- data[[landmark_col]] == landmark
    di <- data[index, , drop = FALSE]
    ri <- risk[index, , drop = FALSE]
    horizon <- landmark + w
    if (nrow(di) == 0L || nrow(di) != nrow(ri))
      stop("Risk and outcome rows do not align within landmark.")
    if (any(di[[entry_col]] < landmark))
      stop("Analysis entry cannot precede its landmark.")
    if (any(di[[entry_col]] >= di[[time_col]]) ||
        any(di[[time_col]] > horizon))
      stop("Invalid landmark counting-process interval.")

    cw <- lt_censoring_weights(
      entry = di[[entry_col]], exit = di[[time_col]],
      status = di[[status_col]], horizon = horizon
    )
    observed_event <- di[[status_col]] > 0L & di[[time_col]] <= horizon
    cases <- observed_event & di[[status_col]] == cause
    competing_controls <- observed_event & di[[status_col]] != cause
    horizon_controls <- di[[time_col]] >= horizon & di[[status_col]] == 0L
    case_weight <- as.numeric(cases) * cw$inverse_at_exit
    control_weight <- as.numeric(competing_controls) * cw$inverse_at_exit +
      as.numeric(horizon_controls) * cw$inverse_at_horizon

    do.call(rbind, lapply(seq_along(models), function(j) {
      prediction <- ri[, j]
      if (any(!is.finite(prediction)) || any(prediction < 0) ||
          any(prediction > 1))
        stop("Predicted risks must be finite and lie in [0, 1].")
      outcome <- as.numeric(cases)
      brier_contribution <- numeric(nrow(di))
      brier_contribution[observed_event] <-
        (outcome[observed_event] - prediction[observed_event])^2 *
        cw$inverse_at_exit[observed_event]
      brier_contribution[horizon_controls] <-
        prediction[horizon_controls]^2 *
        cw$inverse_at_horizon[horizon_controls]
      brier <- sum(brier_contribution) / nrow(di)
      auc <- weighted_dynamic_auc(prediction, case_weight, control_weight)
      if (!is.finite(brier) || brier < 0)
        stop("Internal error: invalid Brier score.")
      data.frame(
        landmark = landmark,
        model = models[j],
        n = nrow(di),
        cases = sum(cases),
        controls = sum(competing_controls | horizon_controls),
        censorings_before_horizon = sum(cw$censor_event),
        AUC = auc,
        Brier = brier,
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  stopifnot(nrow(out) == length(landmarks) * length(models))
  out
}


#' Score predictions with left truncation
#'
#' Compute cumulative/dynamic AUC and cause-specific Brier scores using a
#' counting-process Kaplan--Meier estimate of censoring. For each subject, the
#' censoring probability is conditional on that subject's analysis entry.
#'
#' Analytic influence-function standard errors are deliberately not provided:
#' their extension to staggered entry is not established here. If `B >= 2`,
#' subject-level nonparametric bootstrap standard errors and normal intervals
#' are returned, with the same resample shared across landmarks.
#'
#' @param risk Numeric vector, matrix, or data frame of predicted risks. Matrix
#'   columns are treated as models.
#' @param data Stacked landmark data, one row per subject and landmark.
#' @param time_col,status_col,entry_col,landmark_col Column names for exit,
#'   event status, counting-process entry, and landmark.
#' @param id_col Subject identifier column.
#' @param cause Cause of interest; censoring must be coded zero.
#' @param w Prediction-window width.
#' @param B Number of subject-level bootstrap replicates. Values below two
#'   request point estimates only.
#' @param seed Explicit bootstrap seed. Required when `B >= 2`.
#' @param conf.int Confidence level for bootstrap normal intervals.
#'
#' @return An object of class `LMScoreLT` with landmark-specific and summary
#'   scores. Inference is marked as `"none"` or `"subject bootstrap"`.
#' @export
score_left_truncated <- function(risk, data, time_col, status_col, entry_col,
                                 landmark_col = "LM", id_col, cause = 1L, w,
                                 B = 0L, seed = NULL, conf.int = 0.95) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  required <- c(time_col, status_col, entry_col, landmark_col, id_col)
  if (any(!required %in% names(data)))
    stop(paste("Missing required columns:",
               paste(setdiff(required, names(data)), collapse = ", ")))
  if (anyNA(data[, required, drop = FALSE]))
    stop("Scoring columns cannot contain missing values.")
  if (length(cause) != 1L || !is.finite(cause) || cause <= 0)
    stop("cause must be one positive event code.")
  if (length(w) != 1L || !is.finite(w) || w <= 0)
    stop("w must be one positive finite number.")
  if (length(B) != 1L || !is.finite(B) || B < 0 || B != as.integer(B))
    stop("B must be a non-negative integer.")
  if (anyDuplicated(data[, c(id_col, landmark_col), drop = FALSE]))
    stop("There must be at most one row per subject and landmark.")

  if (is.vector(risk) && !is.list(risk)) {
    risk <- matrix(risk, ncol = 1L,
                   dimnames = list(NULL, "model"))
  } else {
    risk <- as.matrix(risk)
    if (is.null(colnames(risk)))
      colnames(risk) <- paste0("model", seq_len(ncol(risk)))
  }
  storage.mode(risk) <- "double"
  if (nrow(risk) != nrow(data) || ncol(risk) == 0L)
    stop("risk must have one row per data row and at least one model column.")

  point <- score_lt_once(risk, data, time_col, status_col, entry_col,
                         landmark_col, cause, w)
  summary_point <- stats::aggregate(
    point[c("AUC", "Brier")], by = list(model = point$model), FUN = mean
  )
  inference <- "none"
  bootstrap_success <- 0L
  bootstrap_scores <- NULL

  if (B >= 2L) {
    if (is.null(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != as.integer(seed))
      stop("An explicit integer seed is required when B >= 2.")
    set.seed(as.integer(seed))
    ids <- unique(data[[id_col]])
    bootstrap_scores <- lapply(seq_len(B), function(b) {
      sampled <- sample(ids, length(ids), replace = TRUE)
      pieces <- lapply(seq_along(sampled), function(k) {
        index <- which(data[[id_col]] == sampled[k])
        if (length(index) == 0L) stop("Internal bootstrap ID mismatch.")
        list(data = data[index, , drop = FALSE], risk = risk[index, , drop = FALSE])
      })
      db <- do.call(rbind, lapply(pieces, `[[`, "data"))
      rb <- do.call(rbind, lapply(pieces, `[[`, "risk"))
      piece_rows <- vapply(pieces, function(piece) nrow(piece$data), integer(1L))
      db[[id_col]] <- rep(seq_along(pieces), piece_rows)
      tryCatch({
        sb <- score_lt_once(rb, db, time_col, status_col, entry_col,
                            landmark_col, cause, w)
        stats::aggregate(
          sb[c("AUC", "Brier")], by = list(model = sb$model), FUN = mean
        )
      }, error = function(e) NULL)
    })
    ok <- !vapply(bootstrap_scores, is.null, logical(1L))
    bootstrap_success <- sum(ok)
    if (bootstrap_success < B)
      warning(B - bootstrap_success,
              " bootstrap replicates failed because a landmark lacked cases or controls.")
    if (bootstrap_success < 2L)
      stop("Fewer than two valid bootstrap replicates; inference is undefined.")
    bootstrap_scores <- do.call(rbind, lapply(which(ok), function(b) {
      transform(bootstrap_scores[[b]], bootstrap = b)
    }))
    alpha <- 1 - conf.int
    z <- stats::qnorm(1 - alpha / 2)
    bootstrap_summary <- do.call(rbind, lapply(unique(summary_point$model), function(m) {
      bm <- bootstrap_scores[bootstrap_scores$model == m, , drop = FALSE]
      p <- summary_point[summary_point$model == m, , drop = FALSE]
      se_auc <- stats::sd(bm$AUC)
      se_brier <- stats::sd(bm$Brier)
      transform(
        p,
        AUC_se = se_auc,
        AUC_lower = pmax(0, AUC - z * se_auc),
        AUC_upper = pmin(1, AUC + z * se_auc),
        Brier_se = se_brier,
        Brier_lower = pmax(0, Brier - z * se_brier),
        Brier_upper = pmin(1, Brier + z * se_brier)
      )
    }))
    summary_point <- bootstrap_summary
    inference <- "subject bootstrap"
  }

  out <- list(
    score = point,
    summary = summary_point,
    inference = inference,
    analytic_influence_function = "PENDING: not derived for staggered entry",
    B_requested = as.integer(B),
    B_success = bootstrap_success,
    seed = if (B >= 2L) as.integer(seed) else NA_integer_,
    bootstrap = bootstrap_scores
  )
  class(out) <- "LMScoreLT"
  out
}
