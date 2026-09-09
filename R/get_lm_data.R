#' Build a landmark dataset
#'
#' @param data Data frame from which to construct landmark super dataset
#' @param outcome A list with items time and status, containing character
#'   strings identifying the names of time and status variables, respectively,
#'   of the survival outcome
#' @param lm The value of the landmark time point at which to construct the
#'   landmark dataset.
#' @param horizon Scalar, the value of the prediction window (ie predict risk
#'   within time w landmark points)
#' @param covs A list with items fixed and varying, containing character strings
#'   specifying column names in the data containing time-fixed and time-varying
#'   covariates, respectively.
#' @param format Character string specifying whether the original data are in
#'   wide (default) or in long format.
#' @param id Character string specifying the column name in data containing the
#'   subject id.
#' @param rtime Character string specifying the column name in data containing
#'   the (running) time variable associated with the time-varying variables;
#'   only needed if format = "long".
#' @param left.open Boolean (default = FALSE), indicating if the intervals for the
#'   time-varying covariates are open on the left (and closed on the right) or
#'   vice-versa.
#' @param split.data List of data split according to ID. Allows for faster
#'   computation.
#' @param entry Optional character string naming subject-specific observation
#'   entry time. Required when `entry_mode` is not `"shared"`.
#' @param entry_mode Observation-entry rule. `"shared"` preserves the original
#'   behavior, `"strict"` requires entry no later than the landmark, and
#'   `"delayed"` permits entry during the prediction window using a
#'   subject-specific counting-process start.
#'
#' @details This function is based on the `cutLM()` implementation from the
#'   archived `dynpred` package, with changes.
#'   The original function was authored by Hein Putter.
#' @references
#'   - van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#'     Clinical Survival Analysis. Chapman & Hall.
#'   - The dynpred package
#'     (https://cran.r-project.org/web/packages/dynpred/index.html),
#'     in particular, the code for cutLM.
#' @return A landmark dataset.
#'
#' @examples
#' \dontrun{
#' data(relapse)
#' outcome <- list(time = "Time", status = "event")
#' covars <- list(fixed = c("age.at.time.0", "male", "stage", "bmi"),
#'                varying = c("treatment"))
#' lm12 <- get_lm_data(relapse, outcome, lm = 12, horizon = 60, covs = covars,
#'                     format = "long", id = "ID", rtime = "T_txgiven")
#' head(lm12)
#' }
#'
#' @seealso [dynamicLM::stack_data()]
#' @export
get_lm_data <- function(data, outcome, lm, horizon, covs,
                        format = c("wide", "long"), id, rtime,
                        left.open = FALSE, split.data, entry = NULL,
                        entry_mode = c("shared", "strict", "delayed")) {
  format <- match.arg(format)
  entry_mode <- match.arg(entry_mode)
  if (entry_mode != "shared") {
    if (is.null(entry) || length(entry) != 1L || !is.character(entry))
      stop("argument 'entry' must name one column when entry_mode is not 'shared'")
    if (!(entry %in% names(data)))
      stop(paste("Entry column", entry, "is not in the data."))
    if (!is.numeric(data[[entry]]) || any(!is.finite(data[[entry]])))
      stop("Observation entry times must be finite numeric values.")
  }
  if (format == "wide") {
    lmdata <- data
    assessment_time <- if (entry_mode == "delayed") {
      pmax(lm, lmdata[[entry]])
    } else {
      rep(lm, nrow(lmdata))
    }
    if (!is.null(covs$varying)) {
      for (col in covs$varying)
        lmdata[[col]] <- 1 - as.numeric(lmdata[[col]] > assessment_time)
    }

  } else {
    if (missing(id))
      stop("argument 'id' should be specified for long format data")
    if (missing(rtime))
      stop("argument 'rtime' should be specified for long format data")
    lookup <- FALSE
    if (missing(split.data)) {
      data <- data[order(data[[id]], data[[rtime]]), ]
      ids <- unique(data[[id]])
      n <- length(ids)
      lookup <- TRUE
    } else {
      n <- length(split.data)
    }

    lmdata <- lapply(1:n, function(i) {
      if (lookup) {
        wh <- which(data[[id]] == ids[i])
        di <- data[wh, ]
      } else {
        di <- split.data[[i]]
      }
      t.fups <- di[[rtime]]

      assessment_time <- lm
      if (entry_mode != "shared") {
        subject_entries <- unique(di[[entry]])
        if (length(subject_entries) != 1L || !is.finite(subject_entries))
          stop("Observation entry must be one finite value per subject.")
        if (entry_mode == "delayed")
          assessment_time <- max(lm, subject_entries)
      }

      idx <- findInterval(assessment_time, c(t.fups, Inf), left.open = left.open)

      if (idx != 0) {
        return(di[idx, ])
      } else {
        out <- di[1, ]
        if (!is.null(covs$varying)) {
          out[, covs$varying] <- NA
          out[, rtime] <- NA
        }
        return(out)
      }
    })
    lmdata <- do.call(rbind, lmdata)
  }

  if (entry_mode == "shared") {
    keep <- lmdata[[outcome$time]] > lm
    analysis_entry <- rep(lm, nrow(lmdata))
  } else if (entry_mode == "strict") {
    keep <- lmdata[[entry]] <= lm & lmdata[[outcome$time]] > lm
    analysis_entry <- rep(lm, nrow(lmdata))
  } else {
    analysis_entry <- pmax(lm, lmdata[[entry]])
    keep <- lmdata[[entry]] < horizon & lmdata[[outcome$time]] > analysis_entry
  }
  lmdata <- lmdata[keep, ]
  analysis_entry <- analysis_entry[keep]
  if (nrow(lmdata) == 0) {
    warning("Landmark dataset for lm = ", lm, " could not be constructed as no individuals are alive after this point.")
    return(NULL)
  }

  lmdata[outcome$status] <- lmdata[[outcome$status]] *
    as.numeric(lmdata[[outcome$time]] <= horizon)
  lmdata[outcome$time] <- pmin(as.vector(lmdata[[outcome$time]]), horizon)
  lmdata$LM <- lm
  if (entry_mode == "delayed") {
    lmdata$.LM_entry <- analysis_entry
    if (any(lmdata$.LM_entry >= lmdata[[outcome$time]]))
      stop("Internal error: non-positive observation interval after entry filtering.")
  }
  if (format == "long")
    cols <- match(c(id, outcome$time, outcome$status, "LM", covs$fixed,
                    covs$varying, rtime,
                    if (entry_mode != "shared") entry,
                    if (entry_mode == "delayed") ".LM_entry"),
                  names(lmdata))
  else cols <- match(c(outcome$time, outcome$status, "LM", covs$fixed,
                       covs$varying,
                       if (entry_mode != "shared") entry,
                       if (entry_mode == "delayed") ".LM_entry"),
                     names(lmdata))
  cols <- unique(cols[!is.na(cols)])
  return(lmdata[, cols])
}
