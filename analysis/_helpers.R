# Shared helpers for reproducible analyses and number-level provenance.

project_path <- function(...) here::here(...)

assert_required_columns <- function(data, columns, label = deparse(substitute(data))) {
  stopifnot(is.data.frame(data), is.character(columns), length(columns) > 0L)
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(label, " is missing required columns: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_single_seed <- function(seed) {
  stopifnot(length(seed) == 1L, !is.na(seed), is.finite(seed))
  invisible(as.integer(seed))
}

git_commit <- function() {
  commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE)
  stopifnot(length(commit) == 1L, nchar(commit) == 40L)
  commit
}

init_trace <- function(reset = FALSE) {
  path <- project_path("results", "trace.csv")
  columns <- c(
    "trace_id", "report", "label", "value", "units", "script", "function",
    "seed", "runtime_seconds", "timestamp_utc", "source_type", "source",
    "git_commit", "notes"
  )
  if (reset || !file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(setNames(as.data.frame(matrix(nrow = 0L,
                                                   ncol = length(columns))),
                              columns), path, row.names = FALSE, na = "")
  }
  current <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(identical(names(current), columns))
  invisible(path)
}

trace_values <- function(report, labels, values, units, script, function_name,
                         seed = NA_integer_, runtime_seconds, source_type = "computed",
                         source = "", notes = "") {
  path <- init_trace()
  stopifnot(
    length(report) == 1L, nzchar(report),
    length(script) == 1L, nzchar(script),
    length(function_name) == 1L, nzchar(function_name),
    length(runtime_seconds) == 1L, is.finite(runtime_seconds),
    length(labels) == length(values),
    length(units) %in% c(1L, length(values)),
    length(source_type) %in% c(1L, length(values)),
    all(nzchar(labels)),
    all(is.finite(as.numeric(values)))
  )
  if (length(units) == 1L) units <- rep(units, length(values))
  if (length(source_type) == 1L) source_type <- rep(source_type, length(values))
  if (length(source) == 1L) source <- rep(source, length(values))
  if (length(notes) == 1L) notes <- rep(notes, length(values))
  if (is.character(seed)) {
    stopifnot(length(seed) == 1L, nzchar(seed))
    seed_text <- seed
  } else if (is.na(seed)) {
    seed_text <- "NA (deterministic)"
  } else {
    seed_text <- as.character(assert_single_seed(seed))
  }

  existing <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  stamp <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  keys <- paste(report, labels, sep = "::")
  rows <- data.frame(
    trace_id = keys,
    report = report,
    label = labels,
    value = format(values, digits = 17L, scientific = FALSE, trim = TRUE),
    units = units,
    script = script,
    `function` = function_name,
    seed = seed_text,
    runtime_seconds = format(runtime_seconds, digits = 17L,
                             scientific = FALSE, trim = TRUE),
    timestamp_utc = stamp,
    source_type = source_type,
    source = source,
    git_commit = git_commit(),
    notes = notes,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  existing <- existing[!existing$trace_id %in% keys, , drop = FALSE]
  out <- rbind(existing, rows)
  stopifnot(!anyDuplicated(out$trace_id))
  utils::write.csv(out, path, row.names = FALSE, na = "")
  invisible(rows)
}

drop_trace_values <- function(report, labels) {
  path <- init_trace()
  stopifnot(length(report) == 1L, nzchar(report), all(nzchar(labels)))
  existing <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  keys <- paste(report, labels, sep = "::")
  out <- existing[!existing$trace_id %in% keys, , drop = FALSE]
  stopifnot(!anyDuplicated(out$trace_id))
  utils::write.csv(out, path, row.names = FALSE, na = "")
  invisible(nrow(existing) - nrow(out))
}

timed <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, runtime_seconds = proc.time()[["elapsed"]] - start)
}

capture_condition <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "captured_error")),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}
