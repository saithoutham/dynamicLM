#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

trace_path <- here::here("results", "trace.csv")
report_paths <- sort(list.files(here::here("reports"), pattern = "\\.md$",
                                full.names = TRUE))
stopifnot(file.exists(trace_path), length(report_paths) > 0L)
trace <- utils::read.csv(trace_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("trace_id", "report", "label", "value", "units", "script",
              "function", "seed", "runtime_seconds", "timestamp_utc",
              "source_type", "source", "git_commit", "notes")
stopifnot(identical(names(trace), required), !anyDuplicated(trace$trace_id),
          all(nzchar(trace$report)), all(nzchar(trace$label)),
          all(nzchar(trace$script)), all(nzchar(trace$`function`)),
          all(nzchar(trace$seed)), all(is.finite(as.numeric(trace$value))),
          all(is.finite(as.numeric(trace$runtime_seconds))),
          all(nchar(trace$git_commit) == 40L))

extract_code_spans <- function(text) {
  matches <- gregexpr("`[^`]+`", text, perl = TRUE)[[1L]]
  if (identical(matches[1L], -1L)) return(character())
  substring(regmatches(text, list(matches))[[1L]], 2L,
            nchar(regmatches(text, list(matches))[[1L]]) - 1L)
}

references <- do.call(rbind, lapply(report_paths, function(path) {
  report <- tools::file_path_sans_ext(basename(path))
  spans <- extract_code_spans(paste(readLines(path, warn = FALSE),
                                    collapse = "\n"))
  full <- spans[grepl("^[0-9]{2}_[a-z_]+::[A-Za-z0-9_.-]+$", spans)]
  shorthand <- spans[grepl("^::[A-Za-z0-9_.-]+$", spans)]
  data.frame(
    report_file = report,
    trace_id = c(full, if (length(shorthand) > 0L)
      paste0(report, shorthand) else character()),
    stringsAsFactors = FALSE
  )
}))
references <- unique(references)
references$found <- references$trace_id %in% trace$trace_id
utils::write.csv(references, here::here("results", "trace_audit.csv"),
                 row.names = FALSE)
if (any(!references$found)) {
  stop("Missing trace references: ",
       paste(references$trace_id[!references$found], collapse = ", "),
       call. = FALSE)
}
stopifnot(nrow(references) > 0L, all(references$found))
