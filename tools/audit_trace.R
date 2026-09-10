#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

trace_path <- here::here("results", "trace.csv")
report_paths <- sort(list.files(here::here("reports"), pattern = "\\.md$",
                                full.names = TRUE))
document_paths <- report_paths
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

extract_code_spans <- function(lines) {
  unlist(lapply(lines, function(line) {
    matches <- gregexpr("`[^`]+`", line, perl = TRUE)[[1L]]
    if (identical(matches[1L], -1L)) return(character())
    matched <- regmatches(line, list(matches))[[1L]]
    substring(matched, 2L, nchar(matched) - 1L)
  }), use.names = FALSE)
}

extract_trace_blocks <- function(lines) {
  text <- paste(lines, collapse = "\n")
  matches <- gregexpr("\\[TRACE:[^]]*\\]", text, perl = TRUE)[[1L]]
  if (identical(matches[1L], -1L)) return(character())
  blocks <- regmatches(text, list(matches))[[1L]]
  substring(blocks, nchar("[TRACE:") + 1L, nchar(blocks) - 1L)
}

references <- do.call(rbind, lapply(document_paths, function(path) {
  report <- tools::file_path_sans_ext(basename(path))
  trace_ids <- character()
  blocks <- extract_trace_blocks(readLines(path, warn = FALSE))
  for (block in blocks) {
    namespace <- NULL
    spans <- extract_code_spans(strsplit(block, "\n", fixed = TRUE)[[1L]])
    for (span in spans) {
      if (grepl("^[0-9]{2}_[a-z_]+::[A-Za-z0-9_.-]+$", span)) {
        namespace <- sub("::.*$", "", span)
        trace_ids <- c(trace_ids, span)
      } else if (grepl("^::[A-Za-z0-9_.-]+$", span)) {
        if (is.null(namespace)) {
          stop("Trace shorthand has no preceding full namespace in ", path,
               ": ", span, call. = FALSE)
        }
        trace_ids <- c(trace_ids, paste0(namespace, span))
      }
    }
  }
  data.frame(
    report_file = rep(report, length(trace_ids)),
    trace_id = trace_ids,
    stringsAsFactors = FALSE
  )
}))
references <- unique(references)
reference_matches <- vapply(
  references$trace_id,
  function(trace_id) sum(trace$trace_id == trace_id),
  integer(1L)
)
references$found <- reference_matches == 1L
utils::write.csv(references, here::here("results", "trace_audit.csv"),
                 row.names = FALSE)
if (any(!references$found)) {
  stop("Missing trace references: ",
       paste(references$trace_id[!references$found], collapse = ", "),
       call. = FALSE)
}
stopifnot(nrow(references) > 0L, all(references$found))
cat(
  "Trace audit passed:", nrow(references), "unique references across",
  length(document_paths), "documents; every reference resolves to one trace row.\n"
)
