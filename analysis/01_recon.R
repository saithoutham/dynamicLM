#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(devtools)
  library(here)
  library(survival)
})
source(here::here("analysis", "_helpers.R"))
devtools::load_all(here::here(), quiet = TRUE)

get_dataset <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  env <- as.environment("package:survival")
  stopifnot(exists(name, envir = env, inherits = FALSE))
  data <- get(name, envir = env, inherits = FALSE)
  stopifnot(is.data.frame(data), nrow(data) > 0L, ncol(data) > 0L)
  data
}

age_summary <- function(dataset, entry_age) {
  stopifnot(length(entry_age) > 0L, all(is.finite(entry_age)))
  quantiles <- unname(stats::quantile(entry_age, c(0, 0.25, 0.5, 0.75, 1),
                                      na.rm = FALSE, type = 7))
  data.frame(
    dataset = dataset,
    minimum = quantiles[1L],
    q1 = quantiles[2L],
    median = quantiles[3L],
    mean = mean(entry_age),
    q3 = quantiles[4L],
    maximum = quantiles[5L],
    sd = stats::sd(entry_age),
    stringsAsFactors = FALSE
  )
}

event_rows <- function(dataset, event) {
  stopifnot(length(event) > 0L, !anyNA(event))
  counts <- table(event, useNA = "no")
  stopifnot(sum(counts) == length(event))
  data.frame(
    dataset = dataset,
    event_code = names(counts),
    count = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

analysis_runtime <- timed({
  nafld1 <- get_dataset("nafld1")
  nafld2 <- get_dataset("nafld2")
  nafld3 <- get_dataset("nafld3")
  pbcseq <- get_dataset("pbcseq")
  mgus2 <- get_dataset("mgus2")
  flchain <- get_dataset("flchain")
  rotterdam <- get_dataset("rotterdam")
  colon <- get_dataset("colon")

  stopifnot(!anyDuplicated(nafld1$id))
  pbc_first <- pbcseq[!duplicated(pbcseq$id), ]
  stopifnot(all(pbc_first$day == ave(pbcseq$day, pbcseq$id, FUN = min)[!duplicated(pbcseq$id)]))
  stopifnot(!anyDuplicated(mgus2$id), !anyDuplicated(rotterdam$pid))
  stopifnot(all(mgus2$ptime[mgus2$pstat == 1L] <= mgus2$futime[mgus2$pstat == 1L]))

  colon_subject <- colon[!duplicated(colon$id), ]
  stopifnot(nrow(colon) == 2L * nrow(colon_subject))
  stopifnot(all(table(colon$id) == 2L), setequal(unique(colon$etype), c(1, 2)))

  mgus_event <- ifelse(mgus2$pstat == 1L, 1L,
                       ifelse(mgus2$death == 1L, 2L, 0L))
  rotterdam_event <- ifelse(
    rotterdam$recur == 1L &
      (rotterdam$death == 0L | rotterdam$rtime <= rotterdam$dtime),
    1L,
    ifelse(rotterdam$death == 1L, 2L, 0L)
  )

  inventory <- rbindlist(list(
    data.frame(dataset = "nafld1", rows = nrow(nafld1), columns = ncol(nafld1),
               unique_ids = length(unique(nafld1$id)),
               variables = paste(names(nafld1), collapse = " | "),
               entry_age_definition = "age at study entry",
               event_definition = "status: death vs censoring"),
    data.frame(dataset = "nafld2", rows = nrow(nafld2), columns = ncol(nafld2),
               unique_ids = length(unique(nafld2$id)),
               variables = paste(names(nafld2), collapse = " | "),
               entry_age_definition = "join to nafld1; age + earliest observed days / 365.25",
               event_definition = "laboratory measurements; no endpoint column"),
    data.frame(dataset = "nafld3", rows = nrow(nafld3), columns = ncol(nafld3),
               unique_ids = length(unique(nafld3$id)),
               variables = paste(names(nafld3), collapse = " | "),
               entry_age_definition = "join to nafld1; age at study entry",
               event_definition = "named event with days since index"),
    data.frame(dataset = "pbcseq", rows = nrow(pbcseq), columns = ncol(pbcseq),
               unique_ids = length(unique(pbcseq$id)),
               variables = paste(names(pbcseq), collapse = " | "),
               entry_age_definition = "age at enrollment",
               event_definition = "status: censoring, transplant, death"),
    data.frame(dataset = "mgus2", rows = nrow(mgus2), columns = ncol(mgus2),
               unique_ids = length(unique(mgus2$id)),
               variables = paste(names(mgus2), collapse = " | "),
               entry_age_definition = "age at MGUS diagnosis",
               event_definition = "first progression, death without prior progression, or censoring"),
    data.frame(dataset = "flchain", rows = nrow(flchain), columns = ncol(flchain),
               unique_ids = nrow(flchain),
               variables = paste(names(flchain), collapse = " | "),
               entry_age_definition = "age when blood sample obtained",
               event_definition = "death vs censoring"),
    data.frame(dataset = "rotterdam", rows = nrow(rotterdam), columns = ncol(rotterdam),
               unique_ids = length(unique(rotterdam$pid)),
               variables = paste(names(rotterdam), collapse = " | "),
               entry_age_definition = "age at surgery",
               event_definition = "first recurrence, death without prior recurrence, or censoring"),
    data.frame(dataset = "colon", rows = nrow(colon), columns = ncol(colon),
               unique_ids = length(unique(colon$id)),
               variables = paste(names(colon), collapse = " | "),
               entry_age_definition = "age at trial baseline",
               event_definition = "two endpoint rows per subject: recurrence and death")
  ), fill = TRUE)
  stopifnot(nrow(inventory) == 8L, !anyDuplicated(inventory$dataset))

  entry_ages <- rbindlist(list(
    age_summary("nafld1", nafld1$age),
    age_summary("pbcseq", pbc_first$age),
    age_summary("mgus2", mgus2$age),
    age_summary("flchain", flchain$age),
    age_summary("rotterdam", rotterdam$age),
    age_summary("colon", colon_subject$age)
  ))
  stopifnot(nrow(entry_ages) == 6L, !anyDuplicated(entry_ages$dataset))

  events <- rbindlist(list(
    event_rows("nafld1", nafld1$status),
    event_rows("pbcseq", pbc_first$status),
    event_rows("mgus2_first_event", mgus_event),
    event_rows("flchain", flchain$death),
    event_rows("rotterdam_first_event", rotterdam_event),
    event_rows("colon_endpoint_status", paste0("etype", colon$etype, "_status", colon$status))
  ))
  stopifnot(nrow(events) > 0L)

  toy <- data.frame(
    id = c(1L, 2L), entry = c(50, 60), exit = c(70, 75),
    event = c(1L, 1L), x = c(0, 1)
  )
  toy_outcome <- list(time = "exit", status = "event")
  toy_covs <- list(fixed = "x", varying = NULL)
  toy_landmark <- 55
  toy_stack <- stack_data(
    toy, toy_outcome, lms = toy_landmark, w = 10, covs = toy_covs,
    format = "long", id = "id", rtime = "entry"
  )
  correct_ids <- toy$id[toy$entry <= toy_landmark & toy_landmark < toy$exit]
  pre_entry_ids <- toy_stack$data$id[toy_stack$data$entry > toy_stack$data$LM]
  stopifnot(nrow(toy_stack$data) == 2L, identical(correct_ids, 1L),
            identical(pre_entry_ids, 2L))

  wide_error <- capture_condition(stack_data(
    toy[c("id", "exit", "event", "x")], toy_outcome,
    lms = toy_landmark, w = 10, covs = toy_covs, format = "wide", id = "id"
  ))
  stopifnot(inherits(wide_error$value, "captured_error"))

  fake_scores <- data.table(
    tLM = c(55, 60), model = "model", AUC = c(0.6, 0.7)
  )
  fake_iid <- data.table(
    ID = c(1L, 2L, 2L, 3L), tLM = c(55, 55, 60, 60),
    model = "model", IF.AUC = c(-0.1, 0.1, -0.2, 0.2)
  )
  summary_error <- capture_condition(dynamicLM:::summary_metric(
    metric = "AUC", df_t = fake_scores, df_c = NULL, df_iid = fake_iid,
    conf_int = 0.95, object = list(), id_col = "ID", B = 1L, se.fit = TRUE
  ))
  stopifnot(inherits(summary_error$value, "captured_error"))

  list(
    inventory = inventory,
    entry_ages = entry_ages,
    events = events,
    toy_stack = toy_stack$data,
    correct_ids = correct_ids,
    pre_entry_ids = pre_entry_ids,
    wide_error = as.character(wide_error$value),
    summary_error = as.character(summary_error$value)
  )
})

out <- analysis_runtime$value
utils::write.csv(out$inventory, project_path("results", "recon_dataset_inventory.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(out$entry_ages, project_path("results", "recon_entry_age_summary.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(out$events, project_path("results", "recon_event_counts.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(out$toy_stack, project_path("results", "recon_minimal_stack.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(
  data.frame(case = c("wide_format_missing_rtime", "staggered_entry_summary_iid"),
             error = c(out$wide_error, out$summary_error)),
  project_path("results", "recon_errors.csv"), row.names = FALSE, na = ""
)

inventory_labels <- unlist(lapply(seq_len(nrow(out$inventory)), function(i) {
  prefix <- out$inventory$dataset[i]
  paste0(prefix, c("_rows", "_columns", "_unique_ids"))
}))
inventory_values <- unlist(lapply(seq_len(nrow(out$inventory)), function(i) {
  unname(as.numeric(out$inventory[i, c("rows", "columns", "unique_ids")]))
}))
trace_values(
  report = "01_recon",
  labels = inventory_labels,
  values = inventory_values,
  units = rep(c("rows", "columns", "subjects"), nrow(out$inventory)),
  script = "analysis/01_recon.R",
  function_name = "base::nrow/base::ncol/base::unique",
  runtime_seconds = analysis_runtime$runtime_seconds,
  notes = "Installed survival-package data; see results/recon_dataset_inventory.csv."
)

age_fields <- c("minimum", "q1", "median", "mean", "q3", "maximum", "sd")
age_labels <- unlist(lapply(seq_len(nrow(out$entry_ages)), function(i) {
  paste0(out$entry_ages$dataset[i], "_entry_age_", age_fields)
}))
age_values <- unlist(lapply(seq_len(nrow(out$entry_ages)), function(i) {
  unname(as.numeric(as.data.frame(out$entry_ages)[i, age_fields]))
}))
trace_values(
  report = "01_recon",
  labels = age_labels,
  values = age_values,
  units = "years",
  script = "analysis/01_recon.R",
  function_name = "age_summary",
  runtime_seconds = analysis_runtime$runtime_seconds,
  notes = "Entry-age distribution; see results/recon_entry_age_summary.csv."
)

trace_values(
  report = "01_recon",
  labels = paste0(out$events$dataset, "_event_", out$events$event_code),
  values = out$events$count,
  units = "records",
  script = "analysis/01_recon.R",
  function_name = "event_rows",
  runtime_seconds = analysis_runtime$runtime_seconds,
  notes = "Observed event-code counts; see results/recon_event_counts.csv."
)

trace_values(
  report = "01_recon",
  labels = c("toy_current_rows", "toy_correct_rows", "toy_preentry_rows",
             "toy_landmark", "toy_window"),
  values = c(nrow(out$toy_stack), length(out$correct_ids),
             length(out$pre_entry_ids), 55, 10),
  units = c("rows", "rows", "rows", "years", "years"),
  script = "analysis/01_recon.R",
  function_name = "dynamicLM::stack_data",
  runtime_seconds = analysis_runtime$runtime_seconds,
  notes = "Deterministic minimal reproducer; see results/recon_minimal_stack.csv."
)

stopifnot(
  file.exists(project_path("results", "recon_dataset_inventory.csv")),
  file.exists(project_path("results", "recon_entry_age_summary.csv")),
  file.exists(project_path("results", "recon_event_counts.csv")),
  file.exists(project_path("results", "recon_errors.csv"))
)
