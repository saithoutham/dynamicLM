#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(devtools)
  library(haven)
  library(here)
  library(survival)
  library(survey)
})
source(here::here("analysis", "_helpers.R"))
source(here::here("analysis", "_age_data.R"))
source(here::here("analysis", "_simulation_engine.R"))
devtools::load_all(here::here(), quiet = TRUE)

analysis_seed <- 20260908L
assert_single_seed(analysis_seed)
set.seed(analysis_seed)

prediction_window <- 10
cv_folds <- 5L
bootstrap_replicates <- as.integer(Sys.getenv("APPLICATION_BOOTSTRAPS", "1000"))
stopifnot(bootstrap_replicates >= 100L)

condition_log <- data.table(
  component = character(), severity = character(), condition = character()
)
record_condition <- function(component, severity, condition) {
  condition_log <<- rbind(
    condition_log,
    data.table(component = component, severity = severity,
               condition = as.character(condition))
  )
  invisible(NULL)
}

download_public_file <- function(url, filename) {
  directory <- project_path("data", "public")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(directory, filename)
  start <- proc.time()[["elapsed"]]
  result <- capture_condition({
    if (!file.exists(destination)) {
      status <- utils::download.file(url, destination, mode = "wb", quiet = TRUE)
      stopifnot(identical(status, 0L))
    }
    stopifnot(file.exists(destination), file.info(destination)$size > 0)
    destination
  })
  runtime <- proc.time()[["elapsed"]] - start
  if (inherits(result$value, "captured_error")) {
    record_condition(filename, "error", as.character(result$value))
    return(data.table(
      file = filename, url = url, success = 0L, bytes = 0,
      md5 = "", runtime_seconds = runtime, local_path = ""
    ))
  }
  for (warning in result$warnings)
    record_condition(filename, "warning", warning)
  data.table(
    file = filename, url = url, success = 1L,
    bytes = as.numeric(file.info(result$value)$size),
    md5 = unname(tools::md5sum(result$value)), runtime_seconds = runtime,
    local_path = result$value
  )
}

prepare_nhanes <- function() {
  downloads <- rbindlist(list(
    download_public_file(
      paste0("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/",
             "linked_mortality/NHANES_1999_2000_MORT_2019_PUBLIC.dat"),
      "NHANES_1999_2000_MORT_2019_PUBLIC.dat"
    ),
    download_public_file(
      paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/",
             "DataFiles/DEMO.XPT"),
      "NHANES_1999_2000_DEMO.XPT"
    ),
    download_public_file(
      paste0("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/",
             "linked_mortality/NHIS_2018_MORT_2019_PUBLIC.dat"),
      "NHIS_2018_MORT_2019_PUBLIC.dat"
    )
  ))
  if (any(downloads$success == 0L))
    return(list(data = NULL, downloads = downloads))

  lmf_path <- downloads[file == "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
                        local_path]
  lines <- readLines(lmf_path, warn = FALSE)
  # The CDC fixed-width records can omit a final blank position; substr()
  # correctly returns the available 46--48 field in both 46- and 48-column
  # records.
  stopifnot(length(lines) > 0L, all(nchar(lines) >= 46L))
  number <- function(first, last)
    suppressWarnings(as.integer(substr(lines, first, last)))
  mortality <- data.table(
    SEQN = number(1L, 6L), ELIGSTAT = number(15L, 15L),
    MORTSTAT = number(16L, 16L), UCOD_LEADING = number(17L, 19L),
    PERMTH_EXM = number(46L, 48L)
  )
  demo_path <- downloads[file == "NHANES_1999_2000_DEMO.XPT", local_path]
  demo <- as.data.table(haven::read_xpt(demo_path))
  required <- c("SEQN", "RIDAGEYR", "RIAGENDR", "WTMEC2YR", "SDMVPSU",
                "SDMVSTRA")
  stopifnot(all(required %in% names(demo)), !anyDuplicated(demo$SEQN),
            !anyDuplicated(mortality$SEQN))
  merged <- merge(demo[, ..required], mortality, by = "SEQN", all = FALSE)
  merged <- merged[
    ELIGSTAT == 1L & RIDAGEYR >= 18 & is.finite(PERMTH_EXM) &
      PERMTH_EXM > 0 & WTMEC2YR > 0 & complete.cases(merged[, ..required])
  ]
  merged[, `:=`(
    id = as.integer(SEQN), entry_age = as.numeric(RIDAGEYR),
    exit_age = as.numeric(RIDAGEYR) + PERMTH_EXM / 12,
    event = as.integer(MORTSTAT == 1L), male = as.integer(RIAGENDR == 1L)
  )]
  stopifnot(nrow(merged) > 0L, !anyDuplicated(merged$id),
            all(merged$entry_age < merged$exit_age),
            all(merged$event %in% 0:1))
  list(data = as.data.frame(merged), downloads = downloads)
}

prepare_mgus2 <- function() {
  data <- as.data.table(survival::mgus2)
  data <- data[complete.cases(data[, .(id, age, sex, hgb, creat, mspike,
                                      ptime, pstat, futime, death)])]
  data[, `:=`(
    entry_age = as.numeric(age),
    progression_age = fifelse(pstat == 1L, age + ptime / 12, Inf),
    death_age = fifelse(death == 1L, age + futime / 12, Inf),
    followup_age = age + futime / 12,
    male = as.integer(sex == "M")
  )]
  data[, exit_age := pmin(progression_age, followup_age)]
  data[, status := fifelse(pstat == 1L & progression_age <= followup_age, 1L,
                           fifelse(death == 1L, 2L, 0L))]
  data[, `:=`(
    hgb_z = as.numeric(scale(hgb)), creat_z = as.numeric(scale(creat)),
    mspike_z = as.numeric(scale(mspike))
  )]
  derived <- c("id", "entry_age", "exit_age", "status", "male", "hgb_z",
               "creat_z", "mspike_z")
  stopifnot(!anyDuplicated(data$id), !anyNA(data[, ..derived]),
            all(data$entry_age < data$exit_age), all(data$status %in% 0:2),
            all(data[pstat == 1L, ptime <= futime]))
  as.data.frame(data)
}

prepare_flchain <- function() {
  data <- as.data.table(survival::flchain)
  data[, id := seq_len(.N)]
  data <- data[complete.cases(data[, .(id, age, sex, kappa, lambda,
                                      creatinine, futime, death)]) &
                 kappa >= 0 & lambda >= 0 & creatinine > 0 & futime > 0]
  data[, `:=`(
    entry_age = as.numeric(age), exit_age = age + futime / 365.25,
    event = as.integer(death == 1L), male = as.integer(sex == "M"),
    log_kappa_z = as.numeric(scale(log1p(kappa))),
    log_lambda_z = as.numeric(scale(log1p(lambda))),
    creat_z = as.numeric(scale(log(creatinine)))
  )]
  derived <- c("id", "entry_age", "exit_age", "event", "male",
               "log_kappa_z", "log_lambda_z", "creat_z")
  stopifnot(!anyDuplicated(data$id), !anyNA(data[, ..derived]),
            all(data$entry_age < data$exit_age), all(data$event %in% 0:1))
  as.data.frame(data)
}

make_stack <- function(data, outcome, landmarks, w, fixed, varying = NULL,
                       format = "wide", rtime = NULL,
                       mode = c("strict", "delayed", "shared")) {
  mode <- match.arg(mode)
  arguments <- list(
    data = data, outcome = outcome, lms = landmarks, w = w,
    covs = list(fixed = fixed, varying = varying), format = format,
    id = "id", entry_mode = mode
  )
  if (mode != "shared") arguments$entry <- "entry_age"
  if (format == "long") arguments$rtime <- rtime
  stack <- do.call(dynamicLM::stack_data, arguments)
  stopifnot(inherits(stack, "LMdataframe"),
            !anyDuplicated(stack$data[, c("id", "LM")]),
            all(stack$data[[stack$entry_col]] < stack$data[[outcome$time]]))
  if (mode == "strict")
    stopifnot(all(stack$data$entry_age <= stack$data$LM))
  if (mode == "delayed")
    stopifnot(all(stack$data$.LM_entry ==
                    pmax(stack$data$entry_age, stack$data$LM)))
  stack
}

fit_cause_specific <- function(stack, dataset, status_col, cause,
                               base_terms, interaction_terms,
                               landmark_center, landmark_scale = 10,
                               survey_design = FALSE) {
  data <- stack$data
  data$analysis_start <- data[[stack$entry_col]]
  data$analysis_event <- as.integer(data[[status_col]] == cause)
  data$LM_scaled <- (data$LM - landmark_center) / landmark_scale
  interaction_names <- paste0(interaction_terms, "_LM_scaled")
  for (i in seq_along(interaction_terms))
    data[[interaction_names[i]]] <- data[[interaction_terms[i]]] * data$LM_scaled
  required <- unique(c("analysis_start", stack$outcome$time, "analysis_event",
                       "id", base_terms, interaction_names,
                       if (survey_design) c("SDMVPSU", "SDMVSTRA", "WTMEC2YR")))
  keep <- complete.cases(data[, required, drop = FALSE]) &
    apply(data[, required, drop = FALSE], 1L,
          function(row) all(is.finite(as.numeric(row))))
  data <- data[keep, , drop = FALSE]
  stopifnot(nrow(data) > 0L, all(data$analysis_start < data[[stack$outcome$time]]))
  right <- paste(c(base_terms, interaction_names), collapse = " + ")
  formula <- stats::as.formula(paste0(
    "Surv(analysis_start, ", stack$outcome$time, ", analysis_event) ~ ",
    right, if (!survey_design) " + cluster(id)" else ""
  ))
  component <- paste(dataset, stack$entry_mode, cause,
                     if (survey_design) "survey" else "coxph", sep = ":")
  captured <- capture_condition({
    if (survey_design) {
      old_lonely <- getOption("survey.lonely.psu")
      on.exit(options(survey.lonely.psu = old_lonely), add = TRUE)
      options(survey.lonely.psu = "adjust")
      design <- survey::svydesign(
        id = ~SDMVPSU + id, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
        nest = TRUE, data = data
      )
      survey::svycoxph(formula, design = design)
    } else {
      survival::coxph(formula, data = data, ties = "breslow", x = TRUE)
    }
  })
  for (warning in captured$warnings)
    record_condition(component, "warning", warning)
  if (inherits(captured$value, "captured_error")) {
    record_condition(component, "error", as.character(captured$value))
    return(NULL)
  }
  fit <- captured$value
  estimates <- stats::coef(fit)
  standard_errors <- sqrt(diag(stats::vcov(fit)))
  stopifnot(length(estimates) > 0L, all(is.finite(estimates)),
            all(is.finite(standard_errors)),
            identical(names(estimates), names(standard_errors)))
  data.table(
    dataset = dataset, entry_mode = stack$entry_mode,
    estimator = if (survey_design) "svycoxph" else "coxph_cluster",
    cause = cause, observations = nrow(data),
    unique_subjects = uniqueN(data$id), events = sum(data$analysis_event),
    term = names(estimates), estimate = unname(estimates),
    standard_error = unname(standard_errors)
  )
}

risk_set_table <- function(stack, dataset, status_col, causes) {
  rbindlist(lapply(sort(unique(stack$data$LM)), function(landmark) {
    data <- stack$data[stack$data$LM == landmark, , drop = FALSE]
    rows <- lapply(causes, function(cause) data.table(
      dataset = dataset, entry_mode = stack$entry_mode, landmark = landmark,
      cause = cause, rows = nrow(data), unique_subjects = uniqueN(data$id),
      events = sum(data[[status_col]] == cause)
    ))
    rbindlist(rows)
  }))
}

add_visit_features <- function(stack, long_data) {
  data <- as.data.table(stack$data)
  visits <- split(long_data$visit_age, long_data$id)
  data[, staleness := LM - visit_age]
  data[, visit_count := mapply(
    function(subject, landmark) {
      ages <- visits[[as.character(subject)]]
      sum(ages <= landmark)
    }, id, LM
  )]
  data[, visit_exposure := pmax(LM - entry_age, 0) + 0.25]
  data[, entry_age_z := as.numeric(scale(entry_age))]
  stopifnot(all(is.finite(data$staleness)), all(data$staleness >= 0),
            all(data$visit_count >= 1L), all(data$visit_exposure > 0),
            !anyNA(data))
  data
}

safe_paired_bootstrap <- function(stack, risk_a, risk_b, B, seed, w) {
  stopifnot(nrow(stack) == length(risk_a), nrow(stack) == length(risk_b),
            B >= 2L, seed == as.integer(seed))
  point_a <- summary_auc_fast(stack, risk_a, w)
  point_b <- summary_auc_fast(stack, risk_b, w)
  ids <- unique(stack$id)
  row_index <- split(seq_len(nrow(stack)), stack$id)
  set.seed(seed)
  errors <- character()
  boot <- vapply(seq_len(B), function(iteration) {
    sampled <- sample(ids, length(ids), replace = TRUE)
    index <- unlist(row_index[as.character(sampled)], use.names = FALSE)
    boot_stack <- stack[index, , drop = FALSE]
    boot_stack$id <- rep(seq_along(sampled),
                         vapply(row_index[as.character(sampled)], length,
                                integer(1L)))
    tryCatch(
      summary_auc_fast(boot_stack, risk_a[index], w) -
        summary_auc_fast(boot_stack, risk_b[index], w),
      error = function(error) {
        errors <<- c(errors, conditionMessage(error))
        NA_real_
      }
    )
  }, numeric(1L))
  successful <- boot[is.finite(boot)]
  if (length(successful) < ceiling(0.90 * B))
    stop("Fewer than 90% of bootstrap resamples were evaluable: ",
         length(successful), "/", B, "; ",
         paste(unique(errors), collapse = " | "), call. = FALSE)
  list(
    estimate_a = point_a, estimate_b = point_b,
    delta = point_a - point_b, standard_error = stats::sd(successful),
    successful = length(successful), failed = sum(!is.finite(boot)),
    errors = errors
  )
}

run_irregular_oof <- function(stack, long_data, dataset, base_terms,
                              interaction_terms, seed) {
  required <- unique(c("id", "LM", "visit_age", "entry_age", "event",
                       stack$outcome$time, base_terms, interaction_terms))
  keep <- complete.cases(stack$data[, required, drop = FALSE]) &
    apply(stack$data[, required, drop = FALSE], 1L,
          function(row) all(is.finite(as.numeric(row))))
  stack$data <- stack$data[keep, , drop = FALSE]
  stopifnot(nrow(stack$data) > 0L,
            !anyDuplicated(stack$data[, c("id", "LM")]))
  data <- add_visit_features(stack, long_data)
  data[, LM_scaled := (LM - mean(sort(unique(LM)))) / 10]
  data[, staleness_z := as.numeric(scale(staleness))]
  interaction_names <- paste0(interaction_terms, "_LM_scaled")
  for (i in seq_along(interaction_terms))
    data[[interaction_names[i]]] <- data[[interaction_terms[i]]] * data$LM_scaled
  eval_stack <- data.frame(
    id = data$id, landmark = data$LM, start = data$LM,
    stop = data[[stack$outcome$time]], event = data$event
  )
  stopifnot(!anyDuplicated(eval_stack[, c("id", "landmark")]),
            all(eval_stack$start < eval_stack$stop))

  ids <- unique(data$id)
  set.seed(seed)
  folds <- setNames(sample(rep(seq_len(cv_folds), length.out = length(ids))), ids)
  stopifnot(length(folds) == length(ids), all(table(folds) > 0))
  predictions <- matrix(NA_real_, nrow(data), 3L,
                        dimnames = list(NULL, c("base", "staleness", "iiw")))
  weight_caps <- numeric(cv_folds)
  for (fold in seq_len(cv_folds)) {
    test <- which(folds[as.character(data$id)] == fold)
    train <- setdiff(seq_len(nrow(data)), test)
    stopifnot(length(intersect(data$id[train], data$id[test])) == 0L)
    train_data <- as.data.frame(data[train])
    test_data <- as.data.frame(data[test])

    intensity <- stats::glm(
      visit_count ~ male + entry_age_z + LM_scaled,
      offset = log(visit_exposure), family = stats::poisson(),
      data = train_data
    )
    predicted_rate <- as.numeric(stats::predict(
      intensity, newdata = train_data, type = "response"
    )) / train_data$visit_exposure
    stopifnot(all(is.finite(predicted_rate)), all(predicted_rate > 0))
    inverse_rate <- 1 / predicted_rate
    cap <- unname(stats::quantile(inverse_rate, 0.99, type = 8))
    train_data$iiw <- pmin(inverse_rate, cap)
    train_data$iiw <- train_data$iiw / mean(train_data$iiw)
    weight_caps[fold] <- cap
    stopifnot(all(is.finite(train_data$iiw)), all(train_data$iiw > 0),
              abs(mean(train_data$iiw) - 1) < 1e-12)

    rhs_base <- paste(c(base_terms, interaction_names), collapse = " + ")
    formulas <- list(
      base = stats::as.formula(paste0(
        "Surv(LM, ", stack$outcome$time, ", event) ~ ",
        rhs_base, " + cluster(id)"
      )),
      staleness = stats::as.formula(paste0(
        "Surv(LM, ", stack$outcome$time, ", event) ~ ",
        rhs_base, " + staleness_z + cluster(id)"
      )),
      iiw = stats::as.formula(paste0(
        "Surv(LM, ", stack$outcome$time, ", event) ~ ",
        rhs_base, " + cluster(id)"
      ))
    )
    for (variant in names(formulas)) {
      captured <- capture_condition(survival::coxph(
        formulas[[variant]], data = train_data, ties = "breslow",
        weights = if (variant == "iiw") train_data$iiw else NULL,
        x = TRUE
      ))
      for (warning in captured$warnings)
        record_condition(paste(dataset, "fold", fold, variant, sep = ":"),
                         "warning", warning)
      if (inherits(captured$value, "captured_error"))
        stop(as.character(captured$value), call. = FALSE)
      predictions[test, variant] <- as.numeric(stats::predict(
        captured$value, newdata = test_data, type = "lp", reference = "zero"
      ))
    }
  }
  stopifnot(all(is.finite(predictions)))
  aucs <- vapply(colnames(predictions), function(variant)
    summary_auc_fast(eval_stack, predictions[, variant], stack$w), numeric(1L))
  comparisons <- rbindlist(lapply(c("staleness", "iiw"), function(variant) {
    bootstrap_seed <- seed + if (variant == "staleness") 10000L else 20000L
    result <- safe_paired_bootstrap(
      eval_stack, predictions[, variant], predictions[, "base"],
      bootstrap_replicates, bootstrap_seed, stack$w
    )
    if (result$failed > 0L)
      record_condition(
        paste(dataset, variant, "bootstrap", sep = ":"), "warning",
        paste(result$failed, "resamples failed:",
              paste(unique(result$errors), collapse = " | "))
      )
    z <- if (result$standard_error > 0)
      result$delta / result$standard_error else NA_real_
    data.table(
      dataset = dataset, variant = variant, base_auc = result$estimate_b,
      variant_auc = result$estimate_a, delta_auc = result$delta,
      bootstrap_se = result$standard_error,
      p_value = if (is.finite(z)) 2 * stats::pnorm(-abs(z)) else NA_real_,
      bootstrap_successful = result$successful,
      bootstrap_failed = result$failed, bootstrap_seed = bootstrap_seed,
      folds = cv_folds, subjects = length(ids), rows = nrow(data),
      events = sum(data$event), median_staleness = stats::median(data$staleness),
      p90_staleness = unname(stats::quantile(data$staleness, 0.90, type = 8)),
      mean_weight_cap = mean(weight_caps)
    )
  }))
  list(summary = comparisons,
       seed_manifest = data.table(dataset = dataset, id = ids,
                                  fold = unname(folds[as.character(ids)]),
                                  fold_seed = seed))
}

trace_table_numeric <- function(report, table, prefix_columns, value_columns,
                                script, function_name, seed, runtime) {
  for (row_index in seq_len(nrow(table))) {
    prefix <- paste(unlist(table[row_index, ..prefix_columns]), collapse = "_")
    prefix <- gsub("[^A-Za-z0-9]+", "_", tolower(prefix))
    for (column in value_columns) {
      value <- as.numeric(table[[column]][row_index])
      if (!is.finite(value)) next
      trace_values(
        report = report, labels = paste0(prefix, "_", column), values = value,
        units = column, script = script, function_name = function_name,
        seed = seed, runtime_seconds = runtime
      )
    }
  }
}

run <- timed({
  pbc <- make_pbc_age()
  pbc$long$log_bili_z <- as.numeric(scale(pbc$long$log_bili))
  pbc$long$albumin_z <- as.numeric(scale(pbc$long$albumin_locf))
  nafld <- make_nafld_age(include_preindex_labs = FALSE)
  nhanes <- prepare_nhanes()
  mgus <- prepare_mgus2()
  flchain <- prepare_flchain()

  stack_specs <- list()
  for (mode in c("strict", "delayed")) {
    stack_specs[[paste("nafld", mode, sep = "_")]] <- make_stack(
      nafld$long, list(time = "exit_age", status = "event"), c(50, 60, 70),
      prediction_window, c("entry_age", "male", "bmi_z"), "hdl_z",
      "long", "visit_age", mode
    )
    stack_specs[[paste("pbcseq_age", mode, sep = "_")]] <- make_stack(
      pbc$long, list(time = "exit_age", status = "event"), c(40, 50, 60),
      5, c("entry_age", "male", "stage0"),
      c("log_bili_z", "albumin_z"), "long", "visit_age", mode
    )
    stack_specs[[paste("mgus2", mode, sep = "_")]] <- make_stack(
      mgus, list(time = "exit_age", status = "status"), c(60, 70, 80, 90),
      5, c("entry_age", "male", "hgb_z", "creat_z", "mspike_z"),
      mode = mode
    )
    stack_specs[[paste("flchain", mode, sep = "_")]] <- make_stack(
      flchain, list(time = "exit_age", status = "event"), c(60, 70, 80, 90),
      5, c("entry_age", "male", "log_kappa_z", "log_lambda_z", "creat_z"),
      mode = mode
    )
    if (!is.null(nhanes$data))
      stack_specs[[paste("nhanes_1999_2000", mode, sep = "_")]] <- make_stack(
        nhanes$data, list(time = "exit_age", status = "event"),
        c(50, 60, 70, 80), prediction_window,
        c("entry_age", "male", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"),
        mode = mode
      )
  }

  pbc_enrollment <- transform(
    pbc$long, exit_enrollment = exit_age - entry_age,
    visit_enrollment = visit_age - entry_age
  )
  stack_specs$pbcseq_enrollment_shared <- make_stack(
    pbc_enrollment, list(time = "exit_enrollment", status = "event"),
    c(0, 2, 4), 5, c("entry_age", "male", "stage0"),
    c("log_bili_z", "albumin_z"), "long", "visit_enrollment", "shared"
  )
  enrollment_counts <- aggregate(id ~ LM,
                                 stack_specs$pbcseq_enrollment_shared$data,
                                 length)$id
  stopifnot(all(diff(enrollment_counts) <= 0))

  # Minimal reproducer for the discarded 10-year PBC evaluation design.  The
  # age-40 risk set has too few horizon controls for stable subject bootstrap
  # resampling, even though the point estimate itself is defined.
  pbc_w10 <- make_stack(
    pbc$long, list(time = "exit_age", status = "event"), c(40, 50, 60),
    prediction_window, c("entry_age", "male", "stage0"),
    c("log_bili_z", "albumin_z"), "long", "visit_age", "strict"
  )
  pbc_w10_data <- as.data.table(pbc_w10$data)
  pbc_w10_data <- pbc_w10_data[complete.cases(
    pbc_w10_data[, .(id, LM, exit_age, event, log_bili_z, albumin_z, visit_age)]
  )]
  pbc_w10_eval <- data.frame(
    id = pbc_w10_data$id, landmark = pbc_w10_data$LM,
    start = pbc_w10_data$LM, stop = pbc_w10_data$exit_age,
    event = pbc_w10_data$event
  )
  pbc_w10_design <- pbc_w10_data[, .(
    rows = .N, events = sum(event),
    horizon_controls = sum(event == 0L & exit_age >= LM + prediction_window)
  ), by = LM]
  pbc_w10_failure <- capture_condition(safe_paired_bootstrap(
    pbc_w10_eval, pbc_w10_data$log_bili_z, pbc_w10_data$log_bili_z,
    100L, analysis_seed + 303L, prediction_window
  ))
  stopifnot(inherits(pbc_w10_failure$value, "captured_error"))
  pbc_w10_failure_text <- as.character(pbc_w10_failure$value)
  matched_counts <- regmatches(
    pbc_w10_failure_text,
    regexec("evaluable: ([0-9]+)/([0-9]+)", pbc_w10_failure_text)
  )[[1L]]
  stopifnot(length(matched_counts) == 3L)
  pbc_w10_bootstrap_successful <- as.integer(matched_counts[2L])
  pbc_w10_bootstrap_attempted <- as.integer(matched_counts[3L])
  record_condition("pbcseq_age_w10_bootstrap_minimal", "error",
                   pbc_w10_failure_text)

  risk_sets <- rbindlist(list(
    risk_set_table(stack_specs$nafld_strict, "nafld", "event", 1L),
    risk_set_table(stack_specs$nafld_delayed, "nafld", "event", 1L),
    risk_set_table(stack_specs$pbcseq_age_strict, "pbcseq_age", "event", 1L),
    risk_set_table(stack_specs$pbcseq_age_delayed, "pbcseq_age", "event", 1L),
    risk_set_table(stack_specs$pbcseq_enrollment_shared, "pbcseq_enrollment",
                   "event", 1L),
    risk_set_table(stack_specs$mgus2_strict, "mgus2", "status", 1:2),
    risk_set_table(stack_specs$mgus2_delayed, "mgus2", "status", 1:2),
    risk_set_table(stack_specs$flchain_strict, "flchain", "event", 1L),
    risk_set_table(stack_specs$flchain_delayed, "flchain", "event", 1L),
    if (!is.null(nhanes$data)) risk_set_table(
      stack_specs$nhanes_1999_2000_strict, "nhanes_1999_2000", "event", 1L
    ),
    if (!is.null(nhanes$data)) risk_set_table(
      stack_specs$nhanes_1999_2000_delayed, "nhanes_1999_2000", "event", 1L
    )
  ), fill = TRUE)
  stopifnot(nrow(risk_sets) > 0L, all(risk_sets$rows == risk_sets$unique_subjects),
            all(risk_sets$events <= risk_sets$rows))

  models <- list()
  add_model <- function(name, ...)
    models[[name]] <<- fit_cause_specific(...)
  for (mode in c("strict", "delayed")) {
    add_model(paste("nafld", mode), stack_specs[[paste("nafld", mode, sep = "_")]],
              "nafld", "event", 1L, c("male", "bmi_z", "hdl_z"),
              c("bmi_z", "hdl_z"), 60)
    add_model(paste("pbcseq_age", mode),
              stack_specs[[paste("pbcseq_age", mode, sep = "_")]],
              "pbcseq_age", "event", 1L,
              c("male", "stage0", "log_bili_z", "albumin_z"),
              c("log_bili_z", "albumin_z"), 50)
    for (cause in 1:2)
      add_model(paste("mgus2", mode, cause),
                stack_specs[[paste("mgus2", mode, sep = "_")]],
                "mgus2", "status", cause,
                c("male", "hgb_z", "creat_z", "mspike_z"),
                c("hgb_z", "creat_z", "mspike_z"), 75)
    add_model(paste("flchain", mode),
              stack_specs[[paste("flchain", mode, sep = "_")]],
              "flchain", "event", 1L,
              c("male", "log_kappa_z", "log_lambda_z", "creat_z"),
              c("log_kappa_z", "log_lambda_z", "creat_z"), 75)
    if (!is.null(nhanes$data)) {
      add_model(paste("nhanes", mode, "unweighted"),
                stack_specs[[paste("nhanes_1999_2000", mode, sep = "_")]],
                "nhanes_1999_2000", "event", 1L, "male", "male", 65)
      add_model(paste("nhanes", mode, "survey"),
                stack_specs[[paste("nhanes_1999_2000", mode, sep = "_")]],
                "nhanes_1999_2000", "event", 1L, "male", "male", 65,
                survey_design = TRUE)
    }
  }
  add_model("pbcseq_enrollment_shared",
            stack_specs$pbcseq_enrollment_shared, "pbcseq_enrollment", "event",
            1L, c("male", "stage0", "log_bili_z", "albumin_z"),
            c("log_bili_z", "albumin_z"), 2)
  model_summary <- rbindlist(models, fill = TRUE)
  stopifnot(nrow(model_summary) > 0L, all(is.finite(model_summary$estimate)),
            all(is.finite(model_summary$standard_error)))

  irregular <- list(
    run_irregular_oof(
      stack_specs$nafld_strict, nafld$long, "nafld",
      c("male", "bmi_z", "hdl_z"), c("bmi_z", "hdl_z"),
      analysis_seed + 101L
    ),
    run_irregular_oof(
      stack_specs$pbcseq_age_strict, pbc$long, "pbcseq",
      c("male", "stage0", "log_bili_z", "albumin_z"),
      c("log_bili_z", "albumin_z"), analysis_seed + 202L
    )
  )
  irregular_summary <- rbindlist(lapply(irregular, `[[`, "summary"))
  irregular_seeds <- rbindlist(lapply(irregular, `[[`, "seed_manifest"))
  stopifnot(nrow(irregular_summary) == 4L,
            all(irregular_summary$bootstrap_successful +
                  irregular_summary$bootstrap_failed == bootstrap_replicates),
            !anyDuplicated(irregular_seeds[, .(dataset, id)]))

  list(
    downloads = nhanes$downloads, risk_sets = risk_sets,
    model_summary = model_summary, irregular_summary = irregular_summary,
    irregular_seeds = irregular_seeds,
    pbc_w10_design = pbc_w10_design,
    pbc_w10_bootstrap_successful = pbc_w10_bootstrap_successful,
    pbc_w10_bootstrap_attempted = pbc_w10_bootstrap_attempted,
    nhanes_rows = if (is.null(nhanes$data)) 0L else nrow(nhanes$data),
    nhanes_events = if (is.null(nhanes$data)) 0L else sum(nhanes$data$event)
  )
})

out <- run$value
fwrite(out$downloads, project_path("results", "application_downloads.csv"))
fwrite(out$risk_sets, project_path("results", "application_risk_sets.csv"))
fwrite(out$model_summary, project_path("results", "application_models.csv"))
fwrite(out$irregular_summary,
       project_path("results", "application_irregular_visits.csv"))
fwrite(out$irregular_seeds,
       project_path("results", "application_seed_manifest.csv"))
fwrite(out$pbc_w10_design,
       project_path("results", "application_pbc_w10_design_failure.csv"))
fwrite(condition_log, project_path("results", "application_conditions.csv"))
fwrite(data.table(runtime_seconds = run$runtime_seconds),
       project_path("results", "application_runtime.csv"))

trace_values(
  report = "06_application",
  labels = c("analysis_seed", "runtime_seconds", "bootstrap_replicates",
             "cv_folds", "condition_count", "nhanes_analysis_rows",
             "nhanes_analysis_events", "main_prediction_window_years",
             "pbc_application_window_years", "mgus_flchain_window_years",
             "landmark_interaction_scale_years",
             "pbc_w10_bootstrap_attempted", "pbc_w10_bootstrap_successful",
             "pbc_w10_bootstrap_minimum_required"),
  values = c(analysis_seed, run$runtime_seconds, bootstrap_replicates, cv_folds,
             nrow(condition_log), out$nhanes_rows, out$nhanes_events,
             prediction_window, 5, 5, 10, out$pbc_w10_bootstrap_attempted,
             out$pbc_w10_bootstrap_successful,
             ceiling(0.90 * out$pbc_w10_bootstrap_attempted)),
  units = c("seed", "seconds", "resamples", "folds", "conditions", "subjects",
            "events", "years", "years", "years", "years", "resamples",
            "resamples", "resamples"),
  script = "analysis/06_application.R", function_name = "run",
  seed = analysis_seed, runtime_seconds = run$runtime_seconds
)
trace_table_numeric(
  "06_application", out$downloads, "file", c("success", "bytes"),
  "analysis/06_application.R", "download_public_file", analysis_seed,
  run$runtime_seconds
)
trace_table_numeric(
  "06_application", out$pbc_w10_design, "LM",
  c("rows", "events", "horizon_controls"),
  "analysis/06_application.R", "safe_paired_bootstrap/minimal_reproducer",
  analysis_seed + 303L, run$runtime_seconds
)
trace_table_numeric(
  "06_application", out$risk_sets,
  c("dataset", "entry_mode", "landmark", "cause"),
  c("rows", "events"), "analysis/06_application.R", "risk_set_table",
  analysis_seed, run$runtime_seconds
)
trace_table_numeric(
  "06_application", out$model_summary,
  c("dataset", "entry_mode", "estimator", "cause", "term"),
  c("observations", "unique_subjects", "events", "estimate", "standard_error"),
  "analysis/06_application.R", "fit_cause_specific", analysis_seed,
  run$runtime_seconds
)
trace_table_numeric(
  "06_application", out$irregular_summary, c("dataset", "variant"),
  c("base_auc", "variant_auc", "delta_auc", "bootstrap_se", "p_value",
    "bootstrap_successful", "bootstrap_failed", "folds", "subjects", "rows",
    "events", "median_staleness", "p90_staleness", "mean_weight_cap"),
  "analysis/06_application.R", "run_irregular_oof",
  "results/application_seed_manifest.csv plus row-specific bootstrap_seed",
  run$runtime_seconds
)

stopifnot(file.exists(project_path("results", "application_models.csv")),
          file.exists(project_path("results", "application_irregular_visits.csv")),
          file.exists(project_path("results", "application_downloads.csv")))
