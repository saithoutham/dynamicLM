#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(parallel)
})
source(here::here("analysis", "_helpers.R"))
source(here::here("analysis", "_simulation_engine.R"))

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L, arguments %in% c("--verify", "--write"))
mode <- sub("^--", "", arguments)
cores <- as.integer(Sys.getenv(
  "REGENERATE_CORES", as.character(min(8L, parallel::detectCores()))
))
stopifnot(cores >= 1L)

truth_beta <- log(1.5)
lambda0 <- 0.01
entry_mean <- 45
prediction_window <- 10
nominal_coverage <- 0.95
methods <- c("naive", "strict", "delayed")
z_value <- stats::qnorm(1 - (1 - nominal_coverage) / 2)

identical_files <- function(first, second) {
  stopifnot(file.exists(first), file.exists(second))
  first_size <- file.info(first)$size
  second_size <- file.info(second)$size
  if (!identical(first_size, second_size)) return(FALSE)
  identical(
    readBin(first, what = "raw", n = first_size),
    readBin(second, what = "raw", n = second_size)
  )
}

regenerate_phase4 <- function() {
  grid <- fread(project_path("results", "simulation_grid.csv"))
  manifest <- fread(project_path("results", "simulation_seed_manifest.csv"))
  setorder(grid, cell_id)
  setorder(manifest, cell_id, replicate)
  stopifnot(
    nrow(grid) == 108L,
    nrow(manifest) == 108000L,
    !anyDuplicated(manifest[, .(cell_id, replicate)]),
    manifest[, all(.N == 1000L), by = cell_id]$V1
  )

  run_cell <- function(i) {
    row <- grid[i]
    cell_seeds <- manifest[cell_id == row$cell_id]
    landmarks <- 50 + c(0, row$landmark_interval,
                        2 * row$landmark_interval)
    output <- vector("list", nrow(cell_seeds) * length(methods))
    cursor <- 1L
    for (replicate_id in seq_len(nrow(cell_seeds))) {
      replicate_seed <- cell_seeds$seed[replicate_id]
      stopifnot(cell_seeds$replicate[replicate_id] == replicate_id)
      set.seed(replicate_seed)
      cohort <- generate_left_truncated_cohort(
        row$n, truth_beta, lambda0, entry_mean, row$entry_sd, row$study_end
      )
      cohort <- apply_independent_censoring(
        cohort, row$censor_rate, row$study_end
      )
      observed_entry_sd <- stats::sd(cohort$entry)
      observed_censoring <- mean(cohort$censored_before_end)
      stopifnot(is.finite(observed_entry_sd), is.finite(observed_censoring))
      for (method in methods) {
        stack <- make_simulation_stack(
          cohort, landmarks, prediction_window, method
        )
        fitted <- tryCatch(
          fit_fast_cluster_cox(stack),
          error = function(error)
            structure(conditionMessage(error), class = "regeneration_fit_error")
        )
        common <- data.table(
          cell_id = row$cell_id, replicate = replicate_id,
          seed = replicate_seed, method = method
        )
        if (inherits(fitted, "regeneration_fit_error")) {
          output[[cursor]] <- cbind(common, data.table(
            estimate = NA_real_, estimated_se = NA_real_,
            covered = NA_integer_, reject_null = NA_integer_,
            rows = nrow(stack), events = sum(stack$event),
            subjects = uniqueN(stack$id),
            observed_entry_sd = observed_entry_sd,
            observed_censoring = observed_censoring,
            converged = 0L, error = as.character(fitted)
          ))
        } else {
          output[[cursor]] <- cbind(common, data.table(
            estimate = fitted$estimate, estimated_se = fitted$robust_se,
            covered = as.integer(abs(fitted$estimate - truth_beta) <=
                                   z_value * fitted$robust_se),
            reject_null = as.integer(
              abs(fitted$estimate / fitted$robust_se) > z_value
            ),
            rows = fitted$rows, events = fitted$events,
            subjects = fitted$subjects,
            observed_entry_sd = observed_entry_sd,
            observed_censoring = observed_censoring,
            converged = 1L, error = ""
          ))
        }
        cursor <- cursor + 1L
      }
    }
    result <- rbindlist(output)
    stopifnot(nrow(result) == 3000L, all(result$cell_id == row$cell_id))
    result
  }

  runs <- mclapply(
    seq_len(nrow(grid)), run_cell, mc.cores = cores,
    mc.preschedule = TRUE, mc.set.seed = FALSE
  )
  if (any(vapply(runs, inherits, logical(1L), "try-error")))
    stop("Phase 4 raw regeneration worker failed.")
  output <- rbindlist(runs)
  stopifnot(
    nrow(output) == 324000L,
    output[, all(.N == 1000L), by = .(cell_id, method)]$V1,
    all(output$converged == 1L)
  )
  output
}

regenerate_phase6 <- function() {
  manifest <- fread(
    project_path("results", "informative_entry_seed_manifest.csv")
  )
  summary <- fread(project_path("results", "informative_entry_summary.csv"))
  calibration <- fread(
    project_path("results", "informative_entry_calibration.csv")
  )
  design_columns <- c(
    "phase", "regime", "cell_id", "gamma", "theta", "delta",
    "landmark_interval", "censoring_target", "n", "entry_sd"
  )
  design <- unique(summary[, ..design_columns])
  calibration_keys <- c(
    "phase", "gamma", "theta", "delta", "landmark_interval",
    "censoring_target", "entry_sd"
  )
  calibration <- unique(calibration[, c(
    calibration_keys, "study_end", "censor_rate"
  ), with = FALSE])
  design <- merge(
    design, calibration, by = calibration_keys, all.x = TRUE, sort = FALSE
  )
  design[, phase_order := match(phase, c("focused", "confirmation", "frailty"))]
  setorder(design, phase_order, cell_id)
  design[, phase_order := NULL]
  manifest[, phase_order := match(
    phase, c("focused", "confirmation", "frailty")
  )]
  setorder(manifest, phase_order, cell_id, replicate)
  manifest[, phase_order := NULL]
  stopifnot(
    nrow(design) == 132L,
    nrow(manifest) == 132000L,
    !anyNA(design),
    !anyDuplicated(manifest[, .(phase, cell_id, replicate)]),
    manifest[, all(.N == 1000L), by = .(phase, cell_id)]$V1
  )

  run_cell <- function(i) {
    row <- design[i]
    cell_seeds <- manifest[phase == row$phase & cell_id == row$cell_id]
    landmarks <- 50 + c(0, row$landmark_interval,
                        2 * row$landmark_interval)
    output <- vector("list", nrow(cell_seeds) * length(methods))
    cursor <- 1L
    for (replicate_id in seq_len(nrow(cell_seeds))) {
      replicate_seed <- cell_seeds$seed[replicate_id]
      stopifnot(cell_seeds$replicate[replicate_id] == replicate_id)
      set.seed(replicate_seed)
      cohort <- generate_left_truncated_cohort(
        row$n, truth_beta, lambda0, entry_mean, row$entry_sd, row$study_end,
        gamma = row$gamma, theta = row$theta, delta = row$delta
      )
      cohort <- apply_independent_censoring(
        cohort, row$censor_rate, row$study_end
      )
      entry_x_correlation <- stats::cor(cohort$entry, cohort$x)
      observed_censoring <- mean(cohort$censored_before_end)
      naive_stack <- make_simulation_stack(
        cohort, landmarks, prediction_window, "naive"
      )
      preentry_fractions <- vapply(landmarks, function(landmark) {
        rows <- naive_stack$landmark == landmark
        mean(cohort$entry[match(naive_stack$id[rows], cohort$id)] > landmark)
      }, numeric(1L))
      stopifnot(
        is.finite(entry_x_correlation), is.finite(observed_censoring),
        all(is.finite(preentry_fractions)),
        all(preentry_fractions >= 0 & preentry_fractions <= 1)
      )
      stacks <- list(
        naive = naive_stack,
        strict = make_simulation_stack(
          cohort, landmarks, prediction_window, "strict"
        ),
        delayed = make_simulation_stack(
          cohort, landmarks, prediction_window, "delayed"
        )
      )
      for (method in methods) {
        fitted <- tryCatch(
          fit_fast_cluster_cox(stacks[[method]]),
          error = function(error)
            structure(conditionMessage(error), class = "regeneration_fit_error")
        )
        common <- data.table(
          phase = row$phase, regime = row$regime, cell_id = row$cell_id,
          replicate = replicate_id, seed = replicate_seed, method = method,
          gamma = row$gamma, theta = row$theta, delta = row$delta,
          landmark_interval = row$landmark_interval,
          censoring_target = row$censoring_target, n = row$n,
          entry_sd = row$entry_sd,
          entry_x_correlation = entry_x_correlation,
          observed_censoring = observed_censoring,
          preentry_fraction_lm1 = preentry_fractions[1L],
          preentry_fraction_lm2 = preentry_fractions[2L],
          preentry_fraction_lm3 = preentry_fractions[3L]
        )
        if (inherits(fitted, "regeneration_fit_error")) {
          output[[cursor]] <- cbind(common, data.table(
            estimate = NA_real_, estimated_se = NA_real_,
            covered = NA_integer_, reject_null = NA_integer_,
            rows = nrow(stacks[[method]]),
            events = sum(stacks[[method]]$event),
            subjects = uniqueN(stacks[[method]]$id),
            converged = 0L, error = as.character(fitted)
          ))
        } else {
          output[[cursor]] <- cbind(common, data.table(
            estimate = fitted$estimate, estimated_se = fitted$robust_se,
            covered = as.integer(abs(fitted$estimate - truth_beta) <=
                                   z_value * fitted$robust_se),
            reject_null = as.integer(
              abs(fitted$estimate / fitted$robust_se) > z_value
            ),
            rows = fitted$rows, events = fitted$events,
            subjects = fitted$subjects, converged = 1L, error = ""
          ))
        }
        cursor <- cursor + 1L
      }
    }
    result <- rbindlist(output)
    stopifnot(nrow(result) == 3000L,
              all(result$phase == row$phase),
              all(result$cell_id == row$cell_id))
    result
  }

  runs <- mclapply(
    seq_len(nrow(design)), run_cell, mc.cores = cores,
    mc.preschedule = TRUE, mc.set.seed = FALSE
  )
  if (any(vapply(runs, inherits, logical(1L), "try-error")))
    stop("Phase 6 raw regeneration worker failed.")
  output <- rbindlist(runs)
  stopifnot(
    nrow(output) == 396000L,
    output[, all(.N == 1000L), by = .(phase, cell_id, method)]$V1,
    all(output$converged == 1L)
  )
  output
}

targets <- data.table(
  target = c("simulation_raw", "informative_entry_raw"),
  path = c(
    project_path("results", "simulation_raw.csv"),
    project_path("results", "informative_entry_raw.csv")
  )
)
if (mode == "verify") stopifnot(all(file.exists(targets$path)))
if (mode == "write" && any(file.exists(targets$path))) {
  stop("Refusing to overwrite an existing raw result. Use --verify first or ",
       "move the existing target aside.", call. = FALSE)
}

started <- proc.time()[["elapsed"]]
verification <- vector("list", nrow(targets))
for (i in seq_len(nrow(targets))) {
  generated <- if (targets$target[i] == "simulation_raw") {
    regenerate_phase4()
  } else {
    regenerate_phase6()
  }
  output_path <- if (mode == "verify") {
    tempfile(pattern = paste0(targets$target[i], "_"), fileext = ".csv")
  } else {
    targets$path[i]
  }
  fwrite(generated, output_path)
  byte_identical <- if (mode == "verify") {
    identical_files(output_path, targets$path[i])
  } else {
    NA
  }
  verification[[i]] <- data.table(
    target = targets$target[i],
    mode = mode,
    rows = nrow(generated),
    bytes = file.info(output_path)$size,
    byte_identical = if (mode == "verify") as.integer(byte_identical) else NA_integer_,
    md5 = unname(tools::md5sum(output_path))
  )
  if (mode == "verify" && !byte_identical) {
    stop("Byte-identical regeneration failed for ", targets$target[i],
         "; original file retained.", call. = FALSE)
  }
}
verification <- rbindlist(verification)
runtime_seconds <- proc.time()[["elapsed"]] - started
verification[, runtime_seconds := runtime_seconds]
fwrite(
  verification,
  project_path("results", "check", "raw_regeneration_verification.csv")
)

trace_values(
  report = "08_hygiene",
  labels = c(
    "regeneration_targets", "regeneration_cores", "regeneration_runtime_seconds"
  ),
  values = c(nrow(targets), cores, runtime_seconds),
  units = c("files", "cores", "seconds"),
  script = "analysis/regenerate_raw.R",
  function_name = paste0("manifest regeneration --", mode),
  seed = "results/simulation_seed_manifest.csv and results/informative_entry_seed_manifest.csv",
  runtime_seconds = runtime_seconds
)
for (i in seq_len(nrow(verification))) {
  row <- verification[i]
  fields <- c("rows", "bytes")
  values <- c(row$rows, row$bytes)
  if (mode == "verify") {
    fields <- c(fields, "byte_identical")
    values <- c(values, row$byte_identical)
  }
  trace_values(
    report = "08_hygiene",
    labels = paste0(row$target, "_", fields),
    values = values, units = fields,
    script = "analysis/regenerate_raw.R",
    function_name = paste0("manifest regeneration --", mode),
    seed = "seed manifests",
    runtime_seconds = runtime_seconds
  )
}
