#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(parallel)
})
source(here::here("analysis", "_helpers.R"))
source(here::here("analysis", "_simulation_engine.R"))

analysis_seed <- 20260909L
assert_single_seed(analysis_seed)
replicates <- as.integer(Sys.getenv("PHASE6_REPLICATES", "1000"))
regression_replicates <- as.integer(Sys.getenv(
  "PHASE6_REGRESSION_REPLICATES", "1000"
))
allow_small <- identical(Sys.getenv("SIM_ALLOW_SMALL", "0"), "1")
skip_confirmation <- identical(Sys.getenv("PHASE6_SKIP_CONFIRMATION", "0"), "1")
if (!allow_small) {
  stopifnot(replicates >= 1000L, regression_replicates == 1000L,
            !skip_confirmation)
}
stopifnot(replicates >= 1L, regression_replicates >= 1L,
          regression_replicates <= 1000L)
cores <- as.integer(Sys.getenv(
  "SIM_CORES", as.character(min(8L, parallel::detectCores()))
))
stopifnot(cores >= 1L)

truth_beta <- log(1.5)
lambda0 <- 0.01
entry_mean <- 45
prediction_window <- 10
nominal_coverage <- 0.95
landmark_intervals <- c(2, 4, 6)
censoring_targets <- c(0, 0.15, 0.30, 0.50)
sample_sizes <- c(500L, 750L, 1500L)
entry_sds <- c(0.25, 4, 10)
gammas <- c(0, -1, -2, -4)
methods <- c("naive", "strict", "delayed")
pilot_n <- 50000L

tolerance_changes <- data.table(
  phase = c("Phase 4", "Phase 4", "Phase 4", "Phase 6"),
  component = c("fresh censoring-calibration validation",
                "weighted AUC numerical bounds",
                "fast Cox versus coxph validation",
                "cross-process AUC truth-table reproduction"),
  original_value = c(1e-5, 0, 1e-10, 0),
  new_value = c(0.01, sqrt(.Machine$double.eps), 1e-8, 1e-12),
  failure_that_prompted_change = c(
    paste0("A separately drawn validation sample could differ from the target ",
           "by sampling error; the exact original failed value was not retained ",
           "and is PENDING."),
    paste0("Floating-point accumulation could put AUC infinitesimally outside ",
           "[0,1]; the exact original failed value was not retained and is PENDING."),
    paste0("Maximum absolute coefficient difference was ",
           format(1.1709017e-09, scientific = TRUE),
           ", exceeding the original threshold."),
    paste0("Recomputing the deterministic AUC truth table in a fresh process ",
           "differed from the saved table by at most 1.113554e-13; exact ",
           "identity failed.")
  ),
  justification = c(
    paste0("The root is calibrated internally; an independent finite pilot is ",
           "reported rather than required to meet root-finding precision."),
    paste0("The relaxed bound is machine-scale and values are clamped only after ",
           "the finite bound assertion."),
    paste0("The observed discrepancy remained far below Monte Carlo resolution ",
           "and both coefficients and robust standard errors were checked."),
    paste0("The threshold is seven orders of magnitude smaller than the ",
           "predeclared 1e-5 fine-versus-coarse quadrature tolerance and only ",
           "accommodates cross-process floating-point accumulation.")
  ),
  phase6_change = c(0L, 0L, 0L, 1L)
)
stopifnot(nrow(tolerance_changes) == 4L,
          all(tolerance_changes$new_value >= tolerance_changes$original_value),
          all(nzchar(tolerance_changes$failure_that_prompted_change)),
          all(nzchar(tolerance_changes$justification)))
fwrite(tolerance_changes, project_path("results", "tolerance_changes.csv"))

phase4_exact_regression <- function() {
  phase4_grid <- fread(project_path("results", "simulation_grid.csv"))
  phase4_seeds <- fread(project_path("results", "simulation_seed_manifest.csv"))
  phase4_seeds <- phase4_seeds[replicate <= regression_replicates]
  stopifnot(nrow(phase4_grid) == 108L,
            nrow(phase4_seeds) == 108L * regression_replicates,
            !anyDuplicated(phase4_seeds[, .(cell_id, replicate)]))

  old_source <- system2(
    "git", c("show", "f6bc7ec:analysis/_simulation_engine.R"),
    stdout = TRUE, stderr = TRUE
  )
  stopifnot(length(old_source) > 0L)
  old_engine <- new.env(parent = globalenv())
  eval(parse(text = old_source), envir = old_engine)

  check_cell <- function(i) {
    design <- phase4_grid[i]
    seeds <- phase4_seeds[cell_id == design$cell_id]
    exact_cohort <- exact_censored <- exact_stack <- exact_fit <- TRUE
    for (replicate_id in seq_len(nrow(seeds))) {
      seed <- seeds$seed[replicate_id]
      set.seed(seed)
      old_uncensored <- old_engine$generate_left_truncated_cohort(
        design$n, truth_beta, lambda0, entry_mean, design$entry_sd,
        design$study_end
      )
      old_cohort <- old_engine$apply_independent_censoring(
        old_uncensored, design$censor_rate, design$study_end
      )
      set.seed(seed)
      new_uncensored <- generate_left_truncated_cohort(
        design$n, truth_beta, lambda0, entry_mean, design$entry_sd,
        design$study_end, gamma = 0, theta = 0, delta = 0
      )
      new_cohort <- apply_independent_censoring(
        new_uncensored, design$censor_rate, design$study_end
      )
      exact_cohort <- exact_cohort && identical(old_uncensored, new_uncensored)
      if (!exact_cohort)
        stop("Regime A cohort mismatch at cell ", design$cell_id,
             ", replicate ", replicate_id, ".", call. = FALSE)

      exact_censored <- exact_censored && identical(old_cohort, new_cohort)
      if (!exact_censored)
        stop("Regime A censoring mismatch at cell ", design$cell_id,
             ", replicate ", replicate_id, ".", call. = FALSE)

      landmarks <- 50 + c(0, design$landmark_interval,
                           2 * design$landmark_interval)
      for (method in methods) {
        old_stack <- old_engine$make_simulation_stack(
          old_cohort, landmarks, prediction_window, method
        )
        new_stack <- make_simulation_stack(
          new_cohort, landmarks, prediction_window, method
        )
        exact_stack <- exact_stack && identical(old_stack, new_stack)
        if (!exact_stack)
          stop("Regime A stack mismatch at cell ", design$cell_id,
               ", replicate ", replicate_id, ", method ", method, ".",
               call. = FALSE)
        old_fit <- old_engine$fit_fast_cluster_cox(old_stack)
        new_fit <- fit_fast_cluster_cox(new_stack)
        exact_fit <- exact_fit && identical(old_fit, new_fit)
        if (!exact_fit)
          stop("Regime A fit mismatch at cell ", design$cell_id,
               ", replicate ", replicate_id, ", method ", method, ".",
               call. = FALSE)
      }
    }
    data.table(
      cell_id = design$cell_id, replicates_checked = nrow(seeds),
      methods_checked = length(methods), exact_cohort = as.integer(exact_cohort),
      exact_censored_cohort = as.integer(exact_censored),
      exact_stack = as.integer(exact_stack), exact_fit = as.integer(exact_fit)
    )
  }

  started <- proc.time()[["elapsed"]]
  checks <- mclapply(
    seq_len(nrow(phase4_grid)), check_cell, mc.cores = cores,
    mc.preschedule = TRUE, mc.set.seed = FALSE
  )
  runtime <- proc.time()[["elapsed"]] - started
  if (any(vapply(checks, inherits, logical(1L), "try-error")))
    stop("At least one Regime A exact-regression worker failed.")
  checks <- rbindlist(checks)
  stopifnot(nrow(checks) == 108L,
            all(checks$replicates_checked == regression_replicates),
            all(checks$methods_checked == length(methods)),
            all(checks$exact_cohort == 1L),
            all(checks$exact_censored_cohort == 1L),
            all(checks$exact_stack == 1L), all(checks$exact_fit == 1L))
  checks[, runtime_seconds := runtime]
  checks
}

make_focused_design <- function() {
  design <- CJ(gamma = gammas, entry_sd = entry_sds)
  design[, `:=`(
    phase = "focused", regime = fifelse(gamma == 0, "A", "B"),
    theta = 0, delta = 0, landmark_interval = 4,
    censoring_target = 0.30, n = 750L
  )]
  setcolorder(design, c("phase", "regime", "gamma", "theta", "delta",
                        "landmark_interval", "censoring_target", "n", "entry_sd"))
  design[, cell_id := .I]
  design
}

make_confirmation_design <- function(selected_gamma) {
  stopifnot(length(selected_gamma) == 1L, selected_gamma %in% gammas)
  design <- as.data.table(expand.grid(
    landmark_interval = landmark_intervals,
    censoring_target = censoring_targets,
    n = sample_sizes, entry_sd = entry_sds,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ))
  design[, `:=`(phase = "confirmation", regime = "B",
                gamma = selected_gamma, theta = 0, delta = 0)]
  setcolorder(design, c("phase", "regime", "gamma", "theta", "delta",
                        "landmark_interval", "censoring_target", "n", "entry_sd"))
  design[, cell_id := .I]
  stopifnot(nrow(design) == 108L)
  design
}

make_frailty_design <- function() {
  design <- CJ(theta = c(0, 0.5), delta = c(0, -2), entry_sd = entry_sds)
  design[, `:=`(
    phase = "frailty", regime = "C", gamma = 0,
    landmark_interval = 4, censoring_target = 0.30, n = 750L
  )]
  setcolorder(design, c("phase", "regime", "gamma", "theta", "delta",
                        "landmark_interval", "censoring_target", "n", "entry_sd"))
  design[, cell_id := .I]
  stopifnot(nrow(design) == 12L)
  design
}

calibrate_design <- function(design, seed_base) {
  keys <- c("phase", "gamma", "theta", "delta", "landmark_interval",
            "censoring_target", "entry_sd")
  calibration <- unique(design[, ..keys])
  setorderv(calibration, keys)
  calibration[, calibration_id := .I]
  rows <- lapply(seq_len(nrow(calibration)), function(i) {
    row <- calibration[i]
    seed <- as.integer(seed_base + i)
    set.seed(seed)
    landmarks <- 50 + c(0, row$landmark_interval,
                        2 * row$landmark_interval)
    study_end <- max(landmarks) + prediction_window
    pilot <- generate_left_truncated_cohort(
      pilot_n, truth_beta, lambda0, entry_mean, row$entry_sd, study_end,
      gamma = row$gamma, theta = row$theta, delta = row$delta
    )
    if (row$censoring_target == 0) {
      rate <- 0
      achieved <- 0
    } else {
      rate <- calibrate_censor_rate(row$censoring_target, pilot, study_end)
      validation_uniform <- stats::runif(nrow(pilot))
      censor_age <- pilot$entry - log(validation_uniform) / rate
      achieved <- mean(censor_age < pilot$event_age & censor_age < study_end)
    }
    stopifnot(is.finite(rate), rate >= 0, is.finite(achieved),
              achieved >= 0, achieved <= 1)
    data.table(
      row, study_end = study_end, censor_rate = rate,
      pilot_achieved_censoring = achieved,
      calibration_difference = achieved - row$censoring_target,
      pilot_n = nrow(pilot), calibration_seed = seed
    )
  })
  out <- rbindlist(rows)
  stopifnot(nrow(out) == nrow(calibration),
            !anyDuplicated(out[, ..keys]))
  out
}

attach_calibration <- function(design, calibration) {
  keys <- c("phase", "gamma", "theta", "delta", "landmark_interval",
            "censoring_target", "entry_sd")
  out <- merge(
    design,
    calibration[, c(keys, "study_end", "censor_rate"), with = FALSE],
    by = keys, all.x = TRUE, sort = FALSE
  )
  setorder(out, cell_id)
  stopifnot(nrow(out) == nrow(design), !anyNA(out),
            all(out$study_end == 50 + 2 * out$landmark_interval +
                  prediction_window))
  out
}

simulate_design <- function(design, seed_base) {
  seed_manifest <- design[, .(replicate = seq_len(replicates)), by = cell_id]
  seed_manifest[, seed := as.integer(seed_base + cell_id * 2000L + replicate)]
  seed_manifest[, phase := unique(design$phase)]
  stopifnot(!anyDuplicated(seed_manifest$seed), all(seed_manifest$seed > 0L))

  run_cell <- function(i) {
    row <- design[i]
    cell_seeds <- seed_manifest[cell_id == row$cell_id]
    landmarks <- 50 + c(0, row$landmark_interval,
                        2 * row$landmark_interval)
    started <- proc.time()[["elapsed"]]
    output <- vector("list", replicates * length(methods))
    cursor <- 1L
    for (replicate_id in seq_len(replicates)) {
      replicate_seed <- cell_seeds$seed[replicate_id]
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
      stopifnot(is.finite(entry_x_correlation),
                all(is.finite(preentry_fractions)),
                all(preentry_fractions >= 0 & preentry_fractions <= 1))
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
            structure(conditionMessage(error), class = "phase6_fit_error")
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
        if (inherits(fitted, "phase6_fit_error")) {
          output[[cursor]] <- cbind(common, data.table(
            estimate = NA_real_, estimated_se = NA_real_,
            covered = NA_integer_, reject_null = NA_integer_,
            rows = nrow(stacks[[method]]), events = sum(stacks[[method]]$event),
            subjects = uniqueN(stacks[[method]]$id), converged = 0L,
            error = as.character(fitted)
          ))
        } else {
          z <- stats::qnorm(1 - (1 - nominal_coverage) / 2)
          output[[cursor]] <- cbind(common, data.table(
            estimate = fitted$estimate, estimated_se = fitted$robust_se,
            covered = as.integer(abs(fitted$estimate - truth_beta) <=
                                   z * fitted$robust_se),
            reject_null = as.integer(abs(fitted$estimate / fitted$robust_se) > z),
            rows = fitted$rows, events = fitted$events,
            subjects = fitted$subjects, converged = 1L, error = ""
          ))
        }
        cursor <- cursor + 1L
      }
    }
    result <- rbindlist(output)
    stopifnot(nrow(result) == replicates * length(methods),
              all(result$cell_id == row$cell_id),
              all(table(result$method) == replicates))
    list(result = result,
         runtime_seconds = proc.time()[["elapsed"]] - started)
  }

  started <- proc.time()[["elapsed"]]
  runs <- mclapply(
    seq_len(nrow(design)), run_cell, mc.cores = cores,
    mc.preschedule = TRUE, mc.set.seed = FALSE
  )
  wall_runtime <- proc.time()[["elapsed"]] - started
  if (any(vapply(runs, inherits, logical(1L), "try-error")))
    stop("At least one informative-entry simulation worker failed.")
  raw <- rbindlist(lapply(runs, `[[`, "result"))
  runtime <- data.table(
    phase = unique(design$phase), cell_id = design$cell_id,
    runtime_seconds = vapply(runs, `[[`, numeric(1L), "runtime_seconds")
  )
  stopifnot(nrow(raw) == nrow(design) * replicates * length(methods),
            raw[, all(.N == replicates), by = .(cell_id, method)]$V1)
  list(raw = raw, seeds = seed_manifest, runtime = runtime,
       wall_runtime = wall_runtime)
}

summarize_simulation <- function(raw) {
  successful <- raw[converged == 1L]
  grouping <- c("phase", "regime", "cell_id", "method", "gamma", "theta",
                "delta", "landmark_interval", "censoring_target", "n",
                "entry_sd")
  summary <- successful[, .(
    successful_replicates = .N,
    bias = mean(estimate - truth_beta),
    empirical_se = stats::sd(estimate),
    mean_estimated_se = mean(estimated_se),
    se_ratio = mean(estimated_se) / stats::sd(estimate),
    coverage = mean(covered),
    coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / .N),
    power = mean(reject_null),
    power_mcse = sqrt(mean(reject_null) * (1 - mean(reject_null)) / .N),
    mean_rows = mean(rows), mean_events = mean(events),
    mean_subjects = mean(subjects),
    mean_entry_x_correlation = mean(entry_x_correlation),
    mean_observed_censoring = mean(observed_censoring),
    mean_preentry_fraction_lm1 = mean(preentry_fraction_lm1),
    mean_preentry_fraction_lm2 = mean(preentry_fraction_lm2),
    mean_preentry_fraction_lm3 = mean(preentry_fraction_lm3)
  ), by = grouping]
  failures <- raw[, .(failed_replicates = sum(converged == 0L)), by = grouping]
  summary <- merge(summary, failures, by = grouping, all.x = TRUE)
  setorderv(summary, c("phase", "cell_id", "method"))
  numeric_fields <- setdiff(
    names(summary), c(grouping, "successful_replicates", "failed_replicates")
  )
  stopifnot(all(summary$successful_replicates + summary$failed_replicates ==
                  replicates),
            all(is.finite(as.matrix(summary[, ..numeric_fields]))))
  summary
}

regression <- phase4_exact_regression()
fwrite(regression, project_path("results", "regime_a_regression_check.csv"))

focused_design <- make_focused_design()
focused_calibration <- calibrate_design(focused_design, 1700000000L)
focused_design <- attach_calibration(focused_design, focused_calibration)
focused <- simulate_design(focused_design, 1200000000L)
focused_summary <- summarize_simulation(focused$raw)

separation <- dcast(
  focused_summary, gamma + entry_sd ~ method, value.var = "bias"
)
stopifnot(all(c("naive", "strict", "delayed") %in% names(separation)))
separation[, absolute_naive_strict_separation := abs(naive - strict)]
gamma_separation <- separation[, .(
  mean_absolute_naive_strict_separation =
    mean(absolute_naive_strict_separation),
  maximum_absolute_naive_strict_separation =
    max(absolute_naive_strict_separation)
), by = gamma]
setorder(gamma_separation, -mean_absolute_naive_strict_separation, gamma)
selected_gamma <- gamma_separation$gamma[1L]
stopifnot(length(selected_gamma) == 1L, selected_gamma %in% gammas)

if (skip_confirmation) {
  confirmation <- list(raw = focused$raw[0], seeds = focused$seeds[0],
                       runtime = focused$runtime[0], wall_runtime = 0)
  confirmation_calibration <- focused_calibration[0]
  confirmation_summary <- focused_summary[0]
} else {
  confirmation_design <- make_confirmation_design(selected_gamma)
  confirmation_calibration <- calibrate_design(
    confirmation_design, 1710000000L
  )
  confirmation_design <- attach_calibration(
    confirmation_design, confirmation_calibration
  )
  confirmation <- simulate_design(confirmation_design, 1300000000L)
  confirmation_summary <- summarize_simulation(confirmation$raw)
}

frailty_design <- make_frailty_design()
frailty_calibration <- calibrate_design(frailty_design, 1720000000L)
frailty_design <- attach_calibration(frailty_design, frailty_calibration)
frailty <- simulate_design(frailty_design, 1600000000L)
frailty_summary <- summarize_simulation(frailty$raw)

raw <- rbindlist(list(focused$raw, confirmation$raw, frailty$raw), fill = TRUE)
summary <- rbindlist(
  list(focused_summary, confirmation_summary, frailty_summary), fill = TRUE
)
seeds <- rbindlist(
  list(focused$seeds, confirmation$seeds, frailty$seeds), fill = TRUE
)
calibration <- rbindlist(
  list(focused_calibration, confirmation_calibration, frailty_calibration),
  fill = TRUE
)
cell_runtime <- rbindlist(
  list(focused$runtime, confirmation$runtime, frailty$runtime), fill = TRUE
)
failures <- raw[converged == 0L,
                .(failure_count = .N), by = .(phase, cell_id, method, error)]
stopifnot(!anyDuplicated(seeds[, .(phase, cell_id, replicate)]),
          nrow(summary) == uniqueN(raw[, .(phase, cell_id, method)]),
          all(calibration$pilot_n == pilot_n))

fwrite(raw, project_path("results", "informative_entry_raw.csv"))
fwrite(summary, project_path("results", "informative_entry_summary.csv"))
fwrite(seeds, project_path("results", "informative_entry_seed_manifest.csv"))
fwrite(calibration,
       project_path("results", "informative_entry_calibration.csv"))
fwrite(cell_runtime,
       project_path("results", "informative_entry_cell_runtime.csv"))
fwrite(failures, project_path("results", "informative_entry_failures.csv"))
fwrite(gamma_separation,
       project_path("results", "informative_entry_gamma_selection.csv"))

total_wall_runtime <- regression$runtime_seconds[1L] +
  focused$wall_runtime + confirmation$wall_runtime + frailty$wall_runtime
trace_values(
  report = "07_informative_entry",
  labels = c(
    "analysis_seed", "truth_beta", "hazard_ratio", "lambda0", "entry_mean",
    "prediction_window", "nominal_coverage", "replicates_per_cell",
    "regression_replicates_per_cell", "regression_cells",
    "regression_exact_cohort_cells", "regression_exact_fit_cells",
    "regression_runtime_seconds", "focused_cells", "confirmation_cells",
    "frailty_cells", "methods", "total_fits", "fit_failures",
    "calibration_pilot_n", "selected_gamma",
    "selected_gamma_mean_absolute_separation", "focused_runtime_seconds",
    "confirmation_runtime_seconds", "frailty_runtime_seconds",
    "total_runtime_seconds", "parallel_cores", "phase6_tolerance_changes"
  ),
  values = c(
    analysis_seed, truth_beta, exp(truth_beta), lambda0, entry_mean,
    prediction_window, nominal_coverage, replicates, regression_replicates,
    nrow(regression), sum(regression$exact_cohort), sum(regression$exact_fit),
    regression$runtime_seconds[1L], nrow(focused_design),
    if (skip_confirmation) 0 else nrow(confirmation_design),
    nrow(frailty_design), length(methods), nrow(raw), nrow(raw) -
      sum(raw$converged), pilot_n, selected_gamma,
    gamma_separation[gamma == selected_gamma,
                     mean_absolute_naive_strict_separation],
    focused$wall_runtime, confirmation$wall_runtime, frailty$wall_runtime,
    total_wall_runtime, cores, sum(tolerance_changes$phase6_change)
  ),
  units = c(
    "seed", "log hazard ratio", "hazard ratio", "hazard per year", "years",
    "years", "proportion", "replicates", "replicates", "cells", "cells",
    "cells", "seconds", "cells", "cells", "cells", "methods", "fits",
    "fits", "subjects", "years per x unit", "log hazard ratio", "seconds",
    "seconds", "seconds", "seconds", "cores", "changes"
  ),
  script = "analysis/07_informative_entry.R",
  function_name = "phase4_exact_regression/simulate_design",
  seed = "results/informative_entry_seed_manifest.csv and Phase 4 manifest",
  runtime_seconds = total_wall_runtime
)

for (i in seq_len(nrow(tolerance_changes))) {
  trace_values(
    report = "07_informative_entry",
    labels = paste0("tolerance_", i, c("_original", "_new")),
    values = c(tolerance_changes$original_value[i],
               tolerance_changes$new_value[i]),
    units = "tolerance", script = "analysis/07_informative_entry.R",
    function_name = "tolerance audit", seed = NA_integer_,
    runtime_seconds = total_wall_runtime
  )
}

summary_fields <- c(
  "successful_replicates", "failed_replicates", "bias", "empirical_se",
  "mean_estimated_se", "se_ratio", "coverage", "coverage_mcse", "power",
  "power_mcse", "mean_rows", "mean_events", "mean_subjects",
  "mean_entry_x_correlation", "mean_observed_censoring",
  "mean_preentry_fraction_lm1", "mean_preentry_fraction_lm2",
  "mean_preentry_fraction_lm3"
)
for (i in seq_len(nrow(summary))) {
  row <- summary[i]
  prefix <- paste0(
    row$phase, "_cell_", sprintf("%03d", row$cell_id), "_", row$method, "_"
  )
  phase_seeds <- seeds[phase == row$phase & cell_id == row$cell_id]$seed
  runtime <- cell_runtime[
    phase == row$phase & cell_id == row$cell_id, runtime_seconds
  ]
  stopifnot(length(phase_seeds) == replicates, length(runtime) == 1L)
  trace_values(
    report = "07_informative_entry",
    labels = paste0(prefix, summary_fields),
    values = as.numeric(row[, ..summary_fields]),
    units = summary_fields,
    script = "analysis/07_informative_entry.R",
    function_name = "fit_fast_cluster_cox/cell aggregation",
    seed = paste0(min(phase_seeds), "-", max(phase_seeds)),
    runtime_seconds = runtime
  )
}

stopifnot(
  file.exists(project_path("results", "regime_a_regression_check.csv")),
  file.exists(project_path("results", "informative_entry_summary.csv")),
  file.exists(project_path("results", "informative_entry_seed_manifest.csv")),
  file.exists(project_path("results", "tolerance_changes.csv"))
)

write_phase6_pooled_summaries <- function(raw, runtime_seconds,
                                          generating_beta) {
  stopifnot(is.data.table(raw), is.finite(runtime_seconds),
            length(generating_beta) == 1L, is.finite(generating_beta),
            all(c("phase", "method", "estimate", "converged",
                  "covered", "reject_null", "theta", "delta") %in% names(raw)))
  pooled <- rbindlist(list(
    raw[phase == "confirmation" & converged == 1L, .(
      replicates = .N,
      bias = mean(estimate - generating_beta),
      coverage = mean(covered),
      coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / .N),
      power = mean(reject_null),
      power_mcse = sqrt(mean(reject_null) * (1 - mean(reject_null)) / .N)
    ), by = .(phase, method)][, `:=`(theta = NA_real_, delta = NA_real_)],
    raw[phase == "frailty" & converged == 1L, .(
      replicates = .N,
      bias = mean(estimate - generating_beta),
      coverage = mean(covered),
      coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / .N),
      power = mean(reject_null),
      power_mcse = sqrt(mean(reject_null) * (1 - mean(reject_null)) / .N)
    ), by = .(phase, theta, delta, method)]
  ), fill = TRUE)
  stopifnot(nrow(pooled) == 15L,
            pooled[phase == "confirmation", all(replicates == 108000L)],
            pooled[phase == "frailty", all(replicates == 3000L)],
            all(is.finite(pooled$bias)), all(is.finite(pooled$coverage)))
  fwrite(pooled,
         project_path("results", "informative_entry_pooled_summary.csv"))

  fields <- c("replicates", "bias", "coverage", "coverage_mcse",
              "power", "power_mcse")
  for (i in seq_len(nrow(pooled))) {
    row <- pooled[i]
    stratum <- if (row$phase == "confirmation") {
      paste0("pooled_confirmation_", row$method, "_")
    } else {
      paste0("pooled_frailty_theta_", gsub("\\.", "p", row$theta),
             "_delta_", gsub("-", "m", row$delta), "_", row$method, "_")
    }
    trace_values(
      report = "07_informative_entry",
      labels = paste0(stratum, fields),
      values = as.numeric(row[, ..fields]), units = fields,
      script = "analysis/07_informative_entry.R",
      function_name = "write_phase6_pooled_summaries",
      seed = "results/informative_entry_seed_manifest.csv",
      runtime_seconds = runtime_seconds
    )
  }
  invisible(pooled)
}

write_phase6_pooled_summaries(raw, total_wall_runtime, truth_beta)

trace_phase6_design_axes <- function(runtime_seconds) {
  axis_labels <- c(
    paste0("focused_gamma_", seq_along(gammas)),
    paste0("entry_sd_", seq_along(entry_sds)),
    paste0("frailty_theta_", seq_along(c(0, 0.5))),
    paste0("frailty_delta_", seq_along(c(0, -2))),
    "anchor_landmark_interval", "anchor_censoring_target", "anchor_sample_size"
  )
  axis_values <- c(gammas, entry_sds, c(0, 0.5), c(0, -2), 4, 0.30, 750)
  stopifnot(length(axis_labels) == length(axis_values),
            all(is.finite(axis_values)))
  trace_values(
    report = "07_informative_entry", labels = axis_labels,
    values = axis_values, units = c(
      rep("years per x unit", length(gammas)),
      rep("years", length(entry_sds)), rep("log hazard coefficient", 2L),
      rep("years per frailty unit", 2L), "years", "proportion", "subjects"
    ),
    script = "analysis/07_informative_entry.R",
    function_name = "trace_phase6_design_axes", seed = NA_integer_,
    runtime_seconds = runtime_seconds
  )
}

trace_phase6_design_axes(total_wall_runtime)
