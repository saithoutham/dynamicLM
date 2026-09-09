#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})
source(here::here("analysis", "_helpers.R"))

analysis_seed <- NA_integer_
truth_beta <- log(1.5)
confidence_level <- 0.95
z_value <- stats::qnorm(1 - (1 - confidence_level) / 2)
methods <- c("naive", "strict", "delayed")

started <- proc.time()[["elapsed"]]
raw <- fread(project_path("results", "informative_entry_raw.csv"))
assert_required_columns(
  raw,
  c("phase", "cell_id", "replicate", "seed", "method", "theta", "delta",
    "entry_sd", "estimate", "converged", "error")
)
frailty <- raw[phase == "frailty"]
stopifnot(
  nrow(frailty) == 36000L,
  all(frailty$converged == 1L),
  all(is.na(frailty$error) | frailty$error == ""),
  setequal(unique(frailty$method), methods),
  setequal(unique(frailty$theta), c(0, 0.5)),
  setequal(unique(frailty$delta), c(0, -2)),
  setequal(unique(frailty$entry_sd), c(0.25, 4, 10)),
  !anyDuplicated(frailty[, .(cell_id, replicate, method)])
)

arm_summary <- frailty[, .(
  arm_replicates = .N,
  mean_estimate = mean(estimate),
  bias = mean(estimate - truth_beta),
  estimate_variance = stats::var(estimate),
  arm_mcse = stats::sd(estimate) / sqrt(.N),
  minimum_seed = min(seed),
  maximum_seed = max(seed)
), by = .(method, theta, delta)]
stopifnot(
  nrow(arm_summary) == length(methods) * 4L,
  all(arm_summary$arm_replicates == 3000L),
  all(is.finite(as.matrix(arm_summary[, .(
    mean_estimate, bias, estimate_variance, arm_mcse
  )]))),
  all(arm_summary$estimate_variance > 0),
  all(arm_summary$arm_mcse > 0)
)

# Arms use different cell IDs and disjoint seed ranges. Methods within an arm
# share cohort seeds, but every contrast below is computed separately by method
# across independent theta/delta arms.
arm_seeds <- unique(frailty[, .(theta, delta, seed)])
stopifnot(nrow(arm_seeds) == 12000L)
seed_overlap <- arm_seeds[, .N, by = seed][N > 1L]
stopifnot(nrow(seed_overlap) == 0L)

contrast_specification <- rbindlist(list(
  data.table(
    contrast = "delta_at_theta_0",
    theta = c(0, 0), delta = c(-2, 0), weight = c(1, -1)
  ),
  data.table(
    contrast = "delta_at_theta_0p5",
    theta = c(0.5, 0.5), delta = c(-2, 0), weight = c(1, -1)
  ),
  data.table(
    contrast = "theta_at_delta_0",
    theta = c(0.5, 0), delta = c(0, 0), weight = c(1, -1)
  ),
  data.table(
    contrast = "theta_at_delta_m2",
    theta = c(0.5, 0), delta = c(-2, -2), weight = c(1, -1)
  ),
  data.table(
    contrast = "theta_by_delta_interaction",
    theta = c(0.5, 0.5, 0, 0), delta = c(-2, 0, -2, 0),
    weight = c(1, -1, -1, 1)
  )
))
stopifnot(
  nrow(contrast_specification) == 12L,
  contrast_specification[, all(sum(weight) == 0), by = contrast]$V1,
  uniqueN(contrast_specification$contrast) == 5L
)

compute_contrast <- function(selected_method, contrast_name) {
  specification <- contrast_specification[contrast == contrast_name]
  arms <- merge(
    specification,
    arm_summary[method == selected_method],
    by = c("theta", "delta"), all.x = TRUE, sort = FALSE
  )
  stopifnot(
    nrow(arms) == nrow(specification),
    !anyNA(arms),
    all(arms$method == selected_method),
    all(arms$arm_replicates == 3000L)
  )
  estimate <- sum(arms$weight * arms$bias)
  monte_carlo_variance <- sum(
    arms$weight^2 * arms$estimate_variance / arms$arm_replicates
  )
  monte_carlo_se <- sqrt(monte_carlo_variance)
  stopifnot(is.finite(estimate), is.finite(monte_carlo_se),
            monte_carlo_se > 0)
  data.table(
    method = selected_method,
    contrast = contrast_name,
    contrast_estimate = estimate,
    monte_carlo_se = monte_carlo_se,
    ci_level = confidence_level,
    ci_lower = estimate - z_value * monte_carlo_se,
    ci_upper = estimate + z_value * monte_carlo_se,
    interval_includes_zero = as.integer(
      estimate - z_value * monte_carlo_se <= 0 &
        estimate + z_value * monte_carlo_se >= 0
    ),
    arms_independent = 1L,
    arm_replicates = unique(arms$arm_replicates),
    seed_overlap_count = nrow(seed_overlap)
  )
}

decomposition <- rbindlist(lapply(methods, function(method) {
  rbindlist(lapply(
    unique(contrast_specification$contrast),
    function(contrast_name) compute_contrast(method, contrast_name)
  ))
}))
setorder(decomposition, method, contrast)
stopifnot(
  nrow(decomposition) == length(methods) * 5L,
  !anyNA(decomposition),
  all(decomposition$ci_lower < decomposition$contrast_estimate),
  all(decomposition$ci_upper > decomposition$contrast_estimate),
  decomposition[grepl("^delta_|interaction$", contrast),
                all(interval_includes_zero == 1L)],
  decomposition[grepl("^theta_at_", contrast),
                all(interval_includes_zero == 0L)],
  all(decomposition$arms_independent == 1L),
  all(decomposition$seed_overlap_count == 0L)
)

fwrite(decomposition,
       project_path("results", "frailty_decomposition.csv"))
runtime_seconds <- proc.time()[["elapsed"]] - started

trace_values(
  report = "08_frailty_decomposition",
  labels = c(
    "truth_beta", "confidence_level", "raw_frailty_rows", "methods",
    "arm_replicates", "contrasts_per_method", "seed_overlap_count",
    "runtime_seconds"
  ),
  values = c(
    truth_beta, confidence_level, nrow(frailty), length(methods),
    unique(decomposition$arm_replicates),
    uniqueN(decomposition$contrast), nrow(seed_overlap), runtime_seconds
  ),
  units = c(
    "log hazard ratio", "proportion", "rows", "methods", "replicates",
    "contrasts", "seeds", "seconds"
  ),
  script = "analysis/08_frailty_decomposition.R",
  function_name = "replicate-level independent-arm decomposition",
  seed = analysis_seed,
  runtime_seconds = runtime_seconds
)

trace_fields <- c(
  "contrast_estimate", "monte_carlo_se", "ci_level", "ci_lower", "ci_upper",
  "interval_includes_zero", "arms_independent", "arm_replicates",
  "seed_overlap_count"
)
for (i in seq_len(nrow(decomposition))) {
  row <- decomposition[i]
  prefix <- paste0(row$method, "_", row$contrast, "_")
  trace_values(
    report = "08_frailty_decomposition",
    labels = paste0(prefix, trace_fields),
    values = as.numeric(row[, ..trace_fields]),
    units = trace_fields,
    script = "analysis/08_frailty_decomposition.R",
    function_name = "compute_contrast",
    seed = "independent arms; results/informative_entry_seed_manifest.csv",
    runtime_seconds = runtime_seconds
  )
}

stopifnot(
  file.exists(project_path("results", "frailty_decomposition.csv")),
  file.info(project_path("results", "frailty_decomposition.csv"))$size > 0L
)
