generate_entry_age <- function(n, mean_age, sd_age, lower, upper) {
  stopifnot(n > 0L, sd_age > 0, lower < upper)
  out <- numeric()
  while (length(out) < n) {
    draw <- stats::rnorm(max(100L, 2L * (n - length(out))), mean_age, sd_age)
    out <- c(out, draw[draw >= lower & draw < upper])
  }
  out[seq_len(n)]
}


generate_left_truncated_cohort <- function(n, beta, lambda0, entry_mean,
                                           entry_sd, study_end) {
  stopifnot(n > 0L, is.finite(beta), lambda0 > 0, entry_sd > 0,
            study_end > entry_mean)
  accepted <- list()
  accepted_n <- 0L
  while (accepted_n < n) {
    batch_n <- max(500L, 2L * (n - accepted_n))
    entry <- generate_entry_age(batch_n, entry_mean, entry_sd,
                                lower = 25, upper = study_end - 1e-6)
    x <- stats::rnorm(batch_n)
    event_age <- stats::rexp(batch_n, rate = lambda0 * exp(beta * x))
    keep <- event_age > entry
    if (any(keep)) {
      accepted[[length(accepted) + 1L]] <- data.frame(
        entry = entry[keep], x = x[keep], event_age = event_age[keep]
      )
      accepted_n <- accepted_n + sum(keep)
    }
  }
  out <- do.call(rbind, accepted)
  out <- out[seq_len(n), , drop = FALSE]
  out$id <- seq_len(n)
  out <- out[, c("id", "entry", "x", "event_age")]
  stopifnot(nrow(out) == n, !anyNA(out), all(out$event_age > out$entry),
            all(out$entry < study_end), !anyDuplicated(out$id))
  out
}


apply_independent_censoring <- function(data, censor_rate, study_end) {
  stopifnot(is.data.frame(data), censor_rate >= 0, is.finite(censor_rate),
            all(data$entry < study_end), all(data$event_age > data$entry))
  censor_age <- if (censor_rate == 0) {
    rep(Inf, nrow(data))
  } else {
    data$entry + stats::rexp(nrow(data), rate = censor_rate)
  }
  data$exit <- pmin(data$event_age, censor_age, study_end)
  data$status <- as.integer(data$event_age <= censor_age &
                              data$event_age <= study_end)
  data$censored_before_end <- as.integer(censor_age < data$event_age &
                                           censor_age < study_end)
  stopifnot(all(data$exit > data$entry), all(data$status %in% 0:1),
            all(data$censored_before_end %in% 0:1),
            all(data$status + data$censored_before_end <= 1L))
  data
}


make_simulation_stack <- function(data, landmarks, w,
                                  method = c("naive", "strict", "delayed")) {
  method <- match.arg(method)
  rows <- lapply(landmarks, function(landmark) {
    horizon <- landmark + w
    if (method == "naive") {
      keep <- data$exit > landmark
      start <- rep(landmark, sum(keep))
    } else if (method == "strict") {
      keep <- data$entry <= landmark & data$exit > landmark
      start <- rep(landmark, sum(keep))
    } else {
      keep <- data$entry < horizon & data$exit > pmax(landmark, data$entry)
      start <- pmax(landmark, data$entry[keep])
    }
    out <- data[keep, c("id", "x"), drop = FALSE]
    out$landmark <- landmark
    out$start <- start
    out$stop <- pmin(data$exit[keep], horizon)
    out$event <- as.integer(data$status[keep] == 1L & data$exit[keep] <= horizon)
    stopifnot(nrow(out) == length(start), all(out$start < out$stop),
              all(out$event %in% 0:1))
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  stopifnot(nrow(out) > 0L, !anyNA(out),
            !anyDuplicated(out[, c("id", "landmark")]))
  out
}


cox_risk_sums <- function(beta, start, stop, x, event_times) {
  eta <- beta * x
  if (any(!is.finite(eta)) || max(eta) > 700)
    stop("Non-finite or overflowing linear predictor.")
  weight <- exp(eta)
  starts <- order(start)
  stops <- order(stop)
  cumulative <- function(values, ordering, counts) {
    sums <- cumsum(values[ordering])
    ifelse(counts == 0L, 0, sums[pmax(counts, 1L)])
  }
  start_count <- findInterval(event_times, start[starts], left.open = TRUE)
  stop_count <- findInterval(event_times, stop[stops], left.open = TRUE)
  s0 <- cumulative(weight, starts, start_count) -
    cumulative(weight, stops, stop_count)
  s1 <- cumulative(weight * x, starts, start_count) -
    cumulative(weight * x, stops, stop_count)
  s2 <- cumulative(weight * x^2, starts, start_count) -
    cumulative(weight * x^2, stops, stop_count)
  if (any(!is.finite(c(s0, s1, s2))) || any(s0 <= 0))
    stop("Invalid Cox risk-set sums.")
  list(s0 = s0, s1 = s1, s2 = s2, weight = weight)
}


fit_fast_cluster_cox <- function(data, tolerance = 1e-9, max_iterations = 40L) {
  stopifnot(all(c("id", "start", "stop", "event", "x") %in% names(data)),
            all(data$start < data$stop), all(data$event %in% 0:1),
            all(is.finite(data$x)))
  event_rows <- which(data$event == 1L)
  if (length(event_rows) < 2L)
    stop("Fewer than two event rows.")
  event_times <- sort(unique(data$stop[event_rows]))
  event_count <- tabulate(match(data$stop[event_rows], event_times),
                          nbins = length(event_times))
  event_x <- as.numeric(rowsum(data$x[event_rows],
                               match(data$stop[event_rows], event_times),
                               reorder = TRUE))
  stopifnot(length(event_x) == length(event_times), sum(event_count) == length(event_rows))

  beta <- 0
  converged <- FALSE
  for (iteration in seq_len(max_iterations)) {
    sums <- cox_risk_sums(beta, data$start, data$stop, data$x, event_times)
    score <- sum(event_x - event_count * sums$s1 / sums$s0)
    information <- sum(event_count *
                         (sums$s2 / sums$s0 - (sums$s1 / sums$s0)^2))
    if (!is.finite(score) || !is.finite(information) || information <= 0)
      stop("Non-positive or non-finite Cox information.")
    step <- score / information
    beta <- beta + step
    if (!is.finite(beta)) stop("Non-finite Cox coefficient.")
    if (abs(step) < tolerance) {
      converged <- TRUE
      break
    }
  }
  if (!converged) stop("Fast Cox Newton iteration did not converge.")

  sums <- cox_risk_sums(beta, data$start, data$stop, data$x, event_times)
  information <- sum(event_count *
                       (sums$s2 / sums$s0 - (sums$s1 / sums$s0)^2))
  jump <- event_count / sums$s0
  jump_xbar <- jump * sums$s1 / sums$s0
  cumulative_jump <- cumsum(jump)
  cumulative_jump_xbar <- cumsum(jump_xbar)
  jump_through <- function(value, cumulative) {
    index <- findInterval(value, event_times)
    ifelse(index == 0L, 0, cumulative[pmax(index, 1L)])
  }
  integrated_jump <- jump_through(data$stop, cumulative_jump) -
    jump_through(data$start, cumulative_jump)
  integrated_jump_xbar <- jump_through(data$stop, cumulative_jump_xbar) -
    jump_through(data$start, cumulative_jump_xbar)
  event_xbar <- numeric(nrow(data))
  event_xbar[event_rows] <-
    (sums$s1 / sums$s0)[match(data$stop[event_rows], event_times)]
  row_score <- data$event * (data$x - event_xbar) -
    sums$weight * (data$x * integrated_jump - integrated_jump_xbar)
  cluster_score <- as.numeric(rowsum(row_score, data$id, reorder = FALSE))
  stopifnot(is.finite(sum(row_score)), abs(sum(row_score)) < 1e-5,
            all(is.finite(cluster_score)))
  robust_variance <- sum(cluster_score^2) / information^2
  se <- sqrt(robust_variance)
  if (!is.finite(se) || se <= 0) stop("Invalid cluster-robust standard error.")
  list(
    estimate = beta,
    robust_se = se,
    events = length(event_rows),
    rows = nrow(data),
    subjects = length(unique(data$id)),
    iterations = iteration,
    score = sum(row_score)
  )
}


calibrate_censor_rate <- function(target, pilot, study_end) {
  stopifnot(target >= 0, target < 1, nrow(pilot) > 0L)
  if (target == 0) return(0)
  u <- stats::runif(nrow(pilot))
  censor_fraction <- function(log_rate) {
    rate <- exp(log_rate)
    censor_age <- pilot$entry - log(u) / rate
    mean(censor_age < pilot$event_age & censor_age < study_end)
  }
  root <- stats::uniroot(function(log_rate) censor_fraction(log_rate) - target,
                         interval = c(log(1e-6), log(10)), tol = 1e-10)
  rate <- exp(root$root)
  stopifnot(is.finite(rate), rate > 0,
            abs(censor_fraction(log(rate)) - target) < 1e-5)
  rate
}


numerical_true_auc <- function(landmark, w, beta, lambda0,
                               grid_step = 0.0025, limit = 6) {
  x <- seq(-limit, limit, by = grid_step)
  source_mass <- stats::dnorm(x) * grid_step
  survival_to_landmark <- exp(-lambda0 * landmark * exp(beta * x))
  interval_survival <- exp(-lambda0 * w * exp(beta * x))
  case_mass <- source_mass * survival_to_landmark * (1 - interval_survival)
  control_mass <- source_mass * survival_to_landmark * interval_survival
  control_below <- c(0, utils::head(cumsum(control_mass), -1L))
  numerator <- sum(case_mass * (control_below + 0.5 * control_mass))
  auc <- numerator / (sum(case_mass) * sum(control_mass))
  stopifnot(is.finite(auc), auc >= 0, auc <= 1)
  auc
}


lt_censoring_weights_fast <- function(start, stop, event, horizon) {
  stopifnot(length(start) == length(stop), length(stop) == length(event),
            all(start < stop), all(stop <= horizon), all(event %in% 0:1))
  censor_rows <- which(event == 0L & stop < horizon)
  if (length(censor_rows) == 0L) {
    return(list(exit = rep(1, length(start)), horizon = rep(1, length(start))))
  }
  censor_times <- sort(unique(stop[censor_rows]))
  censor_count <- tabulate(match(stop[censor_rows], censor_times),
                           nbins = length(censor_times))
  ordered_start <- sort(start)
  ordered_stop <- sort(stop)
  risk_count <- findInterval(censor_times, ordered_start, left.open = TRUE) -
    findInterval(censor_times, ordered_stop, left.open = TRUE)
  stopifnot(all(risk_count > 0L), all(censor_count <= risk_count))
  survival_after <- cumprod(1 - censor_count / risk_count)
  if (any(survival_after <= 0) || any(!is.finite(survival_after)))
    stop("Fast censoring product-limit estimate reached zero.")
  step_g <- function(at, before = FALSE) {
    index <- findInterval(at, censor_times, left.open = before)
    ifelse(index == 0L, 1, survival_after[pmax(index, 1L)])
  }
  g_start <- step_g(start)
  g_exit_before <- step_g(stop, before = TRUE)
  g_horizon <- step_g(rep(horizon, length(start)))
  inverse_exit <- g_start / g_exit_before
  inverse_horizon <- g_start / g_horizon
  stopifnot(all(is.finite(c(inverse_exit, inverse_horizon))),
            all(c(inverse_exit, inverse_horizon) >= 1))
  list(exit = inverse_exit, horizon = inverse_horizon)
}


weighted_auc_fast <- function(risk, case_weight, control_weight) {
  stopifnot(length(risk) == length(case_weight),
            length(risk) == length(control_weight),
            all(is.finite(c(risk, case_weight, control_weight))),
            all(case_weight >= 0), all(control_weight >= 0))
  total_case <- sum(case_weight)
  total_control <- sum(control_weight)
  if (total_case <= 0 || total_control <= 0)
    stop("AUC requires cases and controls.")
  ordering <- order(risk)
  sorted_risk <- risk[ordering]
  sorted_case <- case_weight[ordering]
  sorted_control <- control_weight[ordering]
  groups <- match(sorted_risk, unique(sorted_risk))
  group_case <- as.numeric(rowsum(sorted_case, groups, reorder = FALSE))
  group_control <- as.numeric(rowsum(sorted_control, groups, reorder = FALSE))
  control_below <- c(0, utils::head(cumsum(group_control), -1L))
  auc <- sum(group_case * (control_below + 0.5 * group_control)) /
    (total_case * total_control)
  tolerance <- sqrt(.Machine$double.eps)
  stopifnot(is.finite(auc), auc >= -tolerance, auc <= 1 + tolerance)
  pmin(1, pmax(0, auc))
}


summary_auc_fast <- function(stack, risk, w) {
  stopifnot(nrow(stack) == length(risk),
            all(c("landmark", "start", "stop", "event") %in% names(stack)))
  values <- vapply(sort(unique(stack$landmark)), function(landmark) {
    index <- stack$landmark == landmark
    horizon <- landmark + w
    weights <- lt_censoring_weights_fast(
      stack$start[index], stack$stop[index], stack$event[index], horizon
    )
    cases <- stack$event[index] == 1L
    controls <- stack$event[index] == 0L & stack$stop[index] >= horizon
    weighted_auc_fast(
      risk[index], as.numeric(cases) * weights$exit,
      as.numeric(controls) * weights$horizon
    )
  }, numeric(1L))
  stopifnot(all(is.finite(values)), all(values >= 0 & values <= 1))
  mean(values)
}


bootstrap_summary_auc_fast <- function(stack, risk, w, B, seed) {
  stopifnot(B >= 2L, B == as.integer(B), seed == as.integer(seed),
            !anyDuplicated(stack[, c("id", "landmark")]))
  point <- summary_auc_fast(stack, risk, w)
  ids <- unique(stack$id)
  row_index <- split(seq_len(nrow(stack)), stack$id)
  set.seed(as.integer(seed))
  bootstrap <- vapply(seq_len(B), function(b) {
    sampled <- sample(ids, length(ids), replace = TRUE)
    index <- unlist(row_index[as.character(sampled)], use.names = FALSE)
    boot_stack <- stack[index, , drop = FALSE]
    boot_stack$id <- rep(seq_along(sampled),
                         vapply(row_index[as.character(sampled)], length,
                                integer(1L)))
    summary_auc_fast(boot_stack, risk[index], w)
  }, numeric(1L))
  se <- stats::sd(bootstrap)
  stopifnot(is.finite(point), is.finite(se), se > 0)
  list(estimate = point, se = se, bootstrap = bootstrap)
}


paired_bootstrap_auc_delta_fast <- function(stack_a, risk_a, stack_b, risk_b,
                                            source_ids, w, B, seed) {
  stopifnot(nrow(stack_a) == length(risk_a), nrow(stack_b) == length(risk_b),
            B >= 2L, B == as.integer(B), seed == as.integer(seed),
            all(stack_a$id %in% source_ids), all(stack_b$id %in% source_ids))
  point_a <- summary_auc_fast(stack_a, risk_a, w)
  point_b <- summary_auc_fast(stack_b, risk_b, w)
  index_a <- split(seq_len(nrow(stack_a)), stack_a$id)
  index_b <- split(seq_len(nrow(stack_b)), stack_b$id)
  set.seed(as.integer(seed))
  bootstrap <- vapply(seq_len(B), function(iteration) {
    sampled <- sample(source_ids, length(source_ids), replace = TRUE)
    score_boot <- function(stack, risk, row_index) {
      present <- as.character(sampled) %in% names(row_index)
      selected_ids <- sampled[present]
      index <- unlist(row_index[as.character(selected_ids)], use.names = FALSE)
      boot_stack <- stack[index, , drop = FALSE]
      boot_stack$id <- rep(seq_along(selected_ids),
                           vapply(row_index[as.character(selected_ids)], length,
                                  integer(1L)))
      summary_auc_fast(boot_stack, risk[index], w)
    }
    score_boot(stack_a, risk_a, index_a) - score_boot(stack_b, risk_b, index_b)
  }, numeric(1L))
  se <- stats::sd(bootstrap)
  stopifnot(is.finite(se), se >= 0, all(is.finite(bootstrap)))
  list(estimate_a = point_a, estimate_b = point_b,
       delta = point_a - point_b, se = se, bootstrap = bootstrap)
}
