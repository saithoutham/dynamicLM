test_that("strict and delayed entry create the intended risk sets", {
  d <- data.frame(
    id = c(1, 1, 2, 3, 4),
    visit = c(50, 55, 56, 62, 40),
    entry = c(50, 50, 56, 62, 40),
    exit = c(70, 70, 65, 66, 54),
    event = c(1, 1, 1, 0, 1),
    x = c(1, 2, 3, 4, 5)
  )
  outcome <- list(time = "exit", status = "event")
  covs <- list(fixed = "entry", varying = "x")

  shared <- stack_data(d, outcome, 55, 10, covs, "long",
                       id = "id", rtime = "visit")
  strict <- stack_data(d, outcome, 55, 10, covs, "long",
                       id = "id", rtime = "visit", entry = "entry",
                       entry_mode = "strict")
  delayed <- stack_data(d, outcome, 55, 10, covs, "long",
                        id = "id", rtime = "visit", entry = "entry",
                        entry_mode = "delayed")

  expect_identical(shared$data$id, c(1, 2, 3))
  expect_identical(strict$data$id, 1)
  expect_identical(delayed$data$id, c(1, 2, 3))
  expect_identical(delayed$data$.LM_entry, c(55, 56, 62))
  expect_identical(delayed$data$x, c(2, 3, 4))
  expect_identical(shared$entry_col, "LM")
  expect_identical(strict$entry_col, "LM")
  expect_identical(delayed$entry_col, ".LM_entry")
  expect_true(all(strict$data$entry <= strict$data$LM))
  expect_true(all(delayed$data$.LM_entry < delayed$data$exit))
})


test_that("wide shared-origin calls no longer require an irrelevant rtime", {
  d <- data.frame(id = 1:3, exit = c(2, 3, 4), event = c(1, 0, 1), x = 1:3)
  out <- stack_data(
    d, list(time = "exit", status = "event"), lms = 0, w = 2,
    covs = list(fixed = "x", varying = NULL), format = "wide", id = "id"
  )
  expect_s3_class(out, "LMdataframe")
  expect_equal(nrow(out$data), 3)
})


test_that("delayed entry propagates into penalized survival responses", {
  d <- data.frame(
    id = 1:4, visit = c(50, 56, 62, 52), entry = c(50, 56, 62, 52),
    exit = c(70, 65, 66, 63), event = c(1, 1, 0, 1), x = c(1, 3, 4, 2)
  )
  lmdata <- stack_data(
    d, list(time = "exit", status = "event"), 55, 10,
    list(fixed = c("entry", "x"), varying = NULL), "long",
    id = "id", rtime = "visit", entry = "entry", entry_mode = "delayed"
  )
  checked <- check_penlm_inputs(lmdata, y = "x")
  response <- checked$y[[1]]
  expect_equal(unname(response[, "start"]), lmdata$data$.LM_entry)
  expect_true(all(response[, "start"] < response[, "stop"]))
})


test_that("model formula cannot silently discard delayed entry", {
  d <- data.frame(
    id = rep(1:8, each = 1), visit = c(50, 51, 52, 56, 57, 58, 59, 60),
    entry = c(50, 51, 52, 56, 57, 58, 59, 60),
    exit = c(61, 62, 63, 64, 65, 66, 67, 68),
    event = c(1, 0, 1, 0, 1, 0, 1, 0), x = seq_len(8)
  )
  lmdata <- stack_data(
    d, list(time = "exit", status = "event"), 55, 10,
    list(fixed = c("entry", "x"), varying = NULL), "long",
    id = "id", rtime = "visit", entry = "entry", entry_mode = "delayed"
  )
  expect_error(
    dynamic_lm(lmdata, survival::Surv(LM, exit, event) ~ x + cluster(id),
               "coxph"),
    "must use `.LM_entry`"
  )
})


test_that("left-truncated IPCW conditions on subject entry", {
  entry <- c(0, 0, 3, 3)
  exit <- c(2, 5, 4, 5)
  status <- c(0, 0, 1, 0)
  weights <- lt_censoring_weights(entry, exit, status, horizon = 5)

  expect_equal(weights$survival_at_horizon, rep(0.5, 4))
  expect_equal(weights$inverse_at_horizon, c(2, 2, 1, 1))
  expect_equal(weights$inverse_at_exit, c(1, 2, 1, 1))

  data <- data.frame(
    id = 1:4, LM = 0, analysis_entry = entry,
    observed_exit = exit, event = status
  )
  scored <- score_left_truncated(
    risk = c(0.1, 0.2, 0.9, 0.3), data = data,
    time_col = "observed_exit", status_col = "event",
    entry_col = "analysis_entry", id_col = "id", cause = 1,
    w = 5
  )
  expect_equal(scored$score$AUC, 1)
  expect_equal(scored$score$Brier, 0.045)
  expect_identical(scored$inference, "none")
  expect_match(scored$analytic_influence_function, "PENDING")
})


test_that("formula objects stored in symbols are evaluated", {
  d <- data.frame(
    id = 1:8, exit = 6:13, event = rep(c(1, 0), 4), x = seq_len(8)
  )
  lmdata <- stack_data(
    d, list(time = "exit", status = "event"), 5, 5,
    list(fixed = "x", varying = NULL), "wide", id = "id"
  )
  lmdata$func_covars <- list(identity)
  lmdata$func_lms <- list(identity)
  lmdata$lm_covs <- "x"
  formula_object <- Surv(LM, exit, event) ~ x + cluster(id)
  expect_s3_class(
    suppressWarnings(dynamic_lm(lmdata, formula_object, "coxph")),
    "LMcoxph"
  )
})
