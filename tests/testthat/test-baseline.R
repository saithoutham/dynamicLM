testthat::test_that("PBC competing-risks tutorial baseline remains stable", {
  expected_path <- here::here("results", "baseline_expected.rds")
  helper_path <- here::here("analysis", "_baseline_pipeline.R")
  testthat::skip_if_not(file.exists(expected_path) && file.exists(helper_path),
                        "project-level baseline artifacts are excluded from package builds")
  source(here::here("analysis", "_helpers.R"), local = TRUE)
  source(here::here("analysis", "_baseline_pipeline.R"), local = TRUE)
  testthat::expect_true(file.exists(expected_path))
  expected <- readRDS(expected_path)

  observed <- run_pbc_baseline(seed = expected$seed)

  testthat::expect_equal(observed$pbc_rows, expected$pbc_rows)
  testthat::expect_equal(observed$pbc_ids, expected$pbc_ids)
  testthat::expect_equal(observed$lm_rows, expected$lm_rows)
  testthat::expect_equal(observed$risk_counts, expected$risk_counts)
  testthat::expect_equal(observed$event_counts, expected$event_counts)
  testthat::expect_equal(observed$path_lambda_counts, expected$path_lambda_counts)
  testthat::expect_equal(observed$cv_lambda_min, expected$cv_lambda_min,
                         tolerance = 1e-7)
  testthat::expect_equal(observed$cv_lambda_1se, expected$cv_lambda_1se,
                         tolerance = 1e-7)
  testthat::expect_equal(observed$unpenalized_coefficients,
                         expected$unpenalized_coefficients, tolerance = 1e-7)
  testthat::expect_equal(observed$prediction_ranges,
                         expected$prediction_ranges, tolerance = 1e-7)
  testthat::expect_equal(observed$score_tables, expected$score_tables,
                         tolerance = 1e-7)
  testthat::expect_equal(observed$conditions, expected$conditions)
})
