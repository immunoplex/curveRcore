# =============================================================================
# test-comparison.R — Tests for compare_calibrations and related utilities
#
# Uses mock calibration_result objects to test the comparison functions
# in curveRcore's result_class.R without depending on any fitting package.
# =============================================================================


# ---- Helper: build mock calibration_results ----

make_mock_result <- function(method = "frequentist", model = "logistic4",
                             n_grid = 10, n_samples = 5, seed = 42) {
  set.seed(seed)
  x <- seq(0, 4, length.out = n_grid)

  grid <- data.frame(
    log10_concentration = x,
    concentration       = 10^x,
    x_fit               = x,
    predicted_response  = logistic4(x, a = 1.5, b = 1.0, c = 2.0, d = 4.5) +
      rnorm(n_grid, 0, 0.02),
    ci_lower            = logistic4(x, a = 1.5, b = 1.0, c = 2.0, d = 4.5) - 0.1,
    ci_upper            = logistic4(x, a = 1.5, b = 1.0, c = 2.0, d = 4.5) + 0.1,
    predicted_concentration = x + rnorm(n_grid, 0, 0.01),
    se_concentration    = abs(rnorm(n_grid, 0.05, 0.01)),
    pcov                = runif(n_grid, 10, 100),
    pcov_pass           = rep(TRUE, n_grid),
    stringsAsFactors    = FALSE
  )

  params <- data.frame(
    term      = c("a", "b", "c", "d"),
    estimate  = c(1.5, 1.0, 2.0, 4.5),
    mean      = c(1.5, 1.0, 2.0, 4.5),
    std_error = c(0.05, 0.1, 0.08, 0.06),
    sd        = c(0.05, 0.1, 0.08, 0.06),
    stringsAsFactors = FALSE
  )

  samples <- if (n_samples > 0) {
    data.frame(
      sampleid  = paste0("S", 1:n_samples),
      curve_id  = rep(1L, n_samples),
      predicted_log10_concentration = runif(n_samples, 1, 4),
      final_concentration = 10^runif(n_samples, 1, 4) * 400,
      se_concentration    = runif(n_samples, 0.03, 0.1),
      pcov                = runif(n_samples, 15, 80),
      stringsAsFactors    = FALSE
    )
  } else NULL

  ensemble <- list()
  ensemble[[model]] <- list(
    model_name = model,
    converged  = TRUE,
    parameters = params,
    fit_stats  = list(aic = -10, bic = -8, rss = 0.01)
  )

  meta <- list(
    method = method, package = paste0("curveR", method),
    curve_id = "1",
    response_var = "mfi", independent_var = "concentration",
    is_log_response = TRUE, is_log_independent = TRUE
  )

  new_calibration_result(
    meta = meta, ensemble = ensemble,
    selection = list(best_model_name = model, criterion = "AIC"),
    grid = grid, samples = samples
  )
}


# ============================================================================
# 1. compare_calibrations — grid merge
# ============================================================================

test_that("compare_calibrations merges two grids correctly", {
  freq  <- make_mock_result("frequentist", seed = 1)
  bayes <- make_mock_result("bayesian", seed = 2)

  comp <- compare_calibrations(freq, bayes)

  expect_s3_class(comp, "calibration_comparison")
  expect_true("frequentist_predicted_response" %in% names(comp))
  expect_true("bayesian_predicted_response" %in% names(comp))
  expect_true("log10_concentration" %in% names(comp))
  expect_equal(nrow(comp), 10)
})

test_that("compare_calibrations works with custom labels", {
  a <- make_mock_result("frequentist", seed = 1)
  b <- make_mock_result("bayesian", seed = 2)

  comp <- compare_calibrations(a, b, label_a = "nls", label_b = "stan")

  expect_true("nls_predicted_response" %in% names(comp))
  expect_true("stan_predicted_response" %in% names(comp))
})


# ============================================================================
# 2. compare_parameters
# ============================================================================

test_that("compare_parameters produces paired table", {
  freq  <- make_mock_result("frequentist", seed = 1)
  bayes <- make_mock_result("bayesian", seed = 2)

  pcomp <- compare_parameters(freq, bayes)

  expect_true(is.data.frame(pcomp))
  expect_true("term" %in% names(pcomp))
  expect_true("diff" %in% names(pcomp))
  expect_true("rel_diff" %in% names(pcomp))
  expect_equal(nrow(pcomp), 4)
})


# ============================================================================
# 3. compare_samples
# ============================================================================

test_that("compare_samples merges by sampleid and curve_id", {
  freq  <- make_mock_result("frequentist", n_samples = 5, seed = 1)
  bayes <- make_mock_result("bayesian", n_samples = 5, seed = 2)

  scomp <- compare_samples(freq, bayes)

  expect_true(is.data.frame(scomp))
  expect_equal(nrow(scomp), 5)
  expect_true("frequentist_predicted_log10_concentration" %in% names(scomp))
  expect_true("bayesian_predicted_log10_concentration" %in% names(scomp))
  expect_true("conc_diff" %in% names(scomp))
})

test_that("compare_samples returns NULL when samples missing", {
  freq  <- make_mock_result("frequentist", n_samples = 0, seed = 1)
  bayes <- make_mock_result("bayesian", n_samples = 5, seed = 2)

  expect_warning(scomp <- compare_samples(freq, bayes))
  expect_null(scomp)
})


# ============================================================================
# 4. agreement_metrics
# ============================================================================

test_that("agreement_metrics computes correct values for identical vectors", {
  x <- c(1, 2, 3, 4, 5)
  m <- agreement_metrics(x, x)

  expect_equal(m$n, 5)
  expect_equal(m$bias, 0)
  expect_equal(m$mae, 0)
  expect_equal(m$rmse, 0)
  expect_equal(m$cor, 1)
  expect_equal(m$ccc, 1)
})

test_that("agreement_metrics handles NAs", {
  x <- c(1, 2, NA, 4, 5)
  y <- c(1.1, 2.1, 3.1, NA, 5.1)

  m <- agreement_metrics(x, y, na.rm = TRUE)
  expect_equal(m$n, 3)
  expect_true(m$bias > 0)
})

test_that("agreement_metrics returns NA for too-short vectors", {
  m <- agreement_metrics(1, 1)
  expect_true(is.na(m$cor))
})

test_that("agreement_metrics detects systematic bias", {
  x <- 1:10
  y <- x + 0.5

  m <- agreement_metrics(x, y)
  expect_equal(m$bias, 0.5)
  expect_equal(m$mae, 0.5)
  expect_equal(m$cor, 1)
  expect_true(m$ccc < 1)
})
