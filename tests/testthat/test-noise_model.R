# Additive + proportional response noise (curveRcore 0.8.0).

test_that("constant model returns sigma unchanged", {
  expect_equal(noise_sd(c(-1, 0, 1), "constant", sigma = 0.1), rep(0.1, 3))
  expect_equal(noise_sd(1, NULL, sigma = 0.2), 0.2)
  expect_error(noise_sd(1, "power", sigma = 0.1), "noise_model")
})

test_that("add_prop on the raw scale", {
  s <- noise_sd(c(0, 1, 4), "add_prop", sigma_add = 0.01, cv = 0.05,
                is_log_response = FALSE)
  expect_equal(s, sqrt(0.01^2 + (0.05 * c(0, 1, 4))^2))
})

test_that("add_prop on log10 is the delta-method transform and is finite where log10 response = 0", {
  mu <- log10(c(0.001, 0.01, 1, 4))
  s <- noise_sd(mu, "add_prop", sigma_add = 0.01, cv = 0.05, is_log_response = TRUE)
  raw <- 10^mu
  expect_equal(s, sqrt(0.01^2 + (0.05 * raw)^2) / (raw * log(10)))
  expect_true(all(is.finite(s) & s > 0))       # no collapse at response = 1
  # Low end noisy, high end tends to cv / ln10.
  expect_gt(s[1], s[3])
  expect_equal(noise_sd(log10(1e6), "add_prop", sigma_add = 0.01, cv = 0.05),
               0.05 / log(10), tolerance = 1e-6)
})

test_that("fit_add_prop_noise recovers sigma_add and cv", {
  set.seed(42)
  m_raw <- rep(10^seq(-2.5, 0.6, length.out = 40), 25)
  sd_raw <- sqrt(0.008^2 + (0.06 * m_raw)^2)
  y_raw <- m_raw + rnorm(length(m_raw), 0, sd_raw)
  f_raw <- fit_add_prop_noise(y_raw - m_raw, m_raw, is_log_response = FALSE)
  expect_true(f_raw$converged)
  expect_equal(f_raw$sigma_add, 0.008, tolerance = 0.15)
  expect_equal(f_raw$cv, 0.06, tolerance = 0.15)
  # Same truth seen on the log scale (residuals approximated by the delta method).
  mu <- log10(m_raw)
  r_log <- (y_raw - m_raw) / (m_raw * log(10))
  f_log <- fit_add_prop_noise(r_log, mu, is_log_response = TRUE)
  expect_equal(f_log$cv, 0.06, tolerance = 0.15)
})

test_that("new_study_params defaults to the constant model and validates", {
  expect_equal(new_study_params()$noise_model, "constant")
  expect_equal(new_study_params(noise_model = "add_prop")$noise_model, "add_prop")
  expect_error(new_study_params(noise_model = "x"), "noise_model")
})

test_that("stored estimates use add_prop terms per point (frequentist)", {
  p_const <- list(a = 0.01, b = 0.5, c = 2, d = 4, sigma_resid = 0.05)
  p_addp  <- list(a = 0.01, b = 0.5, c = 2, d = 4, sigma_add = 0.01, cv = 0.05)
  dr <- list(a = rnorm(200, 0.01, 0.001), b = rnorm(200, 0.5, 0.01),
             c = rnorm(200, 2, 0.02), d = rnorm(200, 4, 0.05))
  e1 <- estimate_points_from_stored(c(0.1, 2), model_name = "logistic4",
                                    params = p_const, draws = dr, is_log_response = FALSE)
  e2 <- estimate_points_from_stored(c(0.1, 2), model_name = "logistic4",
                                    params = p_addp, draws = dr, is_log_response = FALSE)
  expect_equal(nrow(e2), 2L)
  expect_true(all(is.finite(e2$se_concentration %||% e2[[grep("^se", names(e2))[1]]])))
  expect_false(isTRUE(all.equal(e1, e2)))
})
