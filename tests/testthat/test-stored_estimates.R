# estimate_points_from_stored(): back-calculation from a stored fit.

th <- c(a = 1, b = 0.4, c = 1, d = 4)          # logistic4, log10 x
x0 <- c(-0.5, 0.5, 1, 1.6, 2.5)
y0 <- logistic4(x0, th[["a"]], th[["b"]], th[["c"]], th[["d"]])
lims <- list(lloq_log10 = 0, uloq_log10 = 2, lower_lod_response = 1.1,
             upper_lod_response = 3.9, lower_lod_log10_conc = -0.8,
             upper_lod_log10_conc = 3, inflect_x = 1)

test_that("frequentist point estimate is the inverse at the stored parameters", {
  out <- estimate_points_from_stored(y0, dilution = 10, model_name = "logistic4",
                                     params = th, limits = lims)
  expect_equal(out$predicted_log10_concentration, x0, tolerance = 1e-8)
  expect_equal(out$final_concentration, 10^x0 * 10, tolerance = 1e-8)
  # no grid, no draws: se unknown -> pcov at the cap, no quantiles
  expect_true(all(is.na(out$se_concentration)))
  expect_equal(out$pcov, rep(150, 5))
  expect_true(all(is.na(out$conc_q_med)))
  expect_equal(attr(out, "se_source"), rep("none", 5))
  expect_equal(out$pcov_gate_class, c("belowLLOQ", "inDynamicRange", "inDynamicRange",
                                      "inDynamicRange", "aboveULOQ"))
  expect_equal(out$censor_class[c(2, 3, 4)], rep("quantified", 3))
  expect_true(all(names(out) %in% point_estimates_contract("controls")$column))
})

test_that("se is read off the stored precision profile when there are no draws", {
  grid <- data.frame(log10_concentration = seq(-1, 3, by = 0.5))
  grid$se_concentration <- 0.02 + 0.01 * (grid$log10_concentration - 1)^2
  out <- estimate_points_from_stored(y0, model_name = "logistic4", params = th,
                                     grid = grid, limits = lims)
  expect_equal(attr(out, "se_source"), rep("profile", 5))
  expect_equal(out$se_concentration,
               stats::approx(grid$log10_concentration, grid$se_concentration, xout = x0)$y,
               tolerance = 1e-8)
  expect_equal(out$pcov, pmin(out$se_concentration * log(10) * 100, 150))
  expect_equal(out$pcov_rmse, out$pcov)
})

test_that("frequentist se with draws reproduces the delta method", {
  Sigma <- diag(c(0.01, 0.002, 0.003, 0.01))
  dimnames(Sigma) <- list(names(th), names(th))
  set.seed(11)
  Z <- matrix(stats::rnorm(4 * 40000), ncol = 4) %*% chol(Sigma)
  P <- sweep(Z, 2, th, "+"); colnames(P) <- names(th)
  draws <- lapply(seq_len(4), function(j) P[, j]); names(draws) <- names(th)
  out <- estimate_points_from_stored(y0[2:4], model_name = "logistic4", params = th,
                                     limits = lims, draws = draws)
  fns <- make_inv_and_grad_fixed("logistic4")
  analytic <- vapply(y0[2:4], function(y) {
    g <- fns$grad(y, th)[names(th)]
    sqrt(as.numeric(t(g) %*% Sigma %*% g))
  }, numeric(1))
  expect_equal(attr(out, "se_source"), rep("draws", 3))
  expect_equal(out$se_concentration, analytic, tolerance = 0.03)
  expect_equal(out$predicted_log10_concentration, x0[2:4], tolerance = 1e-8)
  expect_true(all(out$conc_q_lo < 10^x0[2:4] & 10^x0[2:4] < out$conc_q_hi))
})

test_that("bayesian estimate is the median and sd of the inverted draws", {
  set.seed(3)
  S <- 400
  draws <- list(a = rnorm(S, 1, 0.05), b = rnorm(S, 0.4, 0.01),
                c = rnorm(S, 1, 0.03), d = rnorm(S, 4, 0.05))
  out <- estimate_points_from_stored(y0[3], model_name = "logistic4", params = th,
                                     method = "bayesian", limits = lims, draws = draws)
  xs <- vapply(seq_len(S), function(s)
    inv_logistic4(y0[3], draws$a[s], draws$b[s], draws$c[s], draws$d[s]), numeric(1))
  expect_equal(out$predicted_log10_concentration, stats::median(xs))
  expect_equal(out$se_concentration, stats::sd(xs))
  # measurement noise widens the interval and leaves the caller's RNG untouched
  set.seed(99); before <- .Random.seed
  nz <- estimate_points_from_stored(y0[3], model_name = "logistic4", params = th,
                                    method = "bayesian", limits = lims, draws = draws,
                                    noise = list(sigma_obs = rep(0.05, S), nu = rep(30, S)))
  expect_identical(.Random.seed, before)
  expect_gt(nz$se_concentration, out$se_concentration)
  again <- estimate_points_from_stored(y0[3], model_name = "logistic4", params = th,
                                       method = "bayesian", limits = lims, draws = draws,
                                       noise = list(sigma_obs = rep(0.05, S), nu = rep(30, S)))
  expect_equal(again$se_concentration, nz$se_concentration)
})

test_that("fixed lower asymptote, out-of-range and missing responses", {
  out <- estimate_points_from_stored(c(0.9, NA, y0[3], 4.2), model_name = "logistic4",
                                     params = th[c("b", "c", "d")], fixed_a = 1,
                                     limits = lims)
  expect_true(is.na(out$predicted_log10_concentration[1]))
  expect_equal(out$censor_class[1], "below_lod")
  expect_equal(out$censor_class[2], "no_response")
  expect_equal(out$predicted_log10_concentration[3], 1, tolerance = 1e-8)
  expect_equal(out$censor_class[4], "saturated")
  expect_error(estimate_points_from_stored(1, model_name = "bogus", params = th), "unknown model")
})

test_that("standards get recovery and a bias-aware pcov_rmse", {
  grid <- data.frame(log10_concentration = seq(-1, 3, by = 0.5), se_concentration = 0.03)
  nominal <- x0 + c(0, 0.05, 0, -0.1, 0)
  out <- estimate_points_from_stored(y0, model_name = "logistic4", params = th,
                                     grid = grid, limits = lims, nominal = nominal)
  expect_equal(out$recovery_pct, 100 * 10^(x0 - nominal), tolerance = 1e-8)
  expect_equal(out$pcov_rmse[2], sqrt(0.03^2 + 0.05^2) * log(10) * 100, tolerance = 1e-8)
  expect_true(all(out$pcov_rmse >= out$pcov - 1e-9))
  expect_silent(validate_point_estimates(
    cbind(data.frame(curve_id = 1, method = "frequentist", well = paste0("A", 1:5),
                     dilution = 1:5, included = TRUE, exclusion_reason = "none"), out),
    "standards"))
})
