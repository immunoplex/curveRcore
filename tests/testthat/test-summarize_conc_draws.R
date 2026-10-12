test_that("invert_with_bounds maps out-of-range responses to -Inf / +Inf", {
  a <- 1; b <- 1; c <- 2; d <- 4
  f <- function(y) inv_logistic4_fixed(y, a, b, c, d)
  y <- c(0.5, 1, 2.5, 4, 5, NA)
  x <- invert_with_bounds(f, y, a, d)
  expect_equal(x[1:2], c(-Inf, -Inf))
  expect_equal(x[3], inv_logistic4_fixed(2.5, a, b, c, d))
  expect_equal(x[4:5], c(Inf, Inf))
  expect_true(is.na(x[6]))
})

test_that("invert_with_bounds handles a decreasing curve (a above d)", {
  # negative slope: response falls as concentration rises; a is still the
  # low-concentration asymptote, now the high response.
  a <- 4; b <- 1; c <- 2; d <- 1
  f <- function(y) inv_logistic4_fixed(y, a, b, c, d)
  x <- invert_with_bounds(f, c(5, 0.5), a, d)
  expect_equal(x, c(-Inf, Inf))
})

test_that("invert_with_bounds handles unguarded inverses without NaN leaking", {
  a <- 1; b <- 1; c <- 2; d <- 4
  f <- function(y) inv_gompertz4_fixed(y, a, b, c, d)
  x <- invert_with_bounds(f, c(0.2, 3, 9), a, d)
  expect_equal(x[c(1, 3)], c(-Inf, Inf))
  expect_true(is.finite(x[2]))
})

test_that("invert_with_bounds returns NA for unusable asymptotes", {
  f <- function(y) y
  expect_true(all(is.na(invert_with_bounds(f, 1:3, NA, 4))))
  expect_true(all(is.na(invert_with_bounds(f, 1:3, 2, 2))))
})

test_that("summarize_conc_draws: quantiles, probabilities and censored shares", {
  set.seed(1)
  x <- rbind(rnorm(4000, 1, 0.1),                       # well inside
             c(rep(-Inf, 1000), rnorm(3000, -0.5, 0.2)))  # 25% beyond low asymptote
  s <- summarize_conc_draws(x, dilution = c(10, 2), lloq_log10 = 0, uloq_log10 = 2,
                            decision_cutoff = 5)
  expect_named(s, c("conc_q_lo", "conc_q_med", "conc_q_hi", "p_below_lloq",
                    "p_above_uloq", "p_below_cutoff", "frac_draws_below_a",
                    "frac_draws_above_d"))
  expect_equal(s$conc_q_med[1], 10^1 * 10, tolerance = 0.05)
  expect_lt(s$conc_q_lo[1], s$conc_q_med[1]); expect_gt(s$conc_q_hi[1], s$conc_q_med[1])
  expect_equal(s$p_below_lloq[1], 0)
  expect_equal(s$frac_draws_below_a[2], 0.25)
  expect_equal(s$conc_q_lo[2], 0)                       # censored draws -> 0
  expect_gt(s$p_below_lloq[2], 0.95)
  # cutoff 5 on the final scale = log10(5/2) pre-dilution for row 2
  expect_equal(s$p_below_cutoff[2], mean(x[2, ] < log10(5 / 2)))
})

test_that("summarize_conc_draws: unbounded upper quantile is NA, not Inf", {
  x <- matrix(c(rep(Inf, 100), rnorm(900, 2)), nrow = 1)
  s <- summarize_conc_draws(x, dilution = 1)
  expect_true(is.na(s$conc_q_hi))
  expect_equal(s$frac_draws_above_d, 0.1)
})

test_that("summarize_conc_draws: all draws censored / no draws / bad dilution", {
  s <- summarize_conc_draws(rbind(rep(-Inf, 50), rep(NA_real_, 50), rep(1, 50)),
                            dilution = c(1, 1, NA), lloq_log10 = 0, decision_cutoff = 3)
  expect_equal(s$conc_q_hi[1], 0)
  expect_equal(s$p_below_lloq[1], 1)
  expect_true(all(is.na(unlist(s[2, ]))))
  expect_true(is.na(s$conc_q_med[3]) && is.na(s$p_below_cutoff[3]))
  expect_equal(s$p_below_lloq[3], 0)                    # pre-dilution: no dilution needed
})

test_that("summarize_conc_draws: missing limits and cutoff give NA probabilities", {
  s <- summarize_conc_draws(matrix(rnorm(100), 1), dilution = 1)
  expect_true(is.na(s$p_below_lloq) && is.na(s$p_above_uloq) && is.na(s$p_below_cutoff))
})
