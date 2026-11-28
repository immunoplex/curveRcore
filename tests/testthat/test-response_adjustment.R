# One rule for values <= 0 and blank subtraction (curveRcore 0.8.0).
# Origin: ELISA project 118, wavelength-subtracted absorbance with wells <= 0.

std_df <- function(resp, included = TRUE) {
  data.frame(dilution = 10^(seq_along(resp) - 1), od = resp, included = included)
}
pp_run <- function(std, ..., is_log_response = TRUE) {
  preprocess_standards(std, antigen_settings = list(standard_curve_concentration = 1000),
                       response_variable = "od", independent_variable = "concentration",
                       is_log_response = is_log_response, apply_prozone = FALSE, ...)
}

test_that("all-positive data: no offset, output identical to a plain log10", {
  std <- std_df(c(3, 1.2, 0.3, 0.05, 0.01))
  pp <- pp_run(std)
  expect_equal(pp$derived_stats$response_offset, 0)
  expect_equal(pp$data$od, log10(std$od))
  adj <- compute_response_adjustment(std, "od", TRUE,
                                     samples = data.frame(od = c(0.2, 2)))
  expect_equal(adj$offset, 0)
  expect_equal(adj$n_nonpositive, 0L)
})

test_that("a sample <= 0 sets one offset for standards and samples alike", {
  std <- std_df(c(3, 1.2, 0.3, 0.05, 0.01))
  smp <- data.frame(od = c(-0.004, 0.2, 1.5))
  adj <- compute_response_adjustment(std, "od", TRUE, samples = smp)
  expect_equal(adj$offset, 0.004 + 0.001)
  expect_equal(adj$n_nonpositive, 1L)

  pp <- pp_run(std, response_adjustment = adj)
  expect_equal(pp$data$od, log10(std$od + adj$offset))
  expect_equal(pp$derived_stats$response_offset, adj$offset)
  # Same raw value -> same fitting-scale value, whether standard or sample.
  ys <- transform_sample_response(c(std$od, smp$od), TRUE, pp$response_adjustment)
  expect_equal(ys[seq_len(nrow(std))], pp$data$od)
  # The most negative sample lands exactly delta above zero.
  expect_equal(ys[nrow(std) + 1], log10(0.001))
  # Round trip back to raw.
  expect_equal(untransform_response(ys, TRUE, pp$response_adjustment),
               c(std$od, smp$od))
})

test_that("standards <= 0 are kept (offset), not floored", {
  std <- std_df(c(3, 1.2, 0.3, 0.02, -0.003))
  pp <- pp_run(std)
  expect_equal(pp$derived_stats$response_offset, 0.003 + 0.001)
  expect_equal(pp$data$od, log10(std$od + 0.004))
  expect_true(all(is.finite(pp$data$od)))
})

test_that("blanks never decide the offset", {
  std <- std_df(c(3, 1.2, 0.3, 0.05, 0.01))
  blk <- data.frame(od = c(-0.02, 0.3))
  adj <- compute_response_adjustment(std, "od", TRUE, blank_data = blk)
  expect_equal(adj$offset, 0)
  pp <- pp_run(std, blank_data = blk)
  expect_equal(pp$data$od, log10(std$od))
  expect_true(all(is.finite(pp$blanks$od)))   # blank <= 0 gets the display floor
})

test_that("masked standards and controls do not decide the offset", {
  std <- std_df(c(3, 1.2, 0.3, 0.05, -0.01), included = c(TRUE, TRUE, TRUE, TRUE, FALSE))
  ctl <- data.frame(od = c(-0.5, 0.4), included = c(FALSE, TRUE))
  adj <- compute_response_adjustment(std, "od", TRUE, controls = ctl)
  expect_equal(adj$offset, 0)
  ctl2 <- data.frame(od = c(-0.002, 0.4))
  expect_equal(compute_response_adjustment(std, "od", TRUE, controls = ctl2)$offset, 0.003)
})

test_that("linear response: no offset, negatives fitted as-is", {
  std <- std_df(c(3, 1.2, 0.3, 0.02, -0.003))
  adj <- compute_response_adjustment(std, "od", FALSE, samples = data.frame(od = -1))
  expect_equal(adj$offset, 0)
  pp <- pp_run(std, is_log_response = FALSE)
  expect_equal(pp$data$od, std$od)
  expect_equal(transform_sample_response(c(-1, 2), FALSE, adj), c(-1, 2))
})

test_that("blank subtraction: no floor, arithmetic centre when a blank <= 0, samples subtracted too", {
  std <- std_df(c(3, 1.2, 0.3, 0.05, 0.01))
  blk <- data.frame(od = c(-0.01, 0.03))
  adj <- compute_response_adjustment(std, "od", TRUE, samples = data.frame(od = 0.5),
                                     blank_data = blk, blank_option = "subtracted")
  expect_equal(adj$blank_centre, 0.01)                 # arithmetic mean
  expect_equal(adj$blank_value, 0.01)
  expect_equal(adj$n_nonpositive, 1L)                 # 0.01 - 0.01 = 0 counts as <= 0
  expect_equal(adj$offset, 0.001)

  pp <- pp_run(std, blank_data = blk, blank_option = "subtracted",
               response_adjustment = adj)
  expect_equal(pp$data$od, log10(std$od - 0.01 + 0.001))   # no floor at 1
  expect_equal(transform_sample_response(0.5, TRUE, pp$response_adjustment),
               log10(0.5 - 0.01 + 0.001))

  # All blanks > 0: geometric mean, as before.
  adj_g <- compute_response_adjustment(std, "od", TRUE, blank_data = data.frame(od = c(0.01, 0.04)),
                                       blank_option = "subtracted")
  expect_equal(adj_g$blank_centre, 0.02)
})

test_that("NULL adjustment keeps the legacy sample transform", {
  expect_equal(transform_sample_response(c(-1, 100), TRUE, NULL), log10(c(1e-6, 100)))
})
