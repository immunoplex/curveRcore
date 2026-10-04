# =============================================================================
# test-tidy_extractors.R -- tidy_samples()/tidy_grid() (and friends) on
# calibration_result_multiplate, across plate counts.
#
# Regression coverage for a bug where every *_multiplate tidy accessor called
# `do.call(.cr_rbind_fill, parts)` instead of `.cr_rbind_fill(parts)`.
# `.cr_rbind_fill()` takes a single `dfs` (list of data frames) argument, so
# `do.call(f, parts)` -- which spreads `parts` as *separate* positional
# arguments -- raised "unused arguments" for any multiplate result with more
# than one non-NULL plate, and silently misbehaved (iterating over a single
# data frame's columns instead of a one-element list) for exactly one plate.
# =============================================================================

make_plate <- function(cid, with_samples = TRUE) {
  new_calibration_result(
    meta = list(method = "frequentist", package = "test", curve_id = cid,
               response_var = "mfi", independent_var = "concentration",
               is_log_response = TRUE, is_log_independent = TRUE),
    grid = data.frame(
      predicted_concentration = seq(-1, 1, length.out = 3),
      se_concentration        = seq(0.3, 0.1, length.out = 3),
      pcov                    = seq(70, 20, length.out = 3),
      pcov_pass               = c(FALSE, TRUE, TRUE)
    ),
    samples = if (with_samples) {
      data.frame(
        sampleid                = paste0(cid, "_S", 1:2),
        predicted_concentration = c(-1, 0),
        se_concentration        = c(0.3, 0.2),
        pcov                    = c(70, 46),
        pcov_pass               = c(FALSE, TRUE)
      )
    } else NULL
  )
}

test_that("tidy_samples/tidy_grid row-bind correctly across 2+ plates", {
  mp <- new_calibration_result_multiplate(
    meta = list(method = "frequentist", package = "test",
               curve_ids = c("p1", "p2", "p3"), is_log_independent = TRUE),
    plates = list(p1 = make_plate("p1"), p2 = make_plate("p2"),
                 p3 = make_plate("p3"))
  )

  s <- tidy_samples(mp)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 6L)  # 3 plates x 2 samples
  expect_setequal(unique(s$curve_id), c("p1", "p2", "p3"))
  expect_true(all(c("predicted_concentration", "se_concentration") %in% names(s)))

  g <- tidy_grid(mp)
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 9L)  # 3 plates x 3 grid points
  expect_setequal(unique(g$curve_id), c("p1", "p2", "p3"))
})

test_that("tidy_samples/tidy_grid row-bind correctly for exactly 1 plate", {
  # The do.call() misuse silently "worked" (no error) for a single plate by
  # passing the lone data frame itself -- instead of a one-element list
  # containing it -- as `.cr_rbind_fill()`'s `dfs` argument, which then
  # iterated over the data frame's *columns*. Assert the real, intended
  # shape comes back.
  mp1 <- new_calibration_result_multiplate(
    meta = list(method = "frequentist", package = "test",
               curve_ids = "p1", is_log_independent = TRUE),
    plates = list(p1 = make_plate("p1"))
  )

  s <- tidy_samples(mp1)
  expect_equal(nrow(s), 2L)
  expect_identical(unique(s$curve_id), "p1")
  expect_true(is.numeric(s$predicted_concentration))

  g <- tidy_grid(mp1)
  expect_equal(nrow(g), 3L)
  expect_identical(unique(g$curve_id), "p1")
})

test_that("tidy_hyperparam/tidy_fit_diag row-bind correctly across plates", {
  plate_with_pop <- function(cid) {
    cr <- make_plate(cid)
    cr$population <- list(
      params = data.frame(term = "sigma_obs", estimate = 1, std_error = 0.1,
                          q_lo = 0.8, q_med = 1, q_hi = 1.2),
      fit_diag = list(fit_seconds = 1.2, converged = TRUE)
    )
    cr
  }
  mp <- new_calibration_result_multiplate(
    meta = list(method = "frequentist", package = "test",
               curve_ids = c("p1", "p2"), is_log_independent = TRUE),
    plates = list(p1 = plate_with_pop("p1"), p2 = plate_with_pop("p2"))
  )

  hp <- tidy_hyperparam(mp)
  expect_equal(nrow(hp), 2L)
  expect_setequal(hp$curve_id, c("p1", "p2"))

  fd <- tidy_fit_diag(mp)
  expect_equal(nrow(fd), 2L)
  expect_setequal(fd$curve_id, c("p1", "p2"))
})
