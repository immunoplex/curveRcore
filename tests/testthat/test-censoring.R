# Increasing curve on the fitting scale: a = 1 (low), d = 4 (high).
# LOD responses 1.2 / 3.8 map to log10 conc -0.5 / 2.5; LLOQ 0, ULOQ 2.
lim <- list(lloq_log10 = 0, uloq_log10 = 2, lod_response = 1.2, ulod_response = 3.8,
            lod_log10 = -0.5, ulod_log10 = 2.5, a = 1, d = 4)
cc <- function(resp, x, dil = 1, ...) {
  args <- utils::modifyList(lim, list(...))
  do.call(classify_censoring, c(list(response = resp, x_log10 = x, dilution = dil), args))
}

test_that("each class gets the right type and dilution-corrected bounds", {
  r <- cc(resp = c(1.1, 1.5, 2.5, 3.6, 3.9, NA),
          x    = c(-0.8, -0.2, 1.0, 2.2, 2.7, NA), dil = 10)
  expect_equal(r$censor_class,
               c("below_lod", "lod_to_lloq", "quantified", "above_uloq", "saturated", "no_response"))
  expect_equal(r$censor_type, c("left", "interval", "none", "right", "right", NA))
  expect_equal(r$censor_upper[1], 10^-0.5 * 10)
  expect_equal(c(r$censor_lower[2], r$censor_upper[2]), c(10^-0.5, 10^0) * 10)
  expect_true(is.na(r$censor_lower[3]) && is.na(r$censor_upper[3]))
  expect_equal(r$censor_lower[4:5], rep(10^2 * 10, 2))
  expect_true(all(is.na(r$censor_upper[4:5])))
})

test_that("a non-invertible point estimate is classified from the response", {
  r <- cc(resp = c(0.9, 4.2), x = c(NaN, NaN))
  expect_equal(r$censor_class, c("below_lod", "saturated"))
  r2 <- cc(resp = c(0.9, 4.2), x = c(NaN, NaN), lod_response = NA, lod_log10 = NA,
           ulod_response = NA, ulod_log10 = NA)
  expect_equal(r2$censor_class, c("below_lloq", "saturated"))
})

test_that("missing LOD collapses the low band to below_lloq", {
  r <- cc(resp = 1.5, x = -0.2, lod_response = NA, lod_log10 = NA)
  expect_equal(r$censor_class, "below_lloq")
  expect_equal(r$censor_upper, 1)
})

test_that("LOD at or above LLOQ leaves no lod_to_lloq band", {
  r <- cc(resp = 1.5, x = -0.2, lod_response = 1.4, lod_log10 = 0.3)
  expect_equal(r$censor_class, "below_lloq")
})

test_that("a curve with no LOQs gives NA except for LOD-based classes", {
  r <- cc(resp = c(1.1, 2.5), x = c(-0.8, 1), lloq_log10 = NA, uloq_log10 = NA)
  expect_equal(r$censor_class, c("below_lod", NA))
})

test_that("decreasing curves use log10 LOD comparisons, not response", {
  r <- cc(resp = c(3.5, 2.0), x = c(-0.8, 1), a = 4, d = 1,
          lod_response = NA, ulod_response = NA)
  expect_equal(r$censor_class, c("below_lod", "quantified"))
})

test_that("fixed-a fits (no a estimate) still use response-scale LOD checks", {
  r <- cc(resp = 1.1, x = NaN, a = NA)
  expect_equal(r$censor_class, "below_lod")
})

test_that("bad dilution keeps the class but gives NA bounds", {
  r <- cc(resp = 1.1, x = -0.8, dil = NA)
  expect_equal(r$censor_class, "below_lod")
  expect_true(is.na(r$censor_upper))
})

fake_cr <- function(curve_id, best = "logistic4") {
  list(
    meta = list(curve_id = curve_id, method = "frequentist"),
    selection = list(best_model_name = best),
    ensemble = list(logistic4 = list(
      eligibility = list(lloq = 0, uloq = 2),
      parameters = data.frame(term = c("a", "b", "c", "d"), estimate = c(1, 1, 1, 4),
                              std_error = 0.1))),
    detection_limits = list(lods = list(
      lower_lod_response = 1.2, upper_lod_response = 3.8,
      lower_lod_log10_conc = -0.5, upper_lod_log10_conc = 2.5)),
    samples = data.frame(sampleid = c("s1", "s2", "s3"),
                         observed_response_fit = c(1.1, 2.5, 3.9),
                         predicted_log10_concentration = c(-0.8, 1, 2.7),
                         dilution = c(2, 2, 2)))
}

test_that("classify_censoring_multiplate adds the four columns per curve", {
  mp <- structure(list(plates = list(`1` = fake_cr(1), `2` = fake_cr(2, best = NA))),
                  class = "calibration_result_multiplate")
  out <- classify_censoring_multiplate(mp)
  s1 <- out$plates[["1"]]$samples
  expect_equal(s1$censor_class, c("below_lod", "quantified", "saturated"))
  expect_equal(s1$censor_upper[1], 10^-0.5 * 2)
  expect_equal(names(s1)[1:3], c("sampleid", "observed_response_fit",
                                 "predicted_log10_concentration"))
  s2 <- out$plates[["2"]]$samples
  expect_true(all(c("censor_class", "censor_type", "censor_lower", "censor_upper") %in% names(s2)))
  expect_true(all(is.na(s2$censor_class)))
  # idempotent: re-running does not duplicate columns
  again <- classify_censoring_multiplate(out)
  expect_equal(ncol(again$plates[["1"]]$samples), ncol(s1))
})
