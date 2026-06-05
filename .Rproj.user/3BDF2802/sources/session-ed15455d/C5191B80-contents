# =============================================================================
# test-detection_limits.R
#
# Unit tests for curveRcore::detection_limits.R
# Tests all three tiers:
#   Tier 1: enrich_grid_with_d2y()
#   Tier 2: compute_shape_loq_from_grid()
#   Tier 3: compute_detection_limits(), compute_detection_limits_multiplate()
#
# Uses a synthetic calibration_result built from known 4PL parameters so
# that all outputs can be verified against analytical expectations.
# =============================================================================

# ── Helper: build a minimal calibration_result for testing ───────────────────

.make_test_cr <- function(method = "frequentist") {

  # Known 4PL on the log10 scale: a = 1.2, b = 0.6, c = 1.5, d = 4.3
  # These are realistic log10(MFI) values for a Luminex assay
  params <- c(a = 1.2, b = 0.6, c = 1.5, d = 4.3)

  # Build a 200-point grid
  x_grid <- seq(-1.5, 3.0, length.out = 200)

  # Forward model: 4PL
  y_grid <- params["a"] + (params["d"] - params["a"]) /
    (1 + exp(-(x_grid - params["c"]) / params["b"]))

  grid <- data.frame(
    log10_concentration     = x_grid,
    concentration           = 10^x_grid,
    x_fit                   = x_grid,
    predicted_response      = y_grid,
    ci_lower                = y_grid - 0.05,
    ci_upper                = y_grid + 0.05,
    predicted_concentration = x_grid,
    se_concentration        = rep(0.02, 200),
    pcov                    = rep(5, 200),
    pcov_rmse               = rep(5, 200),
    pcov_pass               = rep(TRUE, 200),
    stringsAsFactors        = FALSE
  )

  if (method == "frequentist") {
    params_df <- data.frame(
      term      = c("a", "b", "c", "d"),
      estimate  = unname(params),
      std_error = c(0.05, 0.03, 0.08, 0.04),
      stringsAsFactors = FALSE
    )
  } else {
    params_df <- data.frame(
      term  = c("a", "b", "c", "d"),
      mean  = unname(params),
      sd    = c(0.05, 0.03, 0.08, 0.04),
      q2.5  = unname(params) - 1.96 * c(0.05, 0.03, 0.08, 0.04),
      q50   = unname(params),
      q97.5 = unname(params) + 1.96 * c(0.05, 0.03, 0.08, 0.04),
      stringsAsFactors = FALSE
    )
  }

  ensemble <- list(
    logistic4 = list(
      model_name  = "logistic4",
      converged   = TRUE,
      parameters  = params_df,
      fit_stats   = list(aic = 10),
      raw_fit     = NULL,
      grid        = grid,
      eligibility = list(
        eligible            = TRUE,
        gates               = data.frame(gate = "rel_se", passed = TRUE,
                                         detail = "", stringsAsFactors = FALSE),
        dynamic_range_log10 = 3.5,
        lloq                = -1.0,
        uloq                = 2.5
      )
    )
  )

  meta <- list(
    method             = method,
    package            = if (method == "frequentist") "curveRfreq" else "curveRbayes",
    curve_id           = "1",
    response_var       = "mfi",
    independent_var    = "concentration",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    pcov_threshold     = 20
  )

  selection <- list(
    best_model_name = "logistic4",
    criterion       = "AIC+eligibility",
    fallback        = FALSE,
    fallback_reason = "",
    assessments     = list(logistic4 = ensemble$logistic4$eligibility),
    eligible_models = "logistic4"
  )

  cr <- list(
    meta      = meta,
    ensemble  = ensemble,
    selection = selection,
    grid      = grid,
    samples   = NULL
  )
  class(cr) <- c("calibration_result", "list")
  cr
}


# ── Tier 1: enrich_grid_with_d2y() ──────────────────────────────────────────

test_that("enrich_grid_with_d2y adds d2y_dx2 column", {
  cr   <- .make_test_cr()
  grid <- cr$grid

  # Grid should NOT have d2y_dx2 yet

  expect_false("d2y_dx2" %in% names(grid))

  enriched <- enrich_grid_with_d2y(grid, is_log_response = TRUE)

  expect_true("d2y_dx2" %in% names(enriched))
  expect_equal(nrow(enriched), nrow(grid))

  # First and last should be NA (boundary points)
  expect_true(is.na(enriched$d2y_dx2[1]))
  expect_true(is.na(enriched$d2y_dx2[nrow(enriched)]))

  # Interior points should be finite
  interior <- enriched$d2y_dx2[5:(nrow(enriched) - 5)]
  expect_true(all(is.finite(interior)))
})

test_that("d2y_dx2 has correct shape: peak before inflection, valley after", {
  cr       <- .make_test_cr()
  enriched <- enrich_grid_with_d2y(cr$grid, is_log_response = TRUE)

  ok  <- is.finite(enriched$d2y_dx2)
  x   <- enriched$log10_concentration[ok]
  d2y <- enriched$d2y_dx2[ok]

  # The inflection point is at c = 1.5
  c_val <- 1.5

  # Find the global max and min

  idx_max <- which.max(d2y)
  idx_min <- which.min(d2y)

  # Peak (max d2y) should be at x < c
  expect_lt(x[idx_max], c_val)

  # Valley (min d2y) should be at x > c
  expect_gt(x[idx_min], c_val)

  # Peak should be positive, valley should be negative
  expect_gt(d2y[idx_max], 0)
  expect_lt(d2y[idx_min], 0)
})

test_that("enrich_grid_with_d2y handles is_log_response = FALSE", {
  cr   <- .make_test_cr()
  grid <- cr$grid

  # Convert predicted_response to natural scale
  grid$predicted_response <- 10^grid$predicted_response

  enriched <- enrich_grid_with_d2y(grid, is_log_response = FALSE)

  expect_true("d2y_dx2" %in% names(enriched))
  interior <- enriched$d2y_dx2[5:(nrow(enriched) - 5)]
  expect_true(all(is.finite(interior)))
})

test_that("enrich_grid_with_d2y handles short grids gracefully", {
  grid <- data.frame(
    log10_concentration = c(0, 1),
    predicted_response  = c(1.5, 3.0)
  )
  enriched <- enrich_grid_with_d2y(grid, is_log_response = TRUE)
  expect_true(all(is.na(enriched$d2y_dx2)))
})


# ── Tier 2: compute_shape_loq_from_grid() ───────────────────────────────────

test_that("compute_shape_loq_from_grid returns correct structure", {
  cr       <- .make_test_cr()
  enriched <- enrich_grid_with_d2y(cr$grid, is_log_response = TRUE)
  result   <- compute_shape_loq_from_grid(enriched)

  expect_type(result, "list")
  expected_names <- c("shape_lloq_log10", "shape_uloq_log10",
                      "shape_lloq_conc",  "shape_uloq_conc",
                      "shape_lloq_response", "shape_uloq_response")
  expect_true(all(expected_names %in% names(result)))
})

test_that("shape-LLOQ < inflection < shape-ULOQ for symmetric 4PL", {
  cr       <- .make_test_cr()
  enriched <- enrich_grid_with_d2y(cr$grid, is_log_response = TRUE)
  result   <- compute_shape_loq_from_grid(enriched)

  c_val <- 1.5  # inflection on log10 scale

  expect_lt(result$shape_lloq_log10, c_val)
  expect_gt(result$shape_uloq_log10, c_val)

  # Natural-scale consistency
  expect_equal(result$shape_lloq_conc, 10^result$shape_lloq_log10,
               tolerance = 1e-10)
  expect_equal(result$shape_uloq_conc, 10^result$shape_uloq_log10,
               tolerance = 1e-10)
})

test_that("shape-LOQ response values are between asymptotes", {
  cr       <- .make_test_cr()
  enriched <- enrich_grid_with_d2y(cr$grid, is_log_response = TRUE)
  result   <- compute_shape_loq_from_grid(enriched)

  a_val <- 1.2;  d_val <- 4.3

  expect_gt(result$shape_lloq_response, a_val)
  expect_lt(result$shape_lloq_response, d_val)
  expect_gt(result$shape_uloq_response, a_val)
  expect_lt(result$shape_uloq_response, d_val)
})

test_that("compute_shape_loq_from_grid returns NAs for grid without d2y_dx2", {
  cr     <- .make_test_cr()
  result <- compute_shape_loq_from_grid(cr$grid)  # no d2y_dx2 column

  expect_true(is.na(result$shape_lloq_log10))
  expect_true(is.na(result$shape_uloq_log10))
})

test_that("shape-LOQ is symmetric around c for the symmetric 4PL", {
  cr       <- .make_test_cr()
  enriched <- enrich_grid_with_d2y(cr$grid, is_log_response = TRUE)
  result   <- compute_shape_loq_from_grid(enriched)

  c_val <- 1.5
  dist_below <- c_val - result$shape_lloq_log10
  dist_above <- result$shape_uloq_log10 - c_val

  # For a symmetric 4PL, shape-LOQs should be roughly equidistant from c
  expect_equal(dist_below, dist_above, tolerance = 0.05)
})


# ── Tier 3: compute_detection_limits() ──────────────────────────────────────

test_that("compute_detection_limits returns cr with $detection_limits", {
  cr     <- .make_test_cr()
  cr_out <- compute_detection_limits(cr, verbose = FALSE)

  expect_s3_class(cr_out, "calibration_result")
  expect_true("detection_limits" %in% names(cr_out))

  dl <- cr_out$detection_limits
  expect_equal(dl$model_name, "logistic4")
  expect_equal(dl$method, "frequentist")
  expect_equal(dl$alpha, 0.05)
})

test_that("LODs are correctly derived from asymptote CIs (frequentist)", {
  cr     <- .make_test_cr()
  cr_out <- compute_detection_limits(cr, alpha = 0.05, verbose = FALSE)

  dl  <- cr_out$detection_limits
  ens <- cr$ensemble$logistic4$parameters

  z   <- qnorm(0.975)
  a   <- ens$estimate[ens$term == "a"]
  se_a <- ens$std_error[ens$term == "a"]
  d   <- ens$estimate[ens$term == "d"]
  se_d <- ens$std_error[ens$term == "d"]

  expected_llod <- a + z * se_a     # upper CI of a
  expected_ulod <- d - z * se_d     # lower CI of d

  expect_equal(dl$lods$lower_lod_response, expected_llod, tolerance = 1e-10)
  expect_equal(dl$lods$upper_lod_response, expected_ulod, tolerance = 1e-10)
})

test_that("LODs are correctly derived from posterior quantiles (Bayesian)", {
  cr     <- .make_test_cr(method = "bayesian")
  cr_out <- compute_detection_limits(cr, alpha = 0.05, verbose = FALSE)

  dl     <- cr_out$detection_limits
  params <- cr$ensemble$logistic4$parameters

  expected_llod <- params$q97.5[params$term == "a"]
  expected_ulod <- params$q2.5[params$term == "d"]

  expect_equal(dl$lods$lower_lod_response, expected_llod, tolerance = 1e-10)
  expect_equal(dl$lods$upper_lod_response, expected_ulod, tolerance = 1e-10)
})

test_that("LOD concentrations are valid inversions of the LOD response", {
  cr     <- .make_test_cr()
  cr_out <- compute_detection_limits(cr, verbose = FALSE)
  dl     <- cr_out$detection_limits

  params <- c(a = 1.2, b = 0.6, c = 1.5, d = 4.3)

  # Re-invert manually
  llod_resp <- dl$lods$lower_lod_response
  manual_x  <- inv_logistic4(llod_resp, params["a"], params["b"],
                             params["c"], params["d"])

  expect_equal(dl$lods$lower_lod_log10_conc, as.numeric(manual_x),
               tolerance = 1e-8)
  expect_equal(dl$lods$lower_lod_conc, 10^as.numeric(manual_x),
               tolerance = 1e-8)
})

test_that("ULOD < LLOD triggers NA for upper LOD", {
  cr <- .make_test_cr()

  # Make std_error of d very large so lower CI of d < upper CI of a
  cr$ensemble$logistic4$parameters$std_error[
    cr$ensemble$logistic4$parameters$term == "d"] <- 5.0

  cr_out <- compute_detection_limits(cr, verbose = FALSE)
  expect_true(is.na(cr_out$detection_limits$lods$upper_lod_response))
})

test_that("RDL uses CI-modified curves (compressed and expanded)", {
  cr     <- .make_test_cr()
  cr_out <- compute_detection_limits(cr, verbose = FALSE)
  dl     <- cr_out$detection_limits

  # MDC and RDL should all be finite for this well-behaved curve
  expect_true(is.finite(dl$mdc_rdl$mdc_lower_log10))
  expect_true(is.finite(dl$mdc_rdl$rdl_lower_log10))

  # RDL lower uses compressed d (lower CI) → shifts LLOQ rightward
  # So rdl_lower should be >= mdc_lower for a well-behaved curve
  # (compressed curve has a smaller gap between a and d)
  expect_gte(dl$mdc_rdl$rdl_lower_log10, dl$mdc_rdl$mdc_lower_log10)
})

test_that("compute_detection_limits handles non-converged model gracefully", {
  cr <- .make_test_cr()
  cr$ensemble$logistic4$converged <- FALSE

  cr_out <- compute_detection_limits(cr, verbose = FALSE)
  dl     <- cr_out$detection_limits

  expect_true(is.na(dl$lods$lower_lod_response))
  expect_true(is.na(dl$mdc_rdl$mdc_lower_log10))
})

test_that("natural-scale concentrations are 10^log10 values", {
  cr     <- .make_test_cr()
  cr_out <- compute_detection_limits(cr, verbose = FALSE)
  dl     <- cr_out$detection_limits

  if (is.finite(dl$lods$lower_lod_log10_conc)) {
    expect_equal(dl$lods$lower_lod_conc,
                 10^dl$lods$lower_lod_log10_conc, tolerance = 1e-10)
  }
  if (is.finite(dl$mdc_rdl$rdl_upper_log10)) {
    expect_equal(dl$mdc_rdl$rdl_upper_conc,
                 10^dl$mdc_rdl$rdl_upper_log10, tolerance = 1e-10)
  }
})


# ── Tier 3: compute_detection_limits_multiplate() ────────────────────────────

test_that("compute_detection_limits_multiplate iterates over all plates", {
  cr1 <- .make_test_cr()
  cr2 <- .make_test_cr()
  cr2$meta$curve_id <- "2"

  mp <- list(
    meta   = list(method = "frequentist"),
    plates = list("1" = cr1, "2" = cr2)
  )
  class(mp) <- c("calibration_result_multiplate", "list")

  mp_out <- compute_detection_limits_multiplate(mp, verbose = FALSE)

  expect_true("detection_limits" %in% names(mp_out$plates[["1"]]))
  expect_true("detection_limits" %in% names(mp_out$plates[["2"]]))
  expect_equal(mp_out$plates[["1"]]$detection_limits$model_name, "logistic4")
  expect_equal(mp_out$plates[["2"]]$detection_limits$model_name, "logistic4")
})
