# ---- point_estimates_contract -------------------------------------------------

test_that("point_estimates_contract has the expected shape and keys", {
  ct <- point_estimates_contract()
  expect_setequal(names(ct), c("table", "column", "type", "key", "nullable",
                               "added_in", "description"))
  c_ <- point_estimates_contract("controls"); s_ <- point_estimates_contract("standards")
  expect_equal(c_$column[c_$key], c("curve_id", "method", "well"))
  expect_equal(s_$column[s_$key], c("curve_id", "method", "well", "dilution"))
  expect_false(anyDuplicated(paste(ct$table, ct$column)) > 0)
  expect_true(all(nzchar(ct$description)))
  # the existing 13 calib_standards columns are untouched; every 0.7.0 addition is nullable
  expect_equal(sum(s_$added_in == "0.4"), 13L)
  expect_true(all(s_$nullable[s_$added_in == "0.7.0"]))
  expect_true("recovery_pct" %in% s_$column)
  expect_false("recovery_pct" %in% c_$column)
  expect_true(all(c("source", "sampleid", "included", "mask_reason") %in% c_$column))
})

test_that("the estimate columns are the same for controls and standards", {
  ce <- point_estimates_contract("controls", estimates_only = TRUE)
  se <- point_estimates_contract("standards", estimates_only = TRUE)
  expect_equal(se$column, c(ce$column, "recovery_pct"))
  expect_equal(se$type[seq_len(nrow(ce))], ce$type)
})

# ---- DDL <-> contract sync ----------------------------------------------------

pe_ddl_path <- function() {
  p <- testthat::test_path("..", "..", "R", "calib_schema.sql")
  if (file.exists(p)) p else NA_character_
}
pe_ddl_columns <- function(sql, table) {
  start <- grep(paste0("CREATE TABLE IF NOT EXISTS madi_results.", table, " \\("), sql)
  stopifnot(length(start) == 1L)
  end <- start + which(grepl("^\\);", sql[(start + 1):length(sql)]))[1]
  body <- sql[(start + 1):(end - 1)]
  body <- body[grepl("^    [a-z_0-9]+ ", body) & !grepl("^    CONSTRAINT", body)]
  sub("^    ([a-z_0-9]+) .*$", "\\1", body)
}
pe_ddl_added <- function(sql, table) {
  start <- grep(paste0("^ALTER TABLE madi_results.", table, "$"), sql)
  stopifnot(length(start) == 1L)
  rest <- sql[(start + 1):length(sql)]
  blk  <- rest[seq_len(which(grepl(";$", rest))[1])]
  blk  <- blk[grepl("ADD COLUMN IF NOT EXISTS", blk)]
  sub("^ *ADD COLUMN IF NOT EXISTS ([a-z_0-9]+) .*$", "\\1", blk)
}

test_that("calib_schema.sql matches point_estimates_contract, column for column and in order", {
  p <- pe_ddl_path()
  testthat::skip_if(is.na(p), "calib_schema.sql not available (installed package)")
  sql <- readLines(p, warn = FALSE)
  expect_equal(pe_ddl_columns(sql, "calib_controls"), point_estimates_contract("controls")$column)
  s_ <- point_estimates_contract("standards")
  expect_equal(pe_ddl_columns(sql, "calib_standards"), s_$column[s_$added_in == "0.4"])
  expect_equal(pe_ddl_added(sql, "calib_standards"), s_$column[s_$added_in == "0.7.0"])
  expect_true(any(grepl("calib_controls_pkey PRIMARY KEY \\(curve_id, method, well\\)", sql)))
  expect_true(any(grepl("calib_controls_curve_id_fkey FOREIGN KEY \\(curve_id\\) REFERENCES madi_results.curve_lookup", sql)))
})

# ---- validate_point_estimates ---------------------------------------------------

ctrl_rows <- function() {
  data.frame(curve_id = c(1, 1), method = "frequentist", well = c("A1", "A2"),
             sampleid = c("A1", "A2"), source = "QC1 (Low)", dilution = 100,
             included = c(TRUE, FALSE), exclusion_reason = c("none", "masked"),
             mask_reason = c(NA, "bubble"),
             predicted_log10_concentration = c(0.5, 0.6), pcov = c(8, 9),
             pcov_gate_class = "inDynamicRange", censor_class = "quantified",
             censor_type = "none", job_id = "j1", stringsAsFactors = FALSE)
}

test_that("validate_point_estimates accepts conforming rows", {
  expect_invisible(validate_point_estimates(ctrl_rows(), "controls"))
  st <- data.frame(curve_id = 1, method = "bayesian", well = "B1", dilution = 2,
                   included = TRUE, exclusion_reason = "none", recovery_pct = 98.5)
  expect_silent(validate_point_estimates(st, "standards"))
})

test_that("validate_point_estimates reports every problem", {
  x <- ctrl_rows()
  x$bogus <- 1
  x$included[1] <- NA
  x$censor_class[2] <- "kinda_low"
  x$recovery_pct <- 100            # standards only
  err <- tryCatch(validate_point_estimates(x, "controls"), error = conditionMessage)
  expect_match(err, "unknown column\\(s\\): bogus, recovery_pct")
  expect_match(err, "required column 'included' has missing values")
  expect_match(err, "'censor_class' holds value\\(s\\) not allowed: kinda_low")
  y <- ctrl_rows(); y$well[2] <- "A1"
  expect_error(validate_point_estimates(y, "controls"), "duplicate primary-key rows")
  z <- ctrl_rows(); z$included <- NULL
  expect_error(validate_point_estimates(z, "controls"), "missing required column\\(s\\): included")
})

# ---- standard_recovery ----------------------------------------------------------

test_that("standard_recovery: log10 x", {
  r <- standard_recovery(predicted = c(1, 1.1, NA, 0.9), nominal = c(1, 1, 1, 1),
                         se = c(0.02, 0.02, 0.02, NA))
  expect_equal(r$recovery_pct, c(100, 100 * 10^0.1, NA, 100 * 10^-0.1))
  # unbiased point: pcov_rmse equals the pcov formula se * ln10 * 100
  expect_equal(r$pcov_rmse[1], 0.02 * log(10) * 100)
  expect_equal(r$pcov_rmse[2], sqrt(0.02^2 + 0.1^2) * log(10) * 100)
  expect_true(is.na(r$pcov_rmse[3]))
  expect_equal(r$pcov_rmse[4], 150)          # estimate without se -> cap
})

test_that("standard_recovery: natural-scale x and the cap", {
  r <- standard_recovery(predicted = c(10, 20), nominal = c(10, 10), se = c(1, 1),
                         is_log_independent = FALSE)
  expect_equal(r$recovery_pct, c(100, 200))
  # unbiased: se/(x ln10) on log10 -> back to CV = se/x*100
  expect_equal(r$pcov_rmse[1], 10)
  expect_equal(r$pcov_rmse[2], sqrt((1 / (20 * log(10)))^2 + log10(2)^2) * log(10) * 100)
  big <- standard_recovery(predicted = 3, nominal = 0, se = 0.01, cv_x_max = 50)
  expect_equal(big$pcov_rmse, 50)
  expect_error(standard_recovery(1:2, 1, 1), "same length")
})

# ---- gating, censoring and recovery over controls and standard estimates ------

pe_fake_cr <- function(curve_id, with_points = TRUE) {
  pts <- function() data.frame(observed_response_fit = c(1.1, 2.5, 3.9),
                               predicted_log10_concentration = c(-0.8, 1, 2.7),
                               predicted_concentration = c(-0.8, 1, 2.7),
                               se_concentration = c(0.2, 0.05, 0.3),
                               pcov = c(60, 10, 70), dilution = c(2, 2, 2))
  cr <- list(
    meta = list(curve_id = curve_id, method = "frequentist", pcov_threshold = 20,
                is_log_independent = TRUE),
    selection = list(best_model_name = "logistic4"),
    ensemble = list(logistic4 = list(
      eligibility = list(lloq = 0, uloq = 2),
      parameters = data.frame(term = c("a", "b", "c", "d"), estimate = c(1, 1, 1, 4),
                              std_error = 0.1))),
    detection_limits = list(lods = list(
      lower_lod_response = 1.2, upper_lod_response = 3.8,
      lower_lod_log10_conc = -0.5, upper_lod_log10_conc = 2.5)),
    samples = cbind(sampleid = c("s1", "s2", "s3"), pts()))
  if (with_points) {
    cr$controls <- cbind(well = c("A1", "A2", "A3"), source = "QC1", pts())
    cr$standard_estimates <- cbind(well = c("B1", "B2", "B3"), pts(),
                                   concentration = c(-0.7, 1.05, 2.5))
  }
  cr
}

test_that("controls and standard estimates are gated and censored like samples", {
  mp <- structure(list(plates = list(`1` = pe_fake_cr(1))),
                  class = "calibration_result_multiplate")
  out <- classify_censoring_multiplate(classify_pcov_gates_multiplate(mp))
  cr <- out$plates[["1"]]
  for (sl in c("controls", "standard_estimates")) {
    expect_equal(cr[[sl]]$censor_class, cr$samples$censor_class, info = sl)
    expect_equal(cr[[sl]]$censor_upper, cr$samples$censor_upper, info = sl)
    expect_equal(cr[[sl]]$pcov_gate_class, cr$samples$pcov_gate_class, info = sl)
    expect_equal(cr[[sl]]$pcov_pass, cr$samples$pcov_pass, info = sl)
  }
  expect_equal(cr$controls$pcov_gate_class, c("belowLLOQ", "inDynamicRange", "aboveULOQ"))
  expect_equal(cr$controls$source, rep("QC1", 3))
})

test_that("a result without the new slots is processed exactly as before", {
  mp <- structure(list(plates = list(`1` = pe_fake_cr(1, with_points = FALSE))),
                  class = "calibration_result_multiplate")
  out <- classify_censoring_multiplate(classify_pcov_gates_multiplate(mp))
  expect_null(out$plates[["1"]]$controls)
  expect_null(out$plates[["1"]]$standard_estimates)
  expect_equal(out$plates[["1"]]$samples$censor_class, c("below_lod", "quantified", "saturated"))
  expect_identical(add_standard_recovery_multiplate(mp), mp)
})

test_that("add_standard_recovery_multiplate fills recovery_pct and a bias-aware pcov_rmse", {
  mp <- structure(list(plates = list(`1` = pe_fake_cr(1))),
                  class = "calibration_result_multiplate")
  out <- add_standard_recovery_multiplate(mp)
  st <- out$plates[["1"]]$standard_estimates
  expect_equal(st$recovery_pct, 100 * 10^(c(-0.8, 1, 2.7) - c(-0.7, 1.05, 2.5)))
  expect_equal(st$pcov_rmse[2], sqrt(0.05^2 + 0.05^2) * log(10) * 100)
  expect_null(out$plates[["1"]]$controls$recovery_pct)
  # single-curve result
  one <- add_standard_recovery_multiplate(pe_fake_cr(7))
  expect_equal(one$standard_estimates$recovery_pct, st$recovery_pct)
})

# ---- control_reference ------------------------------------------------------------

test_that("control_reference: DDL matches the contract and the validator applies its rules", {
  p <- pe_ddl_path()
  testthat::skip_if(is.na(p), "calib_schema.sql not available (installed package)")
  sql <- readLines(p, warn = FALSE)
  expect_equal(pe_ddl_columns(sql, "control_reference"),
               point_estimates_contract("control_reference")$column)
  expect_true(any(grepl("CREATE OR REPLACE VIEW madi_results.calib_controls_qc AS", sql)))

  ok <- data.frame(project_id = 92L, control_source = "QC1 (Low)", antigen = "PT",
                   expected_conc = 12.5, expected_sd = 1.1, basis = "assigned")
  expect_silent(validate_point_estimates(ok, "control_reference"))
  bad <- rbind(ok, ok); bad$basis[2] <- "guess"
  err <- tryCatch(validate_point_estimates(bad, "control_reference"), error = conditionMessage)
  expect_match(err, "'basis' holds value\\(s\\) not allowed: guess")
  expect_match(err, "duplicate primary-key rows")
  expect_error(validate_point_estimates(ok[, -1], "control_reference"),
               "missing required column\\(s\\): project_id")
})
