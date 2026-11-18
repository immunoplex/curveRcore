# =============================================================================
# point_estimates_contract.R -- the data contract for modelled estimates of
# controls (madi_results.calib_controls, new in 0.7.0) and calibrators (the
# estimate columns added to madi_results.calib_standards in 0.7.0).
#
# R/calib_schema.sql holds the DDL; this file holds the same definition in a
# form code can check against (the worker before writing, I-SPI for help pages
# and exports). A test keeps the two in sync. Design: calib_points_SPEC.md.
# =============================================================================

# Estimate columns shared by test samples, controls and standards, in DDL order.
# (name, type, description)
.estimate_cols <- function() {
  e <- function(column, type, description)
    data.frame(column = column, type = type, description = description,
               stringsAsFactors = FALSE)
  rbind(
    e("predicted_concentration", "numeric", "Back-calculated concentration at the dilution the well was run at, on the scale the fit used for x (log10 when is_log_independent)."),
    e("predicted_log10_concentration", "numeric", "Back-calculated concentration, log10, at the dilution the well was run at."),
    e("final_concentration", "numeric", "Back-calculated concentration multiplied by the dilution (natural scale). For a standard: that point's estimate of the stock concentration."),
    e("se_concentration", "numeric", "Standard error of the back-calculated concentration (log10 scale when is_log_independent)."),
    e("pcov", "numeric", "Precision CV (%) of the back-calculated concentration, capped at cv_x_max."),
    e("pcov_rmse", "numeric", "CV (%) of the root-mean-square error. Standards: includes the bias against the nominal value. Samples and controls: equal to pcov (bias unknown)."),
    e("pcov_pass", "boolean", "TRUE if the point lies in the curve's dynamic range [LLOQ, ULOQ] (or passes the pcov threshold when LOQs are undefined)."),
    e("pcov_gate_class", "text", "belowLLOQ, inDynamicRange or aboveULOQ."),
    e("censor_class", "text", "quantified, lod_to_lloq, below_lod, below_lloq, above_uloq, saturated or no_response."),
    e("censor_type", "text", "none, left, interval or right."),
    e("censor_lower", "double precision", "Lower censoring bound, natural scale x dilution; NULL if unbounded."),
    e("censor_upper", "double precision", "Upper censoring bound, natural scale x dilution; NULL if unbounded."),
    e("conc_q_lo", "double precision", "Lower quantile of the concentration draws (x dilution)."),
    e("conc_q_med", "double precision", "Median of the concentration draws (x dilution)."),
    e("conc_q_hi", "double precision", "Upper quantile of the concentration draws (x dilution)."),
    e("p_below_lloq", "double precision", "Share of draws below the LLOQ."),
    e("p_above_uloq", "double precision", "Share of draws above the ULOQ."),
    e("p_below_cutoff", "double precision", "Share of draws below the decision cutoff."),
    e("frac_draws_below_a", "double precision", "Share of draws whose response lies below the lower asymptote (not invertible)."),
    e("frac_draws_above_d", "double precision", "Share of draws whose response lies above the upper asymptote (not invertible).")
  )
}

.point_contract_rows <- function() {
  r <- function(table, column, type, key, nullable, added_in, description)
    data.frame(table = table, column = column, type = type, key = key,
               nullable = nullable, added_in = added_in, description = description,
               stringsAsFactors = FALSE)
  est <- .estimate_cols()
  C <- "controls"; S <- "standards"
  ctrl <- rbind(
    r(C, "curve_id", "bigint", TRUE, FALSE, "0.7.0", "Curve the control was read on (FK to curve_lookup)."),
    r(C, "method", "varchar(20)", TRUE, FALSE, "0.7.0", "Calibration engine: bayesian or frequentist."),
    r(C, "well", "varchar(20)", TRUE, FALSE, "0.7.0", "Plate well of the control."),
    r(C, "sampleid", "varchar(128)", FALSE, TRUE, "0.7.0", "Control sample identifier as imported (often only the well)."),
    r(C, "source", "varchar(64)", FALSE, TRUE, "0.7.0", "Control name, e.g. 'QC1 (Low)'; identifies one control across plates and studies."),
    r(C, "dilution", "numeric", FALSE, TRUE, "0.7.0", "Dilution factor the control was run at."),
    r(C, "included", "boolean", FALSE, FALSE, "0.7.0", "FALSE if the well is masked. Masked controls are still estimated."),
    r(C, "exclusion_reason", "varchar(20)", FALSE, FALSE, "0.7.0", "none or masked."),
    r(C, "mask_reason", "varchar", FALSE, TRUE, "0.7.0", "Why the well was masked, as recorded at import or in the app."),
    r(C, "concentration", "numeric", FALSE, TRUE, "0.7.0", "Expected (nominal) concentration on the fit's x scale. NULL until reference values exist; see calib_points_SPEC.md."),
    r(C, "log10_concentration", "numeric", FALSE, TRUE, "0.7.0", "Expected concentration, log10. NULL until reference values exist."),
    r(C, "response_model", "numeric", FALSE, TRUE, "0.7.0", "Response on the scale the curve was fitted (after any log transform and blank handling)."),
    r(C, "assay_response_raw", "numeric", FALSE, TRUE, "0.7.0", "Response as measured.")
  )
  ctrl <- rbind(ctrl, r(C, est$column, est$type, FALSE, TRUE, "0.7.0", est$description),
                r(C, "job_id", "varchar(64)", FALSE, TRUE, "0.7.0", "Compute job that wrote the row."),
                r(C, "created_at", "timestamptz", FALSE, TRUE, "0.7.0", "Time the row was written."))

  std <- rbind(
    r(S, "curve_id", "bigint", TRUE, FALSE, "0.4", "Curve the standard belongs to."),
    r(S, "method", "varchar(20)", TRUE, FALSE, "0.4", "Calibration engine: bayesian or frequentist."),
    r(S, "well", "varchar(16)", TRUE, FALSE, "0.4", "Plate well of the standard."),
    r(S, "dilution", "numeric", TRUE, FALSE, "0.4", "Dilution factor of the standard point."),
    r(S, "concentration", "numeric", FALSE, TRUE, "0.4", "Nominal concentration at that dilution, natural scale."),
    r(S, "log10_concentration", "numeric", FALSE, TRUE, "0.4", "Nominal concentration, log10 (the x the fit used when is_log_independent)."),
    r(S, "response_model", "numeric", FALSE, TRUE, "0.4", "Response on the scale the curve was fitted."),
    r(S, "assay_response_raw", "numeric", FALSE, TRUE, "0.4", "Response as measured."),
    r(S, "included", "boolean", FALSE, FALSE, "0.4", "TRUE if the point entered the fit. Recovery of included points is optimistic (the curve was fitted to them)."),
    r(S, "exclusion_reason", "varchar(20)", FALSE, FALSE, "0.4", "none or masked."),
    r(S, "job_id", "text", FALSE, TRUE, "0.4", "Compute job that wrote the row."),
    r(S, "created_at", "timestamptz", FALSE, FALSE, "0.4", "Time the row was written."),
    r(S, "mask_reason", "varchar", FALSE, TRUE, "0.4", "Why the point was masked.")
  )
  std <- rbind(std, r(S, est$column, est$type, FALSE, TRUE, "0.7.0", est$description),
               r(S, "recovery_pct", "double precision", FALSE, TRUE, "0.7.0",
                 "100 x estimated / nominal concentration at the point's dilution."))
  rbind(ctrl, std)
}


#' Data contract for control and calibrator estimates
#'
#' Describes `madi_results.calib_controls` (one row per control well, curve and
#' method; new in 0.7.0) and `madi_results.calib_standards` (one row per standard
#' point, curve and method; estimate columns added in 0.7.0): column names, types,
#' keys, the version that added each column, and a short plain-language
#' description. `R/calib_schema.sql` holds the DDL; a package test keeps the two in
#' sync. Design: `calib_points_SPEC.md`.
#'
#' The estimate columns mean the same as in `calib_samples`, with two additions
#' for standards, whose true value is known: `recovery_pct`, and a `pcov_rmse`
#' that includes the bias against the nominal value.
#'
#' @param table `"controls"`, `"standards"` or `"all"` (default).
#' @param estimates_only If `TRUE`, keep only the modelled-estimate columns (those
#'   shared with `calib_samples`, plus `recovery_pct`).
#' @return A data frame with columns `table`, `column`, `type`, `key`, `nullable`,
#'   `added_in`, `description`.
#' @seealso [validate_point_estimates()], [weights_contract()]
#' @examples
#' point_estimates_contract("standards", estimates_only = TRUE)[, c("column", "description")]
#' @export
point_estimates_contract <- function(table = c("all", "controls", "standards"),
                                     estimates_only = FALSE) {
  table <- match.arg(table)
  x <- .point_contract_rows()
  if (table != "all") x <- x[x$table == table, , drop = FALSE]
  if (isTRUE(estimates_only))
    x <- x[x$column %in% c(.estimate_cols()$column, "recovery_pct"), , drop = FALSE]
  rownames(x) <- NULL
  x
}


#' Check control or standard rows against the contract before writing them
#'
#' Verifies that a data frame destined for `calib_controls` or `calib_standards`
#' contains no column the contract does not know, has every key and non-nullable
#' column filled, and that the controlled vocabularies (`exclusion_reason`,
#' `pcov_gate_class`, `censor_class`, `censor_type`) hold only allowed values. Stops
#' with a message listing the problems, so a bad batch is caught before it reaches
#' the database. `created_at` is filled by the database and may be absent.
#'
#' @param x Data frame of rows to be written.
#' @param table `"controls"` or `"standards"`.
#' @return `x`, invisibly, if it conforms.
#' @seealso [point_estimates_contract()]
#' @export
validate_point_estimates <- function(x, table = c("controls", "standards")) {
  table <- match.arg(table)
  ct <- point_estimates_contract(table)
  problems <- character()

  unknown <- setdiff(names(x), ct$column)
  if (length(unknown))
    problems <- c(problems, paste0("unknown column(s): ", paste(unknown, collapse = ", ")))

  need <- setdiff(ct$column[ct$key | !ct$nullable], "created_at")
  miss <- setdiff(need, names(x))
  if (length(miss))
    problems <- c(problems, paste0("missing required column(s): ", paste(miss, collapse = ", ")))
  for (cn in intersect(need, names(x)))
    if (anyNA(x[[cn]]))
      problems <- c(problems, paste0("required column '", cn, "' has missing values"))

  vocab <- list(
    exclusion_reason = c("none", "masked", "prozone", "lod"),
    pcov_gate_class  = c("belowLLOQ", "inDynamicRange", "aboveULOQ"),
    censor_class     = c("quantified", "lod_to_lloq", "below_lod", "below_lloq",
                         "above_uloq", "saturated", "no_response"),
    censor_type      = c("none", "left", "interval", "right"))
  for (cn in intersect(names(vocab), names(x))) {
    v <- unique(as.character(x[[cn]]))
    bad <- setdiff(v[!is.na(v)], vocab[[cn]])
    if (length(bad))
      problems <- c(problems, paste0("column '", cn, "' holds value(s) not allowed: ",
                                     paste(bad, collapse = ", ")))
  }

  if (anyDuplicated(x[intersect(ct$column[ct$key], names(x))]) && !length(miss))
    problems <- c(problems, "duplicate primary-key rows")

  if (length(problems))
    stop("validate_point_estimates(", table, "): ",
         paste(problems, collapse = "; "), call. = FALSE)
  invisible(x)
}


#' Recovery and bias-aware precision of back-calculated standards
#'
#' For calibrators, whose nominal concentration is known, computes the recovery
#' (estimated / nominal, in percent) and the CV of the root-mean-square error,
#' which adds the bias against the nominal value to the standard error:
#' `sqrt(se^2 + bias^2)` on the log10 scale, expressed as a CV in percent the same
#' way `pcov` is, and capped at `cv_x_max`.
#'
#' @param predicted Back-calculated concentration on the fit's x scale
#'   (`predicted_concentration`).
#' @param nominal Nominal concentration on the same scale (the x value the fit
#'   used: log10 when `is_log_independent`).
#' @param se Standard error of `predicted`, same scale (`se_concentration`).
#' @param is_log_independent Logical: is x on the log10 scale?
#' @param cv_x_max Cap for the CV, as for `pcov` (default 150).
#' @return A data frame with `recovery_pct` and `pcov_rmse`, one row per point.
#'   Both are `NA` where the estimate or the nominal value is missing; a point
#'   with an estimate but no finite `se` gets `pcov_rmse = cv_x_max`.
#' @examples
#' standard_recovery(predicted = c(1.02, 2.10), nominal = c(1, 2), se = c(0.03, 0.05))
#' @export
standard_recovery <- function(predicted, nominal, se, is_log_independent = TRUE,
                              cv_x_max = 150) {
  n <- length(predicted)
  if (length(nominal) != n || length(se) != n)
    stop("standard_recovery(): predicted, nominal and se must be the same length.")
  predicted <- suppressWarnings(as.numeric(predicted))
  nominal   <- suppressWarnings(as.numeric(nominal))
  se        <- suppressWarnings(as.numeric(se))

  if (isTRUE(is_log_independent)) {
    pl10 <- predicted; nl10 <- nominal; se10 <- se
  } else {
    ok <- is.finite(predicted) & predicted > 0
    pl10 <- ifelse(ok, log10(pmax(predicted, 1e-300)), NA_real_)
    nl10 <- ifelse(is.finite(nominal) & nominal > 0, log10(pmax(nominal, 1e-300)), NA_real_)
    # delta method: sd(log10 x) = sd(x) / (x ln 10)
    se10 <- ifelse(ok, se / (predicted * log(10)), NA_real_)
  }
  have <- is.finite(pl10) & is.finite(nl10)
  bias <- pl10 - nl10

  recovery <- ifelse(have, 100 * 10^bias, NA_real_)
  rmse_cv  <- ifelse(have & is.finite(se10),
                     pmin(sqrt(se10^2 + bias^2) * log(10) * 100, cv_x_max),
                     ifelse(have, cv_x_max, NA_real_))
  data.frame(recovery_pct = recovery, pcov_rmse = rmse_cv)
}


#' Add recovery and bias-aware pcov_rmse to every curve's standard estimates
#'
#' Call after the engine has back-calculated the standards
#' (`$standard_estimates` on each curve's result). Overwrites `pcov_rmse` and adds
#' `recovery_pct` there, using [standard_recovery()]. The nominal value is the
#' `concentration` column (the x value the fit used). A result without
#' `$standard_estimates` is returned unchanged.
#'
#' @param mp A `calibration_result` or `calibration_result_multiplate`.
#' @param is_log_independent Logical; defaults to the value recorded in the
#'   result's `meta`, else `TRUE`.
#' @return The same class of object.
#' @export
add_standard_recovery_multiplate <- function(mp, is_log_independent = NULL) {
  is_multi <- inherits(mp, "calibration_result_multiplate")
  plates <- if (is_multi) mp$plates
            else stats::setNames(list(mp), as.character(mp$meta$curve_id %||% "1"))
  for (cid in names(plates)) {
    cr <- plates[[cid]]
    se <- cr$standard_estimates
    if (is.null(se) || !is.data.frame(se) || !nrow(se)) next
    lx <- is_log_independent %||% cr$meta$is_log_independent %||%
          (if (is_multi) mp$meta$is_log_independent else NULL) %||% TRUE
    rc <- standard_recovery(
      predicted = se$predicted_concentration %||% rep(NA_real_, nrow(se)),
      nominal   = se$concentration %||% rep(NA_real_, nrow(se)),
      se        = se$se_concentration %||% rep(NA_real_, nrow(se)),
      is_log_independent = lx,
      cv_x_max = cr$meta$cv_x_max %||% 150)
    se$recovery_pct <- rc$recovery_pct
    se$pcov_rmse    <- rc$pcov_rmse
    cr$standard_estimates <- se
    plates[[cid]] <- cr
  }
  if (is_multi) { mp$plates <- plates; mp } else plates[[1]]
}
