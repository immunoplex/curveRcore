# =============================================================================
# response_adjustment.R — one rule for blank subtraction and values <= 0
#
# The response adjustment is a set-level statistic of one fit group (all plates
# of a multiplate group). It is computed ONCE and applied identically to the
# standards, the controls and the test samples, so all three sit on the same
# response axis:
#
#   adjusted = raw - blank_value + offset          (then log10 if log response)
#
# * blank_value — factor x blank centre when blank_option is a subtraction,
#   else 0. Blank subtraction is opt-in and never the default.
# * offset — only on the log-response path. 0 when every included calibrator,
#   control and test-sample response (after any subtraction) is > 0. Otherwise
#   -min + delta (delta = 0.001 response units by default), so the smallest
#   value lands just above zero and every well is kept and estimated.
#   Blanks never decide the offset: they are not calibrators.
# * display_floor — for rows outside the offset set that are still <= 0 after
#   the adjustment (blanks, masked rows, or callers that did not pass every
#   sample): floored at 1% of the smallest included positive adjusted standard
#   response. Identical to the pre-0.8.0 floor whenever offset = 0.
# =============================================================================


# Centre of the included blank responses: geometric mean when every blank is
# > 0 (unchanged behaviour), arithmetic mean when any blank is <= 0 (a geometric
# mean would silently drop those blanks and bias the centre upward).
.blank_centre <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  if (all(x > 0)) geom_mean(x) else mean(x)
}

.blank_factor <- function(blank_option) {
  switch(blank_option,
         subtracted     = 1,
         subtracted_3x  = 3,
         subtracted_10x = 10,
         0)
}


#' Compute the Response Adjustment for One Fit Group
#'
#' Determines the blank-subtraction value and the log-scale offset that
#' [preprocess_standards()] and [transform_sample_response()] apply, so the
#' standards, controls and test samples of a fit group share one response axis.
#'
#' **Offset rule (log response only).** If every included standard, control
#' and test-sample response (after any blank subtraction) is > 0, the offset is
#' 0 and nothing changes. If any is <= 0, `offset = -min + delta` is added to
#' every value, so the smallest lands `delta` above zero. Blanks never decide
#' the offset. The offset is recorded with the fit.
#'
#' @param standards Data frame of standards (raw response).
#' @param response_variable Character. Response column name (same in all frames).
#' @param is_log_response Logical. Offset is only applied on the log path.
#' @param samples,controls Data frames or NULL. Test samples and controls (raw
#'   response). Rows flagged excluded by `include_col` are ignored.
#' @param blank_data Data frame or NULL. Blanks (only for the subtraction value).
#' @param blank_option Character. See [perform_blank_operation()].
#' @param include_col Character. Logical include/mask column; absent = all rows.
#' @param delta Numeric. Margin above zero for the smallest value. Default 0.001
#'   (response units).
#'
#' @return Named list: `blank_value`, `blank_centre`, `offset`, `delta`,
#'   `n_nonpositive` (values <= 0 in the offset set), `min_response`,
#'   `is_log_response`, `blank_option`.
#' @export
compute_response_adjustment <- function(standards, response_variable,
                                        is_log_response,
                                        samples      = NULL,
                                        controls     = NULL,
                                        blank_data   = NULL,
                                        blank_option = "ignored",
                                        include_col  = "included",
                                        delta        = 0.001) {
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta <= 0)
    stop("delta must be a single positive number")

  blank_centre <- NA_real_
  blank_value  <- 0
  if (!is.null(blank_data) && nrow(blank_data) > 0 &&
      response_variable %in% names(blank_data)) {
    binc <- .included_mask(blank_data, include_col)
    blank_centre <- .blank_centre(blank_data[[response_variable]][binc])
  }
  fac <- .blank_factor(blank_option)
  if (fac > 0 && is.finite(blank_centre)) blank_value <- fac * blank_centre

  pick <- function(df) {
    if (is.null(df) || !nrow(df) || !(response_variable %in% names(df)))
      return(numeric(0))
    v <- df[[response_variable]][.included_mask(df, include_col)]
    v[is.finite(v)]
  }
  vals <- c(pick(standards), pick(controls), pick(samples)) - blank_value
  # "included" blank option adds the blank centre as a calibrator point.
  if (identical(blank_option, "included") && is.finite(blank_centre))
    vals <- c(vals, blank_centre)

  # Within rounding error of zero counts as zero (e.g. 0.01 - mean(-0.01, 0.03)
  # is ~1e-18, which would otherwise land at log10 = -17.8).
  tol      <- if (length(vals)) 1e-12 * max(abs(vals)) else 0
  n_nonpos <- sum(vals <= tol)
  min_resp <- if (length(vals)) min(vals) else NA_real_
  offset <- if (isTRUE(is_log_response) && n_nonpos > 0) -min_resp + delta else 0

  list(blank_value     = blank_value,
       blank_centre    = blank_centre,
       offset          = offset,
       delta           = delta,
       n_nonpositive   = as.integer(n_nonpos),
       min_response    = min_resp,
       is_log_response = isTRUE(is_log_response),
       blank_option    = blank_option)
}


#' Put Raw Responses on the Fitting Scale
#'
#' The single transform for test samples and controls, matching what
#' [preprocess_standards()] does to the standards: subtract the blank value
#' (if any), add the offset, then log10 on the log path. Values still <= 0
#' after the adjustment are set to `display_floor`.
#'
#' @param raw Numeric. Raw responses.
#' @param is_log_response Logical.
#' @param adjustment List from [compute_response_adjustment()] (or the
#'   `derived_stats` of [preprocess_standards()]), or NULL. NULL keeps the
#'   legacy behaviour (no subtraction, no offset, floor 1e-6) for callers that
#'   have not been given an adjustment.
#' @param display_floor Numeric or NULL. Floor for values still <= 0; defaults
#'   to `adjustment$response_floor`, else 1e-6.
#'
#' @return Numeric vector on the fitting scale.
#' @export
transform_sample_response <- function(raw, is_log_response, adjustment = NULL,
                                       display_floor = NULL) {
  raw <- suppressWarnings(as.numeric(raw))
  bv  <- adjustment$blank_value %||% 0
  off <- if (isTRUE(is_log_response)) adjustment$offset %||% 0 else 0
  y   <- raw - bv + off
  if (!isTRUE(is_log_response)) return(y)
  fl <- display_floor %||% adjustment$response_floor %||% 1e-6
  if (!is.finite(fl) || fl <= 0) fl <- 1e-6
  bad <- !is.na(y) & y <= 0
  y[bad] <- fl
  log10(y)
}


#' Back-transform Fitted Responses to the Raw Scale
#'
#' Inverse of [transform_sample_response()] for model-scale responses
#' (e.g. predicted curves, asymptotes): `10^y - offset + blank_value` on the
#' log path, `y + blank_value` otherwise.
#'
#' @param y Numeric. Responses on the fitting scale.
#' @param is_log_response Logical.
#' @param adjustment List or NULL (NULL = no subtraction, no offset).
#' @return Numeric vector on the raw response scale.
#' @export
untransform_response <- function(y, is_log_response, adjustment = NULL) {
  bv  <- adjustment$blank_value %||% 0
  off <- if (isTRUE(is_log_response)) adjustment$offset %||% 0 else 0
  if (isTRUE(is_log_response)) 10^y - off + bv else y + bv
}
