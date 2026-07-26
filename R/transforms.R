# =============================================================================
# transforms.R — Data transformation utilities
#
# Upstream transforms applied to raw assay data before model fitting.
# Both curveRfreq and curveRbayes call these identically, ensuring the
# data handed to each fitting engine is on the same scale.
# =============================================================================


# Resolve the logical "included in fit" mask for a frame.
# If `include_col` is absent (or `data` is NULL/empty) every row is treated as
# included, so mask-aware functions are fully backward compatible with callers
# that never supplied the column.
#' @noRd
.included_mask <- function(data, include_col = "included") {
  if (is.null(data) || nrow(data) == 0L) return(logical(0))
  if (!is.null(include_col) && include_col %in% names(data)) {
    inc <- as.logical(data[[include_col]])
    inc[is.na(inc)] <- FALSE          # NA status is treated as masked
    inc
  } else {
    rep(TRUE, nrow(data))
  }
}


#' Compute Concentration from Dilution and Undiluted Standard
#'
#' Converts a `dilution` column to absolute concentration (optionally
#' log10-transformed) and stores it in `independent_variable`.
#'
#' @param data Data frame with a `dilution` column.
#' @param undiluted_sc_concentration Numeric. Concentration of the undiluted
#'   standard (e.g. 10000).
#' @param independent_variable Character. Column name for the output.
#' @param is_log_concentration Logical. Apply log10 after computing
#'   concentration? Default `TRUE`.
#'
#' @return `data` with the concentration column populated.
#' @export
compute_concentration <- function(data,
                                  undiluted_sc_concentration,
                                  independent_variable,
                                  is_log_concentration = TRUE) {
  independent_variable <- unique(independent_variable)
  data[[independent_variable]] <- (1 / data$dilution) * undiluted_sc_concentration

  if (is_log_concentration) {
    data[[independent_variable]] <- log10(data[[independent_variable]])
  }
  data
}


#' Log10-Transform the Assay Response (mask-aware)
#'
#' Applies `log10()` to the response column when `is_log_response` is TRUE.
#' Non-positive values are floored to an adaptive minimum before transform.
#'
#' The adaptive floor is a **set-level statistic**: it is derived from the
#' *included* positive responses only (`data[included & value > 0]`), but the
#' floor and the `log10()` are then applied to **all** rows. This is what lets
#' masked points land on the same response axis as the fitted points without
#' influencing the transform (see [preprocess_standards()]). Pass an explicit
#' `floor_value` to reuse a floor computed elsewhere (e.g. to transform blanks
#' on the same axis as the standards).
#'
#' @param data Data frame.
#' @param response_variable Character. Name of the response column.
#' @param is_log_response Logical. Apply log10? Default `TRUE`.
#' @param include_col Character. Name of the logical column marking rows that
#'   entered the fit (`TRUE`) versus masked rows (`FALSE`). If the column is
#'   absent, every row is treated as included (backward compatible).
#' @param floor_value Numeric or NULL. If supplied, this exact floor is used
#'   for non-positive values instead of deriving one. Used to share a single
#'   floor across the standards and blanks frames.
#' @param floor_method Character. How to derive the floor when `floor_value`
#'   is NULL: `"adaptive"` (default) uses 1\% of the minimum *included*
#'   positive value; `"fixed"` uses 1e-6.
#' @param verbose Logical. Emit messages about floored values.
#'
#' @return `data` with the response column optionally transformed. The floor
#'   actually used is attached as `attr(., "response_floor")`.
#' @export
compute_log_response <- function(data, response_variable,
                                 is_log_response = TRUE,
                                 include_col  = "included",
                                 floor_value  = NULL,
                                 floor_method = "adaptive",
                                 verbose = FALSE) {
  if (!is_log_response) return(data)

  inc      <- .included_mask(data, include_col)
  raw_vals <- data[[response_variable]]

  adaptive_floor <- if (!is.null(floor_value)) {
    floor_value
  } else {
    # Set-level statistic: derived from INCLUDED positives only.
    positive_incl <- raw_vals[inc & is.finite(raw_vals) & raw_vals > 0]
    if (floor_method == "adaptive" && length(positive_incl) > 0) {
      min(positive_incl) * 0.01
    } else {
      1e-6
    }
  }

  # Applied to ALL rows (included and masked alike).
  bad <- is.na(raw_vals) | raw_vals <= 0
  n_floored <- sum(bad, na.rm = TRUE)
  if (n_floored > 0 && verbose) {
    message(sprintf("[compute_log_response] %d/%d values <= 0 floored to %.2e before log10",
                    n_floored, length(raw_vals), adaptive_floor))
  }
  raw_vals[bad] <- adaptive_floor
  data[[response_variable]] <- log10(raw_vals)
  attr(data, "response_floor") <- adaptive_floor
  data
}


# ---- Prozone (hook effect) correction ----

#' Correct the Prozone (Hook) Effect
#'
#' At very high concentrations the measured signal can decrease (hook effect).
#' This function compresses the post-peak delta toward the peak value.
#'
#' The peak (`max_response` and `logc_at_max`) is a **set-level statistic**
#' computed from the *included* points only. The post-peak reflection is then
#' applied to **all** rows relative to that peak, so masked points beyond the
#' hook are dampened onto the same reference as the fitted points. Rows are
#' never dropped (grain is preserved); rows with a missing response or
#' concentration are passed through untouched.
#'
#' @param stdframe Data frame of standard curve data.
#' @param prop_diff Numeric. Dampening factor (e.g. 0.1).
#' @param dil_scale Numeric. Dilution scale factor (e.g. 2).
#' @param response_variable Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#' @param include_col Character. Logical column marking fitted rows. Absent =
#'   all rows included (backward compatible).
#' @param verbose Logical.
#'
#' @return `stdframe` (all rows) with post-peak response values adjusted. The
#'   peak reference is attached as `attr(., "prozone_peak_response")` and
#'   `attr(., "prozone_logc_at_peak")`.
#' @export
correct_prozone <- function(stdframe,
                            prop_diff = 0.1,
                            dil_scale = 2,
                            response_variable  = "mfi",
                            independent_variable = "concentration",
                            include_col = "included",
                            verbose = FALSE) {
  response_variable <- unique(response_variable)

  inc   <- .included_mask(stdframe, include_col)
  resp  <- stdframe[[response_variable]]
  indep <- stdframe[[independent_variable]]

  # Peak reference from INCLUDED, finite points only.
  ref <- inc & is.finite(resp) & is.finite(indep)
  if (!any(ref)) return(stdframe)

  max_response <- max(resp[ref])
  logc_at_max  <- max(indep[ref & resp == max_response])

  if (verbose) message("Peak (from included) = ", max_response, " at x = ", logc_at_max)

  attr(stdframe, "prozone_peak_response") <- max_response
  attr(stdframe, "prozone_logc_at_peak")  <- logc_at_max

  # Reflection applied to ALL finite post-peak rows (included and masked).
  post_peak <- is.finite(indep) & is.finite(resp) & indep > logc_at_max
  if (!any(post_peak)) return(stdframe)

  pp_resp <- resp[post_peak]
  pp_conc <- indep[post_peak]

  stdframe[[response_variable]][post_peak] <-
    max_response + (max_response - pp_resp) * prop_diff /
    ((pp_conc - logc_at_max) * dil_scale)

  stdframe
}


# ---- Blank handling ----

#' Apply a Blank Operation to Standard Curve Data
#'
#' Performs one of five blank-handling strategies:
#' * `"ignored"` — no adjustment (default)
#' * `"included"` — append blank geometric mean as an extra point
#' * `"subtracted"` — subtract geometric mean of blanks
#' * `"subtracted_3x"` — subtract 3× geometric mean
#' * `"subtracted_10x"` — subtract 10× geometric mean
#'
#' After subtraction, values that become non-positive are floored at
#' 0 (linear) or 1 (log scale).
#'
#' @param blank_data Data frame of blank measurements, or NULL.
#' @param data Data frame of standard curve data.
#' @param response_variable Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#' @param is_log_response Logical. Whether the response has been log10-transformed.
#' @param blank_option Character. One of the five options above.
#' @param include_col Character. Logical column marking fitted rows on
#'   `blank_data` (and, for the `"included"` option, on `data`). Absent = all
#'   rows included (backward compatible).
#' @param blank_mean Numeric or NULL. Pre-computed geometric mean of the
#'   *included* blanks. When NULL it is computed here from
#'   `blank_data[included, response_variable]`. Supplying it lets a caller
#'   guarantee the value used matches the one it records in `derived_stats`.
#' @param verbose Logical.
#'
#' @return `data` with the blank operation applied. Blanks are *not* modified
#'   here; [preprocess_standards()] transforms and returns them separately so
#'   they share the standards' response floor.
#' @export
perform_blank_operation <- function(blank_data, data, response_variable,
                                    independent_variable, is_log_response,
                                    blank_option = "ignored",
                                    include_col  = "included",
                                    blank_mean   = NULL,
                                    verbose = FALSE) {
  valid_options <- c("ignored", "included", "subtracted",
                     "subtracted_3x", "subtracted_10x")
  if (!(blank_option %in% valid_options)) {
    stop("blank_option must be one of: ", paste(valid_options, collapse = ", "))
  }

  if (blank_option == "ignored") return(data)

  if (is.null(blank_data) || nrow(blank_data) == 0) {
    warning("Blank data required when blank_option != 'ignored'. Returning data unchanged.")
    return(data)
  }

  if (blank_option == "included") {
    data <- include_blanks_conc(blank_data, data, response_variable,
                                independent_variable, include_col = include_col)
    if (verbose) message("Blank geometric mean (included blanks) added as extra point.")
    return(data)
  }

  # Subtraction variants
  factor <- switch(blank_option,
                   "subtracted"     = 1,
                   "subtracted_3x"  = 3,
                   "subtracted_10x" = 10)

  # Set-level statistic: geometric mean of INCLUDED blanks only.
  if (is.null(blank_mean)) {
    binc       <- .included_mask(blank_data, include_col)
    blank_mean <- geom_mean(blank_data[[response_variable]][binc])
  }
  data[[response_variable]] <- data[[response_variable]] - factor * blank_mean

  floor_val <- if (is_log_response) 1 else 0
  data[[response_variable]] <- pmax(data[[response_variable]], floor_val)

  if (verbose) message("Blank subtraction (x", factor, ") applied.")
  data
}


#' Include Blank Controls as an Extra Standard Curve Point
#'
#' Appends a synthetic row whose response is the geometric mean of the
#' *included* blanks and whose concentration is `log10(2)` below the minimum
#' *included* standard concentration. The appended row is itself a fit point,
#' so its `included` flag is set to `TRUE`.
#'
#' @param blank_data Data frame of blank wells.
#' @param data Data frame of standards.
#' @param response_variable Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#' @param include_col Character. Logical column marking fitted rows on both
#'   frames. Absent = all rows included.
#'
#' @return `data` with one additional row.
#' @keywords internal
include_blanks_conc <- function(blank_data, data, response_variable,
                                independent_variable = "concentration",
                                include_col = "included") {
  binc <- .included_mask(blank_data, include_col)
  sinc <- .included_mask(data, include_col)

  response_blank    <- geom_mean(blank_data[[response_variable]][binc])
  min_concentration <- min(data[[independent_variable]][sinc], na.rm = TRUE)
  conc_blank        <- min_concentration - log10(2)

  new_row <- data[1, , drop = FALSE]
  new_row[1, ] <- NA
  new_row[[response_variable]]    <- response_blank
  new_row[[independent_variable]] <- conc_blank
  if ("dilution" %in% names(new_row))  new_row$dilution  <- NA_real_
  if ("stype"    %in% names(new_row))  new_row$stype     <- "B"
  if ("well"     %in% names(new_row))  new_row$well      <- "blank_mean"
  if ("sampleid" %in% names(new_row))  new_row$sampleid  <- "blank_mean"
  # The synthetic point enters the fit.
  if (include_col %in% names(new_row)) new_row[[include_col]] <- TRUE

  rbind(data, new_row)
}


# ---- Lower asymptote constraint resolution ----

#' Resolve the Fixed Lower Asymptote Value
#'
#' Determines whether the lower asymptote `a` should be fixed based on
#' the constraint method. Returns the fixed value (on the **raw** scale,
#' before any log-transform) or NULL if `a` should be estimated freely.
#'
#' @param l_asy_constraints Named list. Must contain
#'   `l_asy_constraint_method` and, for non-default methods,
#'   `l_asy_min_constraint` and/or `l_asy_max_constraint`.
#'   Typically an `antigen_constraints` object or a plain list.
#'
#' @return Numeric scalar (the fixed value on the raw scale) or NULL
#'   (parameter is free).
#'
#' @details
#' The constraint method controls behaviour:
#' \describe{
#'   \item{`"default"`}{`a` is always free (returns NULL).}
#'   \item{`"user_defined"`}{`a` is fixed at `l_asy_min_constraint`
#'     (which should equal `l_asy_max_constraint`).}
#'   \item{`"range_of_blanks"`}{`a` is fixed at `l_asy_min_constraint`.
#'     The caller is responsible for computing this from blank data
#'     upstream.}
#'   \item{`"geometric_mean_of_blanks"`}{`a` is fixed at
#'     `l_asy_min_constraint`. The caller is responsible for computing
#'     this from blank data upstream.}
#' }
#'
#' @export
resolve_fixed_lower_asymptote <- function(l_asy_constraints) {
  method <- l_asy_constraints$l_asy_constraint_method %||% "default"

  if (method == "default") {
    return(NULL)
  }

  # For all non-default methods, return the constraint value
  l_asy_constraints$l_asy_min_constraint
}


#' Validate a Fixed Lower Asymptote Before Log Transformation
#'
#' Checks that the value is a positive, finite scalar suitable for `log10()`.
#' Returns the value if valid, NULL otherwise.
#'
#' @param fixed_a_result_raw Numeric scalar or NULL.
#' @param verbose Logical.
#' @return `fixed_a_result_raw` if valid; NULL otherwise.
#' @export
validate_fixed_lower_asymptote <- function(fixed_a_result_raw, verbose = FALSE) {
  if (is.null(fixed_a_result_raw)) return(NULL)

  if (!is.numeric(fixed_a_result_raw) || length(fixed_a_result_raw) != 1 ||
      !is.finite(fixed_a_result_raw)) {
    if (verbose) message("[validate_fixed_lower_asymptote] Not a finite scalar; treating as NULL.")
    return(NULL)
  }

  if (fixed_a_result_raw <= 0) {
    if (verbose) message(sprintf("[validate_fixed_lower_asymptote] %.6f <= 0; log10 undefined. Treating as NULL.",
                                 fixed_a_result_raw))
    return(NULL)
  }

  fixed_a_result_raw
}


#' Full Preprocessing Pipeline for Standard Curve Data (mask-aware)
#'
#' Applies concentration computation, prozone correction, blank handling,
#' and optional log10 response transform in the canonical order.
#'
#' @details
#' **Mask-aware contract.** An `include_col` logical column (default
#' `"included"`, `TRUE` = used in the fit, `FALSE` = masked) may be present on
#' `data` and `blank_data`. Every *set-level statistic* — the prozone peak, the
#' blank geometric mean, the adaptive log floor, and the minimum-concentration
#' anchor — is computed from the **included** rows only. The resulting
#' transforms are then applied to **all** rows, so masked points land on the
#' same axes as the fitted points without ever influencing them. If the column
#' is absent, every row is treated as included and the output is identical to
#' the pre-mask behaviour (backward compatible).
#'
#' The function does **not** drop masked rows; downstream fitters are expected
#' to receive only the included subset (e.g. `pp$data[pp$data$included, ]`),
#' which keeps the fit byte-identical to a fit that never saw the masked rows.
#'
#' Both frames retain a pristine `assay_response_raw` column (the response
#' before prozone/blank/log), so callers can persist the raw and model-space
#' responses side by side.
#'
#' **Blanks are never subtracted automatically.** With the default
#' `blank_option = "ignored"` the standard responses are left untouched; the
#' returned `blanks` frame is transformed for display/persistence only and is
#' *not* subtracted from the standards. Subtraction happens **only** when the
#' caller explicitly selects `"subtracted"`, `"subtracted_3x"`, or
#' `"subtracted_10x"` (which subtract 1x/3x/10x the *included*-blank geometric
#' mean), or adds the blank mean as a point via `"included"`. The returned
#' blanks are always the raw and model-space blank responses, never a
#' subtracted quantity.
#'
#' @param data Data frame of standards with a `dilution` column.
#' @param antigen_settings Named list with `standard_curve_concentration`.
#' @param response_variable Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#' @param is_log_response Logical. Log10-transform the response?
#' @param blank_data Data frame of blanks, or NULL.
#' @param blank_option Character. Blank handling method.
#' @param is_log_independent Logical. Log10-transform concentration?
#' @param apply_prozone Logical. Apply prozone correction?
#' @param include_col Character. Name of the logical include/mask column on
#'   `data` and `blank_data`. Absent on a frame = all its rows are included.
#' @param verbose Logical.
#'
#' @return A named list:
#'   \describe{
#'     \item{`data`}{All standard rows, transformed, carrying `include_col`,
#'       the (log10) `concentration`, the model-space response, and
#'       `assay_response_raw`.}
#'     \item{`blanks`}{All blank rows, transformed onto the standards' response
#'       floor, carrying `include_col`, the model-space response, and
#'       `assay_response_raw`; `NULL` when no `blank_data` was supplied.}
#'     \item{`antigen_fit_options`}{Record of the options used.}
#'     \item{`derived_stats`}{The set-level statistics computed from the
#'       included points: `blank_geomean`, `prozone_peak_response`,
#'       `prozone_logc_at_peak`, `response_floor`, `min_included_concentration`.}
#'   }
#' @export
preprocess_standards <- function(data, antigen_settings, response_variable,
                                 independent_variable,
                                 is_log_response,
                                 blank_data         = NULL,
                                 blank_option       = "ignored",
                                 is_log_independent = TRUE,
                                 apply_prozone      = TRUE,
                                 include_col        = "included",
                                 verbose            = FALSE) {

  undiluted_sc <- antigen_settings$standard_curve_concentration

  # Ensure the include column exists on both frames so the returned frames
  # always carry an explicit status (absent -> all included).
  if (!is.null(data) && !(include_col %in% names(data))) {
    data[[include_col]] <- TRUE
  }
  if (!is.null(blank_data) && nrow(blank_data) > 0 &&
      !(include_col %in% names(blank_data))) {
    blank_data[[include_col]] <- TRUE
  }

  # Pristine raw response (before prozone / blank / log), on both frames.
  data[["assay_response_raw"]] <- data[[response_variable]]
  if (!is.null(blank_data) && nrow(blank_data) > 0) {
    blank_data[["assay_response_raw"]] <- blank_data[[response_variable]]
  }

  data <- compute_concentration(data, undiluted_sc, independent_variable,
                                is_log_concentration = is_log_independent)

  prozone_peak <- NA_real_
  prozone_logc <- NA_real_
  if (apply_prozone) {
    data <- correct_prozone(stdframe = data, prop_diff = 0.1, dil_scale = 2,
                            response_variable    = response_variable,
                            independent_variable = independent_variable,
                            include_col          = include_col,
                            verbose              = verbose)
    prozone_peak <- attr(data, "prozone_peak_response") %||% NA_real_
    prozone_logc <- attr(data, "prozone_logc_at_peak")  %||% NA_real_
    attr(data, "prozone_peak_response") <- NULL
    attr(data, "prozone_logc_at_peak")  <- NULL
  }

  # Blank geometric mean from INCLUDED blanks (used by subtraction / included
  # options and recorded in derived_stats so the value can't drift).
  blank_geomean <- NA_real_
  if (!is.null(blank_data) && nrow(blank_data) > 0) {
    binc          <- .included_mask(blank_data, include_col)
    blank_geomean <- geom_mean(blank_data[[response_variable]][binc])
  }

  data <- perform_blank_operation(blank_data = blank_data, data = data,
                                  response_variable    = response_variable,
                                  independent_variable = independent_variable,
                                  is_log_response      = is_log_response,
                                  blank_option         = blank_option,
                                  include_col          = include_col,
                                  blank_mean           = blank_geomean,
                                  verbose              = verbose)

  # Adaptive log floor: 1% of the minimum INCLUDED positive response (post
  # blank op), shared by the standards and the blanks so they land on the same
  # response axis.
  response_floor <- NA_real_
  if (is_log_response) {
    inc  <- .included_mask(data, include_col)
    resp <- data[[response_variable]]
    pos  <- resp[inc & is.finite(resp) & resp > 0]
    response_floor <- if (length(pos) > 0) min(pos) * 0.01 else 1e-6
  }

  data <- compute_log_response(data, response_variable,
                               is_log_response = is_log_response,
                               include_col     = include_col,
                               floor_value     = response_floor,
                               verbose         = verbose)
  attr(data, "response_floor") <- NULL

  # Transform and return the blanks on the SAME floor (blanks are not
  # blank-subtracted; they are only placed on the model-space response axis).
  blanks_out <- NULL
  if (!is.null(blank_data) && nrow(blank_data) > 0) {
    blanks_out <- compute_log_response(blank_data, response_variable,
                                       is_log_response = is_log_response,
                                       include_col     = include_col,
                                       floor_value     = response_floor,
                                       verbose         = FALSE)
    attr(blanks_out, "response_floor") <- NULL
  }

  # Minimum included concentration (the anchor include_blanks_conc uses).
  sinc <- .included_mask(data, include_col)
  cvals <- data[[independent_variable]][sinc]
  cvals <- cvals[is.finite(cvals)]
  min_included_conc <- if (length(cvals)) min(cvals) else NA_real_

  antigen_fit_options <- list(
    is_log_response      = is_log_response,
    blank_option         = blank_option,
    is_log_concentration = is_log_independent,
    apply_prozone        = apply_prozone
  )

  derived_stats <- list(
    blank_geomean             = blank_geomean,
    prozone_peak_response     = prozone_peak,
    prozone_logc_at_peak      = prozone_logc,
    response_floor            = response_floor,
    min_included_concentration = min_included_conc
  )

  list(data = data, blanks = blanks_out,
       antigen_fit_options = antigen_fit_options,
       derived_stats = derived_stats)
}


# =============================================================================
# Adaptive Constraint Profile
# =============================================================================

#' Build an Adaptive Constraint Profile from Observed Data
#'
#' Inspects the response range, dynamic range, and scale to choose
#' appropriate bounds for nonlinear optimisation. The returned profile
#' is consumed by per-model constraint builders in curveRfreq and
#' curveRbayes.
#'
#' Three scale classes are recognised:
#' \describe{
#'   \item{high}{MFI-like (log-max > 2.5 or raw max > 1000)}
#'   \item{medium}{Intermediate signals}
#'   \item{low}{OD/absorbance-like (narrow dynamic range)}
#' }
#' Narrower dynamic ranges receive wider slope and asymmetry bounds
#' to avoid near-singular Jacobians.
#'
#' @param data Data frame. Must contain response and `concentration` cols.
#' @param response_variable Character. Response column name.
#' @param is_log_response Logical. Is the response already log10-transformed?
#' @param antigen_settings List with `l_asy_min_constraint` and
#'   `l_asy_max_constraint`.
#'
#' @return A named list: `y_min`, `y_max`, `dynamic_range`, `conc_range`,
#'   `scale_class`, `slope_min`, `slope_max`, `g_min`, `g_max`,
#'   `conc_pad_frac`, `d_margin_frac`.
#'
#' @export
adaptive_constraint_profile <- function(data,
                                        response_variable,
                                        is_log_response,
                                        antigen_settings) {

  y_vals <- data[[response_variable]]
  y_vals <- y_vals[is.finite(y_vals)]
  y_min  <- min(y_vals)
  y_max  <- max(y_vals)
  dynamic_range <- y_max - y_min

  conc_vals  <- data$concentration[is.finite(data$concentration)]
  conc_range <- diff(range(conc_vals))

  # Classify the response scale
  if (is_log_response) {
    scale_class <- if (y_max > 2.5) "high"
                   else if (y_max > 0.5) "medium"
                   else "low"
  } else {
    scale_class <- if (y_max > 1000) "high"
                   else if (y_max > 10) "medium"
                   else "low"
  }

  # Adapt slope bounds to dynamic range
  slope_max <- switch(scale_class, high = 2.0, medium = 3.0, low = 5.0)
  slope_min <- switch(scale_class, high = 0.1, medium = 0.05, low = 0.01)

  # Asymmetry parameter g
  g_min <- switch(scale_class, high = 0.5, medium = 0.3, low = 0.1)
  g_max <- switch(scale_class, high = 5.0, medium = 7.0, low = 10.0)

  # How far beyond data range to allow inflection point c
  conc_pad_frac <- switch(scale_class, high = 0.5, medium = 0.7, low = 1.0)

  # Upper asymptote d: margin above/below y_max
  d_margin_frac <- switch(scale_class, high = 0.5, medium = 0.3, low = 0.15)

  list(
    y_min         = y_min,
    y_max         = y_max,
    dynamic_range = dynamic_range,
    conc_range    = conc_range,
    scale_class   = scale_class,
    slope_max     = slope_max,
    slope_min     = slope_min,
    g_min         = g_min,
    g_max         = g_max,
    conc_pad_frac = conc_pad_frac,
    d_margin_frac = d_margin_frac
  )
}
