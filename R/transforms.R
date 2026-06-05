# =============================================================================
# transforms.R — Data transformation utilities
#
# Upstream transforms applied to raw assay data before model fitting.
# Both curveRfreq and curveRbayes call these identically, ensuring the
# data handed to each fitting engine is on the same scale.
# =============================================================================


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


#' Log10-Transform the Assay Response
#'
#' Applies `log10()` to the response column when `is_log_response` is TRUE.
#' Non-positive values are floored to an adaptive minimum before transform.
#'
#' @param data Data frame.
#' @param response_variable Character. Name of the response column.
#' @param is_log_response Logical. Apply log10? Default `TRUE`.
#' @param floor_method Character. How to handle non-positive values before
#'   log10: `"adaptive"` (default) uses 1\% of the minimum positive value;
#'   `"fixed"` uses 1e-6.
#' @param verbose Logical. Emit messages about floored values.
#'
#' @return `data` with the response column optionally transformed.
#' @export
compute_log_response <- function(data, response_variable,
                                 is_log_response = TRUE,
                                 floor_method = "adaptive",
                                 verbose = FALSE) {
  if (!is_log_response) return(data)

  raw_vals <- data[[response_variable]]
  positive_vals <- raw_vals[is.finite(raw_vals) & raw_vals > 0]

  adaptive_floor <- if (floor_method == "adaptive" && length(positive_vals) > 0) {
    min(positive_vals) * 0.01
  } else {
    1e-6
  }

  bad <- is.na(raw_vals) | raw_vals <= 0
  n_floored <- sum(bad, na.rm = TRUE)
  if (n_floored > 0 && verbose) {
    message(sprintf("[compute_log_response] %d/%d values <= 0 floored to %.2e before log10",
                    n_floored, length(raw_vals), adaptive_floor))
  }
  data[[response_variable]][bad] <- adaptive_floor
  data[[response_variable]] <- log10(data[[response_variable]])
  data
}


# ---- Prozone (hook effect) correction ----

#' Correct the Prozone (Hook) Effect
#'
#' At very high concentrations the measured signal can decrease (hook effect).
#' This function compresses the post-peak delta toward the peak value.
#'
#' @param stdframe Data frame of standard curve data.
#' @param prop_diff Numeric. Dampening factor (e.g. 0.1).
#' @param dil_scale Numeric. Dilution scale factor (e.g. 2).
#' @param response_variable Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#' @param verbose Logical.
#'
#' @return `stdframe` with post-peak response values adjusted.
#' @export
correct_prozone <- function(stdframe,
                            prop_diff = 0.1,
                            dil_scale = 2,
                            response_variable  = "mfi",
                            independent_variable = "concentration",
                            verbose = FALSE) {
  response_variable <- unique(response_variable)
  stdframe <- stdframe[!is.na(stdframe[[response_variable]]) &
                         !is.na(stdframe[[independent_variable]]), ]

  max_response <- max(stdframe[[response_variable]], na.rm = TRUE)
  logc_at_max  <- max(
    stdframe[stdframe[[response_variable]] == max_response, ][[independent_variable]]
  )

  if (verbose) message("Peak = ", max_response, " at x = ", logc_at_max)

  post_peak <- stdframe[[independent_variable]] > logc_at_max
  if (sum(post_peak) == 0) return(stdframe)

  pp_resp <- stdframe[post_peak, response_variable]
  pp_conc <- stdframe[post_peak, independent_variable]

  stdframe[post_peak, response_variable] <-
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
#' @param verbose Logical.
#'
#' @return `data` with the blank operation applied.
#' @export
perform_blank_operation <- function(blank_data, data, response_variable,
                                    independent_variable, is_log_response,
                                    blank_option = "ignored",
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
                                independent_variable)
    if (verbose) message("Blank geometric mean included as extra point.")
    return(data)
  }

  # Subtraction variants
  factor <- switch(blank_option,
                   "subtracted"     = 1,
                   "subtracted_3x"  = 3,
                   "subtracted_10x" = 10)

  blank_mean <- geom_mean(blank_data[[response_variable]])
  data[[response_variable]] <- data[[response_variable]] - factor * blank_mean

  floor_val <- if (is_log_response) 1 else 0
  data[[response_variable]] <- pmax(data[[response_variable]], floor_val)

  if (verbose) message("Blank subtraction (x", factor, ") applied.")
  data
}


#' Include Blank Controls as an Extra Standard Curve Point
#'
#' Appends a synthetic row whose response is the geometric mean of the
#' blanks and whose concentration is half the minimum standard.
#'
#' @param blank_data Data frame of blank wells.
#' @param data Data frame of standards.
#' @param response_variable Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#'
#' @return `data` with one additional row.
#' @keywords internal
include_blanks_conc <- function(blank_data, data, response_variable,
                                independent_variable = "concentration") {
  response_blank     <- geom_mean(blank_data[[response_variable]])
  min_concentration  <- min(data[[independent_variable]], na.rm = TRUE)
  conc_blank         <- min_concentration - log10(2)

  new_row <- data[1, , drop = FALSE]
  new_row[1, ] <- NA
  new_row[[response_variable]]    <- response_blank
  new_row[[independent_variable]] <- conc_blank
  if ("dilution" %in% names(new_row))  new_row$dilution  <- NA_real_
  if ("stype"    %in% names(new_row))  new_row$stype     <- "B"
  if ("well"     %in% names(new_row))  new_row$well      <- "blank_mean"
  if ("sampleid" %in% names(new_row))  new_row$sampleid  <- "blank_mean"

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


#' Full Preprocessing Pipeline for Standard Curve Data
#'
#' Applies concentration computation, prozone correction, blank handling,
#' and optional log10 response transform in the canonical order.
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
#' @param verbose Logical.
#'
#' @return A list with `data` (preprocessed) and `antigen_fit_options` (a
#'   record of the options used).
#' @export
preprocess_standards <- function(data, antigen_settings, response_variable,
                                 independent_variable,
                                 is_log_response,
                                 blank_data         = NULL,
                                 blank_option       = "ignored",
                                 is_log_independent = TRUE,
                                 apply_prozone      = TRUE,
                                 verbose            = FALSE) {

  undiluted_sc <- antigen_settings$standard_curve_concentration

  data <- compute_concentration(data, undiluted_sc, independent_variable,
                                is_log_concentration = is_log_independent)

  if (apply_prozone) {
    data <- correct_prozone(stdframe = data, prop_diff = 0.1, dil_scale = 2,
                            response_variable    = response_variable,
                            independent_variable = independent_variable,
                            verbose              = verbose)
  }

  data <- perform_blank_operation(blank_data = blank_data, data = data,
                                  response_variable    = response_variable,
                                  independent_variable = independent_variable,
                                  is_log_response      = is_log_response,
                                  blank_option         = blank_option,
                                  verbose              = verbose)

  data <- compute_log_response(data, response_variable,
                               is_log_response = is_log_response,
                               verbose = verbose)

  antigen_fit_options <- list(
    is_log_response      = is_log_response,
    blank_option         = blank_option,
    is_log_concentration = is_log_independent,
    apply_prozone        = apply_prozone
  )

  list(data = data, antigen_fit_options = antigen_fit_options)
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
