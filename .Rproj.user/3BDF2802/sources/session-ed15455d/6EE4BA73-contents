# =============================================================================
# detection_limits.R
#
# Assay detection and quantification metrics for the calibration_result class.
#
# Three tiers of functions:
#
#   Tier 1 — Grid enrichment (called inside predict_grid_* at Step 2):
#     enrich_grid_with_d2y()
#
#   Tier 2 — Shape-LOQ extraction (called at Step 3 / eligibility):
#     compute_shape_loq_from_grid()
#
#   Tier 3 — Post-hoc detection limits (called after Step 4):
#     compute_detection_limits()      — populates $detection_limits
#     compute_detection_limits_multiplate() — iterates over plates
#
# All tier-3 functions can also be called by the user post-hoc on an
# existing calibration_result, making them safe to add without modifying
# curveRfreq or curveRbayes at all during initial integration.
# =============================================================================


# ─── Tier 1: Grid enrichment ────────────────────────────────────────────────

#' Add a d2y_dx2 column to an existing prediction grid
#'
#' Computes \eqn{d^2(\log_{10} y) / d(\log_{10} x)^2} from the
#' \code{log10_concentration} and \code{predicted_response} columns already
#' present on the grid, using non-uniform central differences.
#'
#' This function is designed to be called inside \code{predict_grid_freq()}
#' or \code{predict_grid_bayes()} immediately after the grid data frame is
#' constructed, adding the column in-place with zero additional model
#' evaluations.
#'
#' When \code{is_log_response = TRUE}, the \code{predicted_response} column
#' is already on the \eqn{\log_{10}} scale, so \code{d2y_dx2} is computed
#' directly.  When \code{is_log_response = FALSE}, the column is first
#' \eqn{\log_{10}}-transformed before differencing.
#'
#' @param grid Data frame with columns \code{log10_concentration} and
#'   \code{predicted_response}.
#' @param is_log_response Logical; whether \code{predicted_response} is
#'   already on the \eqn{\log_{10}} scale (default \code{TRUE}).
#'
#' @return The input \code{grid} with an additional column \code{d2y_dx2}.
#'   Boundary points (first and last rows) receive \code{NA}.
#'
#' @export
enrich_grid_with_d2y <- function(grid, is_log_response = TRUE) {

  u <- grid$log10_concentration
  y <- if (is_log_response) {
    grid$predicted_response
  } else {
    log10(pmax(grid$predicted_response, 1e-300))
  }

  n  <- length(u)
  d2 <- rep(NA_real_, n)

  if (n < 3L) {
    grid$d2y_dx2 <- d2
    return(grid)
  }

  for (i in seq(2L, n - 1L)) {
    if (anyNA(c(y[i - 1L], y[i], y[i + 1L]))) next
    h_l <- u[i] - u[i - 1L]
    h_r <- u[i + 1L] - u[i]
    if (h_l < .Machine$double.eps || h_r < .Machine$double.eps) next
    d2[i] <- 2 * ((y[i + 1L] - y[i]) / h_r -
                    (y[i] - y[i - 1L]) / h_l) / (h_l + h_r)
  }

  grid$d2y_dx2 <- d2
  grid
}


# ─── Tier 2: Shape-LOQ from d2y_dx2 on the grid ─────────────────────────────

#' Compute curvature-based (shape) LOQs from an enriched grid
#'
#' Identifies the lower shape-LOQ as the concentration at the global
#' maximum of \code{d2y_dx2} and the upper shape-LOQ as the
#' concentration at the global minimum.  Vertex positions are refined
#' by 3-point parabolic interpolation.
#'
#' Designed to be called inside \code{assess_model_eligibility()} at
#' Step 3, returning a named list that is merged into the
#' \code{$eligibility} slot.
#'
#' @param grid Data frame with columns \code{log10_concentration},
#'   \code{predicted_response}, and \code{d2y_dx2}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{shape_lloq_log10}}{log10 concentration at the
#'       global maximum of d2y_dx2 (lower shape-LOQ).}
#'     \item{\code{shape_uloq_log10}}{log10 concentration at the
#'       global minimum of d2y_dx2 (upper shape-LOQ).}
#'     \item{\code{shape_lloq_conc}}{Natural-scale lower shape-LOQ.}
#'     \item{\code{shape_uloq_conc}}{Natural-scale upper shape-LOQ.}
#'     \item{\code{shape_lloq_response}}{Predicted response at shape-LLOQ.}
#'     \item{\code{shape_uloq_response}}{Predicted response at shape-ULOQ.}
#'   }
#'
#' @export
compute_shape_loq_from_grid <- function(grid) {

  empty <- list(
    shape_lloq_log10    = NA_real_, shape_uloq_log10    = NA_real_,
    shape_lloq_conc     = NA_real_, shape_uloq_conc     = NA_real_,
    shape_lloq_response = NA_real_, shape_uloq_response = NA_real_
  )

  if (is.null(grid$d2y_dx2) || sum(!is.na(grid$d2y_dx2)) < 3L) return(empty)

  ok <- is.finite(grid$log10_concentration) & is.finite(grid$d2y_dx2)
  x  <- grid$log10_concentration[ok]
  y  <- grid$d2y_dx2[ok]
  yr <- grid$predicted_response[ok]

  if (length(y) < 3L) return(empty)

  # Local extrema via sign changes in first differences
  dy      <- diff(y)
  idx_max <- which(dy[-1L] < 0 & dy[-length(dy)] > 0) + 1L
  idx_min <- which(dy[-1L] > 0 & dy[-length(dy)] < 0) + 1L

  # Refine by 3-point parabolic interpolation
  max_verts <- .parabolic_refine(x, y, idx_max)
  min_verts <- .parabolic_refine(x, y, idx_min)

  # Global max → shape-LLOQ; global min → shape-ULOQ
  lloq_log10 <- if (nrow(max_verts) > 0) {
    max_verts$x[which.max(max_verts$y)]
  } else NA_real_

  uloq_log10 <- if (nrow(min_verts) > 0) {
    min_verts$x[which.min(min_verts$y)]
  } else NA_real_

  # Interpolate predicted_response at shape-LOQ positions
  lloq_response <- .interp_response(x, yr, lloq_log10)
  uloq_response <- .interp_response(x, yr, uloq_log10)

  list(
    shape_lloq_log10    = lloq_log10,
    shape_uloq_log10    = uloq_log10,
    shape_lloq_conc     = if (!is.na(lloq_log10)) 10^lloq_log10 else NA_real_,
    shape_uloq_conc     = if (!is.na(uloq_log10)) 10^uloq_log10 else NA_real_,
    shape_lloq_response = lloq_response,
    shape_uloq_response = uloq_response
  )
}


# ─── Tier 3: Post-hoc detection limits ──────────────────────────────────────

#' Compute and attach detection limits to a calibration_result
#'
#' Computes LODs, MDC, and RDL for the best eligible model (or a specified
#' model) and returns the calibration_result with a new
#' \code{$detection_limits} element.
#'
#' This can be called:
#' \itemize{
#'   \item Automatically inside the fitting pipeline (after Step 4), or
#'   \item By the user post-hoc on an existing \code{calibration_result}.
#' }
#'
#' @param cr A \code{calibration_result} object.
#' @param model_name Character; model to use.  If \code{NULL} (default),
#'   uses \code{cr$selection$best_model_name}.
#' @param alpha Numeric; significance level for CI-based LODs (default 0.05).
#' @param verbose Logical; emit diagnostic messages (default \code{FALSE}).
#'
#' @return The input \code{cr} with \code{$detection_limits} populated.
#'
#' @export
compute_detection_limits <- function(cr, model_name = NULL,
                                     alpha = 0.05, verbose = FALSE) {

  stopifnot(inherits(cr, "calibration_result"))

  model_name <- model_name %||% cr$selection$best_model_name
  method     <- cr$meta$method
  ens        <- cr$ensemble[[model_name]]

  if (is.null(ens) || !isTRUE(ens$converged)) {
    cr$detection_limits <- .empty_detection_limits(
      cr$meta$curve_id, model_name, method, alpha)
    return(cr)
  }

  # ── Extract CIs ──
  pci <- .extract_param_ci(ens, method, alpha)

  # ── LODs on the response (fitting) scale ──
  lower_lod_response <- as.numeric(pci$ci_upper["a"])
  upper_lod_response <- as.numeric(pci$ci_lower["d"])

  if (is.na(upper_lod_response) || is.na(lower_lod_response) ||
      upper_lod_response <= lower_lod_response) {
    if (verbose) message(sprintf(
      "[compute_detection_limits] curve_id=%s ULOD <= LLOD; ULOD set to NA",
      cr$meta$curve_id))
    upper_lod_response <- NA_real_
  }

  est <- pci$estimate

  # ── LODs mapped to concentration ──
  lower_lod_log10 <- .safe_invert_model(model_name, lower_lod_response, est)
  upper_lod_log10 <- .safe_invert_model(model_name, upper_lod_response, est)

  # ── MDC: invert point-estimate curve at LOD responses ──
  mdc_lower <- lower_lod_log10   # same as LOD concentration mapping

  mdc_upper <- upper_lod_log10

  # ── RDL: invert CI-modified curves ──
  est_compressed      <- est
  est_compressed["d"] <- as.numeric(pci$ci_lower["d"])
  rdl_lower <- .safe_invert_model(model_name, lower_lod_response, est_compressed)

  est_expanded      <- est
  est_expanded["d"] <- as.numeric(pci$ci_upper["d"])
  rdl_upper <- .safe_invert_model(model_name, upper_lod_response, est_expanded)

  if (verbose) message(sprintf(
    "[detection_limits] curve_id=%s  LLOD_resp=%.3f  ULOD_resp=%s  MDC=[%s, %s]  RDL=[%s, %s]",
    cr$meta$curve_id,
    lower_lod_response,
    .fmt_na(upper_lod_response),
    .fmt_na(mdc_lower), .fmt_na(mdc_upper),
    .fmt_na(rdl_lower), .fmt_na(rdl_upper)))

  cr$detection_limits <- list(
    model_name = model_name,
    method     = method,
    alpha      = alpha,
    lods = list(
      lower_lod_response   = lower_lod_response,
      upper_lod_response   = upper_lod_response,
      lower_lod_log10_conc = lower_lod_log10,
      upper_lod_log10_conc = upper_lod_log10,
      lower_lod_conc       = .safe_pow10(lower_lod_log10),
      upper_lod_conc       = .safe_pow10(upper_lod_log10)
    ),
    mdc_rdl = list(
      mdc_lower_log10 = mdc_lower,
      mdc_upper_log10 = mdc_upper,
      mdc_lower_conc  = .safe_pow10(mdc_lower),
      mdc_upper_conc  = .safe_pow10(mdc_upper),
      rdl_lower_log10 = rdl_lower,
      rdl_upper_log10 = rdl_upper,
      rdl_lower_conc  = .safe_pow10(rdl_lower),
      rdl_upper_conc  = .safe_pow10(rdl_upper)
    )
  )

  cr
}


#' Compute detection limits for all plates in a multiplate result
#'
#' Iterates over \code{$plates} and calls
#' \code{\link{compute_detection_limits}} on each.
#'
#' @param mp A \code{calibration_result_multiplate} object.
#' @param alpha Numeric; significance level (default 0.05).
#' @param verbose Logical; emit diagnostic messages (default \code{FALSE}).
#'
#' @return The input \code{mp} with \code{$detection_limits} populated
#'   on every plate.
#'
#' @export
compute_detection_limits_multiplate <- function(mp, alpha = 0.05,
                                                verbose = FALSE) {
  stopifnot(inherits(mp, "calibration_result_multiplate"))

  for (cid in names(mp$plates)) {
    cr <- mp$plates[[cid]]
    if (!is.null(cr)) {
      mp$plates[[cid]] <- compute_detection_limits(
        cr, alpha = alpha, verbose = verbose)
    }
  }
  mp
}


# ─── Internal helpers ────────────────────────────────────────────────────────

#' Extract parameter CIs from an ensemble entry
#' @keywords internal
.extract_param_ci <- function(ensemble_entry, method, alpha = 0.05) {
  params <- ensemble_entry$parameters
  z      <- stats::qnorm(1 - alpha / 2)

  if (method == "frequentist") {
    est <- stats::setNames(params$estimate, params$term)
    se  <- stats::setNames(params$std_error, params$term)
    list(estimate = est,
         ci_lower = est - z * se,
         ci_upper = est + z * se)
  } else {
    est <- stats::setNames(params$mean, params$term)
    lo_col <- paste0("q", 100 * alpha / 2)
    hi_col <- paste0("q", 100 * (1 - alpha / 2))
    ci_lo <- if (lo_col %in% names(params)) {
      stats::setNames(params[[lo_col]], params$term)
    } else {
      stats::setNames(params$mean - z * params$sd, params$term)
    }
    ci_hi <- if (hi_col %in% names(params)) {
      stats::setNames(params[[hi_col]], params$term)
    } else {
      stats::setNames(params$mean + z * params$sd, params$term)
    }
    list(estimate = est, ci_lower = ci_lo, ci_upper = ci_hi)
  }
}


#' Safely invert a named model at a response value
#' @keywords internal
.safe_invert_model <- function(model_name, y, params) {
  if (is.na(y)) return(NA_real_)
  tryCatch({
    a <- as.numeric(params["a"]); b <- as.numeric(params["b"])
    c_val <- as.numeric(params["c"]); d <- as.numeric(params["d"])
    g <- if ("g" %in% names(params) && !is.na(params["g"])) {
      as.numeric(params["g"])
    } else NULL

    result <- switch(model_name,
                     logistic4    = inv_logistic4(y, a, b, c_val, d),
                     logistic5    = inv_logistic5(y, a, b, c_val, d, g),
                     gompertz4    = inv_gompertz4(y, a, b, c_val, d),
                     loglogistic4 = inv_loglogistic4(y, a, b, c_val, d),
                     loglogistic5 = inv_loglogistic5(y, a, b, c_val, d, g),
                     stop(sprintf("Unknown model: %s", model_name), call. = FALSE)
    )
    as.numeric(result)
  }, error = function(e) NA_real_)
}


#' 3-point parabolic interpolation for a set of extremum indices
#' @keywords internal
.parabolic_refine <- function(x, y, indices) {
  if (length(indices) == 0L)
    return(data.frame(x = numeric(0), y = numeric(0)))

  verts <- lapply(indices, function(idx) {
    xi <- x[(idx - 1L):(idx + 1L)]
    yi <- y[(idx - 1L):(idx + 1L)]

    denom <- (xi[1] - xi[2]) * (xi[1] - xi[3]) * (xi[2] - xi[3])
    if (abs(denom) < .Machine$double.eps * 1e8)
      return(c(x = xi[2], y = yi[2]))

    a_c <- (xi[3]*(yi[2]-yi[1]) + xi[2]*(yi[1]-yi[3]) + xi[1]*(yi[3]-yi[2])) / denom
    b_c <- (xi[3]^2*(yi[1]-yi[2]) + xi[2]^2*(yi[3]-yi[1]) + xi[1]^2*(yi[2]-yi[3])) / denom

    if (abs(a_c) < .Machine$double.eps)
      return(c(x = xi[2], y = yi[2]))

    xv <- -b_c / (2 * a_c)
    xv <- max(min(xv, max(xi)), min(xi))
    yv <- a_c * xv^2 + b_c * xv + (yi[1] - a_c * xi[1]^2 - b_c * xi[1])
    c(x = xv, y = yv)
  })

  data.frame(
    x = vapply(verts, `[[`, numeric(1), "x"),
    y = vapply(verts, `[[`, numeric(1), "y")
  )
}


#' Linear interpolation of predicted_response at a log10 concentration
#' @keywords internal
.interp_response <- function(x_grid, y_grid, x_target) {
  if (is.na(x_target) || length(x_grid) < 2L) return(NA_real_)
  stats::approx(x_grid, y_grid, xout = x_target, rule = 1)$y
}


#' Safe 10^x, NA-propagating
#' @keywords internal
.safe_pow10 <- function(x) {
  if (is.na(x)) NA_real_ else 10^x
}


#' Format a number for messages, NA-safe
#' @keywords internal
.fmt_na <- function(x, digits = 4) {
  if (is.na(x)) "NA" else format(x, digits = digits)
}


#' Empty detection_limits list for non-converged models
#' @keywords internal
.empty_detection_limits <- function(curve_id, model_name, method, alpha) {
  na_lods <- list(
    lower_lod_response = NA_real_, upper_lod_response = NA_real_,
    lower_lod_log10_conc = NA_real_, upper_lod_log10_conc = NA_real_,
    lower_lod_conc = NA_real_, upper_lod_conc = NA_real_)
  na_mdc <- list(
    mdc_lower_log10 = NA_real_, mdc_upper_log10 = NA_real_,
    mdc_lower_conc = NA_real_, mdc_upper_conc = NA_real_,
    rdl_lower_log10 = NA_real_, rdl_upper_log10 = NA_real_,
    rdl_lower_conc = NA_real_, rdl_upper_conc = NA_real_)
  list(model_name = model_name, method = method, alpha = alpha,
       lods = na_lods, mdc_rdl = na_mdc)
}
