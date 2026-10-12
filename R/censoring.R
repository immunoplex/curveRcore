# =============================================================================
# censoring.R -- censoring class and bounds for back-calculated test samples
#
# A sample outside the quantifiable range is not an imprecise point; it is an
# observation known only to lie in an interval. These functions record which
# interval, so downstream analyses can use censored-data methods (Tobit,
# survival-style, Bayesian censored likelihoods) instead of substituting
# LLOQ/2 or down-weighting.
#
# Classes (checked in this order):
#   no_response  response missing / non-finite
#   below_lod    response below the lower LOD response  -> left,     (0, LOD]
#   saturated    response at/above the upper LOD response -> right,  [ULOQ, inf)
#   lod_to_lloq  conc below LLOQ, above LOD             -> interval, [LOD, LLOQ]
#   below_lloq   conc below LLOQ, LOD unavailable       -> left,     (0, LLOQ]
#   above_uloq   conc above ULOQ                        -> right,    [ULOQ, inf)
#   quantified   inside [LLOQ, ULOQ]                    -> none
#   NA           curve has no LLOQ/ULOQ, and no LOD-based class applies
#
# Response-scale LOD checks come first so a sample whose point estimate could
# not be inverted (response beyond an asymptote) is still classified. They are
# applied only to increasing curves (d > a), the orientation
# compute_detection_limits() assumes; for other curves the LOD check uses the
# log10 LOD concentration instead. LLOQ/ULOQ come from
# cr$ensemble[[best]]$eligibility -- the per-curve location, never
# cr$selection$assessments (see pcov_gates.R for why). Bounds are on the
# natural scale and multiplied by the sample's dilution, the same scale as
# final_concentration; NA means unbounded on that side.
# =============================================================================


#' Classify censoring for one curve's test samples
#'
#' @param response Numeric vector: `observed_response_fit` (fitting scale).
#' @param x_log10 Numeric vector: `predicted_log10_concentration` (pre-dilution).
#' @param dilution Numeric vector (or scalar): each sample's dilution factor.
#' @param lloq_log10,uloq_log10 Numeric scalars (log10, pre-dilution) or `NA`.
#' @param lod_response,ulod_response Numeric scalars: `lower_lod_response` and
#'   `upper_lod_response` from [compute_detection_limits()], or `NA`.
#' @param lod_log10,ulod_log10 Numeric scalars: `lower_lod_log10_conc` and
#'   `upper_lod_log10_conc`, or `NA`.
#' @param a,d Numeric scalars: the best model's asymptote estimates (fitting
#'   scale). Used for the curve's orientation and to place a non-invertible
#'   sample on the correct side.
#'
#' @return A data.frame with one row per sample and columns `censor_class`,
#'   `censor_type` (`none`, `left`, `interval`, `right`, or `NA`),
#'   `censor_lower`, `censor_upper`.
#'
#' @seealso [classify_censoring_multiplate()], [summarize_conc_draws()]
#' @export
classify_censoring <- function(response, x_log10, dilution,
                               lloq_log10 = NA_real_, uloq_log10 = NA_real_,
                               lod_response = NA_real_, ulod_response = NA_real_,
                               lod_log10 = NA_real_, ulod_log10 = NA_real_,
                               a = NA_real_, d = NA_real_) {
  n <- length(response)
  if (length(x_log10) != n)
    stop("classify_censoring(): response and x_log10 must be the same length.")
  if (length(dilution) == 1L) dilution <- rep(dilution, n)
  dil <- suppressWarnings(as.numeric(dilution))
  dil[!is.finite(dil) | dil <= 0] <- NA_real_

  fin <- function(v) length(v) == 1L && is.finite(v)
  # A fixed-a fit carries no `a` estimate; the lower LOD response (just above
  # a) is then the best stand-in for the low-concentration asymptote.
  if (!fin(a) && fin(lod_response)) a <- lod_response
  increasing <- fin(a) && fin(d) && d > a
  have_loqs  <- fin(lloq_log10) && fin(uloq_log10)
  lod_conc   <- if (fin(lod_log10))  10^lod_log10  else NA_real_
  ulod_conc  <- if (fin(ulod_log10)) 10^ulod_log10 else NA_real_
  lloq_conc  <- if (fin(lloq_log10)) 10^lloq_log10 else NA_real_
  uloq_conc  <- if (fin(uloq_log10)) 10^uloq_log10 else NA_real_

  cls     <- rep(NA_character_, n)
  has_y   <- is.finite(response)
  has_x   <- is.finite(x_log10)
  open    <- function() is.na(cls) & has_y

  # Which side a sample sits on when its point estimate is not finite.
  low_side <- rep(NA, n)
  if (fin(a) && fin(d)) low_side <- (response - a) / (d - a) < 0.5

  cls[!has_y] <- "no_response"

  # -- below LOD ------------------------------------------------------------
  if (increasing && fin(lod_response)) {
    cls[open() & response < lod_response] <- "below_lod"
  } else if (fin(lod_log10)) {
    cls[open() & has_x & x_log10 < lod_log10] <- "below_lod"
  }

  # -- saturated (at/above the upper LOD) ------------------------------------
  if (increasing && fin(ulod_response)) {
    cls[open() & response >= ulod_response] <- "saturated"
  } else if (fin(ulod_log10)) {
    cls[open() & has_x & x_log10 >= ulod_log10] <- "saturated"
  }

  # -- non-invertible point estimate: side from the response -----------------
  nx <- open() & !has_x & !is.na(low_side)
  cls[nx &  low_side] <- if (is.finite(lod_conc)) "below_lod" else "below_lloq"
  cls[nx & !low_side] <- "saturated"

  # -- position relative to LLOQ / ULOQ --------------------------------------
  if (have_loqs) {
    below <- open() & has_x & x_log10 < lloq_log10
    cls[below] <- if (is.finite(lod_conc) && lod_conc < lloq_conc) "lod_to_lloq" else "below_lloq"
    cls[open() & has_x & x_log10 > uloq_log10] <- "above_uloq"
    cls[open() & has_x] <- "quantified"
  }

  # -- type and natural-scale, dilution-corrected bounds ---------------------
  type <- c(no_response = NA, below_lod = "left", below_lloq = "left",
            lod_to_lloq = "interval", quantified = "none",
            above_uloq = "right", saturated = "right")[cls]
  low_bound  <- if (is.finite(lod_conc)) lod_conc else lloq_conc
  high_bound <- if (is.finite(uloq_conc)) uloq_conc else ulod_conc

  lower <- upper <- rep(NA_real_, n)
  upper[cls %in% "below_lod"]   <- low_bound * dil[cls %in% "below_lod"]
  upper[cls %in% "below_lloq"]  <- lloq_conc * dil[cls %in% "below_lloq"]
  lower[cls %in% "lod_to_lloq"] <- lod_conc  * dil[cls %in% "lod_to_lloq"]
  upper[cls %in% "lod_to_lloq"] <- lloq_conc * dil[cls %in% "lod_to_lloq"]
  lower[cls %in% c("above_uloq", "saturated")] <-
    high_bound * dil[cls %in% c("above_uloq", "saturated")]

  data.frame(censor_class = cls, censor_type = unname(type),
             censor_lower = lower, censor_upper = upper,
             stringsAsFactors = FALSE)
}


#' Add censoring columns to every curve's test samples
#'
#' Call AFTER [compute_detection_limits_multiplate()] and
#' [classify_pcov_gates_multiplate()] -- the first point at which LODs,
#' LLOQ/ULOQ and samples all exist on the result. Adds `censor_class`,
#' `censor_type`, `censor_lower` and `censor_upper` to `cr$samples` for every
#' curve (all `NA` for a curve with no selected model).
#'
#' @param mp A `calibration_result` or `calibration_result_multiplate`.
#' @return The same class of object with the censoring columns added.
#' @seealso [classify_censoring()]
#' @export
classify_censoring_multiplate <- function(mp) {
  is_multi <- inherits(mp, "calibration_result_multiplate")
  plates <- if (is_multi) mp$plates
            else stats::setNames(list(mp), as.character(mp$meta$curve_id))

  cens_names <- c("censor_class", "censor_type", "censor_lower", "censor_upper")
  na_cols <- function(n) data.frame(censor_class = rep(NA_character_, n),
                                    censor_type  = rep(NA_character_, n),
                                    censor_lower = rep(NA_real_, n),
                                    censor_upper = rep(NA_real_, n),
                                    stringsAsFactors = FALSE)

  for (cid in names(plates)) {
    cr <- plates[[cid]]
    if (is.null(cr) || is.null(cr$samples) || !is.data.frame(cr$samples) ||
        !nrow(cr$samples)) next
    s <- cr$samples

    best <- cr$selection$best_model_name %||% NA_character_
    ens  <- if (!is.na(best)) cr$ensemble[[best]] else NULL
    if (is.null(ens)) {
      cr$samples <- cbind(s[setdiff(names(s), cens_names)], na_cols(nrow(s)))
      plates[[cid]] <- cr
      next
    }

    eb   <- tryCatch(ens$eligibility, error = function(e) NULL) %||% list()
    lods <- cr$detection_limits$lods %||% list()
    est  <- tryCatch(.extract_param_ci(ens, cr$meta$method %||% "frequentist")$estimate,
                     error = function(e) NULL)
    par  <- function(nm) if (!is.null(est) && nm %in% names(est)) as.numeric(est[[nm]]) else NA_real_
    a_hat <- if (is.finite(par("a"))) par("a") else (cr$meta$fixed_a %||% NA_real_)

    cc <- classify_censoring(
      response      = s$observed_response_fit %||% rep(NA_real_, nrow(s)),
      x_log10       = s$predicted_log10_concentration %||% rep(NA_real_, nrow(s)),
      dilution      = s$dilution %||% 1,
      lloq_log10    = eb$lloq %||% NA_real_,
      uloq_log10    = eb$uloq %||% NA_real_,
      lod_response  = lods$lower_lod_response %||% NA_real_,
      ulod_response = lods$upper_lod_response %||% NA_real_,
      lod_log10     = lods$lower_lod_log10_conc %||% NA_real_,
      ulod_log10    = lods$upper_lod_log10_conc %||% NA_real_,
      a = a_hat, d = par("d"))

    cr$samples <- cbind(s[setdiff(names(s), names(cc))], cc)
    plates[[cid]] <- cr
  }

  if (is_multi) { mp$plates <- plates; mp } else plates[[1]]
}
