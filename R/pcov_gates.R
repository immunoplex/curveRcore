# =============================================================================
# pcov_gates.R -- Unified LLOQ/ULOQ-aware pcov_pass + pcov_gate_class
#
# Single source of truth for BOTH curveRfreq and curveRbayes. Replaces the
# old per-point/per-sample "pcov < pcov_threshold" test (which is exactly the
# noisy, jagged-with-few-loops thing that caused samples inside the true
# dynamic range to get marked pcov_pass = FALSE) with an approach that:
#
#   - PREFERS LLOQ/ULOQ (the curve's own dynamic-range boundaries, already
#     derived from the pcov profile via assess_model_eligibility() /
#     .find_lloq_uloq()) when BOTH are defined. Classification is then purely
#     positional: is this point's concentration inside [LLOQ, ULOQ]? pcov
#     itself plays NO role in this branch.
#   - FALLS BACK to the original pcov-vs-threshold test only when either LOQ
#     is undefined, and in that fallback additionally classifies FAILING
#     points as belowLLOQ/aboveULOQ by which side of the inflection point
#     they fall on (since we have no LOQ boundary to compare against).
#
# All comparisons happen in log10, curve-native space: predicted_log10_
# concentration / log10_concentration vs lloq_log10 / uloq_log10 / inflect_x.
# This is deliberate, not incidental -- calib_samples.predicted_concentration
# is NOT reliably log10 (its scale depends on the study's is_log_independent
# setting; see predict_samples.R / predict_bayes.R), and final_concentration
# is dilution-corrected while lloq_conc/uloq_conc are not. The *_log10_*
# quantities are the only ones guaranteed to be on the same, unambiguous
# footing regardless of study settings or a sample's own dilution factor.
#
# IMPORTANT -- eligibility lookup: this file reads eligibility ONLY from
#   cr$ensemble[[best]]$eligibility
# never via cr$selection$assessments. Both curveRfreq and curveRbayes attach
# a per-curve, per-model assessment at cr$ensemble[[model_name]]$eligibility
# identically. cr$selection$assessments means something different in each
# engine: for curveRfreq (single curve per fit) it's trivially per-curve, but
# for curveRbayes (pooled/hierarchical multiplate fits) it is the POOLED,
# model-level assessment shared across every curve in the batch -- using it
# here would silently apply one curve's dynamic range to every plate in a
# Bayesian job. cr$ensemble[[best]]$eligibility is the only location
# confirmed correct and identically-shaped in both engines.
# =============================================================================


# =============================================================================
#' Classify pcov gate membership for a set of points on one curve
#'
#' @param x_log10 Numeric vector. log10-scale position: `predicted_log10_
#'   concentration` for test samples, `log10_concentration` for grid points.
#' @param pcov Numeric vector. Percent CV at each `x_log10`. Only used in the
#'   fallback branch (when LLOQ and/or ULOQ are undefined).
#' @param pcov_threshold Numeric scalar. Only used in the fallback branch.
#' @param lloq_log10,uloq_log10 Numeric scalars (log10 scale), or `NA` if
#'   undefined for this curve/model.
#' @param inflect_x Numeric scalar (log10 scale). Only used in the fallback
#'   branch, to decide which side of the curve a failing point falls on.
#'   If `NA` and a point fails, that point's `pcov_gate_class` stays `NA`
#'   (we know it fails, we just can't say which side) while `pcov_pass`
#'   is still correctly `FALSE`.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{pcov_gate_class}{Factor: `belowLLOQ` / `inDynamicRange` /
#'       `aboveULOQ` / `NA`. In the fallback branch, a failing point with no
#'       `inflect_x` available is `NA` here (we know it fails, not which
#'       side).}
#'     \item{pcov_pass}{Logical, computed independently of
#'       `pcov_gate_class` (NOT derived from it): `TRUE`/`FALSE` whenever the
#'       point could be evaluated at all, `NA` only when `x_log10` (or, in
#'       the fallback branch, `pcov`) itself is unusable. A failing point
#'       with unknown side is `pcov_gate_class = NA` but `pcov_pass =
#'       FALSE`, never `NA`.}
#'   }
#'   Row count and order match the input vectors exactly.
#'
#' @export
classify_pcov_gate <- function(x_log10, pcov, pcov_threshold,
                               lloq_log10 = NA_real_, uloq_log10 = NA_real_,
                               inflect_x  = NA_real_) {
  n <- length(x_log10)
  if (length(pcov) != n)
    stop("classify_pcov_gate(): x_log10 and pcov must be the same length.")

  gate_class <- rep(NA_character_, n)
  pcov_pass  <- rep(NA, n)
  valid_x    <- is.finite(x_log10)

  have_lloq <- length(lloq_log10) == 1L && is.finite(lloq_log10)
  have_uloq <- length(uloq_log10) == 1L && is.finite(uloq_log10)

  if (have_lloq && have_uloq) {
    # ---- Primary path: position relative to LLOQ/ULOQ only. pcov unused. --
    below    <- valid_x & x_log10 <  lloq_log10
    in_range <- valid_x & x_log10 >= lloq_log10 & x_log10 <= uloq_log10
    above    <- valid_x & x_log10 >  uloq_log10

    gate_class[below]    <- "belowLLOQ"
    gate_class[in_range] <- "inDynamicRange"
    gate_class[above]    <- "aboveULOQ"
    pcov_pass[valid_x]   <- in_range[valid_x]

  } else {
    # ---- Fallback: original pcov-vs-threshold test -------------------------
    # pcov_pass is determined the moment we can evaluate pcov < threshold --
    # it does NOT depend on whether we can also assign a side (gate_class).
    # A failing point with no inflect_x still has pcov_pass = FALSE; only its
    # gate_class is left NA (we know it fails, just not which side).
    valid_pcov <- is.finite(pcov)
    evaluable  <- valid_x & valid_pcov
    pass       <- evaluable & (pcov < pcov_threshold)
    pcov_pass[evaluable] <- pass[evaluable]

    gate_class[pass] <- "inDynamicRange"

    fail <- evaluable & !pass
    have_inflect <- length(inflect_x) == 1L && is.finite(inflect_x)
    if (have_inflect) {
      gate_class[fail & x_log10 <  inflect_x] <- "belowLLOQ"
      gate_class[fail & x_log10 >= inflect_x] <- "aboveULOQ"
    }
    # else: fail points keep gate_class = NA; pcov_pass[fail] is already FALSE.
  }

  data.frame(
    pcov_gate_class = factor(gate_class,
                             levels = c("belowLLOQ", "inDynamicRange", "aboveULOQ")),
    pcov_pass = pcov_pass,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
#' Reclassify pcov_pass (grid + samples) and add pcov_gate_class (samples),
#' for every curve in a fitted result, using each curve's own LLOQ/ULOQ
#'
#' Call this AFTER fit_calibration_freq()/fit_calibration_bayes() and AFTER
#' compute_detection_limits_multiplate() -- typically the very next line in
#' the worker pipeline. Overwrites `cr$grid$pcov_pass` and
#' `cr$samples$pcov_pass` for every curve, and adds `cr$samples$
#' pcov_gate_class`. `calib_grid` intentionally gets NO gate_class column
#' (grid pcov_pass is recomputed for display symmetry only; it plays no role
#' in deriving LLOQ/ULOQ, which happens upstream from the raw pcov profile
#' before this function ever runs, so there is no circularity).
#'
#' Also computes each curve's inflection point exactly once and caches it at
#' `cr$detection_limits$inflection` (list(x=, y=)), so `flatten_and_save.R`'s
#' diagnostics-row construction can read the SAME cached value instead of
#' recomputing it independently -- guarantees `calib_diagnostics.inflect_x/y`
#' can never drift from whatever this function used for gating.
#'
#' @param mp A `calibration_result` or `calibration_result_multiplate`.
#' @return The same class of object, with `pcov_pass`/`pcov_gate_class`
#'   (re)computed in place.
#' @export
classify_pcov_gates_multiplate <- function(mp) {
  is_multi <- inherits(mp, "calibration_result_multiplate")
  plates <- if (is_multi) mp$plates
            else stats::setNames(list(mp), as.character(mp$meta$curve_id))

  for (cid in names(plates)) {
    cr <- plates[[cid]]
    if (is.null(cr)) next

    best <- cr$selection$best_model_name %||% NA_character_
    if (is.na(best) || is.null(cr$ensemble[[best]])) {
      plates[[cid]] <- cr
      next
    }

    # Per-curve eligibility -- ALWAYS from $ensemble, never $selection$assessments.
    eb <- tryCatch(cr$ensemble[[best]]$eligibility, error = function(e) NULL) %||% list()
    lloq_log10 <- eb$lloq %||% NA_real_
    uloq_log10 <- eb$uloq %||% NA_real_

    pcov_threshold <- cr$meta$pcov_threshold %||% 20

    infl <- tryCatch(
      compute_inflection(best, cr$ensemble[[best]]$parameters),
      error = function(e) list(x = NA_real_, y = NA_real_)
    )
    if (is.null(cr$detection_limits)) cr$detection_limits <- list()
    cr$detection_limits$inflection <- infl

    # ---- grid: pcov_pass only, no gate_class -------------------------------
    if (!is.null(cr$grid) && is.data.frame(cr$grid) && nrow(cr$grid) > 0) {
      g <- classify_pcov_gate(
        x_log10        = cr$grid$log10_concentration,
        pcov           = cr$grid$pcov,
        pcov_threshold = pcov_threshold,
        lloq_log10     = lloq_log10, uloq_log10 = uloq_log10,
        inflect_x      = infl$x %||% NA_real_
      )
      cr$grid$pcov_pass <- g$pcov_pass
    }

    # ---- samples: pcov_pass + pcov_gate_class ------------------------------
    if (!is.null(cr$samples) && is.data.frame(cr$samples) && nrow(cr$samples) > 0) {
      s <- classify_pcov_gate(
        x_log10        = cr$samples$predicted_log10_concentration,
        pcov           = cr$samples$pcov,
        pcov_threshold = pcov_threshold,
        lloq_log10     = lloq_log10, uloq_log10 = uloq_log10,
        inflect_x      = infl$x %||% NA_real_
      )
      cr$samples$pcov_pass       <- s$pcov_pass
      cr$samples$pcov_gate_class <- as.character(s$pcov_gate_class)
    }

    plates[[cid]] <- cr
  }

  if (is_multi) { mp$plates <- plates; mp } else plates[[1]]
}
