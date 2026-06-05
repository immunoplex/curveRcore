# =============================================================================
# eligibility.R -- Quantification-aware model eligibility gating
#
# Provides assess_model_eligibility() and select_best_eligible(), which are
# called by both curveRfreq and curveRbayes to apply identical gates before
# ranking models by AIC or LOO-CV.
#
# The key insight: AIC and LOO-CV measure forward-fit quality. They are
# agnostic to whether a model can reliably back-calculate concentration.
# These functions add quantification-relevant gates upstream of ranking so
# that an unidentified or boundary-constrained model cannot win selection
# purely on forward-fit criteria.
#
# Gate definitions (all gates must pass for a model to be eligible):
#
#   at_bound        -- No parameter estimate is within bound_tol of its
#                     constraint bound. Boundary estimates invalidate AIC
#                     asymptotics and inflate the covariance matrix.
#                     (Frequentist only; skipped when constraints = NULL.)
#
#   vcov_condition  -- The covariance matrix condition number kappa < 1e8.
#                     Near-singular vcov produces unreliable SE propagation.
#                     (Frequentist only; skipped when vcov_matrix = NULL.)
#
#   rel_se          -- All parameters satisfy SE/|estimate| < max_rel_se.
#                     Parameters with SE larger than the estimate are
#                     unidentified.  Uses std_error (freq) or sd (Bayes).
#
#   dynamic_range   -- The LLOQ-to-ULOQ span at pcov_threshold is at least
#                     min_dynamic_range_log10 log10 units (~3-fold at 0.5).
#                     Zero dynamic range means the model cannot quantify
#                     any sample reliably.
#                     (Requires pcov_profile and grid_x to be supplied.)
# =============================================================================


# -- Internal helper: find LOQ boundaries from a pcov profile -----------------

.find_lloq_uloq <- function(x, pcov, threshold) {
  ok <- is.finite(x) & is.finite(pcov)
  if (sum(ok) < 2L) return(list(lloq = NA_real_, uloq = NA_real_))
  x    <- x[ok];    pcov <- pcov[ok]
  ord  <- order(x); x    <- x[ord]; pcov <- pcov[ord]

  below       <- pcov <= threshold
  if (!any(below)) return(list(lloq = NA_real_, uloq = NA_real_))

  transitions <- diff(below)
  enter       <- which(transitions ==  1L)   # FALSE -> TRUE crossings
  leave       <- which(transitions == -1L)   # TRUE  -> FALSE crossings

  lloq <- if (below[1L]) {
    x[1L]
  } else if (length(enter) > 0L) {
    i <- enter[1L]
    x[i] + (threshold - pcov[i]) / (pcov[i + 1L] - pcov[i]) *
      (x[i + 1L] - x[i])
  } else {
    x[1L]
  }

  uloq <- if (below[length(below)]) {
    x[length(x)]
  } else if (length(leave) > 0L) {
    i <- leave[length(leave)]
    x[i] + (threshold - pcov[i]) / (pcov[i + 1L] - pcov[i]) *
      (x[i + 1L] - x[i])
  } else {
    x[length(x)]
  }

  list(lloq = lloq, uloq = uloq)
}


# =============================================================================
#' Assess Model Eligibility for Quantification
#'
#' Applies a set of identifiability and precision gates to determine whether
#' a fitted model is reliable enough to be considered for selection as the
#' calibration model used for sample quantification.  The same gates are
#' applied in both the frequentist and Bayesian frameworks, with
#' framework-specific gates automatically skipped when the required inputs
#' are not available.
#'
#' @param model_name Character.  Name of the model being assessed.
#' @param parameters Data frame.  Must contain `term` and either `estimate`
#'   + `std_error` (frequentist) or `mean` + `sd` (Bayesian).
#' @param constraints Named list with `$lower` and `$upper` named numeric
#'   vectors for the free parameters, or NULL (Bayesian / fixed-a).
#' @param pcov_profile Numeric vector of pcov values (%) on the prediction
#'   grid, or NULL if the grid has not been computed for this model.
#' @param grid_x Numeric vector of log10_concentration values matching
#'   `pcov_profile`, or NULL.
#' @param pcov_threshold Numeric.  The pcov (%) threshold used to define the
#'   LLOQ and ULOQ.  Default 20.
#' @param bound_tol Numeric.  A parameter estimate within this absolute
#'   distance of a constraint bound is treated as "at bound".  Default 1e-4.
#' @param max_rel_se Numeric.  Maximum permitted relative SE
#'   (`SE / |estimate|`) for any parameter.  Default 5.0.
#' @param min_dynamic_range_log10 Numeric.  Minimum required LLOQ-to-ULOQ
#'   span in log10 concentration units.  Default 0.5 (~3-fold).
#' @param vcov_matrix Numeric matrix.  The parameter covariance matrix from
#'   `vcov(fit)`, or NULL.  Used for the condition-number gate.
#'
#' @return A named list:
#'   \describe{
#'     \item{model_name}{Character.}
#'     \item{eligible}{Logical.  TRUE if all applicable gates pass.}
#'     \item{gates}{Data frame with columns `gate`, `passed`, `detail`.}
#'     \item{dynamic_range_log10}{Numeric or NA.}
#'     \item{lloq}{Numeric or NA (log10 scale).}
#'     \item{uloq}{Numeric or NA (log10 scale).}
#'   }
#'
#' @export
assess_model_eligibility <- function(model_name,
                                     parameters,
                                     constraints        = NULL,
                                     pcov_profile       = NULL,
                                     grid_x             = NULL,
                                     pcov_threshold     = 20,
                                     bound_tol          = 1e-4,
                                     max_rel_se         = 5.0,
                                     min_dynamic_range_log10 = 0.5,
                                     vcov_matrix        = NULL) {

  gates <- list()

  # -- Gate 1: at_bound -----------------------------------------------------
  # Only applies when hard constraint bounds are known (frequentist).
  if (!is.null(constraints) &&
      !is.null(constraints$lower) && !is.null(constraints$upper)) {

    lo  <- constraints$lower
    hi  <- constraints$upper

    # Use the estimate column; tolerate either 'estimate' (freq) or 'mean' (Bayes)
    est_col <- if ("estimate" %in% names(parameters)) "estimate" else "mean"
    ests    <- stats::setNames(parameters[[est_col]], parameters$term)

    at_lo <- at_hi <- character(0)
    for (nm in names(lo)) {
      if (!(nm %in% names(ests))) next
      if (is.finite(lo[nm]) && abs(ests[nm] - lo[nm]) <= bound_tol)
        at_lo <- c(at_lo, nm)
      if (is.finite(hi[nm]) && abs(ests[nm] - hi[nm]) <= bound_tol)
        at_hi <- c(at_hi, nm)
    }

    at_bound_params <- union(at_lo, at_hi)
    passed_bound    <- length(at_bound_params) == 0L

    detail_bound <- if (!passed_bound) {
      parts <- character(0)
      if (length(at_lo) > 0)
        parts <- c(parts, paste0(at_lo, " at lower bound"))
      if (length(at_hi) > 0)
        parts <- c(parts, paste0(at_hi, " at upper bound"))
      paste(parts, collapse = "; ")
    } else ""

    gates[["at_bound"]] <- list(passed = passed_bound, detail = detail_bound)
  }

  # -- Gate 2: vcov_condition ------------------------------------------------
  # Only applies when the vcov matrix is supplied (frequentist).
  if (!is.null(vcov_matrix)) {
    kappa_val <- tryCatch({
      # Use reciprocal condition number from base R: kappa() returns large
      # numbers for ill-conditioned matrices.
      k <- base::kappa(vcov_matrix, exact = FALSE)
      if (!is.finite(k)) Inf else k
    }, error = function(e) Inf)

    passed_kappa <- is.finite(kappa_val) && kappa_val < 1e8
    detail_kappa <- if (!passed_kappa)
      sprintf("condition number = %.3e", kappa_val) else ""

    gates[["vcov_condition"]] <- list(passed = passed_kappa,
                                      detail = detail_kappa)
  }

  # -- Gate 3: rel_se --------------------------------------------------------
  # Applies to both frameworks.  Uses std_error (freq) or sd (Bayes).
  se_col  <- if ("std_error" %in% names(parameters)) "std_error" else "sd"
  est_col <- if ("estimate"  %in% names(parameters)) "estimate"  else "mean"

  ests <- parameters[[est_col]]
  ses  <- parameters[[se_col]]
  terms <- parameters$term

  # Compute relative SE; guard against zero estimates
  rel_se <- ifelse(abs(ests) > 1e-12, ses / abs(ests), Inf)

  bad_rel_se <- terms[is.finite(rel_se) & rel_se > max_rel_se]
  passed_rel_se <- length(bad_rel_se) == 0L

  detail_rel_se <- if (!passed_rel_se) {
    idx <- match(bad_rel_se, terms)
    paste(
      sprintf("%s: rel_se=%.2f", bad_rel_se,
              rel_se[idx]),
      collapse = "; "
    )
  } else ""

  gates[["rel_se"]] <- list(passed = passed_rel_se, detail = detail_rel_se)

  # -- Gate 4: dynamic_range -------------------------------------------------
  # Applies when pcov_profile and grid_x are supplied.
  lloq_val <- uloq_val <- dr_val <- NA_real_

  if (!is.null(pcov_profile) && !is.null(grid_x) &&
      length(pcov_profile) > 0 && length(grid_x) == length(pcov_profile)) {

    loq   <- .find_lloq_uloq(grid_x, pcov_profile, pcov_threshold)
    lloq_val <- loq$lloq
    uloq_val <- loq$uloq

    dr_val <- if (is.finite(lloq_val) && is.finite(uloq_val))
      uloq_val - lloq_val else 0

    passed_dr <- is.finite(dr_val) && dr_val >= min_dynamic_range_log10
    detail_dr <- if (!passed_dr)
      sprintf("dynamic range = %.3g log10 (need >= %.3g)",
              dr_val, min_dynamic_range_log10) else
                sprintf("dynamic range = %.3g log10", dr_val)

    gates[["dynamic_range"]] <- list(passed = passed_dr, detail = detail_dr)
  }

  # -- Assemble gate data frame ----------------------------------------------
  gates_df <- do.call(rbind, lapply(names(gates), function(g) {
    data.frame(gate   = g,
               passed = gates[[g]]$passed,
               detail = gates[[g]]$detail,
               stringsAsFactors = FALSE)
  }))

  # A model with no applicable gates is eligible by default (single-model fit)
  eligible <- if (nrow(gates_df) == 0L) TRUE else all(gates_df$passed)

  list(
    model_name           = model_name,
    eligible             = eligible,
    gates                = gates_df,
    dynamic_range_log10  = dr_val,
    lloq                 = lloq_val,
    uloq                 = uloq_val
  )
}


# =============================================================================
#' Select the Best Eligible Model
#'
#' Given per-model eligibility assessments and a named vector of ranking
#' scores, identifies the best eligible model.  If no model passes all gates,
#' falls back to the model with the widest dynamic range among converged models
#' and records `fallback = TRUE` in the selection.
#'
#' This function is the shared selection layer used by both curveRfreq
#' (ranking by AIC) and curveRbayes (ranking by LOO elpd).  The `aic_selection`
#' or `loo_selection` objects computed by the framework-specific functions are
#' passed through unchanged and augmented with eligibility information.
#'
#' @param assessments Named list of [assess_model_eligibility()] outputs,
#'   one entry per model name.
#' @param ranking_scores Named numeric vector.  For AIC: the AIC value (lower
#'   is better, so provide as-is -- this function finds the minimum).  For LOO:
#'   the elpd value (higher is better -- provide with `higher_is_better = TRUE`).
#' @param criterion Character.  Label stored in the returned selection object,
#'   e.g. `"AIC+eligibility"` or `"LOO+eligibility"`.
#' @param higher_is_better Logical.  If TRUE, higher `ranking_scores` are
#'   preferred (use for LOO elpd).  Default FALSE (use for AIC).
#'
#' @return A named list with:
#'   \describe{
#'     \item{best_model_name}{Character, or NA if no model converged.}
#'     \item{criterion}{Character label.}
#'     \item{fallback}{Logical. TRUE if no eligible model existed.}
#'     \item{fallback_reason}{Character. Describes why fallback occurred.}
#'     \item{assessments}{The full list of [assess_model_eligibility()] outputs.}
#'     \item{eligible_models}{Character vector of eligible model names.}
#'   }
#'
#' @export
select_best_eligible <- function(assessments,
                                 ranking_scores,
                                 criterion        = "AIC+eligibility",
                                 higher_is_better = FALSE) {

  if (length(assessments) == 0L) {
    return(list(
      best_model_name = NA_character_,
      criterion       = criterion,
      fallback        = FALSE,
      fallback_reason = "No models assessed",
      assessments     = assessments,
      eligible_models = character(0)
    ))
  }

  # Separate eligible and ineligible
  eligible_names <- Filter(function(nm) {
    isTRUE(assessments[[nm]]$eligible)
  }, names(assessments))

  # -- Primary path: rank eligible models -----------------------------------
  if (length(eligible_names) > 0L) {
    eligible_scores <- ranking_scores[eligible_names]
    eligible_scores <- eligible_scores[!is.na(eligible_scores)]

    if (length(eligible_scores) == 0L) {
      # Eligible models have no score (shouldn't happen -- treat as fallback)
      best_name    <- eligible_names[1L]
      is_fallback  <- FALSE
      fallback_msg <- ""
    } else {
      best_name <- if (higher_is_better)
        names(which.max(eligible_scores))
      else
        names(which.min(eligible_scores))
      is_fallback  <- FALSE
      fallback_msg <- ""
    }

    # -- Fallback path: no eligible model -------------------------------------
  } else {
    is_fallback <- TRUE

    # Build a summary of why no model passed
    gate_issues <- vapply(names(assessments), function(nm) {
      a    <- assessments[[nm]]
      gdf  <- a$gates
      if (is.null(gdf) || nrow(gdf) == 0L) return("")
      failed <- gdf$gate[!gdf$passed]
      if (length(failed) == 0L) return("")
      paste(failed, collapse = ", ")
    }, character(1))

    gate_summary <- paste(
      sprintf("%s [%s]", names(gate_issues), gate_issues),
      collapse = "; "
    )
    fallback_msg <- paste0("No eligible models. Gate failures -- ", gate_summary)

    if (verbose_fallback <- FALSE) message("[select_best_eligible] ", fallback_msg)

    # Select model with the widest dynamic range, then break ties by ranking score
    dr_vals <- vapply(names(assessments), function(nm) {
      v <- assessments[[nm]]$dynamic_range_log10
      if (is.null(v) || !is.finite(v)) 0 else v
    }, numeric(1))

    max_dr <- max(dr_vals)

    if (max_dr > 0) {
      # Among models with (near-)maximum dynamic range, pick by ranking score
      top_models    <- names(dr_vals)[dr_vals >= max_dr - 0.05]
      top_scores    <- ranking_scores[top_models]
      top_scores    <- top_scores[!is.na(top_scores)]
      if (length(top_scores) > 0L) {
        best_name <- if (higher_is_better)
          names(which.max(top_scores))
        else
          names(which.min(top_scores))
      } else {
        best_name <- names(which.max(dr_vals))
      }
    } else {
      # Truly no usable model: pick by ranking score alone
      scores <- ranking_scores[names(assessments)]
      scores <- scores[!is.na(scores)]
      if (length(scores) > 0L) {
        best_name <- if (higher_is_better)
          names(which.max(scores))
        else
          names(which.min(scores))
      } else {
        best_name <- names(assessments)[1L]
      }
    }
  }

  list(
    best_model_name = best_name,
    criterion       = criterion,
    fallback        = is_fallback,
    fallback_reason = fallback_msg,
    assessments     = assessments,
    eligible_models = eligible_names
  )
}
