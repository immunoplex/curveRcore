# =============================================================================
# result_class.R -- Shared calibration_result S3 class
#
# The single most important deliverable of curveRcore. Both curveRfreq and
# curveRbayes return this exact structure, making cross-method comparison
# trivial.
# =============================================================================


#' Construct a Calibration Result Object
#'
#' Creates a validated `calibration_result` S3 object. This is the
#' canonical output of both `fit_calibration_freq()` (curveRfreq) and
#' `fit_calibration_bayes()` (curveRbayes).
#'
#' @param meta Named list of metadata. Required fields:
#'   `method` (`"frequentist"` or `"bayesian"`), `package`, `version`,
#'   `feature`, `antigen`, `plate`, `response_var`, `independent_var`,
#'   `is_log_response`, `is_log_independent`.
#' @param ensemble Named list of model fit results, one per model attempted.
#'   Each element should be a list with `model_name`, `converged`,
#'   `parameters`, `fit_stats`, and optionally `raw_fit`.
#' @param selection Named list describing best-model selection:
#'   `best_model_name`, `criterion`, `weights` (data.frame).
#' @param grid Data frame of grid predictions from the best model.
#' @param samples Data frame of test sample predictions from the best model,
#'   or NULL if no samples were provided.
#' @param diagnostics Named list of diagnostic quantities (inflection point,
#'   LODs, LOQs, etc.), or NULL.
#' @param standards Data frame of standard data used to fit the curve,
#'   or NULL if not provided.
#' @param blanks Data frame of blank data used for QA,
#'   or NULL if not provided.
#'
#' @return An object of class `calibration_result`.
#'
#' @export
new_calibration_result <- function(meta,
                                   ensemble   = list(),
                                   selection  = list(),
                                   grid       = data.frame(),
                                   samples    = NULL,
                                   diagnostics = NULL,
                                   standards = NULL,
                                   blanks = NULL) {

  # Validate required meta fields
  required_meta <- c("method", "package", "curve_id",
                     "response_var", "independent_var",
                     "is_log_response", "is_log_independent")
  missing <- setdiff(required_meta, names(meta))
  if (length(missing) > 0) {
    stop("meta is missing required fields: ", paste(missing, collapse = ", "))
  }

  if (!(meta$method %in% c("frequentist", "bayesian"))) {
    stop("meta$method must be 'frequentist' or 'bayesian'")
  }

  # Defaults
  if (is.null(meta$version))   meta$version   <- "0.1.0"
  if (is.null(meta$timestamp)) meta$timestamp <- Sys.time()

  out <- list(
    meta        = meta,
    ensemble    = ensemble,
    selection   = selection,
    grid        = grid,
    samples     = samples,
    diagnostics = diagnostics,
    standards = standards,
    blanks = blanks
  )
  class(out) <- c("calibration_result", "list")
  out
}


#' @export
print.calibration_result <- function(x, ...) {
  cat(sprintf("-- calibration_result (%s) --\n", x$meta$method))
  cat(sprintf("  Curve ID : %s\n", x$meta$curve_id))
  cat(sprintf("  Package  : %s v%s\n", x$meta$package, x$meta$version))

  n_models <- length(x$ensemble)
  n_conv   <- sum(vapply(x$ensemble, function(e) isTRUE(e$converged), logical(1)))
  cat(sprintf("  Models   : %d attempted, %d converged\n", n_models, n_conv))

  best <- x$selection$best_model_name
  if (!is.null(best) && !is.na(best)) {
    cat(sprintf("  Best     : %s (by %s)\n", best, x$selection$criterion %||% "?"))
  } else {
    cat("  Best     : (none converged)\n")
  }

  cat(sprintf("  Grid     : %d points\n", nrow(x$grid)))
  if (!is.null(x$standards)) {
    cat(sprintf("  Standards: %d points\n", nrow(x$standards)))
  }
  if (!is.null(x$blanks)) {
    cat(sprintf("  Blanks: %d points\n", nrow(x$blanks)))
  }
  if (!is.null(x$samples)) {
    cat(sprintf("  Samples  : %d predicted\n", nrow(x$samples)))
  }
  invisible(x)
}


#' @export
summary.calibration_result <- function(object, ...) {
  cat(sprintf("-- Summary: calibration_result (%s) curve_id=%s --\n\n",
              object$meta$method, object$meta$curve_id))

  # Selection table
  if (!is.null(object$selection$weights) && nrow(object$selection$weights) > 0) {
    cat("Model selection:\n")
    print(object$selection$weights, row.names = FALSE)
    cat("\n")
  }

  # Best model parameters
  best <- object$selection$best_model_name
  if (!is.null(best) && !is.na(best) && best %in% names(object$ensemble)) {
    ens <- object$ensemble[[best]]
    if (!is.null(ens$parameters)) {
      cat(sprintf("Best model (%s) parameters:\n", best))
      print(ens$parameters, row.names = FALSE)
      cat("\n")
    }
  }

  # Diagnostics
  dx <- object$diagnostics
  if (!is.null(dx)) {
    cat("Diagnostics:\n")
    if (!is.null(dx$inflection_point)) {
      cat(sprintf("  Inflection: x = %.4f, y = %.4f\n",
                  dx$inflection_point$x, dx$inflection_point$y))
    }
    if (!is.null(dx$lloq) && !is.na(dx$lloq))
      cat(sprintf("  LLOQ: %.4f\n", dx$lloq))
    if (!is.null(dx$uloq) && !is.na(dx$uloq))
      cat(sprintf("  ULOQ: %.4f\n", dx$uloq))
  }

  invisible(object)
}


# ---- Multi-plate wrapper ----

#' Construct a Multi-Plate Calibration Result
#'
#' Wraps multiple single-plate `calibration_result` objects into a
#' multi-plate container.
#'
#' @param meta Named list. Multi-plate metadata (must include `plates`
#'   character vector).
#' @param plates Named list of `calibration_result` objects, one per plate.
#'
#' @return An object of class `calibration_result_multiplate`.
#'
#' @export
new_calibration_result_multiplate <- function(meta, plates) {
  stopifnot(is.list(plates))

  out <- list(meta = meta, plates = plates)
  class(out) <- c("calibration_result_multiplate", "list")
  out
}


#' @export
print.calibration_result_multiplate <- function(x, ...) {
  cat(sprintf("-- calibration_result_multiplate (%s) --\n", x$meta$method))
  cat(sprintf("  Curves : %d (%s)\n", length(x$plates),
              paste(names(x$plates), collapse = ", ")))
  invisible(x)
}


# ---- Comparison utility ----

#' Compare Two Calibration Results (Grid Predictions)
#'
#' Merges the prediction grids from two `calibration_result` objects
#' (typically one frequentist and one Bayesian) into a single data frame
#' for paired analysis. Columns are prefixed with the method name.
#'
#' @param result_a A `calibration_result`.
#' @param result_b A `calibration_result`.
#' @param label_a Character prefix for result_a columns. Default uses
#'   `result_a$meta$method`.
#' @param label_b Character prefix for result_b columns.
#'
#' @return A data frame with `log10_concentration` plus prefixed columns
#'   from each grid.
#'
#' @export
compare_calibrations <- function(result_a, result_b,
                                 label_a = NULL, label_b = NULL) {
  stopifnot(inherits(result_a, "calibration_result"),
            inherits(result_b, "calibration_result"))

  label_a <- label_a %||% result_a$meta$method
  label_b <- label_b %||% result_b$meta$method

  fg <- result_a$grid
  bg <- result_b$grid

  merge_col <- if ("log10_concentration" %in% names(fg) &&
                   "log10_concentration" %in% names(bg)) {
    "log10_concentration"
  } else if ("x_fit" %in% names(fg) && "x_fit" %in% names(bg)) {
    "x_fit"
  } else {
    stop("Cannot align grids: no common concentration column found.")
  }

  freq_cols  <- setdiff(names(fg), c(merge_col, "concentration", "x_fit"))
  bayes_cols <- setdiff(names(bg), c(merge_col, "concentration", "x_fit"))
  names(fg)[names(fg) %in% freq_cols]  <- paste0(label_a, "_", freq_cols)
  names(bg)[names(bg) %in% bayes_cols] <- paste0(label_b, "_", bayes_cols)

  merged <- merge(fg, bg, by = intersect(names(fg), names(bg)),
                  all = TRUE, sort = TRUE)

  merged[[paste0(label_a, "_model")]] <- result_a$selection$best_model_name
  merged[[paste0(label_b, "_model")]] <- result_b$selection$best_model_name

  class(merged) <- c("calibration_comparison", "data.frame")
  merged
}


#' Compare Parameters Between Two Calibration Results
#'
#' Produces a side-by-side parameter comparison table for the best model
#' from each result.
#'
#' @param result_a,result_b `calibration_result` objects.
#' @param plate_id_b Integer. For Bayesian results with multiple plates,
#'   which plate's parameters to extract. Default 1.
#'
#' @return Data frame with columns: `term`, `a_estimate`, `b_estimate`,
#'   `a_se`, `b_se`, `diff`, `rel_diff`.
#'
#' @export
compare_parameters <- function(result_a, result_b, plate_id_b = 1L) {
  stopifnot(inherits(result_a, "calibration_result"),
            inherits(result_b, "calibration_result"))

  label_a <- result_a$meta$method
  label_b <- result_b$meta$method

  # Extract parameters from result_a
  best_a <- result_a$selection$best_model_name
  ens_a  <- result_a$ensemble[[best_a]]
  pa <- if (!is.null(ens_a$parameters)) {
    if (is.data.frame(ens_a$parameters)) {
      # frequentist: single data frame
      ens_a$parameters
    } else if (is.list(ens_a$parameters) && length(ens_a$parameters) >= 1) {
      # bayesian: list of data frames per plate
      ens_a$parameters[[1]]
    }
  }

  # Extract parameters from result_b
  best_b <- result_b$selection$best_model_name
  ens_b  <- result_b$ensemble[[best_b]]
  pb <- if (!is.null(ens_b$parameters)) {
    if (is.data.frame(ens_b$parameters)) {
      ens_b$parameters
    } else if (is.list(ens_b$parameters)) {
      ens_b$parameters[[plate_id_b]]
    }
  }

  if (is.null(pa) || is.null(pb)) {
    warning("Cannot extract parameters from one or both results")
    return(data.frame(term = character(), stringsAsFactors = FALSE))
  }

  # Normalise column names
  est_col_a <- intersect(c("estimate", "mean"), names(pa))[1]
  se_col_a  <- intersect(c("std_error", "sd"), names(pa))[1]
  est_col_b <- intersect(c("estimate", "mean"), names(pb))[1]
  se_col_b  <- intersect(c("std_error", "sd"), names(pb))[1]

  pa_slim <- data.frame(term = pa$term,
                         a_estimate = pa[[est_col_a]],
                         a_se = if (!is.na(se_col_a)) pa[[se_col_a]] else NA_real_,
                         stringsAsFactors = FALSE)
  pb_slim <- data.frame(term = pb$term,
                         b_estimate = pb[[est_col_b]],
                         b_se = if (!is.na(se_col_b)) pb[[se_col_b]] else NA_real_,
                         stringsAsFactors = FALSE)

  out <- merge(pa_slim, pb_slim, by = "term", all = TRUE)
  out$diff <- out$b_estimate - out$a_estimate
  out$rel_diff <- ifelse(abs(out$a_estimate) > 1e-10,
                         out$diff / abs(out$a_estimate), NA_real_)

  names(out) <- sub("^a_", paste0(label_a, "_"), names(out))
  names(out) <- sub("^b_", paste0(label_b, "_"), names(out))

  out$model_a <- best_a
  out$model_b <- best_b
  names(out)[names(out) == "model_a"] <- paste0(label_a, "_model")
  names(out)[names(out) == "model_b"] <- paste0(label_b, "_model")

  out
}


#' Compare Sample Predictions Between Two Calibration Results
#'
#' Merges sample predictions by sample identifier and computes agreement.
#'
#' @param result_a,result_b `calibration_result` objects.
#' @param by Character vector of columns to merge on. Default
#'   `c("sampleid", "curve_id")`.
#'
#' @return Data frame with paired predictions and agreement columns.
#'
#' @export
compare_samples <- function(result_a, result_b,
                            by = c("sampleid", "curve_id")) {
  stopifnot(inherits(result_a, "calibration_result"),
            inherits(result_b, "calibration_result"))

  sa <- result_a$samples
  sb <- result_b$samples
  if (is.null(sa) || is.null(sb)) {
    warning("One or both results have no sample predictions")
    return(NULL)
  }

  label_a <- result_a$meta$method
  label_b <- result_b$meta$method

  # Identify which merge keys actually exist
  by <- intersect(by, intersect(names(sa), names(sb)))
  if (length(by) == 0) stop("No common merge keys in sample data")

  # Select prediction columns
  pred_cols <- c("predicted_log10_concentration", "final_concentration",
                 "se_concentration", "pcov")
  a_cols <- intersect(pred_cols, names(sa))
  b_cols <- intersect(pred_cols, names(sb))

  sa_slim <- sa[, c(by, a_cols), drop = FALSE]
  sb_slim <- sb[, c(by, b_cols), drop = FALSE]

  names(sa_slim)[names(sa_slim) %in% a_cols] <- paste0(label_a, "_", a_cols)
  names(sb_slim)[names(sb_slim) %in% b_cols] <- paste0(label_b, "_", b_cols)

  out <- merge(sa_slim, sb_slim, by = by, all = TRUE)

  # Compute agreement on log10 concentration
  fc_a <- paste0(label_a, "_predicted_log10_concentration")
  fc_b <- paste0(label_b, "_predicted_log10_concentration")
  if (fc_a %in% names(out) && fc_b %in% names(out)) {
    out$conc_diff <- out[[fc_b]] - out[[fc_a]]
    both_ok <- is.finite(out[[fc_a]]) & is.finite(out[[fc_b]])
    out$conc_abs_diff <- abs(out$conc_diff)
  }

  out
}


#' Compute Agreement Metrics Between Paired Predictions
#'
#' Given two numeric vectors of paired predictions, computes bias,
#' MAE, RMSE, Pearson correlation, and Lin's concordance correlation.
#'
#' @param x,y Numeric vectors (same length). Paired predictions.
#' @param na.rm Logical. Remove NAs. Default TRUE.
#'
#' @return Named list: `n`, `bias`, `mae`, `rmse`, `cor`, `ccc`.
#'
#' @export
agreement_metrics <- function(x, y, na.rm = TRUE) {
  if (na.rm) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
  }
  n <- length(x)
  if (n < 2) return(list(n = n, bias = NA, mae = NA, rmse = NA,
                          cor = NA, ccc = NA))

  d    <- y - x
  bias <- mean(d)
  mae  <- mean(abs(d))
  rmse <- sqrt(mean(d^2))
  r    <- stats::cor(x, y)

  # Lin's concordance correlation coefficient
  mx <- mean(x); my <- mean(y)
  sx <- stats::sd(x); sy <- stats::sd(y)
  ccc <- (2 * r * sx * sy) / (sx^2 + sy^2 + (mx - my)^2)

  list(n = n, bias = bias, mae = mae, rmse = rmse, cor = r, ccc = ccc)
}
