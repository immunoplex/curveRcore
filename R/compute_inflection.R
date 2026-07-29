# =============================================================================
# compute_inflection.R  (curveRcore)
#
# Exact, closed-form inflection point for each forward model, computed from the
# FITTED PARAMETERS rather than searched for on the prediction grid.
#
# WHY THIS EXISTS
#   The worker previously derived the inflection with
#       i <- which.min(abs(grid$d2y_dx2));  x <- grid$log10_concentration[i]
#   i.e. the grid row of minimum |d2y/dx2|. On a grid that starts at
#   grid_min_conc = 1e-4 (log10 = -4), the entire flat lower asymptote has
#   d2y/dx2 ≈ 0 to machine precision — smaller than the finite-difference |d2y|
#   at the true mid-curve inflection — so the search snaps to the LEFT EDGE of
#   the grid (≈ log10(1e-4) = -4). That is the source of the spurious
#   inflect_x ≈ -3.95 for curves whose data live near log10 conc 1.4–2.4.
#
#   The inflection of every model here is available in closed form, and the fit
#   already estimates the parameters, so there is no need to search a grid at
#   all. Setting d2y/dx2 = 0 gives:
#
#       logistic4     x* = c                      y* = (a + d)/2
#       gompertz4     x* = c                      y* = a + (d - a) e^{-1}
#       loglogistic5  x* = c                      y* = a + (d - a)(1 + g)^{-1/g}
#       logistic5     x* = c + b·ln(g)            y* = a + (d - a)(g/(g+1))^g
#       loglogistic4  x* = c (EC50, raw x-scale)  y* = (a + d)/2
#
#   (For g = 1 the logistic5 result collapses to x* = c, matching logistic4.)
#
#   x* is returned on the SAME x-scale the model was fitted on — i.e.
#   log10(concentration) for the log-independent fits this app uses — so it lines
#   up with calib_grid.log10_concentration, calib_standards, and the plot axis.
# =============================================================================


#' Closed-form inflection point of a calibration model
#'
#' @param model_name Character. One of "logistic4", "logistic5", "gompertz4",
#'   "loglogistic4", "loglogistic5".
#' @param params Named numeric vector, or a data frame carrying a `term` column
#'   plus a value column (`mean` for Bayesian, `estimate` for frequentist).
#'   Must supply a, b, c, d (and g for the 5-parameter models).
#'
#' @return A list with:
#'   \describe{
#'     \item{x}{Inflection x on the fitting scale (log10 concentration for
#'       log-independent fits).}
#'     \item{y}{Response at the inflection (fitting/response scale).}
#'     \item{source}{"analytic".}
#'   }
#'   Returns `x = NA, y = NA` if parameters are missing or invalid.
#'
#' @export
compute_inflection <- function(model_name, params) {

  pv <- .as_param_vector(params)
  if (is.null(pv)) return(list(x = NA_real_, y = NA_real_, source = "analytic"))

  a <- pv[["a"]]; b <- pv[["b"]]; c <- pv[["c"]]; d <- pv[["d"]]
  g <- if ("g" %in% names(pv)) pv[["g"]] else NA_real_

  # Require the core four to be finite; b must be non-zero to define a shape.
  if (any(!is.finite(c(a, b, c, d))) || abs(b) < .Machine$double.eps)
    return(list(x = NA_real_, y = NA_real_, source = "analytic"))

  x_star <- switch(model_name,
    logistic4    = c,
    gompertz4    = c,
    loglogistic5 = c,
    loglogistic4 = c,                                   # EC50 on the raw x-scale
    logistic5    = if (is.finite(g) && g > 0) c + b * log(g) else NA_real_,
    NA_real_
  )
  if (!is.finite(x_star))
    return(list(x = NA_real_, y = NA_real_, source = "analytic"))

  # Evaluate the response at x_star with the model's own forward function so y is
  # always internally consistent with x (no hand-transcribed y formulas to drift).
  y_star <- tryCatch(switch(model_name,
    logistic4    = logistic4(x_star, a, b, c, d),
    logistic5    = logistic5(x_star, a, b, c, d, g),
    gompertz4    = gompertz4(x_star, a, b, c, d),
    loglogistic4 = loglogistic4(x_star, a, b, c, d),
    loglogistic5 = loglogistic5(x_star, a, b, c, d, g),
    NA_real_
  ), error = function(e) NA_real_)

  list(x = as.numeric(x_star), y = as.numeric(y_star), source = "analytic")
}


#' Coerce a params vector/data frame to a named numeric vector a/b/c/d(/g)
#' @keywords internal
.as_param_vector <- function(params) {
  if (is.numeric(params) && !is.null(names(params))) {
    v <- params
  } else if (is.data.frame(params) && "term" %in% names(params)) {
    val_col <- if ("mean" %in% names(params)) "mean"
               else if ("estimate" %in% names(params)) "estimate"
               else return(NULL)
    v <- stats::setNames(params[[val_col]], params$term)
  } else {
    return(NULL)
  }
  # Normalise the location term name: some tables use "c", the Stan draws use
  # "c_par". Accept either.
  if (!("c" %in% names(v)) && "c_par" %in% names(v))
    names(v)[names(v) == "c_par"] <- "c"
  if (!all(c("a", "b", "c", "d") %in% names(v))) return(NULL)
  v
}
