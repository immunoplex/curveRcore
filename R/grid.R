# =============================================================================
# grid.R — Prediction grid generation
#
# Creates the concentration grid used for curve prediction, se_concentration
# and pcov profiling. Both curveRfreq and curveRbayes use the same grid.
# =============================================================================


#' Generate a Prediction Grid of Concentrations
#'
#' Builds a data frame with evenly-spaced (on the fitting scale) concentration
#' values spanning the full range from `grid_min_conc` to `grid_max_conc`.
#'
#' @param std_curve_conc Numeric scalar, or a list with
#'   `$standard_curve_concentration` (backwards compatibility with S3
#'   antigen_constraints objects).
#' @param n_grid Integer, or a list with `$n_grid`, `$grid_min_conc`,
#'   `$grid_max_conc` (backwards compatibility with S3 fit_options objects).
#'   Default 200.
#' @param grid_min_conc Numeric. Minimum concentration (raw scale).
#'   Default 1e-4.
#' @param grid_max_conc Numeric or NULL. Maximum concentration (raw scale).
#'   NULL uses `std_curve_conc`.
#' @param is_log_independent Logical. If TRUE, the grid is generated on the
#'   log10 scale.
#'
#' @return A data frame with columns `log10_concentration`, `concentration`,
#'   `x_fit`.
#'
#' @export
generate_prediction_grid <- function(std_curve_conc,
                                     n_grid = 200L,
                                     grid_min_conc = 1e-4,
                                     grid_max_conc = NULL,
                                     is_log_independent = TRUE) {

  # Backwards compatibility: S3 object dispatch
  if (is.list(std_curve_conc) &&
      !is.null(std_curve_conc$standard_curve_concentration)) {
    ac <- std_curve_conc
    fo <- if (is.list(n_grid) && !is.null(n_grid$n_grid)) n_grid else NULL
    return(generate_prediction_grid(
      std_curve_conc = ac$standard_curve_concentration,
      n_grid         = if (!is.null(fo)) fo$n_grid else 200L,
      grid_min_conc  = if (!is.null(fo)) fo$grid_min_conc else 1e-4,
      grid_max_conc  = if (!is.null(fo))
        (fo$grid_max_conc %||% ac$standard_curve_concentration)
      else ac$standard_curve_concentration,
      is_log_independent = is_log_independent
    ))
  }

  max_conc <- grid_max_conc %||% std_curve_conc
  stopifnot(grid_min_conc > 0, max_conc > grid_min_conc, n_grid >= 2)

  if (is_log_independent) {
    log10_conc <- seq(log10(grid_min_conc), log10(max_conc), length.out = n_grid)
    conc       <- 10^log10_conc
    x_fit      <- log10_conc
  } else {
    conc       <- seq(grid_min_conc, max_conc, length.out = n_grid)
    log10_conc <- log10(conc)
    x_fit      <- conc
  }

  data.frame(
    log10_concentration = log10_conc,
    concentration       = conc,
    x_fit               = x_fit,
    stringsAsFactors    = FALSE
  )
}


#' Compute Predicted Response for a Grid
#'
#' Evaluates a model forward function at every grid point.
#'
#' @param grid Data frame from [generate_prediction_grid()].
#' @param model_name Character. One of the five canonical model names.
#' @param params Named numeric vector of model parameters.
#'
#' @return Numeric vector of predicted responses.
#'
#' @export
predict_grid_response <- function(grid, model_name, params) {
  x <- grid$x_fit
  a <- params[["a"]]; b <- params[["b"]]
  c_val <- params[["c"]]; d <- params[["d"]]

  switch(model_name,
         logistic4    = logistic4(x, a, b, c_val, d),
         loglogistic4 = loglogistic4(x, a, b, c_val, d),
         gompertz4    = gompertz4(x, a, b, c_val, d),
         logistic5    = logistic5(x, a, b, c_val, d, params[["g"]]),
         loglogistic5 = loglogistic5(x, a, b, c_val, d, params[["g"]]),
         stop("Unknown model: ", model_name)
  )
}


#' Compute Confidence Interval for Fitted Curve
#'
#' Uses the delta method to compute pointwise confidence intervals.
#'
#' @param grid Data frame from [generate_prediction_grid()].
#' @param model_name Character. Model name.
#' @param fit Fitted model object (must support `coef()` and `vcov()`).
#' @param level Numeric. Confidence level. Default 0.95.
#' @param independent_variable Character. Default `"concentration"`.
#'
#' @return Data frame with columns `yhat`, `ci_lower`, `ci_upper`, `se_y`.
#'
#' @export
compute_curve_ci <- function(grid, model_name, fit, level = 0.95,
                             independent_variable = "concentration") {
  x <- grid$x_fit
  n_pts <- length(x)

  theta <- stats::coef(fit)
  V     <- stats::vcov(fit)
  p     <- length(theta)
  rhs   <- as.list(stats::formula(fit))[[3]]
  z     <- stats::qnorm(1 - (1 - level) / 2)

  yhat <- numeric(n_pts)
  se_y <- numeric(n_pts)

  for (i in seq_len(n_pts)) {
    nd <- stats::setNames(data.frame(x[i]), independent_variable)
    yhat[i] <- as.numeric(stats::predict(fit, newdata = nd))

    grad <- vapply(seq_len(p), function(j) {
      eps <- max(abs(theta[j]) * 1e-6, 1e-8)
      theta_up <- theta
      theta_up[j] <- theta[j] + eps
      env <- c(as.list(theta_up), as.list(nd))
      (as.numeric(eval(rhs, envir = env)) - yhat[i]) / eps
    }, numeric(1))

    se_y[i] <- sqrt(as.numeric(crossprod(grad, V %*% grad)))
  }

  data.frame(
    yhat     = yhat,
    ci_lower = yhat - z * se_y,
    ci_upper = yhat + z * se_y,
    se_y     = se_y
  )
}
