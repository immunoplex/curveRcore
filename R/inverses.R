# =============================================================================
# inverses.R — Analytical inverse (back-calculation) functions
#
# Each inverse solves for x given observed response y and model parameters.
# Two variants per model:
#   inv_<model>(y, a, b, c, d [, g])         — a is free (estimated)
#   inv_<model>_fixed(y, fixed_a, b, c, d [, g]) — a is externally fixed
#
# All functions return NA for y values outside the open interval (a, d).
# =============================================================================


# ---- logistic4 inverse ----

#' Inverse of the 4PL Model
#'
#' Solves analytically for `x`:
#' \deqn{x = c + b\,\log\!\frac{y - a}{d - y}}
#'
#' @param y Numeric vector. Observed response values.
#' @param a,b,c,d Numeric scalars. Model parameters.
#' @param tol Numeric scalar. Buffer from asymptotes. Default `1e-6`.
#'
#' @return Numeric vector of `x` values. `NA` for out-of-range `y`.
#'
#' @family inverse-functions
#' @export
inv_logistic4 <- function(y, a, b, c, d, tol = 1e-6) {
  lo <- min(a, d) + tol
  hi <- max(a, d) - tol
  result <- rep(NA_real_, length(y))
  ok <- !is.na(y) & y > lo & y < hi
  if (any(ok))
    result[ok] <- c + b * log((y[ok] - a) / (d - y[ok]))
  result
}


#' @rdname inv_logistic4
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @keywords internal
#' @export
inv_logistic4_fixed <- function(y, fixed_a, b, c, d) {
  c + b * log((y - fixed_a) / (d - y))
}


# ---- logistic5 inverse ----

#' Inverse of the 5PL Model
#'
#' Solves for `x`:
#' \deqn{x = c - b\,\log\!\left(\left(\frac{d - a}{y - a}\right)^{1/g}
#'   - 1\right)}
#'
#' @param y Numeric vector. Observed response.
#' @param a,b,c,d Numeric scalars. Model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#' @param tol Numeric scalar. Buffer from asymptotes.
#'
#' @return Numeric vector of `x` values. `NA` for out-of-range `y`.
#'
#' @family inverse-functions
#' @export
inv_logistic5 <- function(y, a, b, c, d, g, tol = 1e-6) {
  lo <- min(a, d) + tol
  hi <- max(a, d) - tol
  result <- rep(NA_real_, length(y))
  ok <- !is.na(y) & y > lo & y < hi
  if (any(ok))
    result[ok] <- c - b * log(((d - a) / (y[ok] - a))^(1 / g) - 1)
  result
}


#' @rdname inv_logistic5
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @keywords internal
#' @export
inv_logistic5_fixed <- function(y, fixed_a, b, c, d, g) {
  c - b * log(((d - fixed_a) / (y - fixed_a))^(1 / g) - 1)
}


# ---- loglogistic4 inverse ----

#' Inverse of the loglogistic4 Model
#'
#' Solves for `x` (requires `x > 0`):
#' \deqn{x = \frac{c}{\left(\frac{d - y}{y - a}\right)^{1/b}}}
#'
#' @param y Numeric vector. Observed response.
#' @param a,b,c,d Numeric scalars. Model parameters.
#'
#' @return Numeric vector of `x` values.
#'
#' @family inverse-functions
#' @export
inv_loglogistic4 <- function(y, a, b, c, d) {
  c / ((d - y) / (y - a))^(1 / b)
}


#' @rdname inv_loglogistic4
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @keywords internal
#' @export
inv_loglogistic4_fixed <- function(y, fixed_a, b, c, d) {
  c / ((d - y) / (y - fixed_a))^(1 / b)
}


# ---- loglogistic5 inverse ----

#' Inverse of the loglogistic5 Model
#'
#' Solves for `x`:
#' \deqn{x = c - \frac{1}{b}\,\log\!\frac{
#'   \left(\frac{y - a}{d - a}\right)^{-g} - 1}{g}}
#'
#' @param y Numeric vector. Observed response.
#' @param a,b,c,d Numeric scalars. Model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return Numeric vector of `x` values.
#'
#' @family inverse-functions
#' @export
inv_loglogistic5 <- function(y, a, b, c, d, g) {
  c - (1 / b) * (log(((y - a) / (d - a))^(-g) - 1) - log(g))
}


#' @rdname inv_loglogistic5
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @keywords internal
#' @export
inv_loglogistic5_fixed <- function(y, fixed_a, b, c, d, g) {
  c - (1 / b) * (log(((y - fixed_a) / (d - fixed_a))^(-g) - 1) - log(g))
}


# ---- gompertz4 inverse ----

#' Inverse of the Gompertz Model
#'
#' Solves for `x`:
#' \deqn{x = c - \frac{1}{b}\,\log\!\left(-\log\frac{y - a}{d - a}\right)}
#'
#' @param y Numeric vector. Observed response.
#' @param a,b,c,d Numeric scalars. Model parameters.
#'
#' @return Numeric vector of `x` values.
#'
#' @family inverse-functions
#' @export
inv_gompertz4 <- function(y, a, b, c, d) {
  c - (1 / b) * log(-log((y - a) / (d - a)))
}


#' @rdname inv_gompertz4
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @keywords internal
#' @export
inv_gompertz4_fixed <- function(y, fixed_a, b, c, d) {
  c - (1 / b) * log(-log((y - fixed_a) / (d - fixed_a)))
}
