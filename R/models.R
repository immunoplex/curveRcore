# =============================================================================
# models.R — Five canonical forward model functions for calibration curves
#
# Convention (enforced across curveRcore, curveRfreq, curveRbayes):
#   a  = lower asymptote (baseline), always a < d
#   d  = upper asymptote (saturation), always d > a
#   b  = scale/slope parameter, always b > 0
#   c  = inflection-point location (on whatever x-scale is passed in)
#   g  = asymmetry parameter (5-param models), always g > 0
#   x  = independent variable (externally transformed; typically log10(conc))
#
# All five models are monotonically increasing from a to d when b > 0.
# The external transformation of x (e.g. log10) is done upstream in R;
# the models themselves are scale-agnostic.
# =============================================================================


#' Four-Parameter Logistic (4PL) Forward Function
#'
#' Computes the response for a four-parameter logistic curve:
#' \deqn{y = a + \frac{d - a}{1 + \exp\!\left(-\frac{x - c}{b}\right)}}
#'
#' Symmetric about its inflection point at \eqn{(c, (a+d)/2)}.
#' Always monotonically increasing when \eqn{b > 0} and \eqn{a < d}.
#'
#' @param x Numeric vector. Independent variable (typically log10-concentration).
#' @param a Numeric scalar. Lower asymptote (baseline response).
#' @param b Numeric scalar. Scale parameter (\eqn{b > 0}); controls steepness.
#' @param c Numeric scalar. Inflection-point location on the x-axis.
#' @param d Numeric scalar. Upper asymptote (saturation response).
#'
#' @return Numeric vector of predicted response values, same length as `x`.
#'
#' @family forward-models
#' @export
logistic4 <- function(x, a, b, c, d) {
  a + (d - a) / (1 + exp(-(x - c) / b))
}


#' Five-Parameter Logistic (5PL) Forward Function
#'
#' Extends [logistic4()] with an asymmetry parameter \eqn{g}:
#' \deqn{y = a + \frac{d - a}{\bigl(1 +
#'   \exp\!\bigl(-\frac{x - c}{b}\bigr)\bigr)^g}}
#'
#' When \eqn{g = 1} this reduces to [logistic4()].
#' \eqn{g > 1} skews toward the upper asymptote;
#' \eqn{0 < g < 1} skews toward the lower asymptote.
#'
#' @inheritParams logistic4
#' @param g Numeric scalar. Asymmetry parameter (\eqn{g > 0}).
#'
#' @return Numeric vector of predicted response values.
#'
#' @family forward-models
#' @export
logistic5 <- function(x, a, b, c, d, g) {
  a + (d - a) / (1 + exp(-(x - c) / b))^g
}


#' Four-Parameter Log-Logistic (Dose-Response) Forward Function
#'
#' Classical Hill-equation parameterisation for positive-valued x:
#' \deqn{y = a + \frac{d - a}{1 + (c / x)^b}}
#'
#' Here \eqn{c} is the EC50 (inflection on the concentration scale) and
#' \eqn{b > 0} is the Hill coefficient. Requires \eqn{x > 0} and
#' \eqn{c > 0}.
#'
#' @section Domain restriction:
#' This model requires positive \code{x}. When the independent variable
#' is log10-transformed (and can be negative), use [logistic4()] instead
#' — on the log scale the two models are algebraically equivalent.
#'
#' @inheritParams logistic4
#'
#' @return Numeric vector of predicted response values.
#'
#' @family forward-models
#' @export
loglogistic4 <- function(x, a, b, c, d) {
  a + (d - a) / (1 + (c / x)^b)
}


#' Five-Parameter Generalised Logistic (Richards) Forward Function
#'
#' A five-parameter generalised logistic curve with asymmetry:
#' \deqn{y = a + (d - a)\,\bigl(1 + g\,\exp(-b\,(x - c))\bigr)^{-1/g}}
#'
#' When \eqn{g = 1} this reduces to [logistic4()] (after rescaling \eqn{b}).
#' Works with any real \eqn{x} (including log10-transformed).
#'
#' @inheritParams logistic4
#' @param g Numeric scalar. Asymmetry (Richards) parameter (\eqn{g > 0}).
#'
#' @return Numeric vector of predicted response values.
#'
#' @family forward-models
#' @export
loglogistic5 <- function(x, a, b, c, d, g) {
  a + (d - a) * (1 + g * exp(-b * (x - c)))^(-1 / g)
}


#' Four-Parameter Gompertz Forward Function
#'
#' Computes the response for a Gompertz saturation curve:
#' \deqn{y = a + (d - a)\,\exp\!\bigl(-\exp(-b\,(x - c))\bigr)}
#'
#' The Gompertz is intrinsically asymmetric (skewed toward the upper
#' asymptote) without requiring a fifth parameter.
#'
#' @inheritParams logistic4
#'
#' @return Numeric vector of predicted response values.
#'
#' @family forward-models
#' @export
gompertz4 <- function(x, a, b, c, d) {
  a + (d - a) * exp(-exp(-b * (x - c)))
}
