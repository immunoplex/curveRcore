# =============================================================================
# derivatives.R — dy/dx and d²y/dx² for all five forward models
#
# All first derivatives are positive when b > 0 and a < d (increasing curve).
# Second derivatives are analytical where tractable, numerical (central
# difference) for the 5-parameter models.
# =============================================================================


# ---- logistic4 derivatives ----

#' First Derivative of the 4PL Model
#'
#' \deqn{\frac{dy}{dx} = \frac{(d - a)\,u}{b\,(1 + u)^2}
#'   \quad\text{where } u = \exp\!\left(-\frac{x - c}{b}\right)}
#'
#' @inheritParams logistic4
#' @return Numeric vector of dy/dx values.
#' @family derivatives
#' @export
dydx_logistic4 <- function(x, a, b, c, d) {
  u <- exp(-(x - c) / b)
  (d - a) * u / (b * (1 + u)^2)
}


#' Second Derivative of the 4PL Model
#'
#' \deqn{\frac{d^2y}{dx^2} = -\frac{(d - a)\,u\,(1 - u)}{b^2\,(1 + u)^3}}
#'
#' @inheritParams logistic4
#' @return Numeric vector of d²y/dx² values.
#' @keywords internal
#' @export
d2x_logistic4 <- function(x, a, b, c, d) {
  u <- exp(-(x - c) / b)
  -(d - a) * u * (1 - u) / (b^2 * (1 + u)^3)
}


# ---- logistic5 derivatives ----

#' First Derivative of the 5PL Model
#'
#' \deqn{\frac{dy}{dx} = \frac{g\,(d - a)\,u}{b\,(1 + u)^{g+1}}}
#'
#' @inheritParams logistic5
#' @return Numeric vector of dy/dx values.
#' @family derivatives
#' @export
dydx_logistic5 <- function(x, a, b, c, d, g) {
  u <- exp(-(x - c) / b)
  g * (d - a) * u / (b * (1 + u)^(g + 1))
}


#' Second Derivative of the 5PL Model (Numerical)
#'
#' Central-difference approximation.
#'
#' @inheritParams logistic5
#' @param h Step size. Default `1e-5`.
#' @return Numeric vector of d²y/dx² values.
#' @keywords internal
#' @export
d2x_logistic5 <- function(x, a, b, c, d, g, h = 1e-5) {
  (logistic5(x + h, a, b, c, d, g) - 2 * logistic5(x, a, b, c, d, g) +
     logistic5(x - h, a, b, c, d, g)) / h^2
}


# ---- loglogistic4 derivatives ----

#' First Derivative of the loglogistic4 Model
#'
#' \deqn{\frac{dy}{dx} = \frac{b\,(d - a)\,r}{x\,(1 + r)^2}
#'   \quad\text{where } r = (c/x)^b}
#'
#' @inheritParams loglogistic4
#' @return Numeric vector of dy/dx values.
#' @family derivatives
#' @export
dydx_loglogistic4 <- function(x, a, b, c, d) {
  r <- (c / x)^b
  b * (d - a) * r / (x * (1 + r)^2)
}


#' Second Derivative of the loglogistic4 Model (Numerical)
#'
#' @inheritParams loglogistic4
#' @param h Step size. Default `1e-5`.
#' @return Numeric vector of d²y/dx² values.
#' @keywords internal
#' @export
d2x_loglogistic4 <- function(x, a, b, c, d, h = 1e-5) {
  (loglogistic4(x + h, a, b, c, d) - 2 * loglogistic4(x, a, b, c, d) +
     loglogistic4(x - h, a, b, c, d)) / h^2
}


# ---- loglogistic5 derivatives ----

#' First Derivative of the loglogistic5 Model
#'
#' \deqn{\frac{dy}{dx} = b\,(d - a)\,\exp(-b(x-c))\,
#'   u^{-1/g - 1}
#'   \quad\text{where } u = 1 + g\,\exp(-b(x-c))}
#'
#' @inheritParams loglogistic5
#' @return Numeric vector of dy/dx values.
#' @family derivatives
#' @export
dydx_loglogistic5 <- function(x, a, b, c, d, g) {
  u <- 1 + g * exp(-b * (x - c))
  b * (d - a) * exp(-b * (x - c)) * u^(-1 / g - 1)
}


#' Second Derivative of the loglogistic5 Model (Numerical)
#'
#' @inheritParams loglogistic5
#' @param h Step size. Default `1e-5`.
#' @return Numeric vector of d²y/dx² values.
#' @keywords internal
#' @export
d2x_loglogistic5 <- function(x, a, b, c, d, g, h = 1e-5) {
  (loglogistic5(x + h, a, b, c, d, g) - 2 * loglogistic5(x, a, b, c, d, g) +
     loglogistic5(x - h, a, b, c, d, g)) / h^2
}


# ---- gompertz4 derivatives ----

#' First Derivative of the Gompertz Model
#'
#' \deqn{\frac{dy}{dx} = b\,(d - a)\,v\,\exp(-v)
#'   \quad\text{where } v = \exp(-b(x - c))}
#'
#' @inheritParams gompertz4
#' @return Numeric vector of dy/dx values.
#' @family derivatives
#' @export
dydx_gompertz4 <- function(x, a, b, c, d) {
  v <- exp(-b * (x - c))
  b * (d - a) * v * exp(-v)
}


#' Second Derivative of the Gompertz Model
#'
#' \deqn{\frac{d^2y}{dx^2} = b^2\,(d - a)\,v\,(v - 1)\,\exp(-v)
#'   \quad\text{where } v = \exp(-b(x - c))}
#'
#' @inheritParams gompertz4
#' @return Numeric vector of d²y/dx² values.
#' @keywords internal
#' @export
d2x_gompertz4 <- function(x, a, b, c, d) {
  v <- exp(-b * (x - c))
  b^2 * (d - a) * v * (v - 1) * exp(-v)
}
