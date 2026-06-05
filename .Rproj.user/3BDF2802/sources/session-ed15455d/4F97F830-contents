# =============================================================================
# gradients.R — Analytical gradients of inverse functions for error propagation
#
# Three families:
#   grad_<model>(y, a, b, c, d [, g])  — full gradient, a free
#   grad_inv_<model>_fixed(...)        — gradient with a fixed externally
#   grad_y_<model>_fixed(...)          — ∂x/∂y with a fixed externally
#
# Plus the dispatch closure builder: make_inv_and_grad_fixed()
# =============================================================================


# ============================================================================
# FULL GRADIENT FUNCTIONS (a free) — return grad_theta and grad_y
# ============================================================================

#' Analytical Gradient of the Inverse 4PL
#'
#' Returns \eqn{\partial x / \partial \theta} and \eqn{\partial x / \partial y}.
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Free model parameters.
#' @return List with `grad_theta` (named vector) and scalar `grad_y`.
#' @family gradient-functions
#' @export
grad_logistic4 <- function(y, a, b, c, d) {
  # x = c + b * log((y - a) / (d - y))
  da <- -b / (y - a)
  db <-  log((y - a) / (d - y))
  dc <-  1
  dd <- -b / (d - y)
  dy <-  b * (d - a) / ((y - a) * (d - y))
  list(grad_theta = c(a = da, b = db, c = dc, d = dd), grad_y = dy)
}


#' Analytical Gradient of the Inverse loglogistic4
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Free model parameters.
#' @return List with `grad_theta` and scalar `grad_y`.
#' @family gradient-functions
#' @export
grad_loglogistic4 <- function(y, a, b, c, d) {
  # x = c / Q^(1/b)  where Q = (d - y) / (y - a)
  Q <- (d - y) / (y - a)
  p <- 1 / b
  x_val <- c / Q^p

  da <-  x_val * p * (d - y) / ((y - a) * (d - y))  # simplify below
  # ∂Q/∂a = (d-y)/(y-a)^2,  ∂x/∂a = -c*p * Q^(-p-1) * ∂Q/∂a
  dQda <- (d - y) / (y - a)^2
  da   <- -(c * p) * Q^(-p - 1) * dQda

  db   <- x_val * log(Q) / b^2   # from d/db of Q^(-1/b) = -Q^(-1/b)*log(Q)*(-1/b^2)

  dc   <- 1 / Q^p   # = x_val / c

  dQdd <- -1 / (y - a)
  dd   <- -(c * p) * Q^(-p - 1) * dQdd

  dQdy <- (-(y - a) - (d - y)) / (y - a)^2  # = -(d-a)/(y-a)^2
  dy   <- -(c * p) * Q^(-p - 1) * dQdy

  list(grad_theta = c(a = da, b = db, c = dc, d = dd), grad_y = dy)
}


#' Analytical Gradient of the Inverse Gompertz
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Free model parameters.
#' @return List with `grad_theta` and scalar `grad_y`.
#' @family gradient-functions
#' @export
grad_gompertz4 <- function(y, a, b, c, d) {
  # x = c - (1/b) * log(-log(R))  where R = (y-a)/(d-a)
  R  <- (y - a) / (d - a)
  L  <- -log(R)  # > 0 inside the curve's range
  da <-  1 / (b * L * (y - a))   # via chain rule through R
  db <-  log(L) / b^2
  dc <-  1
  dd <- -1 / (b * L * (d - a))
  dy <-  1 / (b * L * (y - a))   # note: same magnitude as da
  list(grad_theta = c(a = da, b = db, c = dc, d = dd), grad_y = dy)
}


#' Analytical Gradient of the Inverse 5PL
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Free model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#' @return List with `grad_theta` and scalar `grad_y`.
#' @family gradient-functions
#' @export
grad_logistic5 <- function(y, a, b, c, d, g) {
  # x = c - b * log(W)  where W = ((d-a)/(y-a))^(1/g) - 1
  T_val <- (d - a) / (y - a)
  Tg    <- T_val^(1 / g)
  W     <- Tg - 1

  da <-  b * Tg / (g * W * (y - a))          # via ∂T/∂a chain
  db <- -log(W)
  dc <-  1
  dd <- -b * Tg / (g * W * (d - a))          # via ∂T/∂d
  dg <-  b * Tg * log(T_val) / (g^2 * W)
  dy <- -b * (d - a) * Tg / (g * W * (y - a)^2 * T_val)

  # Simplify dy: Tg/T_val = T_val^(1/g - 1)
  dy <- -b * Tg / (g * W * (y - a))

  list(grad_theta = c(a = da, b = db, c = dc, d = dd, g = dg), grad_y = dy)
}


#' Analytical Gradient of the Inverse loglogistic5
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Free model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#' @return List with `grad_theta` and scalar `grad_y`.
#' @family gradient-functions
#' @export
grad_loglogistic5 <- function(y, a, b, c, d, g) {
  # x = c - (1/b) * (log(V) - log(g))  where V = ((y-a)/(d-a))^(-g) - 1
  ratio   <- (y - a) / (d - a)
  ratio_g <- ratio^(-g)
  V       <- ratio_g - 1

  da <-  (1 / (b * V * g)) * ratio^(-g - 1) / (d - a)
  db <-  (log(V) - log(g)) / b^2
  dc <-  1
  dd <-  (1 / (b * V * g)) * ratio^(-g - 1) * (y - a) / (d - a)^2
  dg <-  (1 / (b * g^2)) * (1 - log(g * V))
  dy <-  ratio^(-g - 1) / (b * g * V * (d - a))

  list(grad_theta = c(a = da, b = db, c = dc, d = dd, g = dg), grad_y = dy)
}


# ============================================================================
# FIXED-a GRADIENT FUNCTIONS — partials for free params only (b, c, d [, g])
# ============================================================================

#' @keywords internal
#' @noRd
grad_inv_logistic4_fixed <- function(y, fixed_a, b, c, d) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  db <- log((y - fixed_a) / (d - y))
  dc <- 1.0
  dd <- -b / (d - y)
  c(b = db, c = dc, d = dd)
}


#' @keywords internal
#' @noRd
grad_inv_loglogistic4_fixed <- function(y, fixed_a, b, c, d) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  Q  <- (d - y) / (y - fixed_a)
  p  <- 1 / b
  db <- (c / Q^p) * log(Q) / b^2
  dc <- 1 / Q^p
  dd <- (c * p) * Q^(-p - 1) / (y - fixed_a)
  c(b = db, c = dc, d = dd)
}


#' @keywords internal
#' @noRd
grad_inv_gompertz4_fixed <- function(y, fixed_a, b, c, d) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  L  <- -log((y - fixed_a) / (d - fixed_a))
  db <-  log(L) / b^2
  dc <-  1.0
  dd <- -1.0 / (b * L * (d - fixed_a))
  c(b = db, c = dc, d = dd)
}


#' @keywords internal
#' @noRd
grad_inv_logistic5_fixed <- function(y, fixed_a, b, c, d, g) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c)
  d <- as.numeric(d); g <- as.numeric(g)
  T_val <- (d - fixed_a) / (y - fixed_a)
  Tg    <- T_val^(1 / g)
  W     <- Tg - 1
  db <- -log(W)
  dc <-  1.0
  dd <- -b * Tg / (g * W * (d - fixed_a))
  dg <-  b * Tg * log(T_val) / (g^2 * W)
  c(b = db, c = dc, d = dd, g = dg)
}


#' @keywords internal
#' @noRd
grad_inv_loglogistic5_fixed <- function(y, fixed_a, b, c, d, g) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c)
  d <- as.numeric(d); g <- as.numeric(g)
  ratio   <- (y - fixed_a) / (d - fixed_a)
  ratio_g <- ratio^(-g)
  V       <- ratio_g - 1
  db <-  (log(V) - log(g)) / b^2
  dc <-  1.0
  dd <- -(1 / b) * ratio_g / V * g / (d - fixed_a)
  dg <-  (1 / b) * (ratio_g / V * log(ratio) - 1 / g)
  c(b = db, c = dc, d = dd, g = dg)
}


# ============================================================================
# FIXED-a grad_y FUNCTIONS (∂x/∂y with a externally fixed)
# ============================================================================

#' @keywords internal
#' @noRd
grad_y_logistic4_fixed <- function(y, fixed_a, b, d) {
  b * (d - fixed_a) / ((y - fixed_a) * (d - y))
}


#' @keywords internal
#' @noRd
grad_y_loglogistic4_fixed <- function(y, fixed_a, b, c, d) {
  Q <- (d - y) / (y - fixed_a)
  p <- 1 / b
  c * p * Q^(-p - 1) * (d - fixed_a) / (y - fixed_a)^2
}


#' @keywords internal
#' @noRd
grad_y_gompertz4_fixed <- function(y, fixed_a, b, d) {
  R <- (y - fixed_a) / (d - fixed_a)
  L <- -log(R)
  1.0 / (b * L * (y - fixed_a))
}


#' @keywords internal
#' @noRd
grad_y_logistic5_fixed <- function(y, fixed_a, b, d, g) {
  T_val <- (d - fixed_a) / (y - fixed_a)
  Tg    <- T_val^(1 / g)
  W     <- Tg - 1
  b * Tg / (g * W * (y - fixed_a))
}


#' @keywords internal
#' @noRd
grad_y_loglogistic5_fixed <- function(y, fixed_a, b, d, g) {
  ratio   <- (y - fixed_a) / (d - fixed_a)
  ratio_g <- ratio^(-g)
  V       <- ratio_g - 1
  g / (b * (d - fixed_a)) * ratio^(-g - 1) / V
}


# ============================================================================
# DISPATCH: make_inv_and_grad_fixed
# ============================================================================

#' Build Inverse, Gradient, and grad_y Closures for a Model
#'
#' Returns three closures (`inv`, `grad`, `grad_y`) that evaluate the inverse
#' function, parameter gradient, and response-derivative for a given model.
#' Used by error propagation via the delta method.
#'
#' @param model Character. One of `"logistic4"`, `"logistic5"`,
#'   `"loglogistic4"`, `"loglogistic5"`, `"gompertz4"`.
#' @param fixed_a Numeric scalar or NULL. If non-NULL, `a` is treated as
#'   a known constant and excluded from the parameter gradient.
#'
#' @return A list with closures: `inv(y, p)`, `grad(y, p)`, `grad_y(y, p)`.
#'   Each accepts a response value `y` and named parameter vector `p`
#'   (from `coef(fit)`).
#'
#' @export
make_inv_and_grad_fixed <- function(model, fixed_a = NULL) {

  # ── Branch A: a is a fixed external constant ──
  if (!is.null(fixed_a)) {
    fixed_a <- as.numeric(fixed_a)
    return(switch(model,
                  logistic4 = list(
                    inv    = function(y, p) inv_logistic4_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]]),
                    grad   = function(y, p) grad_inv_logistic4_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]]),
                    grad_y = function(y, p) grad_y_logistic4_fixed(y, fixed_a, p[["b"]], p[["d"]])
                  ),
                  loglogistic4 = list(
                    inv    = function(y, p) inv_loglogistic4_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]]),
                    grad   = function(y, p) grad_inv_loglogistic4_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]]),
                    grad_y = function(y, p) grad_y_loglogistic4_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]])
                  ),
                  gompertz4 = list(
                    inv    = function(y, p) inv_gompertz4_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]]),
                    grad   = function(y, p) grad_inv_gompertz4_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]]),
                    grad_y = function(y, p) grad_y_gompertz4_fixed(y, fixed_a, p[["b"]], p[["d"]])
                  ),
                  logistic5 = list(
                    inv    = function(y, p) inv_logistic5_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]], p[["g"]]),
                    grad   = function(y, p) grad_inv_logistic5_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]], p[["g"]]),
                    grad_y = function(y, p) grad_y_logistic5_fixed(y, fixed_a, p[["b"]], p[["d"]], p[["g"]])
                  ),
                  loglogistic5 = list(
                    inv    = function(y, p) inv_loglogistic5_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]], p[["g"]]),
                    grad   = function(y, p) grad_inv_loglogistic5_fixed(y, fixed_a, p[["b"]], p[["c"]], p[["d"]], p[["g"]]),
                    grad_y = function(y, p) grad_y_loglogistic5_fixed(y, fixed_a, p[["b"]], p[["d"]], p[["g"]])
                  ),
                  stop("Unsupported model: ", model)
    ))
  }

  # ── Branch B: a is FREE ──
  switch(model,
         logistic4 = list(
           inv    = function(y, p) inv_logistic4_fixed(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]]),
           grad   = function(y, p) grad_logistic4(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]])$grad_theta,
           grad_y = function(y, p) grad_logistic4(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]])$grad_y
         ),
         loglogistic4 = list(
           inv    = function(y, p) inv_loglogistic4_fixed(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]]),
           grad   = function(y, p) grad_loglogistic4(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]])$grad_theta,
           grad_y = function(y, p) grad_loglogistic4(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]])$grad_y
         ),
         gompertz4 = list(
           inv    = function(y, p) inv_gompertz4_fixed(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]]),
           grad   = function(y, p) grad_gompertz4(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]])$grad_theta,
           grad_y = function(y, p) grad_gompertz4(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]])$grad_y
         ),
         logistic5 = list(
           inv    = function(y, p) inv_logistic5_fixed(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]], p[["g"]]),
           grad   = function(y, p) grad_logistic5(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]], p[["g"]])$grad_theta,
           grad_y = function(y, p) grad_logistic5(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]], p[["g"]])$grad_y
         ),
         loglogistic5 = list(
           inv    = function(y, p) inv_loglogistic5_fixed(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]], p[["g"]]),
           grad   = function(y, p) grad_loglogistic5(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]], p[["g"]])$grad_theta,
           grad_y = function(y, p) grad_loglogistic5(y, p[["a"]], p[["b"]], p[["c"]], p[["d"]], p[["g"]])$grad_y
         ),
         stop("Unsupported model: ", model)
  )
}
