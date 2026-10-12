# Analytic inverse gradients must match central finite differences of the inverse functions.
num_grad <- function(fun, args, nm, h = 1e-6) {
  a1 <- args; a2 <- args; a1[[nm]] <- a1[[nm]] + h; a2[[nm]] <- a2[[nm]] - h
  (do.call(fun, a1) - do.call(fun, a2)) / (2 * h)
}
check_grad <- function(gfun, ifun, p, ys, fixed = FALSE) {
  for (y in ys) {
    g <- do.call(gfun, c(list(y = y), p))
    gt <- if (is.list(g)) g$grad_theta else g
    for (n in names(gt)) {
      nm <- if (fixed && n == "a") "fixed_a" else n
      if (!nm %in% names(p)) next
      expect_equal(unname(gt[[n]]), num_grad(ifun, c(list(y = y), p), nm), tolerance = 1e-5,
                   label = sprintf("%s d/d%s at y = %g", deparse(substitute(gfun)), n, y))
    }
    if (is.list(g) && !is.null(g$grad_y))
      expect_equal(g$grad_y, num_grad(ifun, c(list(y = y), p), "y"), tolerance = 1e-5,
                   label = sprintf("%s d/dy at y = %g", deparse(substitute(gfun)), y))
  }
}
ys <- c(1.6, 2.6, 3.8)
test_that("full gradients match finite differences", {
  check_grad(grad_logistic4,    inv_logistic4,    list(a = 1, b = 0.6, c = 2, d = 4.2), ys)
  check_grad(grad_logistic5,    inv_logistic5,    list(a = 1, b = 0.6, c = 2, d = 4.2, g = 0.4), ys)
  check_grad(grad_gompertz4,    inv_gompertz4,    list(a = 1, b = 0.6, c = 2, d = 4.2), ys)
  check_grad(grad_loglogistic4, inv_loglogistic4, list(a = 1, b = 1.5, c = 100, d = 4.2), ys)
  check_grad(grad_loglogistic5, inv_loglogistic5, list(a = 1, b = 1.5, c = 2, d = 4.2, g = 0.6), ys)
})
test_that("fixed-a gradients match finite differences", {
  check_grad(grad_inv_logistic4_fixed,    inv_logistic4_fixed,    list(fixed_a = 1, b = 0.6, c = 2, d = 4.2), ys, TRUE)
  check_grad(grad_inv_logistic5_fixed,    inv_logistic5_fixed,    list(fixed_a = 1, b = 0.6, c = 2, d = 4.2, g = 0.4), ys, TRUE)
  check_grad(grad_inv_gompertz4_fixed,    inv_gompertz4_fixed,    list(fixed_a = 1, b = 0.6, c = 2, d = 4.2), ys, TRUE)
  check_grad(grad_inv_loglogistic4_fixed, inv_loglogistic4_fixed, list(fixed_a = 1, b = 1.5, c = 100, d = 4.2), ys, TRUE)
  check_grad(grad_inv_loglogistic5_fixed, inv_loglogistic5_fixed, list(fixed_a = 1, b = 1.5, c = 2, d = 4.2, g = 0.6), ys, TRUE)
})

test_that("grad_y fixed-a variants match finite differences", {
  for (y in ys) {
    expect_equal(grad_y_logistic4_fixed(y, 1, 0.6, 4.2),       num_grad(inv_logistic4_fixed,    list(y = y, fixed_a = 1, b = 0.6, c = 2, d = 4.2), "y"), tolerance = 1e-5)
    expect_equal(grad_y_logistic5_fixed(y, 1, 0.6, 4.2, 0.4),  num_grad(inv_logistic5_fixed,    list(y = y, fixed_a = 1, b = 0.6, c = 2, d = 4.2, g = 0.4), "y"), tolerance = 1e-5)
    expect_equal(grad_y_gompertz4_fixed(y, 1, 0.6, 4.2),       num_grad(inv_gompertz4_fixed,    list(y = y, fixed_a = 1, b = 0.6, c = 2, d = 4.2), "y"), tolerance = 1e-5)
  }
})
test_that("gradients hold at other parameter sets, including g > 1", {
  for (g in c(0.3, 1.7, 3.5)) {
    check_grad(grad_logistic5,    inv_logistic5,    list(a = 0.5, b = 0.9, c = 1.5, d = 3.9, g = g), c(1.0, 2.2, 3.5))
    check_grad(grad_loglogistic5, inv_loglogistic5, list(a = 0.5, b = 1.2, c = 1.5, d = 3.9, g = g), c(1.0, 2.2, 3.5))
  }
  check_grad(grad_gompertz4, inv_gompertz4, list(a = 0.2, b = 1.4, c = 0.8, d = 3.1), c(0.6, 1.6, 2.9))
})
