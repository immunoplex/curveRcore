# Inverse of the 4PL Model

Solves analytically for `x`: \$\$x = c + b\\\log\\\frac{y - a}{d -
y}\$\$

## Usage

``` r
inv_logistic4(y, a, b, c, d, tol = 1e-06)

inv_logistic4_fixed(y, fixed_a, b, c, d)
```

## Arguments

- y:

  Numeric vector. Observed response values.

- a, b, c, d:

  Numeric scalars. Model parameters.

- tol:

  Numeric scalar. Buffer from asymptotes. Default `1e-6`.

- fixed_a:

  Numeric scalar. Externally fixed lower asymptote.

## Value

Numeric vector of `x` values. `NA` for out-of-range `y`.

## See also

Other inverse-functions:
[`inv_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/inv_gompertz4.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic5.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic4.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic5.md)
