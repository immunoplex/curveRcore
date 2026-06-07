# Inverse of the 5PL Model

Solves for `x`: \$\$x = c - b\\\log\\\left(\left(\frac{d - a}{y -
a}\right)^{1/g} - 1\right)\$\$

## Usage

``` r
inv_logistic5(y, a, b, c, d, g, tol = 1e-06)

inv_logistic5_fixed(y, fixed_a, b, c, d, g)
```

## Arguments

- y:

  Numeric vector. Observed response.

- a, b, c, d:

  Numeric scalars. Model parameters.

- g:

  Numeric scalar. Asymmetry parameter.

- tol:

  Numeric scalar. Buffer from asymptotes.

- fixed_a:

  Numeric scalar. Externally fixed lower asymptote.

## Value

Numeric vector of `x` values. `NA` for out-of-range `y`.

## See also

Other inverse-functions:
[`inv_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/inv_gompertz4.md),
[`inv_logistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic4.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic4.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic5.md)
