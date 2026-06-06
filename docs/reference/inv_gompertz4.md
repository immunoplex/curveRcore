# Inverse of the Gompertz Model

Solves for `x`: \$\$x = c - \frac{1}{b}\\\log\\\left(-\log\frac{y -
a}{d - a}\right)\$\$

## Usage

``` r
inv_gompertz4(y, a, b, c, d)

inv_gompertz4_fixed(y, fixed_a, b, c, d)
```

## Arguments

- y:

  Numeric vector. Observed response.

- a, b, c, d:

  Numeric scalars. Model parameters.

- fixed_a:

  Numeric scalar. Externally fixed lower asymptote.

## Value

Numeric vector of `x` values.

## See also

Other inverse-functions:
[`inv_logistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic4.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic5.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic4.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic5.md)
