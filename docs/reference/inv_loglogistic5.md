# Inverse of the loglogistic5 Model

Solves for `x`: \$\$x = c - \frac{1}{b}\\\log\\\frac{ \left(\frac{y -
a}{d - a}\right)^{-g} - 1}{g}\$\$

## Usage

``` r
inv_loglogistic5(y, a, b, c, d, g)

inv_loglogistic5_fixed(y, fixed_a, b, c, d, g)
```

## Arguments

- y:

  Numeric vector. Observed response.

- a, b, c, d:

  Numeric scalars. Model parameters.

- g:

  Numeric scalar. Asymmetry parameter.

- fixed_a:

  Numeric scalar. Externally fixed lower asymptote.

## Value

Numeric vector of `x` values.

## See also

Other inverse-functions:
[`inv_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/inv_gompertz4.md),
[`inv_logistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic4.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic5.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic4.md)
