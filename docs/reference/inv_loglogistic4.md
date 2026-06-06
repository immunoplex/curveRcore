# Inverse of the loglogistic4 Model

Solves for `x` (requires `x > 0`): \$\$x = \frac{c}{\left(\frac{d -
y}{y - a}\right)^{1/b}}\$\$

## Usage

``` r
inv_loglogistic4(y, a, b, c, d)

inv_loglogistic4_fixed(y, fixed_a, b, c, d)
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
[`inv_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/inv_gompertz4.md),
[`inv_logistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic4.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic5.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic5.md)
