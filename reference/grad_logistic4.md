# Analytical Gradient of the Inverse 4PL

Returns \\\partial x / \partial \theta\\ and \\\partial x / \partial
y\\.

## Usage

``` r
grad_logistic4(y, a, b, c, d)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- a, b, c, d:

  Numeric scalars. Free model parameters.

## Value

List with `grad_theta` (named vector) and scalar `grad_y`.

## See also

Other gradient-functions:
[`grad_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/grad_gompertz4.md),
[`grad_logistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic5.md),
[`grad_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic4.md),
[`grad_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic5.md)
