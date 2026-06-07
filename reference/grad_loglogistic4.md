# Analytical Gradient of the Inverse loglogistic4

Analytical Gradient of the Inverse loglogistic4

## Usage

``` r
grad_loglogistic4(y, a, b, c, d)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- a, b, c, d:

  Numeric scalars. Free model parameters.

## Value

List with `grad_theta` and scalar `grad_y`.

## See also

Other gradient-functions:
[`grad_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/grad_gompertz4.md),
[`grad_logistic4()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic4.md),
[`grad_logistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic5.md),
[`grad_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic5.md)
