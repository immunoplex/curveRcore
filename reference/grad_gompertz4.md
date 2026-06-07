# Analytical Gradient of the Inverse Gompertz

Analytical Gradient of the Inverse Gompertz

## Usage

``` r
grad_gompertz4(y, a, b, c, d)
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
[`grad_logistic4()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic4.md),
[`grad_logistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic5.md),
[`grad_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic4.md),
[`grad_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic5.md)
