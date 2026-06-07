# Analytical Gradient of the Inverse 5PL

Analytical Gradient of the Inverse 5PL

## Usage

``` r
grad_logistic5(y, a, b, c, d, g)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- a, b, c, d:

  Numeric scalars. Free model parameters.

- g:

  Numeric scalar. Asymmetry parameter.

## Value

List with `grad_theta` and scalar `grad_y`.

## See also

Other gradient-functions:
[`grad_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/grad_gompertz4.md),
[`grad_logistic4()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic4.md),
[`grad_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic4.md),
[`grad_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic5.md)
