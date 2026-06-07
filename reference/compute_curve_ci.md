# Compute Confidence Interval for Fitted Curve

Uses the delta method to compute pointwise confidence intervals.

## Usage

``` r
compute_curve_ci(
  grid,
  model_name,
  fit,
  level = 0.95,
  independent_variable = "concentration"
)
```

## Arguments

- grid:

  Data frame from
  [`generate_prediction_grid()`](https://immunoplex.github.io/curveRcore/reference/generate_prediction_grid.md).

- model_name:

  Character. Model name.

- fit:

  Fitted model object (must support
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html)).

- level:

  Numeric. Confidence level. Default 0.95.

- independent_variable:

  Character. Default `"concentration"`.

## Value

Data frame with columns `yhat`, `ci_lower`, `ci_upper`, `se_y`.
