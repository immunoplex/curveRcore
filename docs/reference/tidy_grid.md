# Tidy the precision grid from a calibration result

Extracts the best-model `$grid` (precision profile) into a tidy data
frame, attaching `curve_id` for multiplate inputs.

## Usage

``` r
tidy_grid(x, model = NULL, ...)

# S3 method for class 'calibration_result'
tidy_grid(x, model = NULL, ...)

# S3 method for class 'calibration_result_multiplate'
tidy_grid(x, model = NULL, ...)
```

## Arguments

- x:

  A `calibration_result` or `calibration_result_multiplate`.

- model:

  Optional model name. `NULL` (default) uses each plate's selected best
  model (the top-level `$grid`). Otherwise pulls
  `ensemble[[model]]$grid`.

- ...:

  Unused; for method extensibility.

## Value

A data frame of grid rows with columns including `log10_concentration`,
`concentration`, `predicted_concentration`, `se_concentration`, `pcov`,
`pcov_rmse`, `pcov_pass`, `d2y_dx2`, and (for multiplate) `curve_id`.

## See also

[`tidy_samples()`](https://immunoplex.github.io/curveRcore/reference/tidy_samples.md)
