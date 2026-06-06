# Compare Two Calibration Results (Grid Predictions)

Merges the prediction grids from two `calibration_result` objects
(typically one frequentist and one Bayesian) into a single data frame
for paired analysis. Columns are prefixed with the method name.

## Usage

``` r
compare_calibrations(result_a, result_b, label_a = NULL, label_b = NULL)
```

## Arguments

- result_a:

  A `calibration_result`.

- result_b:

  A `calibration_result`.

- label_a:

  Character prefix for result_a columns. Default uses
  `result_a$meta$method`.

- label_b:

  Character prefix for result_b columns.

## Value

A data frame with `log10_concentration` plus prefixed columns from each
grid.
