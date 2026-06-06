# Compare Sample Predictions Between Two Calibration Results

Merges sample predictions by sample identifier and computes agreement.

## Usage

``` r
compare_samples(result_a, result_b, by = c("sampleid", "curve_id"))
```

## Arguments

- result_a, result_b:

  `calibration_result` objects.

- by:

  Character vector of columns to merge on. Default
  `c("sampleid", "curve_id")`.

## Value

Data frame with paired predictions and agreement columns.
