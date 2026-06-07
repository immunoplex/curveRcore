# Compute Agreement Metrics Between Paired Predictions

Given two numeric vectors of paired predictions, computes bias, MAE,
RMSE, Pearson correlation, and Lin's concordance correlation.

## Usage

``` r
agreement_metrics(x, y, na.rm = TRUE)
```

## Arguments

- x, y:

  Numeric vectors (same length). Paired predictions.

- na.rm:

  Logical. Remove NAs. Default TRUE.

## Value

Named list: `n`, `bias`, `mae`, `rmse`, `cor`, `ccc`.
