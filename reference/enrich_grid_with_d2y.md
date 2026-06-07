# Add a d2y_dx2 column to an existing prediction grid

Computes \\d^2(\log\_{10} y) / d(\log\_{10} x)^2\\ from the
`log10_concentration` and `predicted_response` columns already present
on the grid, using non-uniform central differences.

## Usage

``` r
enrich_grid_with_d2y(grid, is_log_response = TRUE)
```

## Arguments

- grid:

  Data frame with columns `log10_concentration` and
  `predicted_response`.

- is_log_response:

  Logical; whether `predicted_response` is already on the \\\log\_{10}\\
  scale (default `TRUE`).

## Value

The input `grid` with an additional column `d2y_dx2`. Boundary points
(first and last rows) receive `NA`.

## Details

This function is designed to be called inside `predict_grid_freq()` or
`predict_grid_bayes()` immediately after the grid data frame is
constructed, adding the column in-place with zero additional model
evaluations.

When `is_log_response = TRUE`, the `predicted_response` column is
already on the \\\log\_{10}\\ scale, so `d2y_dx2` is computed directly.
When `is_log_response = FALSE`, the column is first
\\\log\_{10}\\-transformed before differencing.
