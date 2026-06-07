# Generate a Prediction Grid of Concentrations

Builds a data frame with evenly-spaced (on the fitting scale)
concentration values spanning the full range from `grid_min_conc` to
`grid_max_conc`.

## Usage

``` r
generate_prediction_grid(
  std_curve_conc,
  n_grid = 200L,
  grid_min_conc = 1e-04,
  grid_max_conc = NULL,
  is_log_independent = TRUE
)
```

## Arguments

- std_curve_conc:

  Numeric scalar, or a list with `$standard_curve_concentration`
  (backwards compatibility with S3 antigen_constraints objects).

- n_grid:

  Integer, or a list with `$n_grid`, `$grid_min_conc`, `$grid_max_conc`
  (backwards compatibility with S3 fit_options objects). Default 200.

- grid_min_conc:

  Numeric. Minimum concentration (raw scale). Default 1e-4.

- grid_max_conc:

  Numeric or NULL. Maximum concentration (raw scale). NULL uses
  `std_curve_conc`.

- is_log_independent:

  Logical. If TRUE, the grid is generated on the log10 scale.

## Value

A data frame with columns `log10_concentration`, `concentration`,
`x_fit`.
