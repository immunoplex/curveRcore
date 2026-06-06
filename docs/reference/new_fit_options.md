# Create Model Fitting and Grid Options

Constructs the options that control which models are fit, grid
resolution, and pcov capping.

## Usage

``` r
new_fit_options(
  model_names = available_models(),
  n_grid = 200L,
  cv_x_max = 150,
  grid_min_conc = 1e-04,
  grid_max_conc = NULL
)
```

## Arguments

- model_names:

  Character vector. Which models to fit. Default: all five canonical
  models.

- n_grid:

  Integer. Number of points in the prediction grid. Default 200.

- cv_x_max:

  Numeric. Cap for percent CV of predicted concentration. Default 150.

- grid_min_conc:

  Numeric. Minimum concentration for the grid (on the raw scale).
  Default 1e-4.

- grid_max_conc:

  Numeric or NULL. Maximum concentration. NULL uses the undiluted
  standard concentration from antigen_constraints.

## Value

A named list of class `fit_options`.
