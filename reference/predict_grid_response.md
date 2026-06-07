# Compute Predicted Response for a Grid

Evaluates a model forward function at every grid point.

## Usage

``` r
predict_grid_response(grid, model_name, params)
```

## Arguments

- grid:

  Data frame from
  [`generate_prediction_grid()`](https://immunoplex.github.io/curveRcore/reference/generate_prediction_grid.md).

- model_name:

  Character. One of the five canonical model names.

- params:

  Named numeric vector of model parameters.

## Value

Numeric vector of predicted responses.
