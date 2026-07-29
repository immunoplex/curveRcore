# Closed-form inflection point of a calibration model

Closed-form inflection point of a calibration model

## Usage

``` r
compute_inflection(model_name, params)
```

## Arguments

- model_name:

  Character. One of "logistic4", "logistic5", "gompertz4",
  "loglogistic4", "loglogistic5".

- params:

  Named numeric vector, or a data frame carrying a `term` column plus a
  value column (`mean` for Bayesian, `estimate` for frequentist). Must
  supply a, b, c, d (and g for the 5-parameter models).

## Value

A list with:

- x:

  Inflection x on the fitting scale (log10 concentration for
  log-independent fits).

- y:

  Response at the inflection (fitting/response scale).

- source:

  "analytic".

Returns `x = NA, y = NA` if parameters are missing or invalid.
