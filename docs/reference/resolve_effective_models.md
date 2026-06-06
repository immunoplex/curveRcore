# Resolve Effective Models for a Given Configuration

Filters the requested model list based on whether the independent
variable is log-transformed (which makes loglogistic4 redundant with
logistic4).

## Usage

``` r
resolve_effective_models(fit_options, study_params)
```

## Arguments

- fit_options:

  A `fit_options` object.

- study_params:

  A `study_params` object.

## Value

Character vector of model names to actually fit.
