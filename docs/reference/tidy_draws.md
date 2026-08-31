# Tidy the posterior / sampling draws from a calibration result

Extracts per-parameter draw vectors (posterior draws for Bayesian fits;
asymptotic-MVN or bootstrap samples for frequentist) into a tidy frame
with one row per `(curve_id, term)` and a list-column `draws` holding
the vector, plus `n_draws` and `sample_kind`. Feeds the `calib_draws`
table (the worker packs `draws` into a `float8[]` column). Present only
when the fit was run with `persist_draws = TRUE`; otherwise returns a
zero-row frame.

## Usage

``` r
tidy_draws(x, ...)

# S3 method for class 'calibration_result'
tidy_draws(x, ...)

# S3 method for class 'calibration_result_multiplate'
tidy_draws(x, ...)
```

## Arguments

- x:

  A `calibration_result` or `calibration_result_multiplate`.

- ...:

  Unused; for method extensibility.

## Value

Data frame `curve_id, term, param_scope, sample_kind, n_draws, draws`
(`draws` is a list-column of numeric vectors).

## Details

Draw order is preserved and shared across terms within a fit, so
consumers may column-bind the vectors into a joint posterior (e.g. the
JOB-3 correlation matrix).

## See also

[`tidy_hyperparam()`](https://immunoplex.github.io/curveRcore/reference/tidy_hyperparam.md)
