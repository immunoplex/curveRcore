# Tidy the population / noise parameters from a calibration result

Extracts group-level (pooled) or per-plate population/noise parameters
(`mu_*`, `sigma_*`, `sigma_obs`, `nu`, and — frequentist —
`sigma_resid`) into a tidy data frame keyed by `curve_id`,
`param_scope = "population"`. Feeds the `calib_hyperparam` table.
Returns a zero-row frame when the fit carries no population slot (e.g. a
bare frequentist per-plate fit whose only noise term rides in
`calib_param`).

## Usage

``` r
tidy_hyperparam(x, ...)

# S3 method for class 'calibration_result'
tidy_hyperparam(x, ...)

# S3 method for class 'calibration_result_multiplate'
tidy_hyperparam(x, ...)
```

## Arguments

- x:

  A `calibration_result` or `calibration_result_multiplate`.

- ...:

  Unused; for method extensibility.

## Value

Data frame
`curve_id, term, param_scope, estimate, std_error, q_lo, q_med, q_hi`.

## See also

[`tidy_draws()`](https://immunoplex.github.io/curveRcore/reference/tidy_draws.md),
[`tidy_fit_diag()`](https://immunoplex.github.io/curveRcore/reference/tidy_fit_diag.md)
