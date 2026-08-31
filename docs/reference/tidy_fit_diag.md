# Tidy the per-fit diagnostics from a calibration result

Extracts sampler/optimizer diagnostics into a tidy frame keyed by
`curve_id`. Bayesian columns (`rhat_max`, `ess_bulk_min`,
`ess_tail_min`, `n_divergent`, `pct_divergent`, `max_treedepth_hit`,
`ebfmi_min`) and frequentist columns (`hessian_condition_number`,
`gradient_norm`, `optimizer_code`, `rel_tol_achieved`) coexist; the
irrelevant set is NA for a given engine. Common columns: `fit_seconds`,
`n_iterations`, `converged`, `fit_seed`. Feeds `calib_fit_diag`.

## Usage

``` r
tidy_fit_diag(x, ...)

# S3 method for class 'calibration_result'
tidy_fit_diag(x, ...)

# S3 method for class 'calibration_result_multiplate'
tidy_fit_diag(x, ...)
```

## Arguments

- x:

  A `calibration_result` or `calibration_result_multiplate`.

- ...:

  Unused; for method extensibility.

## Value

One row per `curve_id`.

## See also

[`tidy_hyperparam()`](https://immunoplex.github.io/curveRcore/reference/tidy_hyperparam.md)
