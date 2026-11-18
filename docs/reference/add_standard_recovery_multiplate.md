# Add recovery and bias-aware pcov_rmse to every curve's standard estimates

Call after the engine has back-calculated the standards
(`$standard_estimates` on each curve's result). Overwrites `pcov_rmse`
and adds `recovery_pct` there, using
[`standard_recovery()`](https://immunoplex.github.io/curveRcore/reference/standard_recovery.md).
The nominal value is the `concentration` column (the x value the fit
used). A result without `$standard_estimates` is returned unchanged.

## Usage

``` r
add_standard_recovery_multiplate(mp, is_log_independent = NULL)
```

## Arguments

- mp:

  A `calibration_result` or `calibration_result_multiplate`.

- is_log_independent:

  Logical; defaults to the value recorded in the result's `meta`, else
  `TRUE`.

## Value

The same class of object.
