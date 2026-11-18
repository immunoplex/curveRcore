# Recovery and bias-aware precision of back-calculated standards

For calibrators, whose nominal concentration is known, computes the
recovery (estimated / nominal, in percent) and the CV of the
root-mean-square error, which adds the bias against the nominal value to
the standard error: `sqrt(se^2 + bias^2)` on the log10 scale, expressed
as a CV in percent the same way `pcov` is, and capped at `cv_x_max`.

## Usage

``` r
standard_recovery(
  predicted,
  nominal,
  se,
  is_log_independent = TRUE,
  cv_x_max = 150
)
```

## Arguments

- predicted:

  Back-calculated concentration on the fit's x scale
  (`predicted_concentration`).

- nominal:

  Nominal concentration on the same scale (the x value the fit used:
  log10 when `is_log_independent`).

- se:

  Standard error of `predicted`, same scale (`se_concentration`).

- is_log_independent:

  Logical: is x on the log10 scale?

- cv_x_max:

  Cap for the CV, as for `pcov` (default 150).

## Value

A data frame with `recovery_pct` and `pcov_rmse`, one row per point.
Both are `NA` where the estimate or the nominal value is missing; a
point with an estimate but no finite `se` gets `pcov_rmse = cv_x_max`.

## Examples

``` r
standard_recovery(predicted = c(1.02, 2.10), nominal = c(1, 2), se = c(0.03, 0.05))
#>   recovery_pct pcov_rmse
#> 1     104.7129  8.302089
#> 2     125.8925 25.743684
```
