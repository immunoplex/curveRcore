# Data contract for control and calibrator estimates

Describes `madi_results.calib_controls` (one row per control well, curve
and method; new in 0.7.0) and `madi_results.calib_standards` (one row
per standard point, curve and method; estimate columns added in 0.7.0):
column names, types, keys, the version that added each column, and a
short plain-language description. `R/calib_schema.sql` holds the DDL; a
package test keeps the two in sync. Design: `calib_points_SPEC.md`.

## Usage

``` r
point_estimates_contract(
  table = c("all", "controls", "standards", "control_reference"),
  estimates_only = FALSE
)
```

## Arguments

- table:

  `"controls"`, `"standards"`, `"control_reference"` or `"all"`
  (default).

- estimates_only:

  If `TRUE`, keep only the modelled-estimate columns (those shared with
  `calib_samples`, plus `recovery_pct`).

## Value

A data frame with columns `table`, `column`, `type`, `key`, `nullable`,
`added_in`, `description`.

## Details

The estimate columns mean the same as in `calib_samples`, with two
additions for standards, whose true value is known: `recovery_pct`, and
a `pcov_rmse` that includes the bias against the nominal value.

Also describes `madi_results.control_reference`, the expected values of
controls (owned by I-SPI; joined to the estimates by the
`calib_controls_qc` view).

## See also

[`validate_point_estimates()`](https://immunoplex.github.io/curveRcore/reference/validate_point_estimates.md),
[`weights_contract()`](https://immunoplex.github.io/curveRcore/reference/weights_contract.md)

## Examples

``` r
point_estimates_contract("standards", estimates_only = TRUE)[, c("column", "description")]
#>                           column
#> 1        predicted_concentration
#> 2  predicted_log10_concentration
#> 3            final_concentration
#> 4               se_concentration
#> 5                           pcov
#> 6                      pcov_rmse
#> 7                      pcov_pass
#> 8                pcov_gate_class
#> 9                   censor_class
#> 10                   censor_type
#> 11                  censor_lower
#> 12                  censor_upper
#> 13                     conc_q_lo
#> 14                    conc_q_med
#> 15                     conc_q_hi
#> 16                  p_below_lloq
#> 17                  p_above_uloq
#> 18                p_below_cutoff
#> 19            frac_draws_below_a
#> 20            frac_draws_above_d
#> 21                  recovery_pct
#>                                                                                                                                          description
#> 1                Back-calculated concentration at the dilution the well was run at, on the scale the fit used for x (log10 when is_log_independent).
#> 2                                                                         Back-calculated concentration, log10, at the dilution the well was run at.
#> 3        Back-calculated concentration multiplied by the dilution (natural scale). For a standard: that point's estimate of the stock concentration.
#> 4                                                         Standard error of the back-calculated concentration (log10 scale when is_log_independent).
#> 5                                                                         Precision CV (%) of the back-calculated concentration, capped at cv_x_max.
#> 6  CV (%) of the root-mean-square error. Standards: includes the bias against the nominal value. Samples and controls: equal to pcov (bias unknown).
#> 7                           TRUE if the point lies in the curve's dynamic range [LLOQ, ULOQ] (or passes the pcov threshold when LOQs are undefined).
#> 8                                                                                                            belowLLOQ, inDynamicRange or aboveULOQ.
#> 9                                                              quantified, lod_to_lloq, below_lod, below_lloq, above_uloq, saturated or no_response.
#> 10                                                                                                                    none, left, interval or right.
#> 11                                                                               Lower censoring bound, natural scale x dilution; NULL if unbounded.
#> 12                                                                               Upper censoring bound, natural scale x dilution; NULL if unbounded.
#> 13                                                                                           Lower quantile of the concentration draws (x dilution).
#> 14                                                                                                   Median of the concentration draws (x dilution).
#> 15                                                                                           Upper quantile of the concentration draws (x dilution).
#> 16                                                                                                                    Share of draws below the LLOQ.
#> 17                                                                                                                    Share of draws above the ULOQ.
#> 18                                                                                                         Share of draws below the decision cutoff.
#> 19                                                                    Share of draws whose response lies below the lower asymptote (not invertible).
#> 20                                                                    Share of draws whose response lies above the upper asymptote (not invertible).
#> 21                                                                                  100 x estimated / nominal concentration at the point's dilution.
```
