# Data contract for the precision-weight tables

Describes `madi_results.calib_weights` (one row per sample and method)
and `madi_results.calib_weights_fit` (one row per multiplate group and
method): column names, types, keys, which model each column applies to,
and a short plain-language description. `R/calib_schema.sql` holds the
DDL; this function is the same definition in a form code can check
against, and a package test keeps the two in sync.

## Usage

``` r
weights_contract(table = c("all", "sample", "fit"), model = NULL)
```

## Arguments

- table:

  `"sample"` (`calib_weights`), `"fit"` (`calib_weights_fit`) or `"all"`
  (default).

- model:

  Optional `"multiplicative"` or `"additive"`: keep only columns that
  apply to that model (plus those that apply to both).

## Value

A data frame with columns `table`, `column`, `type`, `key`, `nullable`,
`added_in`, `models`, `description`.

## One weight set per sample

Both tables keep their original primary keys, so the sample hierarchy is
the same as `calib_samples`. There is one set of weights per (sample,
method); the `model` column records which model produced it, and
re-running with a different model replaces the set. A column that does
not apply to the model that wrote the row is `NULL`.

## See also

[`validate_weights_table()`](https://immunoplex.github.io/curveRcore/reference/validate_weights_table.md)

## Examples

``` r
weights_contract("fit", model = "additive")[, c("column", "description")]
#>                 column
#> 1  multiplate_group_id
#> 2               method
#> 3              antigen
#> 4              feature
#> 5          design_cols
#> 6      scale_predictor
#> 7                beta1
#> 8             beta1_lo
#> 9             beta1_hi
#> 10      interpretation
#> 11               n_fit
#> 12               n_eff
#> 13        weight_ratio
#> 14              job_id
#> 15          created_at
#> 16               model
#> 17 weights_pkg_version
#> 18                 tau
#> 19              tau_lo
#> 20              tau_hi
#> 21               kappa
#> 22            kappa_lo
#> 23            kappa_hi
#> 24         reliability
#> 25      reliability_lo
#> 26      reliability_hi
#> 27           fix_beta1
#> 28          n_censored
#> 29     identifiability
#>                                                                                              description
#> 1                                         Group of plates of one antigen/feature/source fitted together.
#> 2                                                                               bayesian or frequentist.
#> 3                                                                                               Antigen.
#> 4                                                               Feature (isotype or functional readout).
#> 5                      Comma-joined design columns defining the saturated cells, e.g. timeperiod,agroup.
#> 6                                             Precision index used: se or pcov. Additive always uses se.
#> 7                 Exponent on se. Multiplicative: estimated. Additive: 1 when fix_beta1, else estimated.
#> 8                                                   Lower 95% credible limit of beta1 (NULL when fixed).
#> 9                                                   Upper 95% credible limit of beta1 (NULL when fixed).
#> 10 Plain-language reading of the fit (multiplicative: precision regime; additive: identifiability note).
#> 11                                                                         Observations used in the fit.
#> 12                                                                 Effective sample size of the weights.
#> 13                                                                            Largest / smallest weight.
#> 14                                                                       Compute job that wrote the row.
#> 15                                                                             Time the row was written.
#> 16                                                     multiplicative (every pre-0.6.0 row) or additive.
#> 17                                                          curveRweights version that produced the fit.
#> 18                                               SD of the true between-subject variation (log10 scale).
#> 19                                                                      Lower 95% credible limit of tau.
#> 20                                                                      Upper 95% credible limit of tau.
#> 21                      Calibration factor on the reported se; 1 means se is exactly the measurement SD.
#> 22                                                                    Lower 95% credible limit of kappa.
#> 23                                                                    Upper 95% credible limit of kappa.
#> 24                                 Overall reliability tau2 / (tau2 + mean v_meas) over quantified rows.
#> 25                                                      Lower 95% credible limit of overall reliability.
#> 26                                                      Upper 95% credible limit of overall reliability.
#> 27                                                   TRUE if beta1 was fixed at 1 rather than estimated.
#> 28                                                        Censored observations entering the likelihood.
#> 29                              Whether the signal/measurement split is set by the data or by the prior.
```
