# Add censoring columns to every curve's test samples

Call AFTER
[`compute_detection_limits_multiplate()`](https://immunoplex.github.io/curveRcore/reference/compute_detection_limits_multiplate.md)
and
[`classify_pcov_gates_multiplate()`](https://immunoplex.github.io/curveRcore/reference/classify_pcov_gates_multiplate.md)
– the first point at which LODs, LLOQ/ULOQ and samples all exist on the
result. Adds `censor_class`, `censor_type`, `censor_lower` and
`censor_upper` to `cr$samples` for every curve (all `NA` for a curve
with no selected model).

## Usage

``` r
classify_censoring_multiplate(mp)
```

## Arguments

- mp:

  A `calibration_result` or `calibration_result_multiplate`.

## Value

The same class of object with the censoring columns added.

## See also

[`classify_censoring()`](https://immunoplex.github.io/curveRcore/reference/classify_censoring.md)
