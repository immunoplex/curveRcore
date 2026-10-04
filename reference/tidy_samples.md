# Tidy the per-sample predictions from a calibration result

Extracts the `$samples` table from a `calibration_result` or
`calibration_result_multiplate` into a single tidy data frame, attaching
`curve_id` for multiplate inputs. This is the canonical, supported way
for downstream packages (e.g. curveRweights) to read sample-level
concentration and precision; they must not reach into the object
internals directly.

## Usage

``` r
tidy_samples(x, ...)

# S3 method for class 'calibration_result'
tidy_samples(x, ...)

# S3 method for class 'calibration_result_multiplate'
tidy_samples(x, ...)
```

## Arguments

- x:

  A `calibration_result` or `calibration_result_multiplate`.

- ...:

  Unused; for method extensibility.

## Value

A data frame of the per-sample predictions. For multiplate input the
rows of every plate are row-bound with a `curve_id` column. Includes the
carried-through original sample columns plus `predicted_concentration`,
`se_concentration`, `pcov`, `pcov_pass`, etc. Returns a zero-row frame
if no samples are present.

## Passthrough contract for study-design columns

Any column present on the `samples` data frame supplied to the fitting
call (`curveRfreq::fit_calibration_freq()`/`_multiplate()`,
`curveRbayes::fit_calibration_bayes()`) survives into `$samples`
verbatim – neither
[`new_calibration_result()`](https://immunoplex.github.io/curveRcore/reference/new_calibration_result.md)/[`new_calibration_result_multiplate()`](https://immunoplex.github.io/curveRcore/reference/new_calibration_result_multiplate.md)
nor the fitters' sample-prediction step filter or whitelist columns;
they only ever *add* columns (`predicted_concentration`,
`se_concentration`, `pcov`, ...) to the frame they were given. This is
the supported mechanism for threading study-design metadata (e.g.
`timeperiod`, `agroup`/cohort arm) through to downstream consumers, most
notably
`curveRweights::as_weight_data(design = c("timeperiod", "agroup"))`'s
saturated cell-means grouping. No curveRcore/curveRfreq/curveRbayes code
change is needed to add a new design column – just ensure it is present
on `samples` before fitting.

## See also

[`tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.md),
[`pcov_from_se()`](https://immunoplex.github.io/curveRcore/reference/pcov_se_conversion.md)
