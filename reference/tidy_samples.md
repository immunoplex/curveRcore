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

## See also

[`tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.md),
[`pcov_from_se()`](https://immunoplex.github.io/curveRcore/reference/pcov_se_conversion.md)
