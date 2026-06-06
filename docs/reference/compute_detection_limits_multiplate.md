# Compute detection limits for all plates in a multiplate result

Iterates over `$plates` and calls
[`compute_detection_limits`](https://immunoplex.github.io/curveRcore/reference/compute_detection_limits.md)
on each.

## Usage

``` r
compute_detection_limits_multiplate(mp, alpha = 0.05, verbose = FALSE)
```

## Arguments

- mp:

  A `calibration_result_multiplate` object.

- alpha:

  Numeric; significance level (default 0.05).

- verbose:

  Logical; emit diagnostic messages (default `FALSE`).

## Value

The input `mp` with `$detection_limits` populated on every plate.
