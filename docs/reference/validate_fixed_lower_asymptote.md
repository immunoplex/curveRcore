# Validate a Fixed Lower Asymptote Before Log Transformation

Checks that the value is a positive, finite scalar suitable for
[`log10()`](https://rdrr.io/r/base/Log.html). Returns the value if
valid, NULL otherwise.

## Usage

``` r
validate_fixed_lower_asymptote(fixed_a_result_raw, verbose = FALSE)
```

## Arguments

- fixed_a_result_raw:

  Numeric scalar or NULL.

- verbose:

  Logical.

## Value

`fixed_a_result_raw` if valid; NULL otherwise.
