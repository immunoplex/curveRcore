# Resolve the Response Column Name

Auto-detects `"mfi"` (bead arrays) or `"absorbance"` (ELISA) from the
column names of `df`. Falls back to `default`.

## Usage

``` r
resolve_response_col(df, default = "mfi")
```

## Arguments

- df:

  Data frame.

- default:

  Character. Fallback column name.

## Value

Character scalar: the detected response column name.
