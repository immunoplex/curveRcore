# Version of the calib schema defined by this package

The version recorded in `madi_results.schema_version` (component
`calib_schema`) when `R/calib_schema.sql` is applied. A database is
compatible with this package's contract when its latest `calib_schema`
row is at least this version.

## Usage

``` r
calib_schema_version()
```

## Value

Character scalar, e.g. `"0.7.0"`.

## Examples

``` r
calib_schema_version()
#> [1] "0.7.0"
```
