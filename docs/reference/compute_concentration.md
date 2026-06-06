# Compute Concentration from Dilution and Undiluted Standard

Converts a `dilution` column to absolute concentration (optionally
log10-transformed) and stores it in `independent_variable`.

## Usage

``` r
compute_concentration(
  data,
  undiluted_sc_concentration,
  independent_variable,
  is_log_concentration = TRUE
)
```

## Arguments

- data:

  Data frame with a `dilution` column.

- undiluted_sc_concentration:

  Numeric. Concentration of the undiluted standard (e.g. 10000).

- independent_variable:

  Character. Column name for the output.

- is_log_concentration:

  Logical. Apply log10 after computing concentration? Default `TRUE`.

## Value

`data` with the concentration column populated.
