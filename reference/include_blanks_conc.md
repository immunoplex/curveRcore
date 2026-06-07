# Include Blank Controls as an Extra Standard Curve Point

Appends a synthetic row whose response is the geometric mean of the
blanks and whose concentration is half the minimum standard.

## Usage

``` r
include_blanks_conc(
  blank_data,
  data,
  response_variable,
  independent_variable = "concentration"
)
```

## Arguments

- blank_data:

  Data frame of blank wells.

- data:

  Data frame of standards.

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

## Value

`data` with one additional row.
