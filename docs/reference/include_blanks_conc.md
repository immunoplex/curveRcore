# Include Blank Controls as an Extra Standard Curve Point

Appends a synthetic row whose response is the geometric mean of the
*included* blanks and whose concentration is `log10(2)` below the
minimum *included* standard concentration. The appended row is itself a
fit point, so its `included` flag is set to `TRUE`.

## Usage

``` r
include_blanks_conc(
  blank_data,
  data,
  response_variable,
  independent_variable = "concentration",
  include_col = "included"
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

- include_col:

  Character. Logical column marking fitted rows on both frames. Absent =
  all rows included.

## Value

`data` with one additional row.
