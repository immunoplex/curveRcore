# Apply a Blank Operation to Standard Curve Data

Performs one of five blank-handling strategies:

- `"ignored"` — no adjustment (default)

- `"included"` — append blank geometric mean as an extra point

- `"subtracted"` — subtract geometric mean of blanks

- `"subtracted_3x"` — subtract 3× geometric mean

- `"subtracted_10x"` — subtract 10× geometric mean

## Usage

``` r
perform_blank_operation(
  blank_data,
  data,
  response_variable,
  independent_variable,
  is_log_response,
  blank_option = "ignored",
  include_col = "included",
  blank_mean = NULL,
  verbose = FALSE
)
```

## Arguments

- blank_data:

  Data frame of blank measurements, or NULL.

- data:

  Data frame of standard curve data.

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

- is_log_response:

  Logical. Whether the response has been log10-transformed.

- blank_option:

  Character. One of the five options above.

- include_col:

  Character. Logical column marking fitted rows on `blank_data` (and,
  for the `"included"` option, on `data`). Absent = all rows included
  (backward compatible).

- blank_mean:

  Numeric or NULL. Pre-computed geometric mean of the *included* blanks.
  When NULL it is computed here from
  `blank_data[included, response_variable]`. Supplying it lets a caller
  guarantee the value used matches the one it records in
  `derived_stats`.

- verbose:

  Logical.

## Value

`data` with the blank operation applied. Blanks are *not* modified here;
[`preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)
transforms and returns them separately so they share the standards'
response floor.

## Details

The blank value is the geometric mean of the included blanks when all
are0, the arithmetic mean when any is \<= 0. Values that become \<= 0
aftersubtraction are kept (no floor); on the log path the offset
fromcompute_response_adjustment() handles them. Subtraction is opt-in
only.
