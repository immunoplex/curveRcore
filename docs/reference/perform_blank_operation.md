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

- verbose:

  Logical.

## Value

`data` with the blank operation applied.

## Details

After subtraction, values that become non-positive are floored at 0
(linear) or 1 (log scale).
