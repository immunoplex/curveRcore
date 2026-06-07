# Log10-Transform the Assay Response

Applies [`log10()`](https://rdrr.io/r/base/Log.html) to the response
column when `is_log_response` is TRUE. Non-positive values are floored
to an adaptive minimum before transform.

## Usage

``` r
compute_log_response(
  data,
  response_variable,
  is_log_response = TRUE,
  floor_method = "adaptive",
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame.

- response_variable:

  Character. Name of the response column.

- is_log_response:

  Logical. Apply log10? Default `TRUE`.

- floor_method:

  Character. How to handle non-positive values before log10:
  `"adaptive"` (default) uses 1\\ `"fixed"` uses 1e-6.

- verbose:

  Logical. Emit messages about floored values.

## Value

`data` with the response column optionally transformed.
