# Log10-Transform the Assay Response (mask-aware)

Applies [`log10()`](https://rdrr.io/r/base/Log.html) to the response
column when `is_log_response` is TRUE. Non-positive values are floored
to an adaptive minimum before transform.

## Usage

``` r
compute_log_response(
  data,
  response_variable,
  is_log_response = TRUE,
  include_col = "included",
  floor_value = NULL,
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

- include_col:

  Character. Name of the logical column marking rows that entered the
  fit (`TRUE`) versus masked rows (`FALSE`). If the column is absent,
  every row is treated as included (backward compatible).

- floor_value:

  Numeric or NULL. If supplied, this exact floor is used for
  non-positive values instead of deriving one. Used to share a single
  floor across the standards and blanks frames.

- floor_method:

  Character. How to derive the floor when `floor_value` is NULL:
  `"adaptive"` (default) uses 1\\ positive value; `"fixed"` uses 1e-6.

- verbose:

  Logical. Emit messages about floored values.

## Value

`data` with the response column optionally transformed. The floor
actually used is attached as `attr(., "response_floor")`.

## Details

The adaptive floor is a **set-level statistic**: it is derived from the
*included* positive responses only (`data[included & value > 0]`), but
the floor and the [`log10()`](https://rdrr.io/r/base/Log.html) are then
applied to **all** rows. This is what lets masked points land on the
same response axis as the fitted points without influencing the
transform (see
[`preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)).
Pass an explicit `floor_value` to reuse a floor computed elsewhere (e.g.
to transform blanks on the same axis as the standards).
