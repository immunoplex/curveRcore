# Compute the Response Adjustment for One Fit Group

Determines the blank-subtraction value and the log-scale offset that
[`preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)
and
[`transform_sample_response()`](https://immunoplex.github.io/curveRcore/reference/transform_sample_response.md)
apply, so the standards, controls and test samples of a fit group share
one response axis.

## Usage

``` r
compute_response_adjustment(
  standards,
  response_variable,
  is_log_response,
  samples = NULL,
  controls = NULL,
  blank_data = NULL,
  blank_option = "ignored",
  include_col = "included",
  delta = 0.001
)
```

## Arguments

- standards:

  Data frame of standards (raw response).

- response_variable:

  Character. Response column name (same in all frames).

- is_log_response:

  Logical. Offset is only applied on the log path.

- samples, controls:

  Data frames or NULL. Test samples and controls (raw response). Rows
  flagged excluded by `include_col` are ignored.

- blank_data:

  Data frame or NULL. Blanks (only for the subtraction value).

- blank_option:

  Character. See
  [`perform_blank_operation()`](https://immunoplex.github.io/curveRcore/reference/perform_blank_operation.md).

- include_col:

  Character. Logical include/mask column; absent = all rows.

- delta:

  Numeric. Margin above zero for the smallest value. Default 0.001
  (response units).

## Value

Named list: `blank_value`, `blank_centre`, `offset`, `delta`,
`n_nonpositive` (values \<= 0 in the offset set), `min_response`,
`is_log_response`, `blank_option`.

## Details

**Offset rule (log response only).** If every included standard, control
and test-sample response (after any blank subtraction) is \> 0, the
offset is 0 and nothing changes. If any is \<= 0,
`offset = -min + delta` is added to every value, so the smallest lands
`delta` above zero. Blanks never decide the offset. The offset is
recorded with the fit.
