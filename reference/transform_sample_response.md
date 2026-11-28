# Put Raw Responses on the Fitting Scale

The single transform for test samples and controls, matching what
[`preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)
does to the standards: subtract the blank value (if any), add the
offset, then log10 on the log path. Values still \<= 0 after the
adjustment are set to `display_floor`.

## Usage

``` r
transform_sample_response(
  raw,
  is_log_response,
  adjustment = NULL,
  display_floor = NULL
)
```

## Arguments

- raw:

  Numeric. Raw responses.

- is_log_response:

  Logical.

- adjustment:

  List from
  [`compute_response_adjustment()`](https://immunoplex.github.io/curveRcore/reference/compute_response_adjustment.md)
  (or the `derived_stats` of
  [`preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)),
  or NULL. NULL keeps the legacy behaviour (no subtraction, no offset,
  floor 1e-6) for callers that have not been given an adjustment.

- display_floor:

  Numeric or NULL. Floor for values still \<= 0; defaults to
  `adjustment$response_floor`, else 1e-6.

## Value

Numeric vector on the fitting scale.
