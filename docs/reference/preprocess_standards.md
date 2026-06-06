# Full Preprocessing Pipeline for Standard Curve Data

Applies concentration computation, prozone correction, blank handling,
and optional log10 response transform in the canonical order.

## Usage

``` r
preprocess_standards(
  data,
  antigen_settings,
  response_variable,
  independent_variable,
  is_log_response,
  blank_data = NULL,
  blank_option = "ignored",
  is_log_independent = TRUE,
  apply_prozone = TRUE,
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame of standards with a `dilution` column.

- antigen_settings:

  Named list with `standard_curve_concentration`.

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

- is_log_response:

  Logical. Log10-transform the response?

- blank_data:

  Data frame of blanks, or NULL.

- blank_option:

  Character. Blank handling method.

- is_log_independent:

  Logical. Log10-transform concentration?

- apply_prozone:

  Logical. Apply prozone correction?

- verbose:

  Logical.

## Value

A list with `data` (preprocessed) and `antigen_fit_options` (a record of
the options used).
