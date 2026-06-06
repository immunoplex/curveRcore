# Build NLS Formulas for Candidate Models

Constructs [`stats::nls()`](https://rdrr.io/r/stats/nls.html)-compatible
formula objects for the requested models, using the curveRcore
parameterisation (b \> 0, always increasing).

## Usage

``` r
build_nls_formulas(
  model_names,
  response_variable,
  independent_variable = "concentration",
  fixed_a = NULL,
  is_log_response = TRUE
)
```

## Arguments

- model_names:

  Character vector. Which models to build formulas for. Must be a subset
  of
  [`available_models()`](https://immunoplex.github.io/curveRcore/reference/available_models.md).

- response_variable:

  Character. Name of the response column in the fitting data frame.

- independent_variable:

  Character. Name of the concentration column in the fitting data frame.
  Default `"concentration"`.

- fixed_a:

  Numeric scalar or NULL. If non-NULL, the lower asymptote is fixed at
  this value (on the **raw** scale before any log transform).

- is_log_response:

  Logical. If TRUE and `fixed_a` is non-NULL, log10-transform the fixed
  value before substitution.

## Value

Named list of formula objects, one per model.

## Details

When `fixed_a` is non-NULL and `is_log_response = TRUE`, the raw value
is validated and log10-transformed before substitution. When `fixed_a`
is non-NULL and `is_log_response = FALSE`, it is substituted directly.
