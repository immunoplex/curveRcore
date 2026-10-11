# Classify censoring for one curve's test samples

Classify censoring for one curve's test samples

## Usage

``` r
classify_censoring(
  response,
  x_log10,
  dilution,
  lloq_log10 = NA_real_,
  uloq_log10 = NA_real_,
  lod_response = NA_real_,
  ulod_response = NA_real_,
  lod_log10 = NA_real_,
  ulod_log10 = NA_real_,
  a = NA_real_,
  d = NA_real_
)
```

## Arguments

- response:

  Numeric vector: `observed_response_fit` (fitting scale).

- x_log10:

  Numeric vector: `predicted_log10_concentration` (pre-dilution).

- dilution:

  Numeric vector (or scalar): each sample's dilution factor.

- lloq_log10, uloq_log10:

  Numeric scalars (log10, pre-dilution) or `NA`.

- lod_response, ulod_response:

  Numeric scalars: `lower_lod_response` and `upper_lod_response` from
  [`compute_detection_limits()`](https://immunoplex.github.io/curveRcore/reference/compute_detection_limits.md),
  or `NA`.

- lod_log10, ulod_log10:

  Numeric scalars: `lower_lod_log10_conc` and `upper_lod_log10_conc`, or
  `NA`.

- a, d:

  Numeric scalars: the best model's asymptote estimates (fitting scale).
  Used for the curve's orientation and to place a non-invertible sample
  on the correct side.

## Value

A data.frame with one row per sample and columns `censor_class`,
`censor_type` (`none`, `left`, `interval`, `right`, or `NA`),
`censor_lower`, `censor_upper`.

## See also

[`classify_censoring_multiplate()`](https://immunoplex.github.io/curveRcore/reference/classify_censoring_multiplate.md),
[`summarize_conc_draws()`](https://immunoplex.github.io/curveRcore/reference/summarize_conc_draws.md)
