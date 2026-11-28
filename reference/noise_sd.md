# Observation-Noise SD on the Fitting Scale

Observation-Noise SD on the Fitting Scale

## Usage

``` r
noise_sd(
  mu,
  noise_model = "constant",
  sigma = NA_real_,
  sigma_add = NA_real_,
  cv = NA_real_,
  is_log_response = TRUE
)
```

## Arguments

- mu:

  Numeric. Mean response on the fitting scale (log10 when
  `is_log_response`).

- noise_model:

  `"constant"` or `"add_prop"`.

- sigma:

  Numeric scalar. Constant-model SD (fitting scale).

- sigma_add, cv:

  Numeric scalars. Additive SD (raw response units) and coefficient of
  variation (fraction) of the `add_prop` model.

- is_log_response:

  Logical.

## Value

Numeric vector, one SD per `mu`, on the fitting scale.
