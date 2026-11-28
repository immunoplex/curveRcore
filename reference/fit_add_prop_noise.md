# Fit the Additive + Proportional Noise Parameters to Residuals

Maximum-likelihood estimate of `sigma_add` and `cv` from residuals and
fitted means of one fit group (Gaussian on the fitting scale). Used by
the frequentist engine to build weights.

## Usage

``` r
fit_add_prop_noise(resid, mu, is_log_response = TRUE, start = NULL)
```

## Arguments

- resid:

  Numeric. Residuals on the fitting scale.

- mu:

  Numeric. Fitted means on the fitting scale.

- is_log_response:

  Logical.

- start:

  Optional named numeric `c(sigma_add, cv)`.

## Value

Named list `sigma_add`, `cv`, `loglik`, `converged`.
