# Assess Model Eligibility for Quantification

Applies a set of identifiability and precision gates to determine
whether a fitted model is reliable enough to be considered for selection
as the calibration model used for sample quantification. The same gates
are applied in both the frequentist and Bayesian frameworks, with
framework-specific gates automatically skipped when the required inputs
are not available.

## Usage

``` r
assess_model_eligibility(
  model_name,
  parameters,
  constraints = NULL,
  pcov_profile = NULL,
  grid_x = NULL,
  pcov_threshold = 20,
  bound_tol = 1e-04,
  max_rel_se = 5,
  min_dynamic_range_log10 = 0.5,
  vcov_matrix = NULL
)
```

## Arguments

- model_name:

  Character. Name of the model being assessed.

- parameters:

  Data frame. Must contain `term` and either `estimate`

  - `std_error` (frequentist) or `mean` + `sd` (Bayesian).

- constraints:

  Named list with `$lower` and `$upper` named numeric vectors for the
  free parameters, or NULL (Bayesian / fixed-a).

- pcov_profile:

  Numeric vector of pcov values (%) on the prediction grid, or NULL if
  the grid has not been computed for this model.

- grid_x:

  Numeric vector of log10_concentration values matching `pcov_profile`,
  or NULL.

- pcov_threshold:

  Numeric. The pcov (%) threshold used to define the LLOQ and ULOQ.
  Default 20.

- bound_tol:

  Numeric. A parameter estimate within this absolute distance of a
  constraint bound is treated as "at bound". Default 1e-4.

- max_rel_se:

  Numeric. Maximum permitted relative SE (`SE / |estimate|`) for any
  parameter. Default 5.0.

- min_dynamic_range_log10:

  Numeric. Minimum required LLOQ-to-ULOQ span in log10 concentration
  units. Default 0.5 (~3-fold).

- vcov_matrix:

  Numeric matrix. The parameter covariance matrix from `vcov(fit)`, or
  NULL. Used for the condition-number gate.

## Value

A named list:

- model_name:

  Character.

- eligible:

  Logical. TRUE if all applicable gates pass.

- gates:

  Data frame with columns `gate`, `passed`, `detail`.

- dynamic_range_log10:

  Numeric or NA.

- lloq:

  Numeric or NA (log10 scale).

- uloq:

  Numeric or NA (log10 scale).
