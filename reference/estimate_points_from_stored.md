# Back-calculate points from a stored fit

Rebuilds, for a set of points (control wells or standard points) read on
one stored curve, every estimate column of the 0.7.0 contract
([`point_estimates_contract()`](https://immunoplex.github.io/curveRcore/reference/point_estimates_contract.md))
from the stored best-model parameters, precision profile, limits and,
when present, parameter draws. Nothing is refitted. Intended for the
one-off backfill of curves fitted before 0.7.0; see
`calib_points_SPEC.md` for how close each column is to the engine's own
value.

## Usage

``` r
estimate_points_from_stored(
  response,
  dilution = 1,
  model_name,
  params,
  method = c("frequentist", "bayesian"),
  grid = NULL,
  limits = list(),
  draws = NULL,
  fixed_a = NULL,
  noise = NULL,
  is_log_independent = TRUE,
  nominal = NULL,
  is_log_response = TRUE,
  cv_x_max = 150,
  pcov_threshold = 20,
  decision_cutoff = NULL,
  seed = 1L
)
```

## Arguments

- response:

  Numeric vector: each point's response on the fitting scale (for
  standards, `calib_standards.response_model`; for controls, the raw
  response transformed as the fit's `is_log_response` says).

- dilution:

  Numeric vector or scalar: each point's dilution factor.

- model_name:

  Best model: `logistic4`, `logistic5`, `loglogistic4`, `loglogistic5`
  or `gompertz4`.

- params:

  Named numeric vector of the best model's parameter estimates (`a`,
  `b`, `c`, `d`, and `g` for five-parameter models). `a` may be absent
  when the lower asymptote was fixed; give it as `fixed_a`.

- method:

  `"frequentist"` or `"bayesian"`.

- grid:

  Optional data frame: the stored precision profile (`calib_grid` rows
  of the curve and method) with `log10_concentration` (or, when
  `is_log_independent = FALSE`, `concentration`) and `se_concentration`.

- limits:

  Named list from the curve's `calib_diagnostics` row: `lloq_log10`,
  `uloq_log10`, `lower_lod_response`, `upper_lod_response`,
  `lower_lod_log10_conc`, `upper_lod_log10_conc`, `inflect_x`. Missing
  entries count as unknown.

- draws:

  Optional named list (or data frame) of parameter draws, one element
  per term (`calib_draws` rows with `param_scope = 'curve'`). A length-1
  element is a fixed parameter and is recycled.

- fixed_a:

  Numeric scalar, the fixed lower asymptote, or `NULL`.

- noise:

  Optional list for Bayesian fits, per draw and aligned with `draws`:
  `sigma_obs`, `nu`, and for the additive + proportional model
  `sigma_add`, `cv` (with `noise_model = "add_prop"`; see
  [`noise_sd()`](https://immunoplex.github.io/curveRcore/reference/noise_sd.md)).
  When given, each draw's response is perturbed by the fitted
  observation noise before inversion, as the engine does with
  `include_measurement_error = TRUE`.

- is_log_independent:

  Logical: is x (concentration) on the log10 scale in the fit?

- nominal:

  Optional numeric vector: the points' nominal concentration on the
  fit's x scale (standards). When given, `recovery_pct` and a bias-aware
  `pcov_rmse` are added via
  [`standard_recovery()`](https://immunoplex.github.io/curveRcore/reference/standard_recovery.md).

- is_log_response:

  Logical: is the response on the log10 scale in the fit? Only used by
  the additive + proportional noise model (frequentist terms
  `sigma_add`, `cv`, or Bayesian `noise$noise_model = "add_prop"`).

- cv_x_max, pcov_threshold:

  As in the engines (defaults 150 and 20).

- decision_cutoff:

  Decision cutoff on the dilution-corrected natural scale, or `NULL`.

- seed:

  Integer seed for the measurement-noise draws (only used with `noise`).
  The caller's random-number state is left unchanged.

## Value

A data frame with one row per point and the estimate columns of
[`point_estimates_contract()`](https://immunoplex.github.io/curveRcore/reference/point_estimates_contract.md)
(`recovery_pct` only when `nominal` is given). Attribute `se_source`
records, per point, where se came from: `"draws"` (Bayesian draws, or
frequentist draws outside the profile), `"profile"` (the stored
precision profile) or `"none"`.

## See also

[`point_estimates_contract()`](https://immunoplex.github.io/curveRcore/reference/point_estimates_contract.md),
[`standard_recovery()`](https://immunoplex.github.io/curveRcore/reference/standard_recovery.md)
