# Construct a Calibration Result Object

Creates a validated `calibration_result` S3 object. This is the
canonical output of both `fit_calibration_freq()` (curveRfreq) and
`fit_calibration_bayes()` (curveRbayes).

## Usage

``` r
new_calibration_result(
  meta,
  ensemble = list(),
  selection = list(),
  grid = data.frame(),
  samples = NULL,
  diagnostics = NULL,
  standards = NULL,
  blanks = NULL,
  population = NULL
)
```

## Arguments

- meta:

  Named list of metadata. Required fields: `method` (`"frequentist"` or
  `"bayesian"`), `package`, `version`, `feature`, `antigen`, `plate`,
  `response_var`, `independent_var`, `is_log_response`,
  `is_log_independent`.

- ensemble:

  Named list of model fit results, one per model attempted. Each element
  should be a list with `model_name`, `converged`, `parameters`,
  `fit_stats`, and optionally `raw_fit`.

- selection:

  Named list describing best-model selection: `best_model_name`,
  `criterion`, `weights` (data.frame).

- grid:

  Data frame of grid predictions from the best model.

- samples:

  Data frame of test sample predictions from the best model, or NULL if
  no samples were provided.

- diagnostics:

  Named list of diagnostic quantities (inflection point, LODs, LOQs,
  etc.), or NULL.

- standards:

  Data frame of standard data used to fit the curve, or NULL if not
  provided.

- blanks:

  Data frame of blank data used for QA, or NULL if not provided.

- population:

  Named list or NULL. Per-plate population/noise parameters for
  per-plate arms (single-plate Bayes: this plate's `sigma_obs`/`nu` and
  hyperprior-informed `mu_*`; frequentist: `sigma_resid`). Same shape as
  the multiplate `population` slot (see
  [`new_calibration_result_multiplate()`](https://immunoplex.github.io/curveRcore/reference/new_calibration_result_multiplate.md)).
  NULL for the pooled multiplate arm, where group params live on the
  multiplate container instead.

## Value

An object of class `calibration_result`.
