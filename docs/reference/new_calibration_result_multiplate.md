# Construct a Multi-Plate Calibration Result

Wraps multiple single-plate `calibration_result` objects into a
multi-plate container.

## Usage

``` r
new_calibration_result_multiplate(meta, plates, population = NULL)
```

## Arguments

- meta:

  Named list. Multi-plate metadata (must include `plates` character
  vector).

- plates:

  Named list of `calibration_result` objects, one per plate.

- population:

  Named list or NULL. Group-level (population/noise) parameters of the
  pooled Bayesian multiplate fit. NULL for per-plate fits (single-plate
  Bayes or frequentist), which carry no pooled level. Expected shape
  (all optional):

  params

  :   Data frame `term, estimate, std_error, q_lo, q_med, q_hi` — the
      group scalars (`mu_*`, `sigma_*`, `sigma_obs`, `nu`, ...).

  draws

  :   Named list `term -> numeric()` of iteration-ordered posterior
      draws (present only when `persist_draws = TRUE`); all terms share
      one draw order so they can be column-bound into a joint posterior.

  fit_diag

  :   Named list of per-fit sampler diagnostics (`rhat_max`,
      `ess_bulk_min`, `n_divergent`, `fit_seconds`, `fit_seed`, ...).

## Value

An object of class `calibration_result_multiplate`.
