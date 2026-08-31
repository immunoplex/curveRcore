# Create Study-Level Fitting Parameters

Constructs the study-wide settings that control data preprocessing.

## Usage

``` r
new_study_params(
  is_log_response = TRUE,
  is_log_independent = TRUE,
  apply_prozone = TRUE,
  blank_option = "ignored",
  persist_draws = FALSE,
  bayes_single_plate = FALSE
)
```

## Arguments

- is_log_response:

  Logical. Log10-transform the assay response? Default TRUE.

- is_log_independent:

  Logical. Log10-transform the concentration? Default TRUE.

- apply_prozone:

  Logical. Apply prozone (hook effect) correction? Default TRUE.

- blank_option:

  Character. Blank handling method. One of `"ignored"`, `"included"`,
  `"subtracted"`, `"subtracted_3x"`, `"subtracted_10x"`. Default
  `"ignored"`.

- persist_draws:

  Logical. If TRUE, the fitting engine persists the full posterior
  (Bayesian) or asymptotic-MVN (frequentist) parameter draws to
  `calib_draws`. Heavy output; default FALSE. Gate only — the light
  `calib_hyperparam`/`calib_fit_diag` tables are always written.

- bayes_single_plate:

  Logical. Bayesian only. If TRUE, each plate is fit independently
  (N_plates = 1, no cross-plate pooling) instead of as one multiplate
  hierarchy. Default FALSE (multiplate pooling). The frequentist engine
  is always per-plate and ignores this flag.

## Value

A named list of class `study_params`.
