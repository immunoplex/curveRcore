# Select the Best Eligible Model

Given per-model eligibility assessments and a named vector of ranking
scores, identifies the best eligible model. If no model passes all
gates, falls back to the model with the widest dynamic range among
converged models and records `fallback = TRUE` in the selection.

## Usage

``` r
select_best_eligible(
  assessments,
  ranking_scores,
  criterion = "AIC+eligibility",
  higher_is_better = FALSE
)
```

## Arguments

- assessments:

  Named list of
  [`assess_model_eligibility()`](https://immunoplex.github.io/curveRcore/reference/assess_model_eligibility.md)
  outputs, one entry per model name.

- ranking_scores:

  Named numeric vector. For AIC: the AIC value (lower is better, so
  provide as-is – this function finds the minimum). For LOO: the elpd
  value (higher is better – provide with `higher_is_better = TRUE`).

- criterion:

  Character. Label stored in the returned selection object, e.g.
  `"AIC+eligibility"` or `"LOO+eligibility"`.

- higher_is_better:

  Logical. If TRUE, higher `ranking_scores` are preferred (use for LOO
  elpd). Default FALSE (use for AIC).

## Value

A named list with:

- best_model_name:

  Character, or NA if no model converged.

- criterion:

  Character label.

- fallback:

  Logical. TRUE if no eligible model existed.

- fallback_reason:

  Character. Describes why fallback occurred.

- assessments:

  The full list of
  [`assess_model_eligibility()`](https://immunoplex.github.io/curveRcore/reference/assess_model_eligibility.md)
  outputs.

- eligible_models:

  Character vector of eligible model names.

## Details

This function is the shared selection layer used by both curveRfreq
(ranking by AIC) and curveRbayes (ranking by LOO elpd). The
`aic_selection` or `loo_selection` objects computed by the
framework-specific functions are passed through unchanged and augmented
with eligibility information.
