```markdown
# curveR Ecosystem Architecture Summary

## As of 2026-06-10

---

## Three-Package Architecture

### curveRcore (shared foundation)

-   **models.R**: Five forward models (logistic4, logistic5, loglogistic4,
    loglogistic5, gompertz4)
-   **inverses.R**: Analytical inverse functions for all five models
-   **derivatives.R**: dy/dx and d²y/dx² for all five models
-   **gradients.R**: Analytical gradients for delta-method error propagation +
    `make_inv_and_grad_fixed()` dispatcher.
    **API:** `make_inv_and_grad_fixed(model, fixed_a = NULL)` returns a list of
    three closures `$inv(y, p)`, `$grad(y, p)`, `$grad_y(y, p)`. `y` is passed
    to the closures at call time, not to the factory. This allows the factory to
    be constructed once per model and called inside a tight loop over grid points
    or samples without rebuilding closures per point.
-   **formulas.R**: `build_nls_formulas()` — generates NLS formula objects from
    model names (used by curveRfreq)
-   **transforms.R**: `preprocess_standards()`, `compute_concentration()`,
    `correct_prozone()`, `compute_log_response()`,
    `resolve_fixed_lower_asymptote()`, `adaptive_constraint_profile()`
-   **grid.R**: `generate_prediction_grid()` (flat args + S3 backwards compat),
    `predict_grid_response()`, `compute_curve_ci()`
-   **eligibility.R** *(new)*: `assess_model_eligibility()`,
    `select_best_eligible()` — shared quantification-aware model selection logic
    called by both curveRfreq and curveRbayes. See **Model Selection Philosophy**
    section for full description.
-   **settings.R**: S3 convenience constructors (`new_antigen_constraints`,
    `new_study_params`, `new_fit_options`) — optional, not required by fitting
    packages
-   **result_class.R**: `new_calibration_result()` (requires curve_id, NOT
    antigen/plate), `new_calibration_result_multiplate()`, comparison utilities
    (`compare_calibrations`, `compare_parameters`, `compare_samples`,
    `agreement_metrics`)
-   **data_helpers.R**: `filter_by_curve_id()` — utility for upstream callers
-   **utils.R**: `%||%`, `available_models()`, `model_params()`, `geom_mean()`
-   **data.R**: `bead_assay_example` (MFI, 6 plates, 2 antigens) +
    `elisa_assay_example` (OD, 6 plates, 1 antigen)

---

### curveRfreq (frequentist NLS fitting)

-   Receives ALREADY-PREPROCESSED stacked data frames with curve_id +
    concentration + response columns
-   No blanks, no lookup table, no antigen/plate/feature metadata — curve_id
    only
-   **fit_calibration.R**: `fit_calibration_freq()` — single curve, flat args.
    Runs the full pipeline: fit ensemble → compute per-model precision grids →
    run eligibility gates → select best eligible model by AIC → build
    calibration_result
-   **fit_multiplate.R**: `fit_calibration_freq_multiplate()` — splits by
    curve_id, returns `calibration_result_multiplate`
-   **constraints.R**: `compute_model_constraints()` — data-adaptive bounds,
    takes `fixed_a` directly (not antigen_settings)
-   **fit_one_plate.R**: `fit_ensemble_nls()` — multi-start Levenberg-Marquardt
    with fallbacks; `summarise_ensemble()` — extracts AIC, BIC, RSS per model
-   **predict_grid.R**: `predict_grid_freq()` — delta method with
    `se_response = sigma_fit`, outputs pcov + pcov_rmse. Called once per
    converged model (for eligibility gating) and the result reused as the
    plate-level `$grid` for the selected best model
-   **predict_samples.R**: `predict_samples_freq()` — same delta method for
    samples, applied to the best-eligible model only
-   **select_aic.R**: `select_best_aic()` (pure AIC ranking, unchanged),
    `extract_best_parameters()`

---

### curveRbayes (Bayesian hierarchical fitting via Stan)

-   Receives ALREADY-PREPROCESSED stacked data frames (same as curveRfreq)
-   Fits ALL curve_ids simultaneously via hierarchical Stan models
-   Returns `calibration_result_multiplate` — one `calibration_result` per
    curve_id (matches curveRfreq structure)
-   **fit_calibration_bayes.R**: `fit_calibration_bayes()` — flat args matching
    curveRfreq. Runs the full pipeline: fit all model families → compute
    per-model CDAN grids → run eligibility gates → select best eligible model
    by LOO-CV → build multiplate result
-   **fit_bayes.R**: `fit_bayes_single()`, `compile_stan_model()`,
    `extract_curve_params()`
-   **priors.R**: `compute_dynamic_priors()` — data-adaptive, supports `fixed_a`
    soft constraint, called per model_family
-   **stan_data.R**: `build_stan_data()` — stacked df → Stan data list with
    curve_id_map
-   **predict_bayes.R**: `predict_grid_bayes()` — **CDAN implementation**
    (O'Malley 2008): injects Student-t noise from posterior sigma_obs/nu before
    inverse. `predict_samples_bayes()` — no noise injection (sample IS the noisy
    observation)
-   **summary_bayes.R**: `summary_table_bayes()` — handles both multiplate (new)
    and single result (legacy)
-   **loo_selection.R**: `compute_loo()`, `compare_models_loo()` — PSIS-LOO +
    stacking weights (pure LOO ranking, unchanged)
-   **inst/stan/**: 5 Stan models (hierarchical, non-centered, Student-t
    likelihood, data-adaptive priors)

---

## Model Selection Philosophy

### The problem with AIC and LOO alone

The purpose of a calibration curve is **quantification**: given an observed
instrument response, back-calculate the analyte concentration with a known
uncertainty. AIC and LOO-CV measure **forward-fit quality** — how well the
model reproduces the observed standard-point responses. They are completely
agnostic to whether the fitted model can reliably invert.

A 5-parameter model can win on AIC by shaving a tiny amount off the residuals
using an extra flexibility parameter, while simultaneously being catastrophic
for inversion because that parameter is near-unidentified, boundary-constrained,
or produces a covariance matrix so ill-conditioned that the delta method returns
inflated standard errors for every grid point. The result: the model is selected
as "best" but produces a precision profile (pcov) that is pinned at the ceiling
(`cv_x_max = 150%`) for the entire concentration range — zero usable dynamic
range.

**AIC and LOO never see the precision profile. They cannot detect this failure.**

### The four-step selection pipeline

Both curveRfreq and curveRbayes now implement a four-step selection pipeline
using shared logic from `curveRcore::eligibility.R`:

```
Step 1  Fit all models
        curveRfreq: fit_ensemble_nls() via nlsLM
        curveRbayes: fit_bayes_single() via Stan/HMC

Step 2  Compute per-model precision profiles
        curveRfreq: predict_grid_freq() called for every converged model
                    using the delta method (O'Connell et al. 1993)
        curveRbayes: predict_grid_bayes() called for every model
                     using CDAN (O'Malley 2008)
        Profiles are stored at ensemble[[model]]$grid within each plate

Step 3  Apply eligibility gates
        curveRcore::assess_model_eligibility() called per model
        A model must pass ALL applicable gates to be eligible:

          Gate            Condition to pass               Active in
          ──────────────  ──────────────────────────────  ──────────────
          at_bound        No parameter estimate within    Frequentist only
                          bound_tol of a constraint       (needs hard bounds)
                          bound

          vcov_condition  Covariance matrix condition     Frequentist only
                          number κ < 1e8                  (no vcov in Bayes)

          rel_se          SE/|estimate| < max_rel_se      Both
                          for all parameters              (uses σ or SD)

          dynamic_range   LLOQ-to-ULOQ span ≥             Both
                          min_dynamic_range_log10         (needs pcov profile)
                          at pcov_threshold (default 20%)

        The Bayesian path runs only rel_se and dynamic_range because the
        prior regularization already mitigates boundary pathologies, and
        there is no vcov matrix (the posterior handles uncertainty).

Step 4  Rank eligible models and select
        curveRcore::select_best_eligible():
          - Among eligible models: rank by AIC (lower) or LOO elpd (higher)
          - If no model is eligible: fall back to the model with the widest
            dynamic range, set selection$fallback = TRUE
          - Attach full audit trail to selection object
```

### Why this matters in practice

On a typical bead-based immunoassay dataset:

-   logistic4 and gompertz4 are usually identifiable: well-conditioned
    covariance matrix, parameters interior to bounds, meaningful dynamic range
-   logistic5 and loglogistic5 often converge with `g` pinned at its lower
    constraint bound (boundary estimate), SE(g) ≈ SE(g)/|g| >> max_rel_se,
    and pcov capped at 150% everywhere — despite having a lower AIC than the
    4-parameter models

Without gating, logistic5 wins AIC selection and every sample gets reported with
pcov = 150% (the cap). With gating, logistic5 is excluded and logistic4 is
selected — producing a meaningful precision profile and a usable dynamic range.

The Bayesian regularization (log-normal prior on g centred at g=1, i.e., the
4PL) partly mitigates this: the posterior for g is pulled away from the boundary
and the dynamic_range gate may still pass. This is why the Bayesian framework
may validly select logistic5 while the frequentist framework correctly rejects
it — they are applying the same gate logic to precision profiles estimated by
different (and differently regularized) methods.

### Tunable gate parameters

All gate thresholds are exposed as top-level arguments on both fitting functions
so that users can tighten or loosen them for specific assay contexts:

| Parameter                   | Default | Meaning                                          |
|-----------------------------|---------|--------------------------------------------------|
| `pcov_threshold`            | 20      | % CV defining LLOQ/ULOQ for dynamic range gate  |
| `min_dynamic_range_log10`   | 0.5     | Min LLOQ-to-ULOQ span (~3-fold) for eligibility |
| `max_rel_se`                | 5.0     | Max SE/\|estimate\| for any parameter            |
| `bound_tol`                 | 1e-4    | Abs. distance from bound counted as "at bound"   |

---

## Key Design Decisions

### Data contract

-   curveRfreq and curveRbayes receive identically preprocessed data from
    curveRcore
-   Stacked data frames with `curve_id`, `concentration`, response column — all
    on fitting scale
-   No blanks in fitting packages — blank handling is upstream in
    curveRcore::preprocess_standards()
-   No antigen/plate/feature metadata — everything identified by curve_id only

### fixed_a handling

-   `l_asy_constraint_method = "default"` → `fixed_a = NULL` → a is free
    (estimated from data)
-   Non-default methods → `fixed_a` is a numeric value computed upstream
-   Frequentist: fixed_a baked into NLS formula as constant
-   Bayesian: tight prior on a (soft constraint) centered on fixed_a

### gradients.R API

`make_inv_and_grad_fixed(model, fixed_a = NULL)` is a **factory**: it builds
and returns three closures for the specified model and fixed_a setting. The
closures accept `(y, p)` where `y` is the observed response and `p` is the
named parameter vector from `coef(fit)`.

```r
# Build once outside the loop
fns <- curveRcore::make_inv_and_grad_fixed("logistic4", fixed_a = NULL)

# Call inside the loop — no per-iteration closure construction
for (i in seq_len(n_grid)) {
  x_hat   <- fns$inv(y_i, theta)
  grad_th <- fns$grad(y_i, theta)   # ∂x/∂θ — used for var_par
  grad_y  <- fns$grad_y(y_i, theta) # ∂x/∂y — used for var_y (observation noise)
}
```

Parameters are extracted by name (`p[["b"]]`, `p[["g"]]`, etc.) inside every
closure, so the order of elements in `p` does not matter. The old API required
`y` as a positional argument to the factory itself, which forced closure
reconstruction at every grid point.

### pcov implementation

-   **Frequentist**: Delta method with `se_response = summary(fit)$sigma`
    (residual SE). Includes both parameter covariance (var_par) and observation
    noise (var_y) terms per O'Connell et al. (1993) eq. 3.3
-   **Bayesian**: CDAN approach (O'Malley 2008): for each posterior draw,
    generate noisy response via Student-t(ν, μ, σ_obs), then back-calculate
    concentration. Posterior SD = pcov
-   Both output `pcov` (CV-based) and `pcov_rmse` (relative RMSE including bias)
-   pcov correlation between methods: ~0.96

### Output structure

-   Both packages return `calibration_result_multiplate` with `$plates` named
    list
-   Each plate is a `calibration_result` with `$grid`, `$samples`, `$ensemble`,
    `$selection`, `$meta`
-   `meta` requires: method, package, curve_id, response_var, independent_var,
    is_log_response, is_log_independent
-   `meta` also carries: `pcov_threshold` (threshold used for LOQ gating)

### Ensemble entry structure

Each `ensemble[[model_name]]` entry within a plate now contains:

```
$model_name   Character
$converged    Logical
$parameters   Data frame — term, estimate/mean, std_error/sd
$fit_stats    List — aic/bic/rss (freq) or n_divergent/ebfmi (Bayes)
$raw_fit      The nlsLM fit object (freq) or fit_bayes_single() output (Bayes)
$grid         Data frame — per-model pcov profile (ALWAYS present for converged
              models; populated during eligibility assessment in Step 2)
$eligibility  List — output of assess_model_eligibility():
                $eligible             Logical
                $gates                Data frame: gate, passed, detail
                $dynamic_range_log10  Numeric
                $lloq                 Numeric (log10 scale)
                $uloq                 Numeric (log10 scale)
```

The plate-level `$grid` is always the best-eligible model's grid (backwards
compatible). The per-model `ensemble[[model]]$grid` enables direct model-form
comparison of precision profiles in the vignette.

### Selection object structure

`calibration_result$selection` now contains:

```
$best_model_name    Character — the eligibility-gated winner
$criterion          "AIC+eligibility" (freq) or "LOO+eligibility" (Bayes)
$fallback           Logical — TRUE if no eligible model existed
$fallback_reason    Character — describes which gates failed on which models
$eligible_models    Character vector — models that passed all gates
$assessments        Named list — assess_model_eligibility() output per model
$weights            Data frame — AIC weights table (freq); from select_best_aic()
$aic_best           Character — AIC winner before gating (freq only)
$loo_comparison     loo_compare() matrix (Bayes only)
$loo_weights        Stacking weights (Bayes only)
```

The `$weights` / `$loo_comparison` tables are always present and unchanged from
the pre-gating selection, so all existing code that reads `selection$weights`
for reporting continues to work without modification.

### Per-model precision profiles

-   For curveRfreq: `predict_grid_freq()` is called for every converged model
    during eligibility assessment (Step 2). The best-model result is reused as
    the plate-level `$grid` — no recomputation.
-   For curveRbayes: CDAN grids are always computed for all models when
    `length(model_names) > 1` (forced internally for eligibility gating).
    Best model uses `n_draws_predict` (default 500) for full resolution.
    Non-best models use `n_draws_ensemble` (default 260) for faster profiles.
-   `compute_all_grids` argument: meaningful only for single-model Bayes fits
    where the user wants the ensemble grid for diagnostic purposes. For
    multi-model fits, grids are always produced regardless of this argument.
-   Provenance recorded in `meta`: `compute_all_grids`, `n_draws_predict`,
    `n_draws_ensemble`, `pcov_threshold`

---

## Current State

### Tests passing

-   curveRcore: **136 tests** ✓
-   curveRfreq: **124 tests** ✓ (up from 96; eligibility gating + per-model grids)
-   curveRbayes: **253 tests** ✓ (up from 108; eligibility gating + per-model grids)

### What changed in this development cycle

1.  **gradients.R API refactored** — `make_inv_and_grad_fixed()` factory
    signature changed; closures now take `(y, p)` instead of `(p)`. All
    call sites in curveRfreq updated. 200× faster per-grid-point evaluation.

2.  **Eligibility gating implemented** — `curveRcore::eligibility.R` added
    with `assess_model_eligibility()` and `select_best_eligible()`. Both
    fitting packages now gate on quantification-relevant criteria before
    AIC/LOO ranking. Fixes the longstanding issue where 5-parameter models
    with boundary-constrained asymmetry parameters won AIC selection despite
    having zero usable dynamic range.

3.  **Per-model precision grids standardised** — every converged ensemble entry
    now carries `$grid` with a full pcov profile, enabling model-form comparison
    of dynamic range in the vignette and downstream analysis.

4.  **Selection object enriched** — `selection$criterion` updated to
    `"AIC+eligibility"` / `"LOO+eligibility"`; full audit trail of gate
    outcomes available in `selection$assessments`.

5.  **Reference fixtures regenerated** — both `reference_freq_curve1.rds` and
    `reference_bayes_curve123.rds` regenerated. The frequentist fixture now
    stores `logistic4` as best model (previously `logistic5` won AIC despite
    being boundary-constrained).

### Known issues / next steps

1.  **Vignette** needs updating to reflect the new model selection logic,
    use `bayes_result$plates[["1"]]` access pattern throughout, and demonstrate
    per-model pcov comparison using `ensemble[[model]]$grid` for both methods
2.  **pcov_rmse** on frequentist grid may have edge cases with NA values at
    concentration extremes (minor; does not affect selection or quantification)
3.  **Architecture summary** should be updated after vignette changes are
    complete to reflect the final narrative structure

### Key files and their locations

-   curveRcore: `C:/Users/d78039e/Documents/R-git/curveRcore/`
-   curveRfreq: `C:/Users/d78039e/Documents/R-git/curveRfreq/`
-   curveRbayes: `C:/Users/d78039e/Documents/R-git/curveRbayes/`

---

## References

-   O'Malley, A.J. (2008). A Bayesian precision profile for measuring the
    quality of immunoassay experiments. *Phil. Trans. R. Soc. A*
    366:2301–2312. — CDAN precision profiles
-   O'Connell, M.A., Belanger, B.A. & Harland, P.D. (1993). Calibration
    and assay development using the four-parameter logistic model.
    *Chemometr. Intell. Lab. Syst.* 20:97–114. — Delta method for 4PL
-   Ekins, R.P. (1983). The precision profile: its use in assay design,
    assessment and quality control. In *Immunoassays for Clinical Chemistry*,
    pp. 76–105. — Precision profile framework
```
