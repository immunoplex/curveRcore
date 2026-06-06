<div style="padding-top: 1.5rem;">

# curveRcore <img src="man/figures/logo_small.png" align="right" height="80" alt="curveRcore logo"/>

<!-- badges: start -->
[![License: AGPL 3.0](https://img.shields.io/badge/License-AGPL-yellow.svg)](https://opensource.org/license/agpl-3-0)
<!-- badges: end -->

**curveRcore** is the shared engine of the curveR suite. It provides the
forward models, preprocessing pipelines, analytical inverses and gradients,
eligibility gating, detection limits, and the `calibration_result` S3 class
that both `curveRfreq` and `curveRbayes` return. All fitting engines consume
and produce the same objects, so switching between frequentist and Bayesian
approaches — or comparing them directly — requires no data wrangling.

</div>

---

## The curveR ecosystem

```
curveRcore          ← you are here
├── curveRfreq      frequentist NLS (Levenberg-Marquardt, AIC selection)
├── curveRbayes     Bayesian (Stan, LOO-CV selection)
└── curveRweights   downstream variance / weighting models
```

Both fitting engines return a `calibration_result` object defined here,
making any downstream code — plots, extractors, weight models — completely
engine-agnostic.

---

## Installation

```r
# Install curveRcore first; the fitting packages depend on it
devtools::install_github("immunoplex/curveRcore")

# Then install whichever fitting engine(s) you need
devtools::install_github("immunoplex/curveRfreq")
devtools::install_github("immunoplex/curveRbayes")
```

---

## What curveRcore provides

### Five canonical forward models

All five models are parameterised consistently: `a` = lower asymptote,
`d` = upper asymptote, `b` = scale/slope (always > 0), `c` = inflection
point, `g` = asymmetry (5-parameter models only). All are monotonically
increasing and scale-agnostic — the `x` argument accepts whatever the
fitting engine passes in (typically log₁₀-concentration).

```r
library(curveRcore)

x <- seq(-2, 2, length.out = 200)   # log10-concentration

# 4PL — symmetric sigmoid
y4 <- logistic4(x, a = 100, b = 0.5, c = 0, d = 50000)

# 5PL — asymmetric sigmoid (g ≠ 1 skews toward one asymptote)
y5 <- logistic5(x, a = 100, b = 0.5, c = 0, d = 50000, g = 1.8)

# Gompertz — intrinsically asymmetric, no fifth parameter needed
yg <- gompertz4(x, a = 100, b = 2.0, c = 0, d = 50000)

# All five available:
available_models()
#> [1] "logistic4"    "loglogistic4" "gompertz4"    "logistic5"    "loglogistic5"
```

### Preprocessing pipeline

`preprocess_standards()` applies four steps in order before any model is
fitted. Both `curveRfreq` and `curveRbayes` call this function identically,
guaranteeing they receive data on the same scale.

```r
result <- preprocess_standards(
  data                 = standards_df,
  antigen_settings     = list(standard_curve_concentration = 30),
  response_variable    = "mfi",
  independent_variable = "concentration",
  is_log_response      = TRUE,
  blank_data           = blanks_df,
  blank_option         = "subtracted",   # "ignored" | "included" | "subtracted"
                                         # "subtracted_3x" | "subtracted_10x"
  is_log_independent   = TRUE,
  apply_prozone        = TRUE
)

preprocessed <- result$data
# $concentration now holds log10-concentration
# $mfi now holds log10-MFI
```

| Step | Function | What it does |
|---|---|---|
| 1 | `compute_concentration()` | Converts `dilution` column to absolute concentration; optionally log₁₀-transforms |
| 2 | `correct_prozone()` | Dampens the hook (prozone) effect at high concentrations |
| 3 | `perform_blank_operation()` | Applies the chosen blank strategy |
| 4 | `compute_log_response()` | log₁₀-transforms the response; floors non-positive values adaptively |

### Settings objects

Three constructor functions build the validated settings objects that
both fitting engines accept:

```r
# Per-antigen constraints
ac <- new_antigen_constraints(
  antigen        = "alpha",
  std_curve_conc = 10000,
  pcov_threshold = 15
)

# Study-level options
sp <- new_study_params(
  is_log_response    = TRUE,
  is_log_independent = TRUE,
  apply_prozone      = TRUE,
  blank_option       = "subtracted"
)

# Model and grid options
fo <- new_fit_options(
  model_names   = c("logistic4", "logistic5", "gompertz4"),
  n_grid        = 200L,
  cv_x_max      = 150,
  grid_min_conc = 1e-4
)
```

### The `calibration_result` object

Every fitting function in the suite returns a `calibration_result`. The
structure is identical regardless of which engine produced it:

```
calibration_result
├── $meta         provenance: method, package, curve_id, timestamps, scale flags
├── $ensemble     one entry per model attempted — parameters, fit_stats,
│                 precision grid, eligibility gates, raw fit object
├── $selection    best_model_name, criterion, AIC/LOO table, fallback status,
│                 and full per-model gate results
├── $grid         precision profile of the best model (200 points):
│                 predicted_response, se_concentration, pcov, d2y_dx2
├── $samples      back-calculated concentrations for test samples, or NULL
└── $diagnostics  inflection point, LLOQ, ULOQ, LODs, MDC, RDL
```

```r
# Print a compact summary
print(cr)
#> -- calibration_result (frequentist) --
#>   Curve ID : plate1_alpha
#>   Package  : curveRfreq v0.1.0
#>   Models   : 3 attempted, 3 converged
#>   Best     : gompertz4 (by AIC+eligibility)
#>   Grid     : 200 points
#>   Samples  : 20 predicted

# Full summary with parameters and diagnostics
summary(cr)
```

### Tidy extractors

```r
# Precision grid for the best model — one row per grid point
tidy_grid(cr)

# Sample back-calculations — one row per replicate
tidy_samples(cr)

# Both dispatch to multiplate objects automatically
tidy_grid(mp)      # row-binds all plates; adds curve_id column
tidy_samples(mp)
```

### Eligibility gating

Before AIC or LOO-CV ranking, every fitted model passes through four gates.
A model that fails any gate cannot be selected as the best model, preventing
an unidentified or boundary-constrained fit from winning purely on
forward-fit criteria.

| Gate | What it checks |
|---|---|
| `at_bound` | No parameter estimate is sitting on a constraint boundary |
| `vcov_condition` | Covariance matrix condition number < 10⁸ |
| `rel_se` | All parameter relative SEs (SE / \|estimate\|) below threshold |
| `dynamic_range` | LLOQ-to-ULOQ span at the pcov threshold is ≥ 0.5 log₁₀ units |

```r
# Inspect gate results for every model on one curve
cr$selection$assessments[["logistic4"]]$gates
#>            gate passed                    detail
#> 1      at_bound  FALSE     b at lower bound
#> 2 vcov_condition   TRUE
#> 3        rel_se   TRUE
#> 4 dynamic_range   TRUE     dynamic range = 2.14 log10
```

### Detection and quantification limits

`compute_detection_limits()` can be called automatically inside the fitting
pipeline or post-hoc on any existing result:

```r
cr <- compute_detection_limits(cr)

cr$detection_limits$lods$lower_lod_conc    # Lower LOD, natural scale
cr$detection_limits$lods$upper_lod_conc    # Upper LOD, natural scale
cr$detection_limits$mdc_rdl$rdl_lower_conc # Reliable detection limit
```

### Cross-engine comparison

```r
# Compare grid predictions from two engines
compare_calibrations(cr_freq, cr_bayes)

# Side-by-side parameter table
compare_parameters(cr_freq, cr_bayes)

# Paired sample predictions with agreement metrics
compare_samples(cr_freq, cr_bayes)

agreement_metrics(
  x = cr_freq$samples$predicted_log10_concentration,
  y = cr_bayes$samples$predicted_log10_concentration
)
#> $n     [1] 120
#> $bias  [1] -0.008
#> $mae   [1]  0.031
#> $rmse  [1]  0.044
#> $cor   [1]  0.998
#> $ccc   [1]  0.997
```

---

## Example datasets

Two synthetic datasets ship with curveRcore for testing and vignettes:

```r
data(bead_assay_example)   # Luminex bead assay: 2 antigens × 3 plates
data(elisa_assay_example)  # ELISA: 1 antigen × 6 plates, mixed curve shapes

str(bead_assay_example, max.level = 1)
#> List of 6
#>  $ standards      : 'data.frame': 60 obs. of 8 variables
#>  $ blanks         : 'data.frame': 24 obs. of 7 variables
#>  $ samples        : 'data.frame': 120 obs. of 13 variables
#>  $ curve_id_lookup: 'data.frame': 6 obs. of 5 variables
#>  $ response_var   : chr "mfi"
#>  $ indep_var      : chr "concentration"
```

---

## Getting help

- Full function reference: `https://immunoplex.github.io/curveRcore/reference/`
- Getting started vignette: `vignette("getting-started", package = "curveRcore")`
- Model family overview: `vignette("model-overview", package = "curveRcore")`
- File a bug: `https://github.com/immunoplex/curveRcore/issues`
