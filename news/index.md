# Changelog

## curveRcore 0.8.0 (development)

### Additive + proportional response noise model

- New
  [`noise_sd()`](https://immunoplex.github.io/curveRcore/reference/noise_sd.md),
  `NOISE_MODELS`,
  [`check_noise_model()`](https://immunoplex.github.io/curveRcore/reference/check_noise_model.md)
  and
  [`fit_add_prop_noise()`](https://immunoplex.github.io/curveRcore/reference/fit_add_prop_noise.md):
  a noise model shared by both engines. `"constant"` (default; existing
  fits unchanged) or `"add_prop"`, sigma^2 = sigma_add^2 + (cv \* mu)^2
  on the raw response, converted to the log10 scale by the delta method.
  Valid on either response scale; replaces the power-of-mean term
  `exp(log_sigma0 + slope * log|mu|)`, which was undefined on a log10
  response (mu crosses 0 where the response is 1) and was never enabled
  by the worker.
- [`new_study_params()`](https://immunoplex.github.io/curveRcore/reference/new_study_params.md)
  gains `noise_model` (default `"constant"`).
- [`estimate_points_from_stored()`](https://immunoplex.github.io/curveRcore/reference/estimate_points_from_stored.md)
  gains `is_log_response`; reads the frequentist terms `sigma_add`, `cv`
  (per-point response SD) and Bayesian `noise$noise_model = "add_prop"`
  with `sigma_add`, `cv` draws. The `use_hetero` / `log_sigma0` path is
  removed.

### One rule for values \<= 0 and blank subtraction

- New
  [`compute_response_adjustment()`](https://immunoplex.github.io/curveRcore/reference/compute_response_adjustment.md):
  one response adjustment per fit group, computed from the included
  calibrators (standard points), controls and test samples — never the
  blanks. On the log-response path, if every value is \> 0 nothing
  changes; if any is \<= 0, the offset `-min + delta` (default
  `delta = 0.001` response units) is added to every value so each well
  is kept and estimated. Replaces the floors that put standards at 1% of
  the smallest positive response and samples at 1e-6.
- New
  [`transform_sample_response()`](https://immunoplex.github.io/curveRcore/reference/transform_sample_response.md)
  and
  [`untransform_response()`](https://immunoplex.github.io/curveRcore/reference/untransform_response.md):
  the single transform (and inverse) for samples and controls, matching
  what
  [`preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)
  does to the standards.
- [`preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)
  gains `response_adjustment` (the group-wide adjustment; computed from
  the standards alone when NULL) and `offset_delta`. It returns
  `response_adjustment`, and `derived_stats` now records `blank_centre`,
  `blank_value`, `response_offset`, `offset_delta` and `n_nonpositive`.
  Data whose responses are all \> 0 give identical results.
- Blank subtraction (opt-in only; the default stays `"ignored"`): values
  that become \<= 0 are no longer floored at 1 (log) / 0 (linear); the
  offset handles them. The blank centre is the geometric mean when every
  included blank is \> 0 (unchanged) and the arithmetic mean when any
  blank is \<= 0 (the geometric mean silently dropped those blanks).

## curveRcore 0.7.1

### Data contract: options and diagnostics of the weights fit

- `calib_weights_fit` gains a nullable `fit_options` column (text
  holding a JSON object): the options the additive precision-weights fit
  actually used (`use_censoring`, `patient_col`, `plate_col`,
  `fix_beta1`, `fix_kappa`), so a stored weight set says how it was
  made. `NULL` for multiplicative rows and for rows written before
  0.7.1. Added to
  [`weights_contract()`](https://immunoplex.github.io/curveRcore/reference/weights_contract.md)
  and `R/calib_schema.sql`.
- `calib_weights_fit` gains six nullable diagnostic columns, filled for
  both weights models so they are judged on the same footing:
  `rhat_max`, `n_divergent` and `ess_min` (did the sampler converge),
  `ident_sd_ratio` (did the data rather than the prior set the key scale
  parameter: log kappa for the additive model, beta1 for the
  multiplicative one) and `resid_sd`, `resid_sd_ratio` (does the
  variance model describe the scaled residuals across the range of se).
  `NULL` for rows written before 0.7.1.
- [`calib_schema_version()`](https://immunoplex.github.io/curveRcore/reference/calib_schema_version.md)
  is now `"0.7.1"`; applying `calib_schema.sql` records it in
  `madi_results.schema_version`. In i-spi-compute this is migration 012.

## curveRcore 0.7.0

### Data contract: estimates and expected values for controls and calibrators

- New table `calib_controls` in `R/calib_schema.sql`: one row per
  control well, curve and method, with the same modelled estimates
  `calib_samples` holds for test samples (back-calculated concentration,
  pcov, precision gate, censoring, draw quantiles and probabilities),
  plus the control’s `source` (its name, e.g. `QC1 (Low)`), `included` /
  `exclusion_reason` / `mask_reason` (masked controls are estimated too)
  and the response columns.
- `calib_standards` gains the same estimate columns, plus `recovery_pct`
  (estimated / nominal, %), as nullable `ADD COLUMN IF NOT EXISTS`
  additions. For standards `pcov_rmse` includes the bias against the
  nominal value. The key and the existing columns are unchanged.
- New table `control_reference` (expected values of controls, with
  acceptance limits, SD, lot, basis and validity period; owned by I-SPI)
  and view `calib_controls_qc` (each control estimate with the value in
  force on its run date: recovery, z-score, within limits).
- New
  [`point_estimates_contract()`](https://immunoplex.github.io/curveRcore/reference/point_estimates_contract.md)
  and
  [`validate_point_estimates()`](https://immunoplex.github.io/curveRcore/reference/validate_point_estimates.md)
  for the three tables; a test keeps the DDL and the contract in sync.
- New
  [`standard_recovery()`](https://immunoplex.github.io/curveRcore/reference/standard_recovery.md)
  and
  [`add_standard_recovery_multiplate()`](https://immunoplex.github.io/curveRcore/reference/add_standard_recovery_multiplate.md).
- New
  [`estimate_points_from_stored()`](https://immunoplex.github.io/curveRcore/reference/estimate_points_from_stored.md):
  back-calculates controls or standards from a stored fit (parameters,
  precision profile, limits, draws) without refitting. Used once to fill
  the new columns for curves fitted before 0.7.0. Checked against
  curveRfreq and curveRbayes (see `calib_points_SPEC.md`).
- [`classify_pcov_gates_multiplate()`](https://immunoplex.github.io/curveRcore/reference/classify_pcov_gates_multiplate.md)
  and
  [`classify_censoring_multiplate()`](https://immunoplex.github.io/curveRcore/reference/classify_censoring_multiplate.md)
  now also process `$controls` and `$standard_estimates` on each curve’s
  result when present. Results without them are handled exactly as
  before.

### Release manifest

- New table `schema_version`; applying `calib_schema.sql` records
  `calib_schema` 0.7.0. New
  [`calib_schema_version()`](https://immunoplex.github.io/curveRcore/reference/calib_schema_version.md).
- `calib_run` gains `manifest` (jsonb): the composite release manifest a
  job ran under. New
  [`release_manifest()`](https://immunoplex.github.io/curveRcore/reference/release_manifest.md)
  (R, curveR packages, CmdStan, calib schema, renv lock hash, plus the
  caller’s component, version and commit) and
  [`manifest_json()`](https://immunoplex.github.io/curveRcore/reference/manifest_json.md).

### Other

- `calib_schema.sql` no longer names a personal database role in its
  comments.
- Design: `calib_points_SPEC.md`.

## curveRcore 0.6.0

### Data contract: precision-weight tables

- `R/calib_schema.sql` now defines `calib_weights` (per sample) and
  `calib_weights_fit` (per multiplate group), transcribed verbatim from
  the live DB. They were created by i-spi-compute migration 004 and were
  missing from the contract.
- Both tables are extended so one table holds the weights from either
  model, with keys unchanged (same sample hierarchy as `calib_samples`):
  a `model` column (default `multiplicative`, the model of every earlier
  row), `weights_pkg_version`, and the fields of the additive
  measurement-error model (per sample: `v_meas`, `tau2`, `reliability`
  with interval, `censor_type`, `bound`; per fit: `tau`, `kappa`,
  `reliability`, `fix_beta1`, `n_censored`, `identifiability`, with
  intervals). All new columns are nullable, so existing rows stay valid.
  The additions are idempotent `ADD COLUMN IF NOT EXISTS` statements.
  Design: `calib_weights_SPEC.md`.
- New
  [`weights_contract()`](https://immunoplex.github.io/curveRcore/reference/weights_contract.md)
  returns the column-level definition (type, key, which model a column
  applies to, description); new
  [`validate_weights_table()`](https://immunoplex.github.io/curveRcore/reference/validate_weights_table.md)
  checks rows against it before they are written. A test keeps the DDL
  and the contract in sync.

## curveRcore 0.5.2

### Documentation

- [`vignette("getting-started")`](https://immunoplex.github.io/curveRcore/articles/getting-started.md)
  gains a **Censoring** section (classes, bounds,
  [`classify_censoring()`](https://immunoplex.github.io/curveRcore/reference/classify_censoring.md),
  [`summarize_conc_draws()`](https://immunoplex.github.io/curveRcore/reference/summarize_conc_draws.md),
  and why censored samples should not be dropped or substituted) and a
  **Precision weights with curveRweights** section (the weight model,
  the workflow, and why weights and censoring are separate tools). The
  ecosystem overview now includes curveRweights.

## curveRcore 0.5.1

### Fixed analytic inverse-curve gradients (delta-method SEs)

- [`grad_logistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic5.md)
  and
  [`grad_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/grad_gompertz4.md):
  dx/da was set equal to (minus) dx/dy, missing the (d - a) dependence
  of the inverse;
  [`grad_logistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic5.md)
  also returned dx/dy with the wrong sign.
- [`grad_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic5.md):
  dx/da, dx/dd and dx/dg corrected;
  [`grad_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic4.md)
  and `grad_inv_loglogistic4_fixed()`: dx/dd sign corrected;
  `grad_inv_loglogistic5_fixed()`: sign of the 1/g term of dx/dg
  corrected.
- These feed the delta-method sample SEs in curveRfreq; NLS 5PL,
  Gompertz and log-logistic SEs (and CV-based reportability, coverage
  and any precision-profile gate) computed with earlier versions are
  affected. 4PL is unaffected.
- New tests check every analytic gradient, including `grad_y_*_fixed()`,
  against finite differences of the inverse functions at several
  parameter sets (including g \> 1). \# curveRcore 0.5.0

### Out-of-range test samples: censoring, asymmetric intervals, threshold probabilities

- New
  [`classify_censoring()`](https://immunoplex.github.io/curveRcore/reference/classify_censoring.md)
  /
  [`classify_censoring_multiplate()`](https://immunoplex.github.io/curveRcore/reference/classify_censoring_multiplate.md):
  give every back-calculated sample a `censor_class` (`quantified`,
  `lod_to_lloq`, `below_lod`, `below_lloq`, `above_uloq`, `saturated`,
  `no_response`), a `censor_type` (`none`/`left`/`interval`/`right`) and
  natural-scale, dilution-corrected `censor_lower`/`censor_upper`
  bounds. Run after
  [`classify_pcov_gates_multiplate()`](https://immunoplex.github.io/curveRcore/reference/classify_pcov_gates_multiplate.md).
  Classification uses the response scale first, so samples whose point
  estimate could not be inverted (response beyond an asymptote) are now
  placed on the correct side instead of `NA`.
- New
  [`summarize_conc_draws()`](https://immunoplex.github.io/curveRcore/reference/summarize_conc_draws.md):
  the shared per-sample summary both engines use for
  `conc_q_lo`/`conc_q_med`/`conc_q_hi` (2.5/50/97.5% predictive
  quantiles, natural scale, dilution-corrected), `p_below_lloq`,
  `p_above_uloq`, `p_below_cutoff` and `frac_draws_below_a`/`_above_d`.
- New
  [`invert_with_bounds()`](https://immunoplex.github.io/curveRcore/reference/invert_with_bounds.md)
  (internal-use export): keeps draws beyond an asymptote as
  `-Inf`/`+Inf` (censored) instead of dropping them.
- `se_concentration` and `pcov` are unchanged (frozen, §3.3). The new
  `frac_draws_*` columns report how many draws the legacy summaries
  drop.
- `test-pcov_gates.R` is now a real testthat suite (it previously ran as
  an empty test).

## curveRcore 0.4.4

### Documented the study-design-column passthrough contract on `tidy_samples()`

- Any column present on the `samples` frame handed to
  `fit_calibration_freq()`/ `_multiplate()` or `fit_calibration_bayes()`
  survives verbatim into `$samples` – neither
  [`new_calibration_result()`](https://immunoplex.github.io/curveRcore/reference/new_calibration_result.md)/`_multiplate()`
  nor the fitters’ sample-prediction step ever filter columns, they only
  add to the frame they were given. This is now an explicit, documented
  contract (a new `@section` on
  [`tidy_samples()`](https://immunoplex.github.io/curveRcore/reference/tidy_samples.md))
  rather than an implicit side effect, since it is the supported
  mechanism for threading study-design metadata (e.g. `timeperiod`,
  `agroup`) through to `curveRweights::as_weight_data(design = ...)`. No
  code change – documentation only.

## curveRcore 0.4.3

### Bug fix: multiplate tidy accessors mis-used `do.call()` for row-binding

- [`tidy_samples()`](https://immunoplex.github.io/curveRcore/reference/tidy_samples.md),
  [`tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.md),
  [`tidy_hyperparam()`](https://immunoplex.github.io/curveRcore/reference/tidy_hyperparam.md),
  and
  [`tidy_fit_diag()`](https://immunoplex.github.io/curveRcore/reference/tidy_fit_diag.md)
  on a `calibration_result_multiplate` all row-bound plates with
  `do.call(.cr_rbind_fill, parts)`. `.cr_rbind_fill()` takes a single
  `dfs` (list of data frames) argument, so
  [`do.call()`](https://rdrr.io/r/base/do.call.html) — which spreads
  `parts` as *separate* positional arguments — raised
  `"unused arguments"` for any multiplate result with 2+ non-NULL
  plates, and silently misbehaved (the lone data frame’s *columns* were
  iterated instead of a one-element list) for exactly one plate.
- Fixed all five call sites to call `.cr_rbind_fill(parts)` directly.
  Added `tests/testthat/test-tidy_extractors.R`, covering both the
  1-plate and 2+-plate cases for every affected accessor.
- No other curveRcore function routes through `.cr_rbind_fill()` or
  these accessors
  ([`compute_detection_limits_multiplate()`](https://immunoplex.github.io/curveRcore/reference/compute_detection_limits_multiplate.md)
  and
  [`classify_pcov_gates_multiplate()`](https://immunoplex.github.io/curveRcore/reference/classify_pcov_gates_multiplate.md)
  loop over `$plates` directly), so this only affects callers of
  [`tidy_samples()`](https://immunoplex.github.io/curveRcore/reference/tidy_samples.md)/[`tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.md)/[`tidy_hyperparam()`](https://immunoplex.github.io/curveRcore/reference/tidy_hyperparam.md)/
  [`tidy_fit_diag()`](https://immunoplex.github.io/curveRcore/reference/tidy_fit_diag.md)
  themselves — notably curveRweights’ `as_weight_data()`.

## curveRcore 0.4.2

- Created a new pcov_gate_class and changed the basis for pcov_pass
  classifcations.

## curveRcore 0.4.1

- new_antigen_constraints() default pcov_threshold is now 20 (was 15),
  matching the fitting-engine default and removing the mismatch that let
  sample-level pcov_pass diverge from the grid/LOQ threshold.

## curveRcore 0.4.0 (2026-07-29)

- **Fixed the inflection-point calculation.** `inflect_x` no longer
  snaps to the flat lower asymptote (~log10(grid_min_conc) = -4). Added
  [`compute_inflection()`](https://immunoplex.github.io/curveRcore/reference/compute_inflection.md),
  which returns the exact closed-form inflection from the fitted
  parameters (x = c for logistic4/gompertz4/loglogistic5; x = c +
  b\*ln(g) for logistic5), with y evaluated through the model’s own
  forward function.
- Hardened the grid-based inflection fallback to locate the steepest
  point (argmax \|dy/dx\|) instead of min \|d2y/dx2\|, making it immune
  to flat grid tails.

## curveRcore 0.3.0

### Mask-aware preprocessing

- [`preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)
  and its helpers
  ([`correct_prozone()`](https://immunoplex.github.io/curveRcore/reference/correct_prozone.md),
  [`perform_blank_operation()`](https://immunoplex.github.io/curveRcore/reference/perform_blank_operation.md),
  [`include_blanks_conc()`](https://immunoplex.github.io/curveRcore/reference/include_blanks_conc.md),
  [`compute_log_response()`](https://immunoplex.github.io/curveRcore/reference/compute_log_response.md))
  are now **mask-aware**. A new `include_col` argument (default
  `"included"`, a logical column; `TRUE` = used in the fit, `FALSE` =
  masked) threads through the pipeline. Every set-level statistic —
  prozone peak, blank geometric mean, adaptive log floor, and
  minimum-concentration anchor — is computed from the **included** rows
  only, while the resulting transforms are applied to **all** rows.
  Masked points therefore land on the same axes as the fitted points
  without ever influencing them.
- **Backward compatible.** When `include_col` is absent from a frame,
  every row is treated as included and the numeric output is identical
  to 0.2.0. This is covered by a regression guard in
  `tests/testthat/test-preprocess-mask.R`.
- [`correct_prozone()`](https://immunoplex.github.io/curveRcore/reference/correct_prozone.md)
  no longer drops rows with missing response/concentration; grain is
  preserved so masked/blank rows survive the pipeline. The peak
  reference is exposed via `attr(., "prozone_peak_response")` /
  `"prozone_logc_at_peak"`.
- [`compute_log_response()`](https://immunoplex.github.io/curveRcore/reference/compute_log_response.md)
  gains a `floor_value` argument so a single adaptive floor (derived
  from included standards) can be shared across the standards and blanks
  frames, keeping both on the same response axis.
- **Blanks are never subtracted automatically.** Documented explicitly:
  with the default `blank_option = "ignored"` the returned `blanks` are
  transformed for display/persistence only and the standards are
  untouched. Subtraction of the included-blank geometric mean occurs
  solely under the explicit `"subtracted"`, `"subtracted_3x"`,
  `"subtracted_10x"` options (and `"included"` adds it as an extra fit
  point). All blank statistics use the *included* blanks only.

### New `preprocess_standards()` return shape

- The return is now a list of `data`, **`blanks`** (all blank rows,
  transformed onto the standards’ floor, with `included` and a pristine
  `assay_response_raw` column), `antigen_fit_options`, and
  **`derived_stats`** (`blank_geomean`, `prozone_peak_response`,
  `prozone_logc_at_peak`, `response_floor`,
  `min_included_concentration`). `$data` is unchanged for existing
  callers and additionally carries the `included` and
  `assay_response_raw` columns. This is the material the batch worker
  persists to `madi_results.calib_standards` / `calib_blanks`.

### Downstream engines unaffected

- curveRfreq, curveRbayes, and curveRweights require **no change**: the
  worker filters to the included subset (`pp$data[pp$data$included, ]`)
  before fitting, so the fitters never see a masked row and the fit is
  byte-identical to before.

## curveRcore 0.2.0

- Initial release.
- Five canonical forward models: `logistic4`, `logistic5`,
  `loglogistic4`, `loglogistic5`, `gompertz4`.
- Shared `calibration_result` S3 class.
- Eligibility gating via
  [`assess_model_eligibility()`](https://immunoplex.github.io/curveRcore/reference/assess_model_eligibility.md)
  and
  [`select_best_eligible()`](https://immunoplex.github.io/curveRcore/reference/select_best_eligible.md).
- Full detection/quantification limit suite (LOD, LLOQ, ULOQ).
