# curveRcore 0.3.0

## Mask-aware preprocessing

* `preprocess_standards()` and its helpers (`correct_prozone()`,
  `perform_blank_operation()`, `include_blanks_conc()`, `compute_log_response()`)
  are now **mask-aware**. A new `include_col` argument (default `"included"`,
  a logical column; `TRUE` = used in the fit, `FALSE` = masked) threads through
  the pipeline. Every set-level statistic — prozone peak, blank geometric mean,
  adaptive log floor, and minimum-concentration anchor — is computed from the
  **included** rows only, while the resulting transforms are applied to **all**
  rows. Masked points therefore land on the same axes as the fitted points
  without ever influencing them.
* **Backward compatible.** When `include_col` is absent from a frame, every row
  is treated as included and the numeric output is identical to 0.2.0. This is
  covered by a regression guard in `tests/testthat/test-preprocess-mask.R`.
* `correct_prozone()` no longer drops rows with missing response/concentration;
  grain is preserved so masked/blank rows survive the pipeline. The peak
  reference is exposed via `attr(., "prozone_peak_response")` /
  `"prozone_logc_at_peak"`.
* `compute_log_response()` gains a `floor_value` argument so a single adaptive
  floor (derived from included standards) can be shared across the standards and
  blanks frames, keeping both on the same response axis.
* **Blanks are never subtracted automatically.** Documented explicitly: with the
  default `blank_option = "ignored"` the returned `blanks` are transformed for
  display/persistence only and the standards are untouched. Subtraction of the
  included-blank geometric mean occurs solely under the explicit `"subtracted"`,
  `"subtracted_3x"`, `"subtracted_10x"` options (and `"included"` adds it as an
  extra fit point). All blank statistics use the *included* blanks only.

## New `preprocess_standards()` return shape

* The return is now a list of `data`, **`blanks`** (all blank rows, transformed
  onto the standards' floor, with `included` and a pristine `assay_response_raw`
  column), `antigen_fit_options`, and **`derived_stats`** (`blank_geomean`,
  `prozone_peak_response`, `prozone_logc_at_peak`, `response_floor`,
  `min_included_concentration`). `$data` is unchanged for existing callers and
  additionally carries the `included` and `assay_response_raw` columns. This is
  the material the batch worker persists to `madi_results.calib_standards` /
  `calib_blanks`.

## Downstream engines unaffected

* curveRfreq, curveRbayes, and curveRweights require **no change**: the worker
  filters to the included subset (`pp$data[pp$data$included, ]`) before fitting,
  so the fitters never see a masked row and the fit is byte-identical to before.

# curveRcore 0.2.0

* Initial release.
* Five canonical forward models: `logistic4`, `logistic5`, `loglogistic4`,
  `loglogistic5`, `gompertz4`.
* Shared `calibration_result` S3 class.
* Eligibility gating via `assess_model_eligibility()` and
  `select_best_eligible()`.
* Full detection/quantification limit suite (LOD, LLOQ, ULOQ).
