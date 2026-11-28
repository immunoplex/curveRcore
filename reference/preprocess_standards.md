# Full Preprocessing Pipeline for Standard Curve Data (mask-aware)

Applies concentration computation, prozone correction, blank handling,
and optional log10 response transform in the canonical order.

## Usage

``` r
preprocess_standards(
  data,
  antigen_settings,
  response_variable,
  independent_variable,
  is_log_response,
  blank_data = NULL,
  blank_option = "ignored",
  is_log_independent = TRUE,
  apply_prozone = TRUE,
  include_col = "included",
  response_adjustment = NULL,
  offset_delta = 0.001,
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame of standards with a `dilution` column.

- antigen_settings:

  Named list with `standard_curve_concentration`.

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

- is_log_response:

  Logical. Log10-transform the response?

- blank_data:

  Data frame of blanks, or NULL.

- blank_option:

  Character. Blank handling method.

- is_log_independent:

  Logical. Log10-transform concentration?

- apply_prozone:

  Logical. Apply prozone correction?

- include_col:

  Character. Name of the logical include/mask column on `data` and
  `blank_data`. Absent on a frame = all its rows are included.

- response_adjustment:

  List from
  [`compute_response_adjustment()`](https://immunoplex.github.io/curveRcore/reference/compute_response_adjustment.md)
  for the whole fit group (standards + controls + samples), or NULL to
  compute it from the standards alone.

- offset_delta:

  Numeric. `delta` used when `response_adjustment` is NULL.

- verbose:

  Logical.

## Value

A named list:

- `data`:

  All standard rows, transformed, carrying `include_col`, the (log10)
  `concentration`, the model-space response, and `assay_response_raw`.

- `blanks`:

  All blank rows, transformed onto the standards' response floor,
  carrying `include_col`, the model-space response, and
  `assay_response_raw`; `NULL` when no `blank_data` was supplied.

- `antigen_fit_options`:

  Record of the options used.

- `derived_stats`:

  The set-level statistics computed from the included points:
  `blank_geomean`, `prozone_peak_response`, `prozone_logc_at_peak`,
  `response_floor`, `min_included_concentration`, `blank_centre`,
  `blank_value`, `response_offset`, `offset_delta`, `n_nonpositive`.

- `response_adjustment`:

  The adjustment used (with `response_floor` added); pass it to
  [`transform_sample_response()`](https://immunoplex.github.io/curveRcore/reference/transform_sample_response.md)
  for samples and controls.

## Details

**Mask-aware contract.** An `include_col` logical column (default
`"included"`, `TRUE` = used in the fit, `FALSE` = masked) may be present
on `data` and `blank_data`. Every *set-level statistic* — the prozone
peak, the blank geometric mean, the adaptive log floor, and the
minimum-concentration anchor — is computed from the **included** rows
only. The resulting transforms are then applied to **all** rows, so
masked points land on the same axes as the fitted points without ever
influencing them. If the column is absent, every row is treated as
included and the output is identical to the pre-mask behaviour (backward
compatible).

The function does **not** drop masked rows; downstream fitters are
expected to receive only the included subset (e.g.
`pp$data[pp$data$included, ]`), which keeps the fit byte-identical to a
fit that never saw the masked rows.

Both frames retain a pristine `assay_response_raw` column (the response
before prozone/blank/log), so callers can persist the raw and
model-space responses side by side.

**Blanks are never subtracted automatically.** With the default
`blank_option = "ignored"` the standard responses are left untouched;
the returned `blanks` frame is transformed for display/persistence only
and is *not* subtracted from the standards. Subtraction happens **only**
when the caller explicitly selects `"subtracted"`, `"subtracted_3x"`, or
`"subtracted_10x"` (which subtract 1x/3x/10x the *included*-blank
geometric mean), or adds the blank mean as a point via `"included"`. The
returned blanks are always the raw and model-space blank responses,
never a subtracted quantity.

**Values \<= 0 (log response).** One rule, shared with the controls and
test samples via `response_adjustment` (see
[`compute_response_adjustment()`](https://immunoplex.github.io/curveRcore/reference/compute_response_adjustment.md)):
if every included calibrator, control and test-sample response is \> 0
nothing changes; otherwise the offset `-min + delta` is added to every
value before log10, so every well is kept. Blanks never decide the
offset. When `response_adjustment` is NULL it is computed from the
standards alone.
