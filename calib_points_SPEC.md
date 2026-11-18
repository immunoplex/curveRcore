# calib_controls / calib_standards: modelled estimates for controls and calibrators (0.7.0)

**Goal.** Store, for every control well and every calibrator (standard-curve point),
the same modelled estimates `calib_samples` holds for test samples: back-calculated
concentration, precision (pcov), precision gate, censoring, and posterior or
predictive quantiles. This lets I-SPI and exports read them directly instead of
recomputing them, and supports replicate linearity, standard-source bridging,
batch and multisite normalization, and long-term control and calibrator tracking
within and across studies.

**Ownership.** The schema is defined here (`R/calib_schema.sql`, mirrored by
`point_estimates_contract()`); curveRfreq and curveRbayes produce the estimates;
curveRcore adds precision gates, censoring and standard recovery; the
i-spi-compute worker writes the rows; I-SPI reads them. A test keeps the DDL and
the contract in sync.

## Decisions

1. **Estimates never change a fit.** Both engines back-calculate points after the
   curve is fitted (`predict_samples_freq()`, `predict_samples_bayes()`). Controls
   and standards go through the same step, so fitted curves, parameters and
   existing `calib_samples` rows are unchanged.
2. **Controls get their own table**, `calib_controls`, one row per
   (curve_id, method, well). A separate table keeps controls out of everything that
   reads `calib_samples`, in particular the precision-weight fits.
3. **A control is estimated through every curve it can be read on**: the same
   matching as test samples (plate, nominal sample dilution, antigen, feature,
   wavelength; every standard source), via the worker's new `control_for_fit` view.
4. **Masked controls are estimated too**, with `included = FALSE`,
   `exclusion_reason = 'masked'` and the `mask_reason`. (Masked test samples are
   still dropped before estimation, unchanged.)
5. **`calib_standards` keeps its key and its existing columns**; the estimate
   columns are added as nullable columns (`ADD COLUMN IF NOT EXISTS`). Rows written
   before 0.7.0 keep NULL estimates until the curve is refitted.
6. **Standards know their true value**, so they carry two things samples cannot:
   * `recovery_pct` = 100 x estimated / nominal concentration (both at the
     dilution the point was run at, i.e. `10^(predicted_log10_concentration -
     log10_concentration)`).
   * a bias-aware `pcov_rmse`: the CV of the root-mean-square error,
     sqrt(se^2 + bias^2), where bias is the log10 difference between the estimate
     and the nominal value. For samples and controls bias is unknown and
     `pcov_rmse = pcov`, as before.
7. **A standard that entered the fit is checked against a curve it helped shape**,
   so its recovery is optimistic. `included` says which points were in the fit;
   masked standards are independent checks. Readers should show this.
8. **For a standard, `final_concentration` = estimate x dilution** is that point's
   estimate of the stock concentration. Across points, sources and plates it is the
   quantity to compare when bridging standard sources.
9. **Expected (nominal) values of controls are not stored with the estimates.**
   `calib_controls.concentration` / `log10_concentration` are present for symmetry
   with `calib_standards` and are NULL for now. Control recovery will come from a
   separate reference-value table joined at read time (see *Later* below), because
   reference values are assigned and revised independently of fitting.
10. **Same vocabularies as `calib_samples`**: `censor_class`, `censor_type`,
    `pcov_gate_class` carry the same CHECK constraints.

## Columns

`calib_controls` (new): curve_id, method, well, sampleid, source, dilution, included,
exclusion_reason, mask_reason, concentration, log10_concentration, response_model,
assay_response_raw, predicted_concentration, predicted_log10_concentration,
final_concentration, se_concentration, pcov, pcov_rmse, pcov_pass, pcov_gate_class,
censor_class, censor_type, censor_lower, censor_upper, conc_q_lo, conc_q_med,
conc_q_hi, p_below_lloq, p_above_uloq, p_below_cutoff, frac_draws_below_a,
frac_draws_above_d, job_id, created_at.
Key (curve_id, method, well); FK curve_id -> curve_lookup ON DELETE CASCADE.

`calib_standards` (added): predicted_concentration, predicted_log10_concentration,
final_concentration, se_concentration, pcov, pcov_rmse, pcov_pass, pcov_gate_class,
censor_class, censor_type, censor_lower, censor_upper, conc_q_lo, conc_q_med,
conc_q_hi, p_below_lloq, p_above_uloq, p_below_cutoff, frac_draws_below_a,
frac_draws_above_d, recovery_pct.

Per-column meaning: `point_estimates_contract()`.

## Result object (engines -> worker)

A fitted `calibration_result` may carry two new data frames beside `$samples`:

* `$controls`: the control rows handed to the engine, with the sample estimate
  columns added.
* `$standard_estimates`: the standard rows (included and masked), with the same
  estimate columns added. `concentration` is the x value the fit used (log10 when
  `is_log_independent`).

`classify_pcov_gates_multiplate()` and `classify_censoring_multiplate()` now
process `$samples`, `$controls` and `$standard_estimates` identically.
`add_standard_recovery_multiplate()` adds `recovery_pct` and the bias-aware
`pcov_rmse` to `$standard_estimates`. A result without the new slots is handled
exactly as before.

## Rollout

1. curveRcore 0.7.0: this contract, DDL, gating/censoring over all point types,
   standard recovery (plus the release-manifest pieces: `schema_version`, run
   provenance).
2. curveRfreq / curveRbayes 0.6.0: `controls =` argument and an
   `estimate_standards` option (defaults leave output unchanged).
3. i-spi-compute migration 011 (same statements as the 0.7.0 DDL additions, plus the
   `control_for_fit` view) and worker 1.1.0. Applied before the new worker starts;
   additive only, no downtime.
4. I-SPI 1.1.0: clone/delete include `calib_controls`; Replicate Linearity and
   exports read stored estimates, falling back to in-app computation for curves
   fitted before 0.7.0.

## Later: reference values for controls

Not part of 0.7.0. A control's expected concentration is assigned by the lab
(certificate, bridging run, or accumulated history), can differ by antigen, feature
and lot, and is revised when a new lot is introduced. Storing it per fitted row
would freeze it at fit time; instead a reference table is joined when reading:

```
madi_results.control_reference (
  reference_id      bigint identity primary key,
  project_id        integer not null,
  study_accession   varchar(15) not null default '__none__',   -- '__none__' = project-wide
  control_source    varchar(64) not null,                      -- calib_controls.source
  control_lot       varchar(64) not null default '__none__',
  antigen           varchar(64) not null,
  feature           varchar(15) not null default '__none__',
  units             text,
  expected_conc     numeric,          -- target, at stock (dilution 1)
  lower_limit       numeric,          -- acceptance range, same scale
  upper_limit       numeric,
  basis             text,             -- certificate | assigned | historical
  valid_from, valid_to date,
  assigned_by text, assigned_at timestamptz, note text,
  unique (project_id, study_accession, control_source, control_lot, antigen, feature, valid_from)
)
```

A view `calib_controls_qc` would join `calib_controls` (via curve_lookup for
project, study, antigen, feature, run date) to the reference in force, giving
`recovery_pct = 100 * final_concentration / expected_conc`, the z-score against
the acceptance range or history, and Westgard-rule flags. Owned by I-SPI (the
fitting engines never need it).
