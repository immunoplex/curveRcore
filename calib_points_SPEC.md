# calib_controls / calib_standards / control_reference: estimates and expected values for controls and calibrators (0.7.0)

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
   with `calib_standards` and stay NULL. Expected values live in `control_reference`
   and are joined at read time by the `calib_controls_qc` view, because they are
   assigned and revised independently of fitting (see *Control reference values*).
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
   standard recovery, `control_reference` + `calib_controls_qc`, the backfill
   function, and the release-manifest pieces (`schema_version`,
   `calib_run.manifest`, `release_manifest()`).
2. curveRfreq / curveRbayes 0.6.0: `controls =` argument and an
   `estimate_standards` option (defaults leave output unchanged).
3. i-spi-compute migration 011 (same statements as the 0.7.0 DDL additions, plus the
   `control_for_fit` view) and worker 1.1.0. Applied before the new worker starts;
   additive only, no downtime.
4. One-off backfill (below) of every curve fitted before 0.7.0.
5. I-SPI 1.1.0: clone/delete include `calib_controls` (and the `control_reference`
   rows of a cloned study); Replicate Linearity and exports read only stored
   estimates. Because of the backfill there is no in-app fallback computation.

## Backfill of curves fitted before 0.7.0

`estimate_points_from_stored()` rebuilds every estimate column for a set of points
read on one stored curve, from what the database already holds; nothing is refitted
and existing fit rows are not touched. A one-off runner in i-spi-compute (with
migration 011) reads, per (curve, method): the best model's `calib_param` rows
(including `sigma_resid` for frequentist fits), `calib_grid`, `calib_diagnostics`,
curve-scope `calib_draws` (and, for Bayesian fits, the population `sigma_obs`, `nu`,
`log_sigma0`, `log_sigma_slope` draws and the noise model recorded in
`calib_run.params`). It then updates the `calib_standards` rows (response:
`response_model`; nominal: `log10_concentration`) and inserts `calib_controls` rows
(response: the raw control response, log-transformed as the fit was) with
`job_id = 'backfill-0.7.0'`.

How close the result is to what the engine would have stored (checked against
curveRfreq 0.5.1 and curveRbayes 0.5.1 on `bead_assay_example`, 20 samples):

| Column | Frequentist | Bayesian |
|---|---|---|
| point estimate | identical (inverse at the estimates) | median of the inverted draws: within the engine's own run-to-run Monte Carlo spread; inverse at the stored estimates when no draws are stored |
| se, pcov | stored profile: within 0.1%; outside the profile, from the draws: within 0.3% | SD of the inverted draws with the fitted noise: within the engine's run-to-run spread; the stored profile when no draws are stored |
| gate, censoring | identical (20/20) | identical (20/20) |
| quantiles, probabilities | same parameter draws with fresh noise draws: equal in distribution | same |

Without stored draws (fits before curveRfreq 0.5.0, or Bayesian fits without
`persist_draws`) the quantile and probability columns stay NULL, as they are in
`calib_samples` for those fits.

## Control reference values

`control_reference` (owned by I-SPI; the engines never read it): one row per
expected value and validity period of a control (`control_source` =
`calib_controls.source`), lot, antigen and feature, within a project and optionally
one study (`'__none__'` = project-wide / every feature). Values are at stock
(dilution 1), the scale of `calib_controls.final_concentration`: `expected_conc`,
`expected_sd`, `lower_limit`, `upper_limit`, `units`, `basis` (certificate,
assigned, historical), `valid_from` (default `-infinity`), `valid_to` (NULL = still
in force), `assigned_by`, `assigned_at`, `note`. A changed value is a new row with a
new `valid_from`, and the old row is closed with `valid_to`, so history is kept.

`calib_controls_qc` joins each control estimate to the value in force on the
plate's acquisition date (`xmap_header`), preferring a study-specific over a
project-wide value and a feature-specific over an every-feature value, and adds
`recovery_pct = 100 * final_concentration / expected_conc`,
`z_score = (final_concentration - expected_conc) / expected_sd` and
`within_limits`. Westgard rules need the sequence of runs and are applied by I-SPI
over this view.

Integration in I-SPI (after 1.1.0): an editor for values with a "suggest from
history" action (mean and SD of stored estimates over chosen runs, saved only when
a person accepts it), Levey-Jennings charts per control, antigen and standard
source within and across studies, normalization factors, and export columns.

## Release manifest

`madi_results.schema_version` records each schema component and version applied
(`calib_schema` = `calib_schema_version()`, written by `calib_schema.sql` itself).
`calib_run.manifest` (jsonb) holds, per job, the manifest the worker ran under:
`release_manifest()` (R, curveR package versions, CmdStan, calib schema version,
renv lock MD5) plus the worker's component name, image version and commit,
serialised with `manifest_json()`. NULL for jobs run before 0.7.0.
