# calib_weights / calib_weights_fit: integrated precision-weight storage (0.6.0)

**Goal.** Store precision weights from more than one weighting model in the existing
two tables, so users download one bundle with one set of weights, and add the
fields the additive measurement-error model needs. Nothing about the sample
hierarchy changes.

**Ownership.** The schema is defined here (`R/calib_schema.sql`, mirrored by
`weights_contract()`); curveRweights produces the values; the i-spi-compute worker
writes them; I-SPI reads them. A test keeps the DDL and `weights_contract()` in sync.

## Decisions
1. **Keys unchanged.** `calib_weights`: (curve_id, method, sampleid, patientid,
   timeperiod, dilution), the same sample hierarchy as `calib_samples`.
   `calib_weights_fit`: (multiplate_group_id, method).
2. **One weight set per (sample, method).** `model` is an ordinary column, not part
   of the key. Re-running with another model replaces the set (the workers already
   delete by key then insert). A side-by-side store can be added later by adding
   `model` to the keys without disturbing readers.
3. **`sigma`, `w`, `w_norm` mean the same in every model**: total SD, inverse-variance
   weight, mean-1 weight under the fitted model. Consumers of those columns are
   model-agnostic.
4. **New columns are nullable**, NULL = does not apply to the model that wrote the
   row, except `model` (NOT NULL DEFAULT 'multiplicative': what every pre-0.6.0 row
   used). Pre-0.6.0 rows therefore remain valid and unchanged.
5. **No CHECK on `model`** (family convention); `validate_weights_table()` enforces it
   before writing. `censor_type` carries the same CHECK as `calib_samples`.
6. **Not added:** a foreign key from `calib_weights` to `calib_samples`. It would
   guarantee identical samples but cascade-delete weights when a curve is
   recalibrated, a behaviour change; decide separately.

## Added columns (0.6.0)
`calib_weights`: model, weights_pkg_version, v_meas, tau2, reliability,
reliability_lo, reliability_hi, censor_type, bound.
`calib_weights_fit`: model, weights_pkg_version, tau, tau_lo, tau_hi, kappa, kappa_lo,
kappa_hi, reliability, reliability_lo, reliability_hi, fix_beta1, n_censored,
identifiability. Per-column meaning: `weights_contract()`.

## Row semantics
* **multiplicative**: sigma = phi * se^beta1. Fills phi, beta1 (+ intervals). tau,
  kappa, reliability, v_meas, tau2, censor_type, bound stay NULL.
* **additive**: sigma^2 = tau2 + v_meas, v_meas = kappa^2 * se^(2 beta1). Fills tau,
  kappa, reliability, fix_beta1, n_censored, identifiability; beta1 = 1 when fixed
  (CI NULL); phi, phi_lo, phi_hi NULL. Censored samples keep their row: censor_type
  and bound set, `w_norm` NULL (represented by the bound, not a weight).

## Applying
`R/calib_schema.sql` is idempotent. On the live DB the CREATEs skip and the
`ALTER TABLE ... ADD COLUMN IF NOT EXISTS` blocks add the columns; on a fresh DB the
result has the identical column order (verified on a scratch PostgreSQL 16: upgraded
vs fresh layouts identical; second apply a no-op). i-spi-compute migration 010 is
the same two ALTER statements.
