-- ============================================================================
-- calib_settings_seed_new_params.sql
-- Register the two new settings-cascade parameters (persist_draws,
-- bayes_single_plate) as SYSTEM-WIDE defaults in madi_results.calib_settings,
-- and describe them in madi_results.calib_settings_meta.
--
-- Both default FALSE. persist_draws gates calib_draws; bayes_single_plate flips
-- Bayesian fitting to per-plate (see the schema-expansion spec). Idempotent.
-- ============================================================================

-- ----------------------------------------------------------------------------
-- STEP 0 — match your conventions first. I do NOT know your exact param_group /
-- param_control_type vocabulary, so copy an existing BOOLEAN param's rows as a
-- template and mirror its group/control_type in the INSERTs below.
--   SELECT * FROM madi_results.calib_settings
--     WHERE param_name IN ('apply_prozone','include_measurement_error','is_log_response')
--       AND project_id = -1;
--   SELECT * FROM madi_results.calib_settings_meta
--     WHERE param_name IN ('apply_prozone','include_measurement_error','is_log_response');
-- Then adjust param_group / param_control_type / param_label wording to match.
-- ----------------------------------------------------------------------------


-- ----------------------------------------------------------------------------
-- STEP 1 — system-wide default rows (project_id = -1, all scope '__none__',
-- which satisfies calib_settings_system_wildcard_chk). calib_settings_id comes
-- from the sequence default; do not supply it.
-- ----------------------------------------------------------------------------
INSERT INTO madi_results.calib_settings
    (param_name, param_data_type, param_boolean_value, param_group)
VALUES
    ('persist_draws',      'boolean', FALSE, 'calibration'),   -- confirmed 2026-08
    ('bayes_single_plate', 'boolean', FALSE, 'calibration')    -- confirmed group
ON CONFLICT (project_id, study_accession, experiment_accession,
             feature, antigen, param_name) DO NOTHING;


-- ----------------------------------------------------------------------------
-- STEP 2 — metadata rows (drive UI/validation). param_control_type = 'switchInput'
-- (confirmed against the live DB's other boolean params).
-- ----------------------------------------------------------------------------
INSERT INTO madi_results.calib_settings_meta
    (param_name, param_group, param_label, param_data_type,
     param_control_type, param_description)
VALUES
    ('persist_draws', 'calibration', 'Persist posterior/sampling draws', 'boolean',
     'switchInput',
     'If TRUE, the fitting engine writes the full posterior (Bayesian) or '
     || 'asymptotic-MVN (frequentist) parameter draws to calib_draws. Heavy '
     || 'output; default FALSE. Does not gate calib_hyperparam / calib_fit_diag.'),
    ('bayes_single_plate', 'calibration', 'Bayesian single-plate (unpooled) fit',
     'boolean', 'switchInput',
     'Bayesian only. If TRUE, each plate is fit independently (N_plates = 1, no '
     || 'cross-plate pooling) instead of one multiplate hierarchy. Default FALSE '
     || '(multiplate pooling). The frequentist engine is always per-plate and '
     || 'ignores this flag.')
ON CONFLICT (param_name) DO NOTHING;


-- ----------------------------------------------------------------------------
-- STEP 3 — verify the resolver surfaces them (should return the two rows).
-- ----------------------------------------------------------------------------
-- SELECT param_name, param_data_type, param_boolean_value, param_group
-- FROM madi_results.calib_settings
-- WHERE param_name IN ('persist_draws','bayes_single_plate') AND project_id = -1;
