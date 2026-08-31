-- ============================================================================
-- calib_schema.sql — canonical DDL for madi_results.calib_*
--
-- Existing tables: transcribed VERBATIM from the live DB (extract_calib_ddl.sql
-- output), reordered into FK-dependency order for fresh-apply safety.
-- New tables (schema expansion): calib_hyperparam, calib_draws, calib_fit_diag —
-- reconciled to the family's real conventions (method text; FK (curve_id,method,
-- model_name)->calib_fit ON DELETE CASCADE, exactly like calib_param/calib_loo;
-- NO method CHECK on children; job_id varchar(64) + created_at like calib_loo).
--
-- PREREQUISITE (external, not created here): madi_results.curve_lookup — the
-- injection/lookup table that calib_fit / calib_diagnostics / calib_grid /
-- calib_samples FK to.
--
-- Idempotent (IF NOT EXISTS). Safe to re-run against the live DB.
-- ============================================================================

CREATE SCHEMA IF NOT EXISTS madi_results;

-- Sequence backing calib_settings.calib_settings_id (authoritative definition
-- from the live DB). The CREATE goes here (the column DEFAULT needs it); the
-- ALTER ... OWNED BY is deferred to after calib_settings is created (it
-- references that column). OWNER TO is environment-specific — see the tail.
CREATE SEQUENCE IF NOT EXISTS madi_results.calib_settings_calib_settings_id_seq
    INCREMENT 1
    START 1
    MINVALUE 1
    MAXVALUE 9223372036854775807
    CACHE 1;


-- ============================================================================
-- Job + fit hub
-- ============================================================================

CREATE TABLE IF NOT EXISTS madi_results.calib_run (
    job_id character varying(64) NOT NULL,
    method text NOT NULL,
    package text,
    version text,
    best_model text,
    params jsonb,
    status text,
    started_at timestamp with time zone DEFAULT now(),
    finished_at timestamp with time zone,
    CONSTRAINT calib_run_pkey PRIMARY KEY (job_id),
    CONSTRAINT calib_run_method_check CHECK ((method = ANY (ARRAY['frequentist'::text, 'bayesian'::text])))
);

CREATE TABLE IF NOT EXISTS madi_results.calib_fit (
    curve_id bigint NOT NULL,
    method text NOT NULL,
    model_name text NOT NULL,
    package text,
    version text,
    converged boolean,
    eligible boolean,
    is_best boolean NOT NULL DEFAULT false,
    is_fallback boolean,
    criterion text,
    score_type text,
    selection_score numeric,
    selection_weight numeric,
    dynamic_range_log10 numeric,
    n_params integer,
    job_id character varying(64),
    created_at timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_fit_pkey PRIMARY KEY (curve_id, method, model_name),
    CONSTRAINT calib_fit_curve_id_fkey FOREIGN KEY (curve_id) REFERENCES madi_results.curve_lookup(curve_id) ON DELETE CASCADE,
    CONSTRAINT calib_fit_method_check CHECK ((method = ANY (ARRAY['frequentist'::text, 'bayesian'::text]))),
    CONSTRAINT calib_fit_score_type_check CHECK ((score_type = ANY (ARRAY['aic'::text, 'loo_elpd'::text])))
);


-- ============================================================================
-- Children FK -> calib_fit(curve_id, method, model_name)
-- ============================================================================

CREATE TABLE IF NOT EXISTS madi_results.calib_param (
    curve_id bigint NOT NULL,
    method text NOT NULL,
    model_name text NOT NULL,
    term text NOT NULL,
    estimate numeric,
    std_error numeric,
    q_lo numeric,
    q_med numeric,
    q_hi numeric,
    job_id text,
    CONSTRAINT calib_param_pkey PRIMARY KEY (curve_id, method, model_name, term),
    CONSTRAINT calib_param_curve_id_method_model_name_fkey FOREIGN KEY (curve_id, method, model_name) REFERENCES madi_results.calib_fit(curve_id, method, model_name) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS madi_results.calib_loo (
    curve_id bigint NOT NULL,
    method text NOT NULL DEFAULT 'bayesian'::text,
    model_name text NOT NULL,
    elpd_loo numeric,
    se_elpd_loo numeric,
    p_loo numeric,
    looic numeric,
    elpd_diff numeric,
    se_diff numeric,
    pareto_k_good integer,
    pareto_k_ok integer,
    pareto_k_bad integer,
    pareto_k_vbad integer,
    job_id character varying(64),
    created_at timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_loo_pkey PRIMARY KEY (curve_id, method, model_name),
    CONSTRAINT calib_loo_curve_id_method_model_name_fkey FOREIGN KEY (curve_id, method, model_name) REFERENCES madi_results.calib_fit(curve_id, method, model_name) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS madi_results.calib_gate (
    curve_id bigint NOT NULL,
    method text NOT NULL,
    model_name text NOT NULL,
    gate text NOT NULL,
    passed boolean,
    detail text,
    job_id text,
    CONSTRAINT calib_gate_pkey PRIMARY KEY (curve_id, method, model_name, gate),
    CONSTRAINT calib_gate_curve_id_method_model_name_fkey FOREIGN KEY (curve_id, method, model_name) REFERENCES madi_results.calib_fit(curve_id, method, model_name) ON DELETE CASCADE
);

-- ---- NEW: calib_hyperparam (population/noise params; always written) --------
-- Sibling of calib_param (per curve/method/model/term). Bayesian-only content
-- (mu_*, sigma_*, sigma_obs, nu), stored per-curve (duplicated across the group).
CREATE TABLE IF NOT EXISTS madi_results.calib_hyperparam (
    curve_id bigint NOT NULL,
    method text NOT NULL,
    model_name text NOT NULL,
    term text NOT NULL,
    param_scope text,
    estimate numeric,
    std_error numeric,
    q_lo numeric,
    q_med numeric,
    q_hi numeric,
    job_id character varying(64),
    created_at timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_hyperparam_pkey PRIMARY KEY (curve_id, method, model_name, term),
    CONSTRAINT calib_hyperparam_curve_id_method_model_name_fkey FOREIGN KEY (curve_id, method, model_name) REFERENCES madi_results.calib_fit(curve_id, method, model_name) ON DELETE CASCADE
);

-- ---- NEW: calib_draws (posterior/sampling draws; gated by persist_draws) ----
CREATE TABLE IF NOT EXISTS madi_results.calib_draws (
    curve_id bigint NOT NULL,
    method text NOT NULL,
    model_name text NOT NULL,
    term text NOT NULL,
    param_scope text,
    sample_kind text,
    draws double precision[] NOT NULL,
    n_draws integer NOT NULL,
    job_id character varying(64),
    created_at timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_draws_pkey PRIMARY KEY (curve_id, method, model_name, term),
    CONSTRAINT calib_draws_curve_id_method_model_name_fkey FOREIGN KEY (curve_id, method, model_name) REFERENCES madi_results.calib_fit(curve_id, method, model_name) ON DELETE CASCADE
);

-- ---- NEW: calib_fit_diag (per-fit sampler/optimizer diagnostics; always) ----
-- Sibling of calib_loo (per curve/method/model). Distinct from calib_diagnostics
-- (the detection-limit/LOQ table). Bayesian + frequentist column banks coexist.
CREATE TABLE IF NOT EXISTS madi_results.calib_fit_diag (
    curve_id bigint NOT NULL,
    method text NOT NULL,
    model_name text NOT NULL,
    fit_seconds numeric,
    n_iterations integer,
    converged boolean,
    fit_seed bigint,
    rhat_max numeric,
    ess_bulk_min numeric,
    ess_tail_min numeric,
    n_divergent integer,
    pct_divergent numeric,
    max_treedepth_hit integer,
    ebfmi_min numeric,
    hessian_condition_number numeric,
    gradient_norm numeric,
    optimizer_code integer,
    rel_tol_achieved numeric,
    job_id character varying(64),
    created_at timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_fit_diag_pkey PRIMARY KEY (curve_id, method, model_name),
    CONSTRAINT calib_fit_diag_curve_id_method_model_name_fkey FOREIGN KEY (curve_id, method, model_name) REFERENCES madi_results.calib_fit(curve_id, method, model_name) ON DELETE CASCADE
);


-- ============================================================================
-- Children FK -> curve_lookup(curve_id)  (per-curve, model-agnostic)
-- ============================================================================

CREATE TABLE IF NOT EXISTS madi_results.calib_diagnostics (
    curve_id bigint NOT NULL,
    method text NOT NULL,
    model_name text,
    lloq_log10 numeric,
    uloq_log10 numeric,
    lloq_conc numeric,
    uloq_conc numeric,
    shape_lloq_log10 numeric,
    shape_uloq_log10 numeric,
    shape_lloq_conc numeric,
    shape_uloq_conc numeric,
    inflect_x numeric,
    inflect_y numeric,
    inflect_x_lower numeric,
    inflect_x_upper numeric,
    lower_lod_response numeric,
    upper_lod_response numeric,
    lower_lod_log10_conc numeric,
    upper_lod_log10_conc numeric,
    lower_lod_conc numeric,
    upper_lod_conc numeric,
    mdc_lower_log10 numeric,
    mdc_upper_log10 numeric,
    mdc_lower_conc numeric,
    mdc_upper_conc numeric,
    rdl_lower_log10 numeric,
    rdl_upper_log10 numeric,
    rdl_lower_conc numeric,
    rdl_upper_conc numeric,
    pcov_threshold numeric,
    cv_x_max numeric,
    alpha numeric,
    job_id character varying(64),
    created_at timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_diagnostics_pkey PRIMARY KEY (curve_id, method),
    CONSTRAINT calib_diagnostics_curve_id_fkey FOREIGN KEY (curve_id) REFERENCES madi_results.curve_lookup(curve_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS madi_results.calib_grid (
    curve_id bigint NOT NULL,
    method text NOT NULL,
    point_index integer NOT NULL,
    model_name text,
    log10_concentration numeric,
    concentration numeric,
    predicted_response numeric,
    ci_lower numeric,
    ci_upper numeric,
    predicted_concentration numeric,
    se_concentration numeric,
    pcov numeric,
    pcov_rmse numeric,
    pcov_pass boolean,
    d2y_dx2 numeric,
    noise_mode text,
    job_id character varying(64),
    created_at timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_grid_pkey PRIMARY KEY (curve_id, method, point_index),
    CONSTRAINT calib_grid_curve_id_fkey FOREIGN KEY (curve_id) REFERENCES madi_results.curve_lookup(curve_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS madi_results.calib_samples (
    curve_id bigint NOT NULL,
    method text NOT NULL,
    sampleid character varying(128) NOT NULL,
    patientid character varying(128) NOT NULL DEFAULT '__none__'::character varying,
    timeperiod character varying(64) NOT NULL DEFAULT '__none__'::character varying,
    dilution character varying(64) NOT NULL DEFAULT '__none__'::character varying,
    predicted_concentration numeric,
    final_concentration numeric,
    se_concentration numeric,
    pcov numeric,
    pcov_rmse numeric,
    pcov_pass boolean,
    job_id character varying(64),
    created_at timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_samples_pkey PRIMARY KEY (curve_id, method, sampleid, patientid, timeperiod, dilution),
    CONSTRAINT calib_samples_curve_id_fkey FOREIGN KEY (curve_id) REFERENCES madi_results.curve_lookup(curve_id) ON DELETE CASCADE
);


-- ============================================================================
-- Children with NO FK (input echoes; keyed by curve/method/well)
-- ============================================================================

CREATE TABLE IF NOT EXISTS madi_results.calib_standards (
    curve_id bigint NOT NULL,
    method character varying(20) NOT NULL,
    well character varying(16) NOT NULL,
    dilution numeric NOT NULL,
    concentration numeric,
    log10_concentration numeric,
    response_model numeric,
    assay_response_raw numeric,
    included boolean NOT NULL,
    exclusion_reason character varying(20) NOT NULL DEFAULT 'none'::character varying,
    job_id text,
    created_at timestamp with time zone NOT NULL DEFAULT now(),
    mask_reason character varying,
    CONSTRAINT calib_standards_pkey PRIMARY KEY (curve_id, method, well, dilution)
);

CREATE TABLE IF NOT EXISTS madi_results.calib_blanks (
    curve_id bigint NOT NULL,
    method character varying(20) NOT NULL,
    well character varying(16) NOT NULL,
    response_model numeric,
    assay_response_raw numeric,
    included boolean NOT NULL,
    exclusion_reason character varying(20) NOT NULL DEFAULT 'none'::character varying,
    job_id text,
    created_at timestamp with time zone NOT NULL DEFAULT now(),
    mask_reason character varying,
    CONSTRAINT calib_blanks_pkey PRIMARY KEY (curve_id, method, well)
);


-- ============================================================================
-- Settings cascade + logs
-- ============================================================================

CREATE TABLE IF NOT EXISTS madi_results.calib_settings (
    calib_settings_id bigint NOT NULL DEFAULT nextval('madi_results.calib_settings_calib_settings_id_seq'::regclass),
    project_id integer NOT NULL DEFAULT '-1'::integer,
    study_accession character varying(15) NOT NULL DEFAULT '__none__'::character varying,
    experiment_accession character varying(15) NOT NULL DEFAULT '__none__'::character varying,
    feature character varying(15) NOT NULL DEFAULT '__none__'::character varying,
    antigen character varying(64) NOT NULL DEFAULT '__none__'::character varying,
    param_name character varying(64) NOT NULL,
    param_data_type character varying(15) NOT NULL,
    param_integer_value integer,
    param_numeric_value numeric,
    param_boolean_value boolean,
    param_character_value text,
    param_user text,
    updated_at timestamp with time zone DEFAULT now(),
    param_group character varying(24),
    CONSTRAINT calib_settings_scope_param_uniq UNIQUE (project_id, study_accession, experiment_accession, feature, antigen, param_name),
    CONSTRAINT calib_settings_pkey PRIMARY KEY (calib_settings_id),
    CONSTRAINT calib_settings_antigen_needs_feature_chk CHECK ((((antigen)::text = '__none__'::text) OR ((feature)::text <> '__none__'::text))),
    CONSTRAINT calib_settings_experiment_needs_study_chk CHECK ((((experiment_accession)::text = '__none__'::text) OR ((study_accession)::text <> '__none__'::text))),
    CONSTRAINT calib_settings_feature_needs_experiment_chk CHECK ((((feature)::text = '__none__'::text) OR ((experiment_accession)::text <> '__none__'::text))),
    CONSTRAINT calib_settings_param_data_type_check CHECK (((param_data_type)::text = ANY ((ARRAY['integer'::character varying, 'numeric'::character varying, 'boolean'::character varying, 'character'::character varying])::text[]))),
    CONSTRAINT calib_settings_system_wildcard_chk CHECK (((project_id <> '-1'::integer) OR (((study_accession)::text = '__none__'::text) AND ((experiment_accession)::text = '__none__'::text) AND ((feature)::text = '__none__'::text) AND ((antigen)::text = '__none__'::text))))
);

CREATE TABLE IF NOT EXISTS madi_results.calib_settings_meta (
    param_name character varying(64) NOT NULL,
    param_group character varying(24),
    param_label character varying(256),
    param_data_type character varying(15) NOT NULL,
    param_control_type character varying(64),
    param_choices_list character varying(256),
    param_char_len numeric,
    updated_at timestamp with time zone DEFAULT now(),
    param_description text,
    CONSTRAINT calib_settings_meta_pkey PRIMARY KEY (param_name),
    CONSTRAINT calib_settings_meta_param_data_type_check CHECK (((param_data_type)::text = ANY ((ARRAY['integer'::character varying, 'numeric'::character varying, 'boolean'::character varying, 'character'::character varying])::text[])))
);

CREATE TABLE IF NOT EXISTS madi_results.calib_settings_migration_log (
    project_id integer,
    study_accession character varying(15),
    experiment_accession text,
    antigen character varying(64),
    model_form_list text,
    standard_curve_concentration numeric,
    l_asy_min_constraint numeric,
    l_asy_max_constraint numeric,
    l_asy_constraint_method character varying(40),
    feature text,
    feature_source text,
    n_feat_in_context bigint
);

CREATE TABLE IF NOT EXISTS madi_results.calib_settings_studycfg_log (
    project_id integer,
    study_accession character varying,
    param_name character varying,
    param_group character varying(24),
    param_data_type text,
    param_integer_value integer,
    param_boolean_value boolean,
    param_character_value text,
    study_raw character varying(15),
    param_raw character varying(50)
);


-- ============================================================================
-- Sequence ownership (deferred: references calib_settings, created above).
-- OWNED BY ties the sequence's lifecycle to the column (dropped with it).
-- ============================================================================
ALTER SEQUENCE madi_results.calib_settings_calib_settings_id_seq
    OWNED BY madi_results.calib_settings.calib_settings_id;

-- OWNER is environment-specific (role 'd78039e' exists on the current DB but may
-- not on a fresh target — uncomment and set the correct role when applying there):
-- ALTER SEQUENCE madi_results.calib_settings_calib_settings_id_seq OWNER TO d78039e;
