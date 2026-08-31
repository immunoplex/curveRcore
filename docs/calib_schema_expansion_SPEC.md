# calib\_\* Schema Expansion — Design Spec (rev. 2)

**Goal:** let the deployed **i-spi-compute-sim** worker expose what a
trustworthy, worker-based SBC needs — posterior **draws** and the
**population/noise parameters** — plus sampler diagnostics, without
breaking existing `calib_*` consumers.

**Rev. 2 note:** revised after seeing the live DDLs for `calib_fit`,
`calib_loo`, `calib_run`. Three corrections vs rev. 1: (a)
`selection_weight` already exists — no new column; (b) model-selection +
degenerate-fit signals are already largely persisted
(`calib_fit.converged/eligible/selection_weight`,
`calib_loo.elpd_diff/ se_diff/pareto_k_*`); (c) the contract’s universal
pattern is **per-curve + FK→calib_fit + ON DELETE CASCADE**, so new
tables adopt that grain (population values duplicated across the group),
not a group key.

**Ownership boundary:** `calib_*` schema is defined/enforced by the
contract in `curveRcore`; settings are read in
`curveRfreq`/`curveRbayes`, which generate the params/draws; the worker
flattens to the DB. Each change is tagged accordingly. `[CONFIRM]` marks
the one item that depends on curveRbayes internals.

------------------------------------------------------------------------

## 0. Deciding principle: conform to the existing contract pattern

Every existing `calib_*` table is keyed
`(curve_id, method, model_name[, term])`, FKs to `calib_fit`, cascades
on delete, and uses `job_id varchar(64)` +
`created_at timestamptz DEFAULT now()`. New tables MUST match this, so:

| Parameter class | Examples | Home | Written |
|----|----|----|----|
| Per-curve terms | a, b, c, d, g | `calib_param` (EXISTING, unchanged) | always |
| Population / noise | sigma_obs, nu, sigma_log_b, mu_a…mu_g, cross-plate SDs | **`calib_hyperparam`** (NEW) | always |
| Posterior draws | draws of any term | **`calib_draws`** (NEW) | gated (`persist_draws`) |
| Sampler diagnostics | rhat, ess, divergences, treedepth, seconds, seed | **`calib_mcmc_diag`** (NEW) | always |

**Grain (CONFIRMED against curveRbayes source — fit_bayes.R /
stan_data.R / priors.R):** - **Per-curve** `a/b/c/d/g` are Stan
`a[p]/b[p]/c_par[p]/d[p]/g[p]`, indexed by plate `p`; one Stan fit
covers the whole group (`N_plates = n_curves`). These are the
`calib_param` terms. (`c_par` is emitted as `c`.) - **Per-group**
`sigma_obs`, `nu`, and the population means (`mu_*`) are single,
non-plate-indexed scalars in the one group fit → group-grain. They are
NOT extracted today (`extract_curve_params` pulls only the plate-indexed
params), which is why they’re absent from `calib_param`. The full
posterior is present in `posterior::as_draws_df(mcmc_fit$draws())`; the
package must add extraction of the non-indexed columns to populate
`calib_hyperparam`/`calib_draws`. - **`sigma_log_b` IS a free parameter
(CORRECTION to config §4).** The `.stan parameters{}` block declares
`real<lower=0> sigma_log_b;` with prior `sigma_log_b ~ normal(0, 0.5)`
in `model{}`. The “0.5” is the **prior scale, not a fixed value** — the
cross-plate SD of log_b is estimated, has a posterior, and CAN be
SBC-ranked directly. No `free_crossplate_sd` switch is needed (already
free). The stage-1 “b non-uniformity ⇒ fixed-pooling” framing must be
revised accordingly (see §6). - Because population/noise params are
group-scalars, storing them per-curve (§1) simply duplicates them across
the group’s curves — cheap, and it keeps the contract’s per-curve +
FK→calib_fit + cascade pattern; extraction runs once per fit and is
written to each curve of the group. - **Term vocabulary (LOCKED from the
`.stan parameters{}` blocks):** - `calib_param` (per-curve, EXISTING):
`a, b, c, d, g` (`b=exp(log_b)`, `g=exp(log_g)`, `c` from `c_par`). -
`calib_hyperparam` (per-group): population hypers
`mu_a, sigma_a, mu_d, sigma_d, mu_log_b, sigma_log_b, mu_c, sigma_c` +
(5PL only) `mu_log_g, sigma_log_g`; noise `sigma_obs, nu, sigma_blank`.
Optionally `log_sigma0, log_sigma_slope` (homoscedastic study →
nuisance; include only if wanted). - Do NOT persist `raw_*`
(non-centered auxiliaries). `log_lik` already feeds `calib_loo`;
`y_pred` is available in `generated quantities` if a per-obs PPC is ever
wanted (not in this diff).

**What is ALREADY covered (no change needed):** - Model selection:
`calib_fit.selection_weight/selection_score/criterion/ score_type/is_best/is_fallback` +
`calib_loo.elpd_diff/se_diff` quantify how strongly one model beats
another per curve — enough to measure “selection drops asymmetry (4PL
over 5PL as g rises)” with zero schema change. - Degenerate-fit signals:
`calib_fit.converged/eligible` + `calib_loo.pareto_k_*` (bad/vbad
Pareto-k = influential points / unreliable fit) already exist. The ONLY
missing diagnostic is sampler-level rhat/ess/divergences →
`calib_mcmc_diag`.

------------------------------------------------------------------------

## 1. `calib_hyperparam` — population/noise params (ALWAYS written)

Mirrors `calib_param` column-for-column; adds a `param_scope` marker so
a population term is self-describing even though it’s stored per-curve.

``` sql
CREATE TABLE IF NOT EXISTS madi_results.calib_hyperparam
(
    curve_id    bigint NOT NULL,
    method      text   NOT NULL,
    model_name  text   NOT NULL,
    term        text   NOT NULL,          -- 'sigma_obs','nu','sigma_log_b','mu_a',...
    param_scope text,                      -- 'population' (informational)
    estimate    numeric,
    std_error   numeric,
    q_lo        numeric,                   -- 2.5%
    q_med       numeric,                   -- 50%
    q_hi        numeric,                   -- 97.5%
    job_id      character varying(64),
    created_at  timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_hyperparam_pkey PRIMARY KEY (curve_id, method, model_name, term),
    CONSTRAINT calib_hyperparam_fkey FOREIGN KEY (curve_id, method, model_name)
        REFERENCES madi_results.calib_fit (curve_id, method, model_name)
        MATCH SIMPLE ON UPDATE NO ACTION ON DELETE CASCADE
);
```

Kept OUT of `calib_param` on purpose: co-mingling population terms with
the per-curve model terms under one PK would let existing `a/b/c/d/g`
consumers pick up `sigma_obs` rows unless every query filters by term. A
sibling table keeps the two classes cleanly separated while reusing the
identical shape + FK + cascade.

------------------------------------------------------------------------

## 2. `calib_draws` — posterior draws (GATED by `persist_draws`, default OFF)

One table, arrays not long-form, covering both per-curve and population
terms.

``` sql
CREATE TABLE IF NOT EXISTS madi_results.calib_draws
(
    curve_id    bigint NOT NULL,
    method      text   NOT NULL,
    model_name  text   NOT NULL,
    term        text   NOT NULL,          -- a,b,c,d,g,sigma_obs,nu,sigma_log_b,mu_*
    param_scope text,                      -- 'curve' | 'population'
    draws       double precision[] NOT NULL,  -- iteration-ordered posterior draws
    n_draws     integer NOT NULL,
    job_id      character varying(64),
    created_at  timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_draws_pkey PRIMARY KEY (curve_id, method, model_name, term),
    CONSTRAINT calib_draws_fkey FOREIGN KEY (curve_id, method, model_name)
        REFERENCES madi_results.calib_fit (curve_id, method, model_name)
        MATCH SIMPLE ON UPDATE NO ACTION ON DELETE CASCADE
);
```

**Joint-draw ordering contract:** within one (curve, model) all `term`
arrays MUST be written in the same iteration order, so the joint
posterior is reconstructed by column-binding arrays. This is what lets
the 5×5 cross-parameter correlation matrix (currently JOB 3’s local
re-fit) be computed straight from the population `mu_*` arrays —
retiring the local capture.

**Sizing / why gated + arrays:** long-form is `n_draws × terms × curves`
(≈ 5 × 3500 × 2250 ≈ 39M rows/stage); arrays are one row per term (≈
11k) and a rank is a single-row read. Population-term draws are
duplicated across the group’s curves; because the table is gated OFF by
default and SBC pins one model on a bounded N, this is acceptable. \[If
storage ever bites, the documented exception is to write population-term
draw rows only for a group-representative curve — but keep summaries
(§1) per-curve for pattern consistency.\]

------------------------------------------------------------------------

## 3. `calib_fit_diag` — fit diagnostics, engine-agnostic (ALWAYS written)

Renamed from `calib_mcmc_diag` to serve BOTH engines (see §10). Common
core + engine-specific nullable columns. Per-fit grain (group;
duplicated per-curve to match the contract pattern). Fills the scorer’s
all-NA `diagnostics` slot and is the instrument for the
converged-but-degenerate fits.

``` sql
CREATE TABLE IF NOT EXISTS madi_results.calib_fit_diag
(
    curve_id          bigint NOT NULL,
    method            text   NOT NULL,     -- 'bayesian' | 'frequentist'
    model_name        text   NOT NULL,
    -- common
    fit_seconds       numeric,
    n_iterations      integer,
    converged         boolean,             -- also in calib_fit; mirrored for convenience
    fit_seed          bigint,              -- sampler seed OR optimizer start seed
    -- Bayesian (NULL for frequentist)
    rhat_max          numeric,
    ess_bulk_min      numeric,
    ess_tail_min      numeric,
    n_divergent       integer,
    pct_divergent     numeric,
    max_treedepth_hit integer,
    ebfmi_min         numeric,
    -- Frequentist (NULL for bayesian) [PROVISIONAL — confirm against curveRfreq]
    hessian_condition_number numeric,      -- ~ identifiability; the freq analog of
                                           --   the degenerate b~20 step-function ridge
    gradient_norm            numeric,
    optimizer_code           integer,      -- convergence/status code
    rel_tol_achieved         numeric,
    job_id            character varying(64),
    created_at        timestamp with time zone DEFAULT now(),
    CONSTRAINT calib_fit_diag_pkey PRIMARY KEY (curve_id, method, model_name),
    CONSTRAINT calib_fit_diag_fkey FOREIGN KEY (curve_id, method, model_name)
        REFERENCES madi_results.calib_fit (curve_id, method, model_name)
        MATCH SIMPLE ON UPDATE NO ACTION ON DELETE CASCADE
);
```

Bayesian columns populate from `fit_bayes_single`’s
`mcmc_fit$diagnostic_summary()` (`num_divergent`, `num_max_treedepth`,
`ebfmi`) + `mcmc_fit$summary()` (`rhat`, `ess_bulk`, `ess_tail`) +
`mcmc_fit$time()` + the `seed`. Frequentist columns are PROVISIONAL
pending curveRfreq source (§10).

------------------------------------------------------------------------

## 4. New settings-cascade parameters

Two new toggles, both added to curveRcore’s
[`new_study_params()`](https://immunoplex.github.io/curveRcore/reference/new_study_params.md)
(`settings.R` — the study-level fitting-options object, next to
`apply_prozone`), to the DB settings cascade `resolve_settings_batch`
reads, AND honored as job-level overrides on the submit params
(`calib_run.params` is `jsonb`, so job params are captured with the run
— the same channel `sampling=3500` uses). Precedence everywhere: **job
param \> DB setting \> `study_params` default.**

- **`persist_draws`** (boolean, default `FALSE`) — gates writing
  `calib_draws` only. `calib_hyperparam`/`calib_fit_diag` are NOT gated
  (always written, light).
- **`bayes_single_plate`** (boolean, default `FALSE`) — see §11.
  Bayesian-only; the frequentist arm is always per-plate and ignores it.
  Default `FALSE` preserves current multiplate-pooling behavior.

Add `is.logical`/length-1 validation to
[`new_study_params()`](https://immunoplex.github.io/curveRcore/reference/new_study_params.md)
and thread both through the resolve path that already surfaces
`is_log_response`, `model_form_list`, `bayes_sampling`, etc.

------------------------------------------------------------------------

## 5. Where each change lands

- **curveRcore (contract + result objects + tidy accessors):**
  - Define `calib_hyperparam`, `calib_draws`, `calib_fit_diag` in the
    schema contract + validation; ship the migration. No change to
    `calib_param`, `calib_fit`, `calib_loo`, `calib_run`.
  - **Add a `population` slot to
    `new_calibration_result_multiplate(meta, plates, population = NULL)`**
    — today it is only `meta` + `plates`, so the group-level
    hyperparams/draws/fit-diag have nowhere to live. `population` holds
    the group’s `calib_hyperparam` summary + (gated) draws + fit-diag.
  - **Add tidy accessors**
    [`tidy_hyperparam()`](https://immunoplex.github.io/curveRcore/reference/tidy_hyperparam.md),
    [`tidy_draws()`](https://immunoplex.github.io/curveRcore/reference/tidy_draws.md),
    [`tidy_fit_diag()`](https://immunoplex.github.io/curveRcore/reference/tidy_fit_diag.md)
    following the
    [`tidy_samples()`](https://immunoplex.github.io/curveRcore/reference/tidy_samples.md)/[`tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.md)
    pattern (per-result, multiplate row-binds via `.cr_rbind_fill`).
    `tidy_hyperparam` reads `x$population` and emits per-curve rows
    (fanned across the group’s `curve_id`s) so it lands in the
    per-curve + FK table.
  - Add `persist_draws` + `bayes_single_plate` to
    [`new_study_params()`](https://immunoplex.github.io/curveRcore/reference/new_study_params.md)
    (§4).
- **curveRbayes:**
  - Extract the non-indexed population/noise draws (§0 term list) → fill
    the result `population` slot → `calib_hyperparam`; extract per-fit
    sampler diagnostics
    (`diagnostic_summary()`+[`summary()`](https://rdrr.io/r/base/summary.html)+[`time()`](https://rdrr.io/r/stats/time.html)+seed)
    → `calib_fit_diag`.
  - When `persist_draws` is true, emit draw arrays → `calib_draws`.
  - When `bayes_single_plate` is true, fit per-plate (N_plates=1 loop)
    instead of the multiplate hierarchy (§11) — no `.stan` change.
- **curveRfreq (only if the freq arm is built, §10):** add the `mvrnorm`
  sampler → `calib_draws` (`sample_kind='mvn'`) and `kappa(vcov)` →
  `calib_fit_diag`; emit `sigma_resid` as a `calib_param` term. No
  hyperparam rows.
- **worker / flatten:** flatten the new tidy outputs into the new tables
  (gate draws); extend `read_calib_results()`’s table list.

------------------------------------------------------------------------

## 6. Read-back contract (what `sbc_via_sim` and JOB 3 will call)

- **SBC rank:** read `calib_draws.draws` for (curve, method, model,
  term); `rank = sum(draws < theta_tilde) + 1`. Population-term ranks
  read from a group-representative curve.
- **Correlation matrix (JOB 3, now worker-native):** read the `mu_*`
  arrays for a representative curve, cbind in term order → joint draws →
  [`cor()`](https://rdrr.io/r/stats/cor.html). No local re-fit.
- **Convergence gate:** join `calib_mcmc_diag` (+ existing
  `calib_fit.converged`, `calib_loo.pareto_k_bad/vbad`) to exclude
  non-converged/divergent fits before ranking — SBC must not rank a fit
  that didn’t converge.
- **`sigma_log_b` (CORRECTED):** it IS a monitored group parameter
  (`sigma_log_b ~ normal(0, 0.5)`), so it appears in `calib_draws` and
  is ranked like any other group scalar. The stage-1 “b non-uniformity ⇒
  *fixed* 0.5 is wrong” story is void — the cross-plate SD is estimated.
  The live SBC questions become: is the half-Normal(0,0.5) *prior* on
  `sigma_log_b` well-calibrated (rank `sigma_log_b`), and are the
  per-plate `b[p]` ranks uniform? Non-uniformity now points at prior
  mis-specification or genuine misfit, not a frozen constant.
- **SBC draws from the model PRIOR, not the anchor.** Note the priors
  are data-adaptive (`compute_dynamic_priors`: `prior_a_mu = y_min`,
  etc.) and `nu ~ gamma(2, 0.1)` (prior mean 20) — whereas the recovery
  grid fixes `nu = 2.15` from the anchor. SBC must draw the whole
  hierarchy from these model priors given a fixed calibrator design; it
  is a different regime from the recovery grid and must not reuse
  anchor-fixed truth.
- **Model-selection result (no schema change):**
  `calib_fit.selection_weight` + `calib_loo.elpd_diff/se_diff` per
  model, grouped by the truth’s g, quantify the asymmetry-dropping
  effect directly.

------------------------------------------------------------------------

## 7. Net effect

- `calib_param`, `calib_fit`, `calib_loo`, `calib_run`: **unchanged.**
- `calib_hyperparam`, `calib_mcmc_diag`: NEW, always written, light,
  per-curve + FK + cascade (contract-consistent).
- `calib_draws`: NEW, gated by `persist_draws` (default off; on for
  sim), arrays.
- Result: exact SBC ranks + joint draws + sampler diagnostics available
  FROM THE WORKER; `sbc_via_sim` and a worker-native JOB 3 become
  possible; normal production pays nothing (draws off, three
  small/among-light tables added); and the model-selection +
  degenerate-fit questions are answerable from columns that already
  exist (`selection_weight`, `elpd_diff`, `converged`, `pareto_k`) plus
  the new `rhat/divergences`.

------------------------------------------------------------------------

## 8. Deferred (NOT in this contract diff)

- **Per-curve back-calc for validation samples** (for `conc_cover_95`):
  an **injection** change in `02d` (inject the calibrator points as
  samples so the worker back-calculates them with posterior CIs), not a
  schema change. Flagged so it isn’t lost.

------------------------------------------------------------------------

## 9. Build / rollout order (contract-first)

curveRcore *defines and enforces* the schema, so the DB change flows
FROM the contract — never hand-edit the DB ahead of it (that drifts, and
enforcement then fails against a DB curveRcore didn’t author). All
changes are additive and the draws gate defaults OFF, so each step is
backward-compatible and safe in isolation.

1.  **curveRcore** — add `calib_hyperparam`, `calib_draws`,
    `calib_fit_diag` to the contract + validation; its migration creates
    the tables (the DB change is an artifact of this step). **In
    parallel:** add `persist_draws` (default false) to the settings
    cascade — a leaf that blocks nothing.
2.  **curveRbayes / curveRfreq** — extract the (now-named) non-indexed
    population/ noise params → `calib_hyperparam`; extract per-fit
    diagnostics → `calib_fit_diag`; read `persist_draws` and, when true,
    emit parameter samples → `calib_draws`.
3.  **worker / flatten** — flatten the new outputs into the new tables
    (gate draws); add the new tables to `read_calib_results()`.
4.  **harness** — `sbc_via_sim` (Bayes) and, later, its frequentist
    analog read back via §6.

Because b needs NO Stan change (§0), there is no `.stan` edit in this
sequence — it is pure plumbing plus one new setting.

------------------------------------------------------------------------

## 10. Frequentist parallel — CONFIRMED against curveRfreq source

**curveRfreq is independent per-plate NLS, not mixed-effects.**
`fit_calibration_freq_multiplate` splits by `curve_id` and fits each
plate separately
([`minpack.lm::nlsLM`](https://rdrr.io/pkg/minpack.lm/man/nlsLM.html),
base `nls` port fallback); there is **no pooling, no hierarchy, no
population parameters.** Uncertainty is propagated by the **delta
method** from `coef(fit)` + `vcov(fit)` (`predict_samples_freq`); there
is no native draws/bootstrap. Residual SE is `summary(fit)$sigma`.

Consequences for each table:

- **`calib_param`** — fully serves frequentist: `estimate = coef(fit)`,
  `std_error = sqrt(diag(vcov))`, `q_lo/q_hi` = Wald CIs. Add one
  per-curve term **`sigma_resid`** (`= summary(fit)$sigma`) — the
  frequentist noise analog. It is PER-CURVE, so it belongs here, NOT in
  `calib_hyperparam`.
- **`calib_hyperparam`** — **Bayesian-only.** The frequentist arm has no
  population params, so it writes zero rows. (The table is unchanged; it
  simply has no frequentist content.)
- **`calib_draws`** — for frequentist this is an **asymptotic-MVN
  sample** (`MASS::mvrnorm(n, coef, vcov)`), `sample_kind = 'mvn'`. The
  ingredients exist; generating the sample is a small ADD to curveRfreq.
  This enables a frequentist SBC-analog (rank the true value in the MVN
  sampling distribution) on the same plumbing and same `persist_draws`
  gate. (A parametric bootstrap is the alternative
  `sample_kind = 'bootstrap'` if preferred over the asymptotic MVN.)
- **`calib_fit_diag`** — frequentist columns CONFIRMED: `converged`,
  `aic`, `bic`, `rss`, `df_resid`, `n_iterations` (mostly already
  computed by `summarise_ensemble` / in `calib_fit`/`calib_loo`), plus
  the proposed **`hessian_condition_number = kappa(vcov(fit))`** — a
  cheap ADD and the direct frequentist analog of the degenerate-fit /
  identifiability signal (a near-singular `vcov` is the `b≈20`
  step-function ridge). `rhat/ess/divergences` are NULL for frequentist.

### The strategic point: curveRfreq per-plate IS the Framing-A comparator

The stage-1 gate is “Bayes **wins** on concordant T1,” a comparison the
study does not yet have a second arm for. curveRfreq’s per-plate
independent fit is literally “the plate is the inferential unit” — the
exact contrast to the Bayesian multiplate pooling (Framing A). So
running the frequentist arm through the same worker
(`method = 'frequentist'`) does double duty: it is both the frequentist
SBC-analog target AND the gate’s missing comparator. The user’s instinct
to get ahead of it is right — it is not a separate workload, it is the
comparator.

**Confound to flag for the study (not the schema):** Bayes-multiplate vs
Freq-per-plate differs in BOTH engine (Bayes vs NLS) AND unit (pooled vs
per-plate). To isolate “pooling helps,” a third arm — Bayes
**single-plate** (same engine, no pooling) — cleanly separates the unit
effect from the engine effect. All three run through the same worker;
the schema serves all three unchanged. Decide before
`04_compute_metrics` which arms the paper reports.

### Two small curveRfreq ADDs (only if the freq arm is built)

1.  **MVN sampler** — `MASS::mvrnorm(n_draws, coef(fit), vcov(fit))` per
    curve, emitted to `calib_draws` with `sample_kind='mvn'`, gated by
    `persist_draws`. (Enables the frequentist SBC-analog.)
2.  **`hessian_condition_number`** — `kappa(vcov(fit))` into
    `calib_fit_diag`. (Identifiability / degenerate-fit signal;
    trivial.)

Everything else the frequentist arm needs already exists in curveRfreq
or the existing `calib_*` columns. No new tables beyond §1–3;
`calib_hyperparam` stays Bayesian-only; the migration still serves both
engines in one pass.

------------------------------------------------------------------------

## 11. Bayes single-plate arm (new method, default OFF)

**Decision:** add a Bayesian per-plate (unpooled) fitting mode, toggled
by `bayes_single_plate` (§4, default `FALSE`). It is the third arm that
isolates the *pooling* effect from the *engine* effect.

**Mechanism (no `.stan` change):** when `bayes_single_plate = TRUE`,
curveRbayes fits **each plate as its own `N_plates = 1` fit** — a loop
over `curve_id`, the same structure curveRfreq already uses — reusing
the existing hierarchical `.stan` models unchanged. With one plate, the
population level is just the plate’s own params under the hyperprior
(the cross-plate SDs `sigma_*` are prior-reverted, as expected when
there is nothing to pool). Nothing in the model or schema changes; only
the fit orchestration branches. Default `FALSE` keeps the current single
multiplate hierarchy per group.

**Grain:** each single-plate fit is naturally per-curve, so
`calib_param`, `calib_hyperparam`, `calib_draws`, `calib_fit_diag` all
populate per-curve with no fan-out (unlike the multiplate arm, where
group params are duplicated across the group). The schema serves it
unchanged.

**The three arms this completes** (all through the same worker, one
schema):

| Arm                | Engine   | Unit           | Toggle                      |
|--------------------|----------|----------------|-----------------------------|
| Bayes multiplate   | Bayesian | pooled (group) | default                     |
| Bayes single-plate | Bayesian | per-plate      | `bayes_single_plate = TRUE` |
| Freq per-plate     | NLS      | per-plate      | `method = 'frequentist'`    |

Clean attributions the paper can now make: - **Bayes-multiplate vs
Bayes-single-plate** → the *pooling* effect, engine held fixed. This is
the clean Framing-A test (“the plate is the wrong inferential unit”). -
**Bayes-single-plate vs Freq-per-plate** → the *engine* effect, unit
held fixed. - **Bayes-multiplate vs Freq-per-plate** → the deployed
head-to-head (what the stage-1 “Bayes wins” gate compares in production
terms).

For SBC, note the same distinction matters: SBC of the **multiplate**
model draws the whole hierarchy (incl. `sigma_*`) from the prior; SBC of
the **single-plate** model draws a single plate’s params under the
hyperprior. Both run through `sbc_via_sim` with `bayes_single_plate` set
accordingly.
