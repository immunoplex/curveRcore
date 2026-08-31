# stage1_sanity — STAGE SUMMARY & Gate Assessment

**Date:** 2026-08-30. **Status:** the production fit→score pipeline is
PROVEN end to end at stage scale (3 cells × 50 reps × 15 plates = 2,250
curves, 150 multiplate groups, `sampling=3500`, fit through
**i-spi-compute-sim**). The stage produced a real, interpretable
recovery result. The formal stage-1 gate (“Bayes *wins* on concordant
T1 + SBC uniform”) is **NOT yet decidable** — the plate-wise comparator
arm and a trustworthy SBC are not built. This document records what is
proven, every bound that got established (and how), and the open
decisions blocking `04_compute_metrics`.

Companion docs: `02_fit_STAGE_COMPLETE.md` (fit-path proof + 10 DB
facts), `BUILD_PLAN_02fit_production.md`.

------------------------------------------------------------------------

## 1. What stage 1 tested

`stage1_sanity` = the first multi-cell, multi-replicate run. Three
T1-concordant scenarios at increasing asymmetry — **s02 symmetric (g≈1),
s05 moderate (g≈2), s08 danger_zone (g≈3)** — each 50 replicates × 15
plates. Its job was to (a) exercise the production injector/scorer at
scale, and (b) give SBC real power.

## 2. The proven pipeline path

    truth RDS (per cell,rep,plate: x, y=log10(mfi)-scale response, mu, true a/d/b/c/g)
      -> 02d inject: mfi = 10^y  (raw MFI);  log10_conc = x + center;
                     experiment_accession = per-cell scope token (one group / cell)
      -> DB assigns curve_id (identity) + multiplate_group_id (trigger, from scope)
      -> standard_for_fit (natural-key join) -> resolve_settings_batch (27 rows/curve)
      -> submit to i-spi-compute-sim with sampling=3500 (chunked 25 groups/job)
      -> worker fits each 15-plate group as ONE multiplate hierarchy (4 stacked models,
         best by LOO) -> persists calib_*
      -> 03 score: invert the fitted FORWARD curve (predicted_response -> log10_conc)
         at each calibrator point; window on shape-LOQ; compare to true_log10 = x+center

Every stage in this chain is now verified against live objects, not
inferred.

## 3. Bounds established (the evolving understanding)

Each item below is a constraint the study must respect; most were
discovered as a bug and fixed. This is the “why the bounds are what they
are” record.

**Scale bounds** 1. **Response is on the model’s log10(mfi) scale.**
`mu_curve` builds `mu∈[a,d]` from the anchor’s *fitted* (log10-mfi)
asymptotes. The worker fits with `is_log_response=TRUE`, so injection
must write **raw MFI = 10^y** or the response is double-logged and the
fit recovers `log10(true a,d)`. *Caught by* `recovery_sanity_check`:
pre-fix `est_a == log10(true_a)` to 3 decimals. 2. **x is
median-centered; the worker re-centers on the calibrator median.**
Inject at `x + center`
(`center = mean(log10(std_curve_conc/template_dils))`) to keep data in
the real assay domain; the worker’s own recentring means the `center`
shift cancels for the fitted `c`. So `true_log10 = x + center`.

**Grouping / infrastructure bounds** 3. **multiplate_group_id is
trigger-assigned from SCOPE, not plateid.** Scope = {project, study,
experiment, nominal_sample_dilution, feature, antigen, source,
wavelength}. All curves sharing the template’s scope fuse into ONE group
— a whole stage collapsed to 2,250-in-1 until we varied
**experiment_accession per (cell,replicate)**. `experiment_accession` is
`varchar(15)`, so the token is compact (`s02T1r26`). *Guard added:*
inject asserts `distinct groups == cells`. 4. **Fits MUST target
i-spi-compute-sim, never prod.** A shared secret set `ISPI_COMPUTE_URL`
to the bare `/i-spi-compute` root; the client only *warned*, and a full
stage ran on the production compute stack. *Fixed:* client auto-corrects
the bare root, **hard-errors** on any non-`-sim` URL, and `submit_stage`
refuses to submit unless `base_url` ends in `/i-spi-compute-sim`. 5.
**`sampling=3500` is locked and must be on the POST.** The probe submit
path sent no params → worker default 1000, confounding LOQ
comparability. *Fixed:* `submit_stage` passes
`sampling_to_api_params()`. 6. **The DB connection dies during long
polls.** A ~36-min poll leaves `conn` idle; preprod drops it; the
freshness guard was the first to touch it and mis-read the dead
connection as “curves missing” (alternating chunk failures). *Fixed:*
guard pings/reopens `conn` before each chunk.

**Metric bounds** 7. **Back-calc must invert the FORWARD curve, not read
`predicted_concentration` off the grid.** The grid’s
`predicted_concentration` is `f⁻¹(f(x))≈x` (a round-trip identity) —
comparing it to true x is trivially ~0. Real recovery inverts
`predicted_response -> log10_conc` at the observed response (`= y`,
since `log10(10^y)=y`). 8. **Recovery is only defined inside the
quantification window.** Eligibility `lloq/uloq` are NA in this build;
the scorer cascades to **shape LOQ** (then RDL, LOD). Points outside are
excluded. 9. **Score `mu`-inversion for the recovery gate; `y`-inversion
is “practical”.** Inverting noisy `y` folds in ν=2.15 obs noise;
inverting the noiseless `mu` isolates whether the *curve* was recovered.
10. **Unmeasurable ≠ mis-measured.** Points the *selected* model’s
response range can’t reach (extrapolated) are UNRECOVERABLE; scoring
them as huge squared errors conflated two things. *Fixed:* RMSE over
measurable points only, with `unrecoverable_rate` reported separately.
11. **Summarize recovery by MEDIAN, not mean.** With ν=2.15 (variance
barely finite) and occasional degenerate fits (see §4c), the mean and
max are tail-dominated and unstable run-to-run.

**Generator bound** 12. **`danger_zone` cell.** `sub("_.*$","",·)` (a
tier idiom) stripped the underscore inside the key `"danger_zone"` →
`"danger"` → atomic-vector subscript error. *Fixed:* name-checked direct
lookup. This had blocked every stage past stage 0.

## 4. Recovery findings (the science)

Numbers are windowed (shape-LOQ), measurable points only.

### (a) The typical curve recovers well, flat across asymmetry

`rmse_mu_log10` **median ≈ 0.10 for all three g levels** (0.102 / 0.103
/ 0.108). This tracks the stage-0 smoke group (~0.07). The production
path recovers the concordant curve.

### (b) Model selection drops asymmetry — monotonic in g (a real result)

`unrecoverable_rate` (in-window points the selected model can’t reach):

| asymmetry   | g   | unrecoverable_rate |
|-------------|-----|--------------------|
| symmetric   | 1   | 0.01               |
| moderate    | 2   | 0.02               |
| danger_zone | 3   | 0.08               |

Mechanism, confirmed on a `danger_zone` cell: the worker’s LOO/stacking
selects a **4-parameter (symmetric) model as best** for g\>1 truth. A
4PL can’t bend to g=3, so `a` recovers but `d` compresses ~0.4 low and
the fitted ceiling sits below the true high-standard responses → those
points extrapolate off the top. This is exactly what the `danger_zone`
tier exists to expose: **under Student-t noise, the deployed selector
discards the 5th (asymmetry) parameter.**

### (c) A heavy right tail in EVERY stratum — degenerate fits (not design)

Excluding unmeasurable points barely moved the mean (danger_zone mean
0.45 vs median 0.11), so the inflation is a subset of
*measurable-but-badly-fit* curves present even at g≈1. Symmetric
`rmse_mu` quantiles:
`0.10 (50%) / 0.20 (75%) / 0.57 (90%) / 1.04 (95%) / max 2.02`.

The worst symmetric cell is a **degenerate step-function fit**:
`b ≈ 18–24` (true ≈ 1.3), `c ≈ -0.89` (true ≈ 0.70), `d = 3.0896`
*identical across all 15 plates*. The population level collapsed to a
near-vertical step — a boundary/ multimodal posterior mode, most
plausibly triggered by a ν=2.15 outlier. The metric is faithful (a step
back-calculates concentration terribly); the **fit** is bad on that
replicate. This is a fitting pathology, independent of asymmetry and of
the scorer.

## 5. Scorer status

`03_score_fits.R` emits the full `.summary_schema()` scoring half per
(scenario, tier, replicate, plate), keyed by `fit_key`, plus a per-curve
`recovery` table. Confirmed columns/scales from the live fit path.
Known, documented NA gaps (worker doesn’t persist them): MCMC
`rhat/ess/divergences`, `fit_seconds`, `param_posterior`
q10/q25/q75/q90, and a true posterior back-calc CI (needs calibrator
points injected as validation samples).

## 6. SBC status — NOT yet trustworthy

The SBC panel currently HALTs on all params. **Do not act on it, and
ignore its “loosen priors / re-run anchor” line** (generic
`evaluate_sbc` advice aimed at the population-parameter anchor-SBC,
which already passed at n=300). Two reasons this panel is
uninformative: - Ranks are a **normal approximation** from `(mean, sd)`
(only 3 quantiles are persisted). `b/c/g` are log-scale/skewed →
non-uniform even under perfect calibration. - Ranks compare against
**5PL truth params while the best model is often 4PL** (§4b) → `c/d/g`
are compared across different models — meaningless by construction. The
`b` non-uniformity does have a *real* contributor (the degenerate
high-`b` fits of §4c), but it’s entangled with the two artifacts above.

A trustworthy SBC needs posterior **draws** (or a proper quantile-CDF
rank) AND a **fixed model**. Until then, windowed conc recovery is the
signal to trust.

## 7. What the gate actually requires vs. what is proven

The stage-1 gate is **“Bayes WINS on concordant T1”** — a *comparison*
against the plate-by-plate fit (Framing A: “the plate is the wrong
inferential unit”). What exists today is the **Bayesian multiplate arm
only** (single-arm recovery). **The plate-wise comparator arm is not
built**, so “wins” cannot be evaluated. Honest status: the production
path is proven correct and has produced a real result, but the gate’s
comparison is pending the comparator.

## 8. Open decisions (blocking 04_compute_metrics)

1.  **Concordance policy for T1.** Keep model selection (score back-calc
    only; `unrecoverable_rate` vs g becomes a headline result), OR force
    concordance by constraining the worker’s `model_form_list` to the
    generating form (then T1 truly tests parameter recovery and SBC of a
    known model). This decision determines whether param-recovery/SBC
    are even in the metric set.
2.  **Comparator arm.** Build the plate-wise fit (per-plate, no
    cross-plate pooling) through the same worker, so “Bayes wins” is
    computable.
3.  **Degenerate fits.** Quantify prevalence and decide handling. First
    check: are the high-`rmse_mu` curves flagged
    `converged_flag == FALSE` (then a real pipeline would reject them —
    legitimate filter), or `TRUE` (silent bad mode — more concerning,
    points at boundary priors / multimodality). Relatedly, is `b`
    hitting a model/prior bound at ~20?
4.  **SBC on real draws + fixed model** before any SBC verdict is
    reported.
5.  **Metric summarisation:** median + quantiles, not mean (ν=2.15
    tails).

## 9. One-line status

Pipeline: **proven correct end to end.** Typical-curve recovery: **good
(~0.10 log10, flat across g).** Two real phenomena surfaced —
**asymmetry dropped by model selection** (monotonic in g) and
**degenerate step-function fits** (a heavy tail in all strata). Gate:
**not yet decidable** (no comparator arm; SBC not trustworthy). These
are findings and design decisions, not remaining pipeline bugs.

------------------------------------------------------------------------

## ERRATUM (added after reading the curveRbayes `.stan` source)

Two claims above, inherited from config §4, are **wrong** and are
corrected here:

- **`sigma_log_b` is NOT a fixed constant.**
  `hierarchical_logistic5.stan` declares it in `parameters{}` with prior
  `sigma_log_b ~ normal(0, 0.5)` — the cross-plate SD of log_b is
  **estimated**, and “0.5” is the prior scale. So the §6 “b
  non-uniformity ⇒ the fixed `sigma_log_b = 0.5` pooling is wrong /
  `free_crossplate_sd` switch” framing is void: the parameter is already
  free and has a posterior. The real questions are whether the *prior*
  on `sigma_log_b` is well-calibrated and whether per-plate `b[p]` ranks
  are uniform — both answerable once SBC runs on real draws.
- **The SBC panel that HALTed on all params was doubly invalid**
  (normal-approx ranks + 4PL-vs-5PL model mismatch) and has since been
  removed from the scorer. It is not evidence of miscalibration.

The recovery findings (§4: ~0.10 median rmse_mu flat across g; asymmetry
dropped by model selection; converged-but-degenerate step-function fits
~16%) are unaffected — they come from back-calc, not from `sigma_log_b`
or the SBC panel.
