# Contributing to the curveR ecosystem

This document governs **all four packages**: `curveRcore`, `curveRfreq`,
`curveRbayes`, and `curveRweights`, plus the `curveR` hub. It is the
single source of truth for the conventions every contribution must
follow. The canonical copy lives in the hub repo
(`github.com/immunoplex/curveR`); each spoke either vendors a copy or
links to it from its own `CONTRIBUTING.md` stub. If a spoke copy and
this one ever disagree, **this one wins**.

Org: `immunoplex`. Hub site: <https://immunoplex.github.io/curveR>.

------------------------------------------------------------------------

## 1. What the ecosystem is

Four packages that share one data contract:

| Package | Role | Stan? |
|----|----|----|
| `curveRcore` | Preprocessing, forward models, inverses, gradients, eligibility gating, **the S3 classes** `calibration_result` / `calibration_result_multiplate`, and shared extractors/helpers | No |
| `curveRfreq` | Frequentist multi-start LM NLS fitting; returns the core classes | No |
| `curveRbayes` | Bayesian hierarchical fitting via Stan; returns the **same** core classes | Yes (Suggests) |
| `curveRweights` | Precision-weight estimation (`sigma_i = phi * se_i^beta1`, `w_i = 1/sigma_i^2`) from the core classes | Yes (Suggests) |
| `curveR` (hub) | Meta-package + aggregated pkgdown site + umbrella vignette | — |

Dependency graph — **read §3 before touching any dependency edge**:

    curveRcore   (anchor — depends on nothing internal)
    curveRfreq    → curveRcore
    curveRbayes   → curveRcore        (+ Suggests cmdstanr)
    curveRweights → curveRcore        (+ Suggests cmdstanr)
    curveR (hub)  → curveRcore, curveRfreq ; Suggests curveRbayes, curveRweights, cmdstanr

------------------------------------------------------------------------

## 2. The data contract (do not break without a coordinated major release)

Both fitting engines return an identical object so downstream code is
engine-agnostic. Any change to the shapes below is a **breaking contract
change** and triggers a coordinated major bump of all four packages
(§5).

    calibration_result_multiplate
    ├── $meta    : method, package, curve_ids, n_curves, response_var,
    │              is_log_response, is_log_independent, timestamp
    └── $plates  : named list keyed by curve_id (character) → calibration_result

    calibration_result
    ├── $meta, $ensemble[[model]], $selection, $detection_limits, $grid, $samples

**`$grid` columns** (per converged model):
`log10_concentration | concentration | x_fit | predicted_response | ci_lower | ci_upper | predicted_concentration | se_concentration | pcov | pcov_rmse | pcov_pass | d2y_dx2`

**`$samples` columns**:
`[original sample columns carried through] | raw_assay_response | observed_response_fit | predicted_log10_concentration | predicted_concentration | final_concentration | se_concentration | pcov | pcov_rmse | pcov_pass`

Two contract facts contributors break most often:

- **`pcov` is a capped percent %CV; `se_concentration` is the uncapped
  modeling-scale SD.** By construction
  `pcov = se_concentration × ln(10) × 100`, then capped at `cv_x_max`
  (default 150). Downstream modeling of residual scale (e.g.
  `curveRweights`) must use **`se_concentration`**, not `pcov` — the cap
  destroys the precision gradient. The single `pcov ↔︎ se_concentration`
  conversion lives in **one** core helper
  ([`pcov_from_se()`](https://immunoplex.github.io/curveRcore/reference/pcov_se_conversion.md)
  /
  [`se_from_pcov()`](https://immunoplex.github.io/curveRcore/reference/pcov_se_conversion.md));
  never reimplement it.
- **Design-group metadata is not in the core object.** Columns like
  `timeperiod` / `cohort_arm` appear in `$samples` **only** if they were
  columns in the `samples` data frame passed to the fitting call (the
  engine carries the original sample columns through). Code that needs
  design groups must require them on the input or join them by
  `sampleid` — it cannot recover them from the calibration result.

------------------------------------------------------------------------

## 3. Frozen architectural rules

These are decided. A PR that violates one will be rejected regardless of
how clever it is.

1.  **Dependency direction is sacred — but core may still evolve.**
    `curveRcore` never depends on anything internal; the hub depends on
    the spokes, never the reverse; a spoke may depend only on
    `curveRcore`. Core *is* allowed to change for general reasons —
    efficiency, or collecting common functions to remove duplication
    across spokes — provided the change is generally useful and adds no
    downward dependency. “Generally useful” is the test: a helper that
    only `curveRweights` could ever want does not belong in core; an
    extractor or conversion every consumer shares does.
2.  **Cross-package accessors live in `curveRcore`.** The canonical
    extractors
    [`tidy_samples()`](https://immunoplex.github.io/curveRcore/reference/tidy_samples.md)
    and
    [`tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.md)
    and the
    [`pcov_from_se()`](https://immunoplex.github.io/curveRcore/reference/pcov_se_conversion.md)
    /
    [`se_from_pcov()`](https://immunoplex.github.io/curveRcore/reference/pcov_se_conversion.md)
    helpers are owned by core. Consumers call them; they do not reach
    into raw object internals or re-derive these themselves.
3.  **`pcov` stays capped; expose the uncapped SD instead.** Do not
    remove the cap from `pcov` — it serves the eligibility gates, plots,
    and `pcov_pass` semantics. The uncapped modeling-scale SD is
    `se_concentration` (and, if the cap is ever found to leak into it,
    an explicit `se_concentration_raw` column). Extending core to
    *guarantee* an uncapped SD is allowed and additive; capping
    `se_concentration` is not.
4.  **Do not mutate the core S3 class from a downstream package.**
    `curveRweights` returns a standalone `precision_weights` object plus
    a `join_weights()` helper; it does not write fields back into
    `calibration_result`. Extending core (rule 1) and
    mutating-core-from-a-spoke are different things: the first is
    sometimes right, the second never is.
5.  **Naming and style.** Function verbs are `fit_*`, `predict_*`,
    `as_*`, `compare_*`; columns are `snake_case`; logical flags mirror
    the core (`is_log_response`, `is_log_independent`, …). Keep output
    column names identical across the frequentist and Bayesian paths so
    downstream code stays engine-agnostic.
6.  **Stan is a `Suggests`, never an `Imports`.** `curveRbayes` and
    `curveRweights` keep `cmdstanr` in `Suggests`; CmdStan is a non-CRAN
    external dependency. Functions that need Stan must check for it at
    entry and fail with an actionable message
    ([`cmdstanr::install_cmdstan()`](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html) +
    the install-docs URL), never a low-level stack trace. Any feature
    gated on Stan must be reflected in the capability matrix in the hub
    `installation` article.

------------------------------------------------------------------------

## 4. The contract is shared — coordinate cross-package changes

Because four packages share one object, some changes cannot be made in a
single repo in isolation:

- **Adding a column to `$grid` / `$samples`** is additive and
  non-breaking *if* both engines emit it and no consumer assumes a fixed
  column set. Land it in `curveRcore` first (constructors + extractors),
  then in both engines, then announce it in the hub `NEWS.md`.
- **Renaming or removing a contract field/column** is breaking. It
  requires a coordinated major release (§5) and an entry in the hub
  `NEWS.md` “Breaking changes”.
- **New core helper that a spoke needs** must be exported from
  `curveRcore` and released (even a minor bump) *before* the spoke that
  depends on it can declare the new `curveRcore (>= x.y.z)` minimum in
  its `DESCRIPTION`.

When in doubt, open an issue on the **hub** repo describing the
cross-package change before writing code in any spoke.

------------------------------------------------------------------------

## 5. Versioning and releases

- **Scheme:** coordinated `major.minor.patch` across all four packages.

- **Major:** bump **all four together** on a breaking contract change (a
  `$grid`/`$samples` shape change, a removed/renamed field, a changed
  default that alters results).

- **Minor / patch:** may diverge per package. Additive core work (new
  extractor, new helper, a new additive column) is a **minor** core bump
  and does not force a major elsewhere.

- **Compatibility matrix:** every release records the exact four package
  versions and notable changes in the hub `NEWS.md`. Template:

  ``` markdown
  # curveR <x.y.z>
  ## Package versions in this release
  - curveRcore   a.b.c
  - curveRfreq   a.b.c
  - curveRbayes  a.b.c
  - curveRweights a.b.c
  ## Notable / breaking changes
  - ...
  ```

- **Release flow:** tag a spoke release → the spoke’s GitHub Actions
  release workflow fires a `repository_dispatch` to `immunoplex/curveR`
  (`<spoke>-release`) → the hub rebuilds and redeploys its aggregated
  pkgdown site. Don’t hand-rebuild the hub site.

- **CRAN / r-universe:** `curveRcore`, `curveRfreq`, and the `curveR`
  meta-package target CRAN; the two Stan packages (`curveRbayes`,
  `curveRweights`) live on r-universe + GitHub because `cmdstanr` is not
  on CRAN. The meta-package’s `Suggests` (not `Imports`) on the Stan
  packages is what keeps a CRAN-only install valid — do not “fix” this
  by promoting them to `Imports`.

------------------------------------------------------------------------

## 6. Local development

Standard R-package workflow (adjust tooling to taste, but keep the
checks):

``` r
# from the package root
devtools::load_all()        # or pkgload::load_all()
devtools::document()        # roxygen2 — regenerate man/ and NAMESPACE
devtools::test()            # testthat
devtools::check()           # R CMD check — must be clean before PR
pkgdown::check_pkgdown()    # reference index complete, no orphaned topics
```

For the two Stan packages you also need CmdStan:

``` r
# one-time, if you intend to run/extend Bayesian or weighting code or knit their vignettes
cmdstanr::install_cmdstan()
cmdstanr::cmdstan_version()  # verify
```

Conventions for code that touches Stan-gated paths: guard tests with
`testthat::skip_if_not(requireNamespace("cmdstanr", quietly = TRUE))`
and a CmdStan availability check, so the suite stays green on machines
without Stan.

------------------------------------------------------------------------

## 7. Tests

- Every new exported function needs tests; every bug fix needs a
  regression test.
- Contract-adjacent changes need a test asserting the `$grid` /
  `$samples` column set, so an accidental shape change fails loudly.
- Cross-package behavior (e.g. `curveRweights` consuming a
  `calibration_result_multiplate`) gets at least one end-to-end smoke
  test built from `bead_assay_example`.
- Keep Stan-dependent tests fast (few iterations / small data) and
  skippable as above.

------------------------------------------------------------------------

## 8. Documentation and vignettes

- Roxygen for every exported object, with cross-package `@seealso` links
  so the pkgdown reference forms a web (e.g. a fitting function points
  at
  [`curveRcore::preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)
  and at the other engine).
- **Spoke vignettes are package-local**: a spoke demonstrates itself
  (and at most its direct `curveRcore` dependency), never the full
  multi-package workflow.
- **The umbrella / methods-comparison vignette lives only in the hub.**
  It is the one place the whole ecosystem runs together; keeping it off
  the spokes is deliberate and drives users to the hub. Do not duplicate
  it into a spoke.
- Stan-backed vignettes use the agreed precompute/caching strategy so CI
  site builds stay tractable (see the hub’s site README / Issue-2 plan).
- Every spoke `README` opens with the three-line ecosystem banner
  pointing at the hub URL.

------------------------------------------------------------------------

## 9. Pull-request checklist

Before opening a PR, confirm:

No dependency edge points the wrong way (§3.1); core gained nothing
spoke-specific.

Shared extractors/helpers were used, not reimplemented (§3.2);
`pcov`/`se` handling routes through the core helper (§2).

If `$grid` / `$samples` shape changed: both engines updated, additive vs
breaking classified, hub `NEWS.md` updated, version bumps planned
(§4–§5).

Naming/style conventions followed (§3.5); freq and Bayes output columns
still match.

Stan stays in `Suggests`; Stan-gated functions fail gracefully;
capability docs updated (§3.6).

[`devtools::document()`](https://devtools.r-lib.org/reference/document.html),
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html),
[`devtools::check()`](https://devtools.r-lib.org/reference/check.html),
[`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html)
all clean.

New/changed exports documented with `@seealso` cross-links; vignette
placement respects §8.

Cross-package changes were discussed in a hub issue first (§4).

------------------------------------------------------------------------

## 10. Reporting issues

- Bugs and feature requests that concern **one package** go on that
  package’s repo.
- Anything touching the **shared contract, the dependency graph,
  versioning, or the hub** goes on the `immunoplex/curveR` hub repo, so
  the discussion is visible to all four packages.
