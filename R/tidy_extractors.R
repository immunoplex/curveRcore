# =============================================================================
#
# tidy_extractors.R -- canonical tidy accessors for calibration results,
# and the single pcov <-> se_concentration conversion. Used by downstream
# packages (curveRweights) instead of reaching into object internals.
#
# ============================================================================

#' Convert between posterior CV (pcov) and the log10-scale concentration SD
#'
#' In the curveR ecosystem the back-calculated concentration is reported on the
#' log10 scale (when `is_log_independent = TRUE`).  `se_concentration` is the
#' delta-method standard deviation of that log10 concentration, and `pcov` is
#' the percent coefficient of variation derived from it and then capped:
#'
#' \deqn{\mathrm{pcov} = \mathrm{se\_concentration} \times \ln(10) \times 100,
#'   \quad \text{then capped at } cv\_x\_max.}
#'
#' `se_concentration` is therefore the *uncapped* modelling-scale SD; `pcov` is
#' a censored percent. Downstream variance/weight models should consume
#' `se_concentration`, not `pcov` (the cap destroys the precision gradient).
#'
#' These helpers are the single canonical implementation of the relationship.
#' Do not reimplement it elsewhere.
#'
#' @param se Numeric vector of `se_concentration` values (log10-scale SD).
#' @param pcov Numeric vector of `pcov` values (percent).
#' @return A numeric vector of the same length.
#' @name pcov_se_conversion
#' @export
pcov_from_se <- function(se) {
  se * log(10) * 100
}

#' @rdname pcov_se_conversion
#' @export
se_from_pcov <- function(pcov) {
  pcov / (log(10) * 100)
}


#' Tidy the per-sample predictions from a calibration result
#'
#' Extracts the `$samples` table from a `calibration_result` or
#' `calibration_result_multiplate` into a single tidy data frame, attaching
#' `curve_id` for multiplate inputs. This is the canonical, supported way for
#' downstream packages (e.g. curveRweights) to read sample-level concentration
#' and precision; they must not reach into the object internals directly.
#'
#' @param x A `calibration_result` or `calibration_result_multiplate`.
#' @param ... Unused; for method extensibility.
#' @return A data frame of the per-sample predictions. For multiplate input the
#'   rows of every plate are row-bound with a `curve_id` column. Includes the
#'   carried-through original sample columns plus `predicted_concentration`,
#'   `se_concentration`, `pcov`, `pcov_pass`, etc. Returns a zero-row frame if
#'   no samples are present.
#' @seealso [tidy_grid()], [pcov_from_se()]
#' @export
tidy_samples <- function(x, ...) UseMethod("tidy_samples")

#' @rdname tidy_samples
#' @export
tidy_samples.calibration_result <- function(x, ...) {
  s <- x$samples
  if (is.null(s) || nrow(s) == 0L) return(s %||% data.frame())
  cid <- x$meta$curve_id %||% NA_character_
  if (!"curve_id" %in% names(s)) s$curve_id <- as.character(cid)
  s
}

#' @rdname tidy_samples
#' @export
tidy_samples.calibration_result_multiplate <- function(x, ...) {
  parts <- lapply(names(x$plates), function(cid) {
    cr <- x$plates[[cid]]
    if (is.null(cr)) return(NULL)
    s <- cr$samples
    if (is.null(s) || nrow(s) == 0L) return(NULL)
    if (!"curve_id" %in% names(s)) s$curve_id <- as.character(cid)
    s
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) return(data.frame())
  do.call(.cr_rbind_fill, parts)
}


#' Tidy the precision grid from a calibration result
#'
#' Extracts the best-model `$grid` (precision profile) into a tidy data frame,
#' attaching `curve_id` for multiplate inputs.
#'
#' @inheritParams tidy_samples
#' @param model Optional model name. `NULL` (default) uses each plate's selected
#'   best model (the top-level `$grid`). Otherwise pulls
#'   `ensemble[[model]]$grid`.
#' @return A data frame of grid rows with columns including
#'   `log10_concentration`, `concentration`, `predicted_concentration`,
#'   `se_concentration`, `pcov`, `pcov_rmse`, `pcov_pass`, `d2y_dx2`, and (for
#'   multiplate) `curve_id`.
#' @seealso [tidy_samples()]
#' @export
tidy_grid <- function(x, model = NULL, ...) UseMethod("tidy_grid")

#' @rdname tidy_grid
#' @export
tidy_grid.calibration_result <- function(x, model = NULL, ...) {
  g <- if (is.null(model)) x$grid else x$ensemble[[model]]$grid
  if (is.null(g) || nrow(g) == 0L) return(data.frame())
  g <- as.data.frame(g)
  cid <- x$meta$curve_id %||% NA_character_
  if (!"curve_id" %in% names(g)) g$curve_id <- as.character(cid)
  g
}

#' @rdname tidy_grid
#' @export
tidy_grid.calibration_result_multiplate <- function(x, model = NULL, ...) {
  parts <- lapply(names(x$plates), function(cid) {
    cr <- x$plates[[cid]]
    if (is.null(cr)) return(NULL)
    g <- if (is.null(model)) cr$grid else cr$ensemble[[model]]$grid
    if (is.null(g) || nrow(g) == 0L) return(NULL)
    g <- as.data.frame(g)
    if (!"curve_id" %in% names(g)) g$curve_id <- as.character(cid)
    g
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) return(data.frame())
  do.call(.cr_rbind_fill, parts)
}


# ---- internal utilities ----------------------------------------------------
# NOTE: curveRcore already defines `%||%` (and may define a row-bind helper).
# Do NOT redefine `%||%` here -- reuse the package's existing one (or base R's,
# available since R 4.4.0). The helper below is namespaced to avoid clobbering
# any existing `rbind_fill`; if curveRcore already has an equivalent, delete
# this and call that instead.

# Row-bind data frames with differing column sets (union of columns).
.cr_rbind_fill <- function(dfs) {
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(d) {
    miss <- setdiff(all_cols, names(d))
    for (m in miss) d[[m]] <- NA
    d[all_cols]
  })
  do.call(rbind, dfs2)
}


# =============================================================================
# Population / draws / fit-diagnostic accessors (schema-expansion: calib_hyper-
# param, calib_draws, calib_fit_diag). Same pattern as tidy_samples/tidy_grid:
# a per-result method + a multiplate method that row-binds with curve_id.
#
# Grain rule: the pooled multiplate arm stores group params ONCE on
# x$population; tidy_hyperparam fans them across the group's curve_ids so they
# land in the per-curve + FK calib_hyperparam table. Per-plate arms
# (single-plate Bayes, frequentist) store per-plate population on each plate's
# own $population, which is emitted with that plate's curve_id (no fan-out).
# =============================================================================

# ---- .population accessor: normalise a $population slot to a tidy frame ------
# Returns a data.frame(term, estimate, std_error, q_lo, q_med, q_hi) or NULL.
.pop_params_df <- function(pop) {
  if (is.null(pop) || is.null(pop$params) || !nrow(pop$params)) return(NULL)
  p <- as.data.frame(pop$params)
  need <- c("term", "estimate", "std_error", "q_lo", "q_med", "q_hi")
  for (m in setdiff(need, names(p))) p[[m]] <- NA
  p[need]
}


#' Tidy the population / noise parameters from a calibration result
#'
#' Extracts group-level (pooled) or per-plate population/noise parameters
#' (`mu_*`, `sigma_*`, `sigma_obs`, `nu`, and — frequentist — `sigma_resid`)
#' into a tidy data frame keyed by `curve_id`, `param_scope = "population"`.
#' Feeds the `calib_hyperparam` table. Returns a zero-row frame when the fit
#' carries no population slot (e.g. a bare frequentist per-plate fit whose only
#' noise term rides in `calib_param`).
#'
#' @inheritParams tidy_samples
#' @return Data frame `curve_id, term, param_scope, estimate, std_error,
#'   q_lo, q_med, q_hi`.
#' @seealso [tidy_draws()], [tidy_fit_diag()]
#' @export
tidy_hyperparam <- function(x, ...) UseMethod("tidy_hyperparam")

#' @rdname tidy_hyperparam
#' @export
tidy_hyperparam.calibration_result <- function(x, ...) {
  df <- .pop_params_df(x$population)
  if (is.null(df)) return(data.frame())
  df$curve_id    <- as.character(x$meta$curve_id %||% NA_character_)
  df$param_scope <- "population"
  df[c("curve_id", "term", "param_scope",
       "estimate", "std_error", "q_lo", "q_med", "q_hi")]
}

#' @rdname tidy_hyperparam
#' @export
tidy_hyperparam.calibration_result_multiplate <- function(x, ...) {
  cids <- names(x$plates)
  # Pooled arm: one group-level $population -> fan across every curve_id.
  grp <- .pop_params_df(x$population)
  if (!is.null(grp)) {
    parts <- lapply(cids, function(cid) {
      d <- grp; d$curve_id <- as.character(cid); d$param_scope <- "population"; d
    })
    out <- do.call(.cr_rbind_fill, parts)
    return(out[c("curve_id", "term", "param_scope",
                 "estimate", "std_error", "q_lo", "q_med", "q_hi")])
  }
  # Per-plate arms: collect each plate's own $population.
  parts <- lapply(cids, function(cid) {
    cr <- x$plates[[cid]]
    if (is.null(cr)) return(NULL)
    d <- tidy_hyperparam(cr)
    if (nrow(d) == 0L) return(NULL)
    d
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) return(data.frame())
  do.call(.cr_rbind_fill, parts)
}


#' Tidy the posterior / sampling draws from a calibration result
#'
#' Extracts per-parameter draw vectors (posterior draws for Bayesian fits;
#' asymptotic-MVN or bootstrap samples for frequentist) into a tidy frame with
#' one row per `(curve_id, term)` and a list-column `draws` holding the vector,
#' plus `n_draws` and `sample_kind`. Feeds the `calib_draws` table (the worker
#' packs `draws` into a `float8[]` column). Present only when the fit was run
#' with `persist_draws = TRUE`; otherwise returns a zero-row frame.
#'
#' Draw order is preserved and shared across terms within a fit, so consumers
#' may column-bind the vectors into a joint posterior (e.g. the JOB-3
#' correlation matrix).
#'
#' @inheritParams tidy_samples
#' @return Data frame `curve_id, term, param_scope, sample_kind, n_draws, draws`
#'   (`draws` is a list-column of numeric vectors).
#' @seealso [tidy_hyperparam()]
#' @export
tidy_draws <- function(x, ...) UseMethod("tidy_draws")

# Build the per-(term) draw rows for one result, given a named list of draw
# vectors, a curve_id, a scope label, and a sample_kind.
.draw_rows <- function(draw_list, curve_id, param_scope, sample_kind) {
  if (is.null(draw_list) || length(draw_list) == 0L) return(NULL)
  terms <- names(draw_list)
  data.frame(
    curve_id    = as.character(curve_id),
    term        = terms,
    param_scope = param_scope,
    sample_kind = sample_kind,
    n_draws     = vapply(draw_list, length, integer(1)),
    draws       = I(unname(draw_list)),   # list-column
    row.names   = NULL, stringsAsFactors = FALSE
  )
}

#' @rdname tidy_draws
#' @export
tidy_draws.calibration_result <- function(x, ...) {
  cid  <- x$meta$curve_id %||% NA_character_
  kind <- if (identical(x$meta$method, "frequentist")) "mvn" else "posterior"
  rows <- list()
  # per-curve model params: from the best model's $draws (if present)
  best <- x$selection$best_model_name %||% NA_character_
  ens  <- if (!is.na(best)) x$ensemble[[best]] else NULL
  if (!is.null(ens) && !is.null(ens$draws))
    rows[[length(rows) + 1L]] <- .draw_rows(ens$draws, cid, "curve", kind)
  # population/noise draws
  if (!is.null(x$population) && !is.null(x$population$draws))
    rows[[length(rows) + 1L]] <- .draw_rows(x$population$draws, cid, "population", kind)
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

#' @rdname tidy_draws
#' @export
tidy_draws.calibration_result_multiplate <- function(x, ...) {
  cids <- names(x$plates)
  parts <- list()
  # pooled group draws -> fan across curve_ids
  if (!is.null(x$population) && !is.null(x$population$draws)) {
    kind <- if (identical(x$meta$method, "frequentist")) "mvn" else "posterior"
    for (cid in cids)
      parts[[length(parts) + 1L]] <- .draw_rows(x$population$draws, cid, "population", kind)
  }
  # per-plate draws (curve params, and per-plate population for per-plate arms)
  for (cid in cids) {
    cr <- x$plates[[cid]]
    if (is.null(cr)) next
    d <- tidy_draws(cr)
    if (nrow(d) > 0L) parts[[length(parts) + 1L]] <- d
  }
  parts <- Filter(function(d) !is.null(d) && nrow(d) > 0L, parts)
  if (length(parts) == 0L) return(data.frame())
  do.call(rbind, parts)
}


#' Tidy the per-fit diagnostics from a calibration result
#'
#' Extracts sampler/optimizer diagnostics into a tidy frame keyed by
#' `curve_id`. Bayesian columns (`rhat_max`, `ess_bulk_min`, `ess_tail_min`,
#' `n_divergent`, `pct_divergent`, `max_treedepth_hit`, `ebfmi_min`) and
#' frequentist columns (`hessian_condition_number`, `gradient_norm`,
#' `optimizer_code`, `rel_tol_achieved`) coexist; the irrelevant set is NA for
#' a given engine. Common columns: `fit_seconds`, `n_iterations`, `converged`,
#' `fit_seed`. Feeds `calib_fit_diag`.
#'
#' @inheritParams tidy_samples
#' @return One row per `curve_id`.
#' @seealso [tidy_hyperparam()]
#' @export
tidy_fit_diag <- function(x, ...) UseMethod("tidy_fit_diag")

.fit_diag_cols <- c(
  "fit_seconds", "n_iterations", "converged", "fit_seed",
  "rhat_max", "ess_bulk_min", "ess_tail_min", "n_divergent",
  "pct_divergent", "max_treedepth_hit", "ebfmi_min",
  "hessian_condition_number", "gradient_norm", "optimizer_code",
  "rel_tol_achieved")

.fit_diag_row <- function(diag, curve_id) {
  row <- data.frame(curve_id = as.character(curve_id), stringsAsFactors = FALSE)
  for (col in .fit_diag_cols)
    row[[col]] <- if (!is.null(diag) && !is.null(diag[[col]])) diag[[col]] else NA
  row
}

#' @rdname tidy_fit_diag
#' @export
tidy_fit_diag.calibration_result <- function(x, ...) {
  # prefer an explicit $population$fit_diag; else the best model's fit_stats
  diag <- x$population$fit_diag
  if (is.null(diag)) {
    best <- x$selection$best_model_name %||% NA_character_
    ens  <- if (!is.na(best)) x$ensemble[[best]] else NULL
    diag <- if (!is.null(ens)) ens$fit_stats else NULL
  }
  .fit_diag_row(diag, x$meta$curve_id %||% NA_character_)
}

#' @rdname tidy_fit_diag
#' @export
tidy_fit_diag.calibration_result_multiplate <- function(x, ...) {
  cids <- names(x$plates)
  grp  <- x$population$fit_diag           # pooled arm: one diag -> fan out
  parts <- lapply(cids, function(cid) {
    cr <- x$plates[[cid]]
    if (is.null(cr)) return(NULL)
    if (!is.null(grp)) .fit_diag_row(grp, cid) else tidy_fit_diag(cr)
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) return(data.frame())
  do.call(.cr_rbind_fill, parts)
}
