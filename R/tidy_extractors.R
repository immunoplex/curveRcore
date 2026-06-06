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
