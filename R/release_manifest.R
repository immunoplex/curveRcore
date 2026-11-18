# =============================================================================
# release_manifest.R -- the calib schema version and the composite release
# manifest a job (or an app) runs under.
#
# A manifest is a small, flat-ish named list: what produced a result. The
# worker stores it in calib_run.manifest (jsonb) for every job; I-SPI shows the
# deployed one and carries it on downloads. Keys are stable; a value that cannot
# be determined is NA (JSON null), never omitted.
# =============================================================================

#' Version of the calib schema defined by this package
#'
#' The version recorded in `madi_results.schema_version` (component
#' `calib_schema`) when `R/calib_schema.sql` is applied. A database is compatible
#' with this package's contract when its latest `calib_schema` row is at least
#' this version.
#'
#' @return Character scalar, e.g. `"0.7.0"`.
#' @examples
#' calib_schema_version()
#' @export
calib_schema_version <- function() "0.7.0"


#' Composite release manifest
#'
#' Collects the versions that determine a result: R, the curveR packages,
#' CmdStan (when cmdstanr and an installation are available), the calib schema
#' version, the hash of an renv lockfile, and whatever the caller adds (the
#' component's own name, version and commit). The worker stores it with every
#' job in `calib_run.manifest`; apps show it and attach it to downloads.
#'
#' @param component Optional name of the component building the manifest, e.g.
#'   `"i-spi-compute-worker"` or `"i-spi"`.
#' @param version Optional version of that component (e.g. its image tag).
#' @param commit Optional commit (git SHA) of that component.
#' @param packages Packages whose installed versions are recorded (`NA` when not
#'   installed).
#' @param cmdstan Logical: try to record the CmdStan version through cmdstanr.
#' @param lockfile Optional path to an `renv.lock`; its MD5 is recorded.
#' @param extra Optional named list of further entries (scalars or named lists).
#'
#' @return A named list with entries `manifest_version`, `created_at` (UTC,
#'   ISO 8601), `component`, `version`, `commit`, `r_version`, `platform`,
#'   `packages` (named list), `cmdstan`, `calib_schema`, `renv_lock_md5`, and the
#'   entries of `extra`. Serialise with [manifest_json()].
#' @seealso [manifest_json()], [calib_schema_version()]
#' @examples
#' m <- release_manifest("example", "1.0.0", cmdstan = FALSE)
#' m$calib_schema
#' cat(manifest_json(m))
#' @export
release_manifest <- function(component = NA_character_, version = NA_character_,
                             commit = NA_character_,
                             packages = c("curveRcore", "curveRfreq", "curveRbayes", "curveRweights"),
                             cmdstan = TRUE, lockfile = NULL, extra = list()) {
  pkg_ver <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                                  error = function(e) NA_character_)
  cmdstan_ver <- NA_character_
  if (isTRUE(cmdstan) && nzchar(system.file(package = "cmdstanr"))) {
    f <- tryCatch(getExportedValue("cmdstanr", "cmdstan_version"), error = function(e) NULL)
    if (!is.null(f))
      cmdstan_ver <- tryCatch(as.character(f(error_on_NA = FALSE)) %||% NA_character_,
                              error = function(e) NA_character_)
    if (!length(cmdstan_ver)) cmdstan_ver <- NA_character_
  }
  lock_md5 <- NA_character_
  if (!is.null(lockfile) && length(lockfile) == 1L && file.exists(lockfile))
    lock_md5 <- unname(tools::md5sum(lockfile))

  m <- list(
    manifest_version = 1L,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    component = as.character(component %||% NA_character_),
    version = as.character(version %||% NA_character_),
    commit = as.character(commit %||% NA_character_),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    packages = stats::setNames(lapply(packages, pkg_ver), packages),
    cmdstan = cmdstan_ver,
    calib_schema = calib_schema_version(),
    renv_lock_md5 = lock_md5)
  if (length(extra)) {
    if (is.null(names(extra)) || any(!nzchar(names(extra))))
      stop("release_manifest(): every element of `extra` must be named.")
    m[names(extra)] <- extra
  }
  m
}


#' Serialise a release manifest to JSON
#'
#' Minimal, dependency-free JSON writer for manifests: named lists become
#' objects, unnamed vectors of length > 1 arrays, `NA`/`NULL` `null`. Suitable
#' for `calib_run.manifest` (jsonb) and JSON sidecars.
#'
#' @param x A manifest from [release_manifest()] (or any list of scalars and
#'   named lists).
#' @param pretty Logical: indent with two spaces.
#' @return Character scalar.
#' @export
manifest_json <- function(x, pretty = FALSE) {
  esc <- function(s) {
    s <- gsub("\\", "\\\\", s, fixed = TRUE)
    s <- gsub("\"", "\\\"", s, fixed = TRUE)
    s <- gsub("\n", "\\n", s, fixed = TRUE)
    s <- gsub("\r", "\\r", s, fixed = TRUE)
    s <- gsub("\t", "\\t", s, fixed = TRUE)
    paste0("\"", s, "\"")
  }
  scalar <- function(v) {
    if (is.null(v) || (length(v) == 1L && is.na(v))) return("null")
    if (is.logical(v)) return(if (v) "true" else "false")
    if (is.numeric(v)) return(if (is.finite(v)) format(v, digits = 15, scientific = FALSE, trim = TRUE) else "null")
    esc(as.character(v))
  }
  ind <- function(d) if (pretty) paste0("\n", strrep("  ", d)) else ""
  sep <- if (pretty) ": " else ":"
  enc <- function(v, d) {
    if (is.list(v) && !is.null(names(v))) {
      if (!length(v)) return("{}")
      items <- vapply(names(v), function(k) paste0(ind(d + 1), esc(k), sep, enc(v[[k]], d + 1)), "")
      return(paste0("{", paste(items, collapse = ","), ind(d), "}"))
    }
    if (is.list(v) || length(v) > 1L) {
      if (!length(v)) return("[]")
      items <- vapply(seq_along(v), function(i) paste0(ind(d + 1), enc(v[[i]], d + 1)), "")
      return(paste0("[", paste(items, collapse = ","), ind(d), "]"))
    }
    if (!length(v)) return("null")
    scalar(v)
  }
  enc(x, 0L)
}
