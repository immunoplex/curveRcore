# =============================================================================
# data_helpers.R — Data subsetting and curve-ID helpers
# =============================================================================


#' Filter a Dataset List by Curve ID
#'
#' Subsets all data-frame elements of a dataset list (standards, blanks,
#' samples, curve_id_lookup) to rows matching a specific `curve_id`.
#' Non-data-frame elements and data frames without a `curve_id` column
#' are passed through unchanged.
#'
#' Also attaches `$whole_standards` and `$curve_id_whole_lookup` as
#' unfiltered copies for downstream cross-plate operations.
#'
#' @param loaded_data A named list as returned by `data(bead_assay_example)`
#'   or `data(elisa_assay_example)`. Must contain data frames with a
#'   `curve_id` column.
#' @param curve_id Scalar. The `curve_id` value to filter on.
#' @param target_names Character vector. Which list elements to filter.
#'   Default: `c("standards", "blanks", "samples", "curve_id_lookup")`.
#' @param verbose Logical. Emit messages about skipped elements.
#'
#' @return A copy of `loaded_data` with the target data frames filtered
#'   to the specified `curve_id`, plus `$whole_standards` and
#'   `$curve_id_whole_lookup` preserving the full unfiltered versions.
#'
#' @export
filter_by_curve_id <- function(loaded_data,
                               curve_id,
                               target_names = c("standards", "blanks",
                                                "samples", "curve_id_lookup"),
                               verbose = FALSE) {
  filtered <- loaded_data

  # Preserve unfiltered copies for cross-plate operations
  filtered$curve_id_whole_lookup <- filtered$curve_id_lookup
  filtered$whole_standards       <- filtered$standards

  filtered[target_names] <- lapply(filtered[target_names], function(df) {
    if (!is.data.frame(df)) {
      if (verbose) message("[filter_by_curve_id] Skipping non-data.frame")
      return(df)
    }
    if (nrow(df) == 0) return(df)
    if (!("curve_id" %in% names(df))) {
      if (verbose) message("[filter_by_curve_id] No curve_id column")
      return(df)
    }
    df[as.character(df$curve_id) == as.character(curve_id), , drop = FALSE]
  })

  filtered
}
