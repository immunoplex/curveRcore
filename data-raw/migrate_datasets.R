#!/usr/bin/env Rscript
# =============================================================================
# migrate_datasets.R
#
# Copies bead_assay_example and elisa_assay_example from curveRfreq into
# curveRcore/data/. Run this ONCE from the curveRcore root after installing
# curveRfreq.
#
# Usage:
#   Rscript data-raw/migrate_datasets.R
#
# Prerequisite:
#   install.packages("curveRfreq")  # or devtools::install_github(...)
# =============================================================================

message("Migrating datasets from curveRfreq -> curveRcore ...")

if (!requireNamespace("curveRfreq", quietly = TRUE)) {
  stop("curveRfreq must be installed. Run:\n",
       "  devtools::install_github('immunoplex/curveRfreq')")
}

# Load datasets from curveRfreq
data("bead_assay_example", package = "curveRfreq", envir = environment())
data("elisa_assay_example", package = "curveRfreq", envir = environment())

# Verify structure
stopifnot(is.list(bead_assay_example))
stopifnot(all(c("standards", "blanks", "samples", "curve_id_lookup",
                "response_var", "indep_var") %in% names(bead_assay_example)))

stopifnot(is.list(elisa_assay_example))
stopifnot(all(c("standards", "blanks", "samples", "curve_id_lookup",
                "response_var", "indep_var") %in% names(elisa_assay_example)))

# Save to curveRcore/data/
save(bead_assay_example, file = "data/bead_assay_example.rda", compress = "xz")
save(elisa_assay_example, file = "data/elisa_assay_example.rda", compress = "xz")

message("Done. Saved:")
message("  data/bead_assay_example.rda  (",
        round(file.size("data/bead_assay_example.rda") / 1024, 1), " KB)")
message("  data/elisa_assay_example.rda (",
        round(file.size("data/elisa_assay_example.rda") / 1024, 1), " KB)")
message("\nRemember to add LazyData: true to DESCRIPTION.")
