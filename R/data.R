# =============================================================================
# data.R — Dataset documentation
# =============================================================================


#' Bead-based immunoassay example dataset
#'
#' A named list of simulated multi-plate bead-based immunoassay data
#' spanning 6 plates across 2 analytes (alpha and beta), measured on
#' 3 replicate plates each. Contains 20 patient samples per plate across
#' 3 timepoints (baseline, month3, month6) and 2 treatment groups.
#'
#' @format A named list with six elements:
#' \describe{
#'   \item{standards}{Data frame (60 rows x 8 cols). Standard curve data
#'     with columns: \code{curve_id}, \code{stype}, \code{sampleid},
#'     \code{well}, \code{dilution}, \code{mfi},
#'     \code{assay_response_variable}, \code{assay_independent_variable}.}
#'   \item{blanks}{Data frame (24 rows x 7 cols). Blank well measurements
#'     (4 per plate).}
#'   \item{samples}{Data frame (120 rows x 13 cols). Patient samples with
#'     columns including \code{timeperiod}, \code{patientid}, \code{agroup},
#'     \code{dilution}, \code{mfi}.}
#'   \item{curve_id_lookup}{Data frame (6 rows x 5 cols). Maps \code{curve_id}
#'     to \code{antigen}, \code{study_accession}, \code{experiment_accession},
#'     \code{plate}.}
#'   \item{response_var}{Character: \code{"mfi"}.}
#'   \item{indep_var}{Character: \code{"concentration"}.}
#' }
#'
#' @source Synthetic data. Standard curves simulated using fitted parameters
#'   from real Luminex data with realistic plate-to-plate variability
#'   and heteroscedastic noise. See \code{data-raw/migrate_datasets.R}.
#'
#' @examples
#' data(bead_assay_example)
#' str(bead_assay_example, max.level = 1)
#' head(bead_assay_example$standards)
"bead_assay_example"


#' ELISA example dataset
#'
#' A named list of simulated multi-plate ELISA data for a single analyte
#' (alpha), spanning 6 plates. Plates 1-3 use a 5-parameter logistic
#' standard curve; plates 4-6 use Gompertz, reflecting realistic
#' between-plate variability in curve shape.
#'
#' @format A named list with six elements:
#' \describe{
#'   \item{standards}{Data frame (60 rows x 8 cols). Standard curve data
#'     with columns: \code{curve_id}, \code{stype}, \code{sampleid},
#'     \code{well}, \code{dilution}, \code{od},
#'     \code{assay_response_variable}, \code{assay_independent_variable}.}
#'   \item{blanks}{Data frame (24 rows x 7 cols). Blank well measurements.}
#'   \item{samples}{Data frame (120 rows x 12 cols). Patient samples at
#'     fixed serum dilution 1:400.}
#'   \item{curve_id_lookup}{Data frame (6 rows x 5 cols). Maps \code{curve_id}
#'     to \code{antigen}, \code{study_accession}, \code{experiment_accession},
#'     \code{plate}.}
#'   \item{response_var}{Character: \code{"od"}.}
#'   \item{indep_var}{Character: \code{"concentration"}.}
#' }
#'
#' @source Synthetic data with biologically plausible parameters and
#'   proportional CV matching typical ELISA plate-reader behaviour.
#'
#' @examples
#' data(elisa_assay_example)
#' str(elisa_assay_example, max.level = 1)
#' head(elisa_assay_example$standards)
"elisa_assay_example"
