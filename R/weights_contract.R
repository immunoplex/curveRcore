# =============================================================================
# weights_contract.R -- the data contract for madi_results.calib_weights and
# calib_weights_fit (precision weights from curveRweights).
#
# R/calib_schema.sql holds the DDL; this file holds the same definition in a
# form code can check against (workers before writing, curveRweights in its
# tests, I-SPI when building help pages). A test keeps the two in sync.
# =============================================================================

.WEIGHTS_MODELS <- c("multiplicative", "additive")

# one row per column: table, column, type, key, nullable, added_in, models, description
.weights_contract_rows <- function() {
  r <- function(table, column, type, key, nullable, added_in, models, description)
    data.frame(table = table, column = column, type = type, key = key,
               nullable = nullable, added_in = added_in, models = models,
               description = description, stringsAsFactors = FALSE)
  S <- "sample"; F <- "fit"
  rbind(
    # ---- calib_weights: per-sample (same sample hierarchy as calib_samples) ----
    r(S, "curve_id", "bigint", TRUE, FALSE, "0.5", "both", "Curve the sample was measured on (FK to curve_lookup)."),
    r(S, "method", "text", TRUE, FALSE, "0.5", "both", "Calibration engine the weights are based on: bayesian or frequentist."),
    r(S, "sampleid", "varchar(128)", TRUE, FALSE, "0.5", "both", "Sample identifier."),
    r(S, "patientid", "varchar(128)", TRUE, FALSE, "0.5", "both", "Subject identifier ('__none__' if absent)."),
    r(S, "timeperiod", "varchar(64)", TRUE, FALSE, "0.5", "both", "Time point ('__none__' if absent)."),
    r(S, "dilution", "varchar(64)", TRUE, FALSE, "0.5", "both", "Sample dilution ('__none__' if absent)."),
    r(S, "se", "double precision", FALSE, TRUE, "0.5", "both", "Measurement SE of the back-calculated concentration (log10 scale), the precision index."),
    r(S, "pcov", "double precision", FALSE, TRUE, "0.5", "both", "Posterior CV from the calibration curve; reference only."),
    r(S, "sigma", "double precision", FALSE, TRUE, "0.5", "both", "Total SD of the sample under the fitted model: multiplicative phi*se^beta1; additive sqrt(tau2 + v_meas)."),
    r(S, "w", "double precision", FALSE, TRUE, "0.5", "both", "Inverse-variance weight 1/sigma^2. Additive weights are bounded above by 1/tau2."),
    r(S, "w_norm", "double precision", FALSE, TRUE, "0.5", "both", "w divided by its mean over the fit's weighted rows (mean 1), for software that treats weights as counts. NULL for censored rows."),
    r(S, "job_id", "varchar(64)", FALSE, TRUE, "0.5", "both", "Compute job that wrote the row."),
    r(S, "created_at", "timestamptz", FALSE, TRUE, "0.5", "both", "Time the row was written."),
    r(S, "agroup", "varchar(40)", FALSE, TRUE, "0.5", "both", "Study arm, carried through from calib_samples."),
    r(S, "model", "text", FALSE, FALSE, "0.6.0", "both", "Model that produced the weights: multiplicative (every pre-0.6.0 row) or additive."),
    r(S, "weights_pkg_version", "varchar(32)", FALSE, TRUE, "0.6.0", "both", "curveRweights version that produced the row."),
    r(S, "v_meas", "double precision", FALSE, TRUE, "0.6.0", "additive", "Measurement-error variance of the sample, kappa^2 * se^(2*beta1)."),
    r(S, "tau2", "double precision", FALSE, TRUE, "0.6.0", "additive", "True between-subject (signal) variance estimated for the group."),
    r(S, "reliability", "double precision", FALSE, TRUE, "0.6.0", "additive", "tau2 / (tau2 + v_meas), in (0, 1]: the share of the sample's variance that is signal."),
    r(S, "reliability_lo", "double precision", FALSE, TRUE, "0.6.0", "additive", "Lower 95% credible limit of reliability."),
    r(S, "reliability_hi", "double precision", FALSE, TRUE, "0.6.0", "additive", "Upper 95% credible limit of reliability."),
    r(S, "censor_type", "text", FALSE, TRUE, "0.6.0", "additive", "none, left, interval or right; as in calib_samples. Censored rows are represented by their bound, not a weight."),
    r(S, "bound", "double precision", FALSE, TRUE, "0.6.0", "additive", "log10 censoring limit used in the likelihood for censored rows; NULL for quantified rows."),
    # ---- calib_weights_fit: one row per (multiplate_group_id, method) ----
    r(F, "multiplate_group_id", "uuid", TRUE, FALSE, "0.5", "both", "Group of plates of one antigen/feature/source fitted together."),
    r(F, "method", "text", TRUE, FALSE, "0.5", "both", "bayesian or frequentist."),
    r(F, "antigen", "varchar(64)", FALSE, TRUE, "0.5", "both", "Antigen."),
    r(F, "feature", "varchar(15)", FALSE, TRUE, "0.5", "both", "Feature (isotype or functional readout)."),
    r(F, "design_cols", "text", FALSE, TRUE, "0.5", "both", "Comma-joined design columns defining the saturated cells, e.g. timeperiod,agroup."),
    r(F, "scale_predictor", "text", FALSE, TRUE, "0.5", "both", "Precision index used: se or pcov. Additive always uses se."),
    r(F, "phi", "double precision", FALSE, TRUE, "0.5", "multiplicative", "Baseline scale of sigma = phi * se^beta1. NULL for additive."),
    r(F, "phi_lo", "double precision", FALSE, TRUE, "0.5", "multiplicative", "Lower 95% credible limit of phi."),
    r(F, "phi_hi", "double precision", FALSE, TRUE, "0.5", "multiplicative", "Upper 95% credible limit of phi."),
    r(F, "beta1", "double precision", FALSE, TRUE, "0.5", "both", "Exponent on se. Multiplicative: estimated. Additive: 1 when fix_beta1, else estimated."),
    r(F, "beta1_lo", "double precision", FALSE, TRUE, "0.5", "both", "Lower 95% credible limit of beta1 (NULL when fixed)."),
    r(F, "beta1_hi", "double precision", FALSE, TRUE, "0.5", "both", "Upper 95% credible limit of beta1 (NULL when fixed)."),
    r(F, "interpretation", "text", FALSE, TRUE, "0.5", "both", "Plain-language reading of the fit (multiplicative: precision regime; additive: identifiability note)."),
    r(F, "n_fit", "integer", FALSE, TRUE, "0.5", "both", "Observations used in the fit."),
    r(F, "n_eff", "double precision", FALSE, TRUE, "0.5", "both", "Effective sample size of the weights."),
    r(F, "weight_ratio", "double precision", FALSE, TRUE, "0.5", "both", "Largest / smallest weight."),
    r(F, "job_id", "varchar(64)", FALSE, TRUE, "0.5", "both", "Compute job that wrote the row."),
    r(F, "created_at", "timestamptz", FALSE, TRUE, "0.5", "both", "Time the row was written."),
    r(F, "model", "text", FALSE, FALSE, "0.6.0", "both", "multiplicative (every pre-0.6.0 row) or additive."),
    r(F, "weights_pkg_version", "varchar(32)", FALSE, TRUE, "0.6.0", "both", "curveRweights version that produced the fit."),
    r(F, "tau", "double precision", FALSE, TRUE, "0.6.0", "additive", "SD of the true between-subject variation (log10 scale)."),
    r(F, "tau_lo", "double precision", FALSE, TRUE, "0.6.0", "additive", "Lower 95% credible limit of tau."),
    r(F, "tau_hi", "double precision", FALSE, TRUE, "0.6.0", "additive", "Upper 95% credible limit of tau."),
    r(F, "kappa", "double precision", FALSE, TRUE, "0.6.0", "additive", "Calibration factor on the reported se; 1 means se is exactly the measurement SD."),
    r(F, "kappa_lo", "double precision", FALSE, TRUE, "0.6.0", "additive", "Lower 95% credible limit of kappa."),
    r(F, "kappa_hi", "double precision", FALSE, TRUE, "0.6.0", "additive", "Upper 95% credible limit of kappa."),
    r(F, "reliability", "double precision", FALSE, TRUE, "0.6.0", "additive", "Overall reliability tau2 / (tau2 + mean v_meas) over quantified rows."),
    r(F, "reliability_lo", "double precision", FALSE, TRUE, "0.6.0", "additive", "Lower 95% credible limit of overall reliability."),
    r(F, "reliability_hi", "double precision", FALSE, TRUE, "0.6.0", "additive", "Upper 95% credible limit of overall reliability."),
    r(F, "fix_beta1", "boolean", FALSE, TRUE, "0.6.0", "additive", "TRUE if beta1 was fixed at 1 rather than estimated."),
    r(F, "n_censored", "integer", FALSE, TRUE, "0.6.0", "additive", "Censored observations entering the likelihood."),
    r(F, "identifiability", "text", FALSE, TRUE, "0.6.0", "additive", "Whether the signal/measurement split is set by the data or by the prior.")
  )
}


#' Data contract for the precision-weight tables
#'
#' Describes `madi_results.calib_weights` (one row per sample and method) and
#' `madi_results.calib_weights_fit` (one row per multiplate group and method):
#' column names, types, keys, which model each column applies to, and a short
#' plain-language description. `R/calib_schema.sql` holds the DDL; this function
#' is the same definition in a form code can check against, and a package test
#' keeps the two in sync.
#'
#' @section One weight set per sample:
#' Both tables keep their original primary keys, so the sample hierarchy is the
#' same as `calib_samples`. There is one set of weights per (sample, method);
#' the `model` column records which model produced it, and re-running with a
#' different model replaces the set. A column that does not apply to the model
#' that wrote the row is `NULL`.
#'
#' @param table `"sample"` (`calib_weights`), `"fit"` (`calib_weights_fit`) or
#'   `"all"` (default).
#' @param model Optional `"multiplicative"` or `"additive"`: keep only columns
#'   that apply to that model (plus those that apply to both).
#' @return A data frame with columns `table`, `column`, `type`, `key`,
#'   `nullable`, `added_in`, `models`, `description`.
#' @seealso [validate_weights_table()]
#' @examples
#' weights_contract("fit", model = "additive")[, c("column", "description")]
#' @export
weights_contract <- function(table = c("all", "sample", "fit"), model = NULL) {
  table <- match.arg(table)
  x <- .weights_contract_rows()
  if (table != "all") x <- x[x$table == table, , drop = FALSE]
  if (!is.null(model)) {
    model <- match.arg(model, .WEIGHTS_MODELS)
    x <- x[x$models %in% c("both", model), , drop = FALSE]
  }
  rownames(x) <- NULL
  x
}


#' Check rows against the weights contract before writing them
#'
#' Verifies that a data frame destined for `calib_weights` or
#' `calib_weights_fit` has a valid `model`, contains no column the contract does
#' not know, has every key and non-nullable column filled, and leaves empty the
#' columns that belong to the other model. Stops with a message listing the
#' problems, so a bad batch is caught before it reaches the database.
#'
#' @param x Data frame of rows to be written.
#' @param table `"sample"` or `"fit"`.
#' @param model `"multiplicative"` or `"additive"`: the model that produced the
#'   rows. If `x` has a `model` column it must agree.
#' @return `x`, invisibly, if it conforms.
#' @seealso [weights_contract()]
#' @export
validate_weights_table <- function(x, table = c("sample", "fit"),
                                   model = c("multiplicative", "additive")) {
  table <- match.arg(table); model <- match.arg(model)
  ct <- weights_contract(table)
  problems <- character()

  unknown <- setdiff(names(x), ct$column)
  if (length(unknown))
    problems <- c(problems, paste0("unknown column(s): ", paste(unknown, collapse = ", ")))

  if ("model" %in% names(x)) {
    bad <- setdiff(unique(as.character(x$model)), model)
    if (length(bad))
      problems <- c(problems, paste0("model column holds ", paste(bad, collapse = ", "),
                                     " but rows are declared ", model))
  }

  need <- ct$column[(ct$key | !ct$nullable) & ct$column != "model"]
  miss <- setdiff(need, names(x))
  if (length(miss))
    problems <- c(problems, paste0("missing required column(s): ", paste(miss, collapse = ", ")))
  for (cn in intersect(need, names(x)))
    if (anyNA(x[[cn]]))
      problems <- c(problems, paste0("required column '", cn, "' has missing values"))

  other <- setdiff(.WEIGHTS_MODELS, model)
  foreign <- ct$column[ct$models == other]
  for (cn in intersect(foreign, names(x)))
    if (!all(is.na(x[[cn]])))
      problems <- c(problems, paste0("column '", cn, "' belongs to the ", other,
                                     " model but is filled in ", model, " rows"))

  if (length(problems))
    stop("validate_weights_table(", table, ", ", model, "): ",
         paste(problems, collapse = "; "), call. = FALSE)
  invisible(x)
}
