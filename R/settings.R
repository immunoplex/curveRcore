# =============================================================================
# settings.R — Settings schemas shared across curveRfreq and curveRbayes
#
# Three constructor/validator pairs:
#   new_antigen_constraints()  — per-antigen constraint metadata
#   new_study_params()         — study-level fitting options
#   new_fit_options()          — model selection and grid options
# =============================================================================


#' Create Antigen-Level Constraint Settings
#'
#' Constructs and validates the antigen-specific constraint metadata that
#' controls lower-asymptote handling, concentration scaling, and pcov
#' thresholds.
#'
#' @param antigen Character. Antigen identifier.
#' @param l_asy_min Numeric. Lower bound for the `a` parameter.
#' @param l_asy_max Numeric. Upper bound for the `a` parameter.
#' @param l_asy_method Character. One of `"default"`, `"user_defined"`,
#'   `"range_of_blanks"`, `"geometric_mean_of_blanks"`.
#' @param std_curve_conc Numeric. Undiluted standard concentration
#'   (e.g. 10000).
#' @param pcov_threshold Numeric. Percent CV acceptance threshold
#'   (e.g. 15 or 20).
#' @param std_error_blank Numeric or NULL. SE of blank wells.
#'
#' @return A named list of class `antigen_constraints`.
#'
#' @export
new_antigen_constraints <- function(antigen,
                                    l_asy_min = 0,
                                    l_asy_max = 0,
                                    l_asy_method = "default",
                                    std_curve_conc = 10000,
                                    pcov_threshold = 15,
                                    std_error_blank = NULL) {

  valid_methods <- c("default", "user_defined", "range_of_blanks",
                     "geometric_mean_of_blanks")
  if (!(l_asy_method %in% valid_methods)) {
    stop("l_asy_method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  stopifnot(is.numeric(std_curve_conc), std_curve_conc > 0)
  stopifnot(is.numeric(pcov_threshold), pcov_threshold > 0)

  out <- list(
    antigen                      = antigen,
    l_asy_min_constraint         = l_asy_min,
    l_asy_max_constraint         = l_asy_max,
    l_asy_constraint_method      = l_asy_method,
    standard_curve_concentration = std_curve_conc,
    pcov_threshold               = pcov_threshold,
    std_error_blank              = std_error_blank
  )
  class(out) <- c("antigen_constraints", "list")
  out
}


#' Create Study-Level Fitting Parameters
#'
#' Constructs the study-wide settings that control data preprocessing.
#'
#' @param is_log_response Logical. Log10-transform the assay response?
#'   Default TRUE.
#' @param is_log_independent Logical. Log10-transform the concentration?
#'   Default TRUE.
#' @param apply_prozone Logical. Apply prozone (hook effect) correction?
#'   Default TRUE.
#' @param blank_option Character. Blank handling method. One of
#'   `"ignored"`, `"included"`, `"subtracted"`, `"subtracted_3x"`,
#'   `"subtracted_10x"`. Default `"ignored"`.
#' @param persist_draws Logical. If TRUE, the fitting engine persists the
#'   full posterior (Bayesian) or asymptotic-MVN (frequentist) parameter
#'   draws to `calib_draws`. Heavy output; default FALSE. Gate only — the
#'   light `calib_hyperparam`/`calib_fit_diag` tables are always written.
#' @param bayes_single_plate Logical. Bayesian only. If TRUE, each plate is
#'   fit independently (N_plates = 1, no cross-plate pooling) instead of as
#'   one multiplate hierarchy. Default FALSE (multiplate pooling). The
#'   frequentist engine is always per-plate and ignores this flag.
#'
#' @return A named list of class `study_params`.
#'
#' @export
new_study_params <- function(is_log_response    = TRUE,
                             is_log_independent = TRUE,
                             apply_prozone      = TRUE,
                             blank_option       = "ignored",
                             persist_draws      = FALSE,
                             bayes_single_plate = FALSE) {

  valid_blanks <- c("ignored", "included", "subtracted",
                    "subtracted_3x", "subtracted_10x")
  if (!(blank_option %in% valid_blanks)) {
    stop("blank_option must be one of: ", paste(valid_blanks, collapse = ", "))
  }
  stopifnot(is.logical(is_log_response), is.logical(is_log_independent),
            is.logical(apply_prozone),
            is.logical(persist_draws),      length(persist_draws) == 1L,
            is.logical(bayes_single_plate), length(bayes_single_plate) == 1L)

  out <- list(
    is_log_response    = is_log_response,
    is_log_independent = is_log_independent,
    apply_prozone      = apply_prozone,
    blank_option       = blank_option,
    persist_draws      = persist_draws,
    bayes_single_plate = bayes_single_plate
  )
  class(out) <- c("study_params", "list")
  out
}


#' Create Model Fitting and Grid Options
#'
#' Constructs the options that control which models are fit, grid
#' resolution, and pcov capping.
#'
#' @param model_names Character vector. Which models to fit. Default:
#'   all five canonical models.
#' @param n_grid Integer. Number of points in the prediction grid.
#'   Default 200.
#' @param cv_x_max Numeric. Cap for percent CV of predicted concentration.
#'   Default 150.
#' @param grid_min_conc Numeric. Minimum concentration for the grid
#'   (on the raw scale). Default 1e-4.
#' @param grid_max_conc Numeric or NULL. Maximum concentration. NULL
#'   uses the undiluted standard concentration from antigen_constraints.
#'
#' @return A named list of class `fit_options`.
#'
#' @export
new_fit_options <- function(model_names = available_models(),
                            n_grid      = 200L,
                            cv_x_max    = 150,
                            grid_min_conc = 1e-4,
                            grid_max_conc = NULL) {

  valid <- available_models()
  bad <- setdiff(model_names, valid)
  if (length(bad) > 0) {
    stop("Unknown model(s): ", paste(bad, collapse = ", "),
         ". Must be among: ", paste(valid, collapse = ", "))
  }
  stopifnot(is.numeric(n_grid), n_grid >= 10)
  stopifnot(is.numeric(cv_x_max), cv_x_max > 0)
  stopifnot(is.numeric(grid_min_conc), grid_min_conc > 0)

  out <- list(
    model_names   = model_names,
    n_grid        = as.integer(n_grid),
    cv_x_max      = cv_x_max,
    grid_min_conc = grid_min_conc,
    grid_max_conc = grid_max_conc
  )
  class(out) <- c("fit_options", "list")
  out
}


#' Resolve Effective Models for a Given Configuration
#'
#' Filters the requested model list based on whether the independent
#' variable is log-transformed (which makes loglogistic4 redundant
#' with logistic4).
#'
#' @param fit_options A `fit_options` object.
#' @param study_params A `study_params` object.
#'
#' @return Character vector of model names to actually fit.
#' @export
resolve_effective_models <- function(fit_options, study_params) {
  models <- fit_options$model_names
  if (study_params$is_log_independent) {
    # On log scale, loglogistic4 ≡ logistic4 (see PARAMETERIZATION.md)
    # Keep only one to avoid redundant fits
    models <- setdiff(models, "loglogistic4")
  }
  models
}


#' Resolve the Response Column Name
#'
#' Auto-detects `"mfi"` (bead arrays) or `"absorbance"` (ELISA) from
#' the column names of `df`. Falls back to `default`.
#'
#' @param df Data frame.
#' @param default Character. Fallback column name.
#'
#' @return Character scalar: the detected response column name.
#' @export
resolve_response_col <- function(df, default = "mfi") {
  candidates <- c("mfi", "absorbance")
  present <- intersect(candidates, names(df))
  if (length(present) >= 1) present[1] else default
}
