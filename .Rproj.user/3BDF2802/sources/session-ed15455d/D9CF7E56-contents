# =============================================================================
# formulas.R — Build NLS formula objects from curveRcore model definitions
#
# Canonical formula generation for all five models. Used by curveRfreq
# and curveRbayes. When fixed_a is supplied, parameter 'a' is baked into
# the formula as a numeric constant.
# =============================================================================


#' Build NLS Formulas for Candidate Models
#'
#' Constructs `stats::nls()`-compatible formula objects for the requested
#' models, using the curveRcore parameterisation (b > 0, always increasing).
#'
#' When `fixed_a` is non-NULL and `is_log_response = TRUE`, the raw value
#' is validated and log10-transformed before substitution. When `fixed_a`
#' is non-NULL and `is_log_response = FALSE`, it is substituted directly.
#'
#' @param model_names Character vector. Which models to build formulas for.
#'   Must be a subset of [available_models()].
#' @param response_variable Character. Name of the response column in the
#'   fitting data frame.
#' @param independent_variable Character. Name of the concentration column
#'   in the fitting data frame. Default `"concentration"`.
#' @param fixed_a Numeric scalar or NULL. If non-NULL, the lower asymptote
#'   is fixed at this value (on the **raw** scale before any log transform).
#' @param is_log_response Logical. If TRUE and `fixed_a` is non-NULL,
#'   log10-transform the fixed value before substitution.
#'
#' @return Named list of formula objects, one per model.
#'
#' @export
build_nls_formulas <- function(model_names,
                               response_variable,
                               independent_variable = "concentration",
                               fixed_a = NULL,
                               is_log_response = TRUE) {

  valid <- available_models()
  bad <- setdiff(model_names, valid)
  if (length(bad) > 0)
    stop("Unknown model(s): ", paste(bad, collapse = ", "))

  # Prepare the 'a' token for formula strings
  a_tok <- "a"
  if (!is.null(fixed_a)) {
    if (is_log_response) {
      val <- validate_fixed_lower_asymptote(fixed_a, verbose = FALSE)
      if (!is.null(val)) {
        a_tok <- sprintf("(%.10f)", log10(val + 1e-6))
      }
    } else {
      a_tok <- sprintf("(%.10f)", fixed_a)
    }
  }

  rv <- response_variable
  x  <- independent_variable

  templates <- list(
    logistic4    = sprintf("%s ~ %s + (d - %s) / I(1 + exp(-((%s - c) / b)))",
                           rv, a_tok, a_tok, x),
    logistic5    = sprintf("%s ~ %s + (d - %s) / I((1 + exp(-((%s - c) / b)))^g)",
                           rv, a_tok, a_tok, x),
    loglogistic4 = sprintf("%s ~ %s + (d - %s) / I(1 + (c / %s)^b)",
                           rv, a_tok, a_tok, x),
    loglogistic5 = sprintf("%s ~ %s + (d - %s) * I((1 + g * exp(-b * (%s - c)))^(-1/g))",
                           rv, a_tok, a_tok, x),
    gompertz4    = sprintf("%s ~ %s + (d - %s) * I(exp(-exp(-b * (%s - c))))",
                           rv, a_tok, a_tok, x)
  )

  formulas <- lapply(templates[model_names], function(s) stats::as.formula(s))
  names(formulas) <- model_names
  formulas
}
