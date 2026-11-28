# =============================================================================
# noise_model.R — response noise models shared by both engines
#
#   constant  sigma_fit = sigma                        (default; existing fits)
#   add_prop  sigma_raw(mu) = sqrt(sigma_add^2 + (cv * mu)^2)
#             raw response:  sigma_fit = sigma_raw(mu)
#             log10 response: sigma_fit = sigma_raw(mu) / (mu * ln 10),
#                             mu = 10^(fitted log10 response)  (delta method)
#
# add_prop describes immunoassay noise as an additive floor (reader / blank
# noise, wavelength correction) plus a part proportional to the signal. It is
# defined on the raw mean, so it is valid on either response scale; it replaces
# the earlier power-of-mean term exp(log_sigma0 + slope * log|mu|), which was
# undefined on a log10 response (mu crosses 0 where the response is 1).
#
# On the log path mu is the response after the offset of
# compute_response_adjustment() (y + c); c is 0 unless some value was <= 0 and
# then only 0.001 above -min, so it does not change the noise materially.
# =============================================================================

#' Response Noise Models
#'
#' The noise models the engines accept: `"constant"` (one SD on the fitting
#' scale; the default) and `"add_prop"` (additive + proportional,
#' \eqn{\sigma^2 = \sigma_{add}^2 + (CV \mu)^2} on the raw response).
#' @export
NOISE_MODELS <- c("constant", "add_prop")

#' Validate a Noise Model Name
#' @param noise_model Character scalar or NULL (NULL = `"constant"`).
#' @return The validated name.
#' @export
check_noise_model <- function(noise_model) {
  nm <- noise_model %||% "constant"
  if (!is.character(nm) || length(nm) != 1L || !(nm %in% NOISE_MODELS))
    stop("noise_model must be one of: ", paste(NOISE_MODELS, collapse = ", "))
  nm
}


#' Observation-Noise SD on the Fitting Scale
#'
#' @param mu Numeric. Mean response on the fitting scale (log10 when
#'   `is_log_response`).
#' @param noise_model `"constant"` or `"add_prop"`.
#' @param sigma Numeric scalar. Constant-model SD (fitting scale).
#' @param sigma_add,cv Numeric scalars. Additive SD (raw response units) and
#'   coefficient of variation (fraction) of the `add_prop` model.
#' @param is_log_response Logical.
#'
#' @return Numeric vector, one SD per `mu`, on the fitting scale.
#' @export
noise_sd <- function(mu, noise_model = "constant", sigma = NA_real_,
                     sigma_add = NA_real_, cv = NA_real_,
                     is_log_response = TRUE) {
  noise_model <- check_noise_model(noise_model)
  if (noise_model == "constant") return(rep(sigma, length(mu)))
  m_raw <- if (isTRUE(is_log_response)) 10^mu else mu
  s_raw <- sqrt(sigma_add^2 + (cv * m_raw)^2)
  if (isTRUE(is_log_response)) s_raw / (m_raw * log(10)) else s_raw
}


#' Fit the Additive + Proportional Noise Parameters to Residuals
#'
#' Maximum-likelihood estimate of `sigma_add` and `cv` from residuals and fitted
#' means of one fit group (Gaussian on the fitting scale). Used by the
#' frequentist engine to build weights.
#'
#' @param resid Numeric. Residuals on the fitting scale.
#' @param mu Numeric. Fitted means on the fitting scale.
#' @param is_log_response Logical.
#' @param start Optional named numeric `c(sigma_add, cv)`.
#'
#' @return Named list `sigma_add`, `cv`, `loglik`, `converged`.
#' @export
fit_add_prop_noise <- function(resid, mu, is_log_response = TRUE, start = NULL) {
  ok <- is.finite(resid) & is.finite(mu)
  resid <- resid[ok]; mu <- mu[ok]
  if (length(resid) < 3L) stop("fit_add_prop_noise(): need at least 3 residuals")
  m_raw <- if (isTRUE(is_log_response)) 10^mu else abs(mu)
  # Starting values: additive part from the low end, CV from the high end.
  s0 <- start %||% {
    raw_res <- if (isTRUE(is_log_response)) resid * m_raw * log(10) else resid
    lo <- m_raw <= stats::quantile(m_raw, 0.25)
    hi <- m_raw >= stats::quantile(m_raw, 0.75)
    c(sigma_add = max(stats::sd(raw_res[lo]), 1e-8, na.rm = TRUE),
      cv        = max(stats::sd(raw_res[hi] / m_raw[hi]), 1e-4, na.rm = TRUE))
  }
  nll <- function(lp) {
    s <- noise_sd(mu, "add_prop", sigma_add = exp(lp[1]), cv = exp(lp[2]),
                  is_log_response = is_log_response)
    if (any(!is.finite(s) | s <= 0)) return(1e100)
    -sum(stats::dnorm(resid, 0, s, log = TRUE))
  }
  o <- stats::optim(log(pmax(s0, 1e-10)), nll, method = "Nelder-Mead",
                    control = list(maxit = 2000, reltol = 1e-10))
  list(sigma_add = unname(exp(o$par[1])), cv = unname(exp(o$par[2])),
       loglik = -o$value, converged = o$convergence == 0L)
}
