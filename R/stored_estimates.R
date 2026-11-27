# =============================================================================
# stored_estimates.R -- back-calculate points (controls, standards) from a fit
# that is already stored, without refitting.
#
# Used once, to fill calib_controls and the 0.7.0 estimate columns of
# calib_standards for curves fitted before 0.7.0, so readers (I-SPI, exports)
# need only the stored values and no second, in-app computation. New fits get
# the same columns from the engines directly.
#
# Inputs are what the database holds for one curve and method:
#   calib_param       best-model parameter estimates
#   calib_grid        the precision profile (se_concentration along x)
#   calib_diagnostics LLOQ/ULOQ, LODs, inflection
#   calib_draws       per-draw parameters, when the fit persisted them
#
# How each estimate is rebuilt (and how close it is to the engine's own value):
#   point estimate  frequentist: inverse at the parameter estimates (identical).
#                   bayesian: median over draws of the inverted response when
#                   draws are stored (the engine's definition), else inverse at
#                   the stored estimates.
#   se / pcov       bayesian with draws: SD of the inverted draws with the
#                   fitted measurement noise (`noise`), the engine's definition.
#                   Otherwise the stored precision profile at the estimate,
#                   whenever the point lies within it (frequentist: the engine's
#                   delta-method value to grid resolution, verified within
#                   0.1%); outside the profile, frequentist delta method with
#                   the draws' covariance plus the residual noise (sigma_resid,
#                   stored in calib_param). Otherwise unknown (pcov at the cap,
#                   as the engines do for a point they cannot invert).
#   gate, censoring classify_pcov_gate() / classify_censoring() with the stored
#                   limits: the same functions the pipeline uses.
#   quantiles       summarize_conc_draws() over the stored draws with the same
#                   noise model (fresh noise draws, so equal in distribution,
#                   not draw for draw); NULL without draws.
# =============================================================================


#' Back-calculate points from a stored fit
#'
#' Rebuilds, for a set of points (control wells or standard points) read on one
#' stored curve, every estimate column of the 0.7.0 contract
#' ([point_estimates_contract()]) from the stored best-model parameters,
#' precision profile, limits and, when present, parameter draws. Nothing is
#' refitted. Intended for the one-off backfill of curves fitted before 0.7.0; see
#' `calib_points_SPEC.md` for how close each column is to the engine's own value.
#'
#' @param response Numeric vector: each point's response on the fitting scale
#'   (for standards, `calib_standards.response_model`; for controls, the raw
#'   response transformed as the fit's `is_log_response` says).
#' @param dilution Numeric vector or scalar: each point's dilution factor.
#' @param model_name Best model: `logistic4`, `logistic5`, `loglogistic4`,
#'   `loglogistic5` or `gompertz4`.
#' @param params Named numeric vector of the best model's parameter estimates
#'   (`a`, `b`, `c`, `d`, and `g` for five-parameter models). `a` may be absent
#'   when the lower asymptote was fixed; give it as `fixed_a`.
#' @param method `"frequentist"` or `"bayesian"`.
#' @param grid Optional data frame: the stored precision profile
#'   (`calib_grid` rows of the curve and method) with `log10_concentration` (or,
#'   when `is_log_independent = FALSE`, `concentration`) and `se_concentration`.
#' @param limits Named list from the curve's `calib_diagnostics` row:
#'   `lloq_log10`, `uloq_log10`, `lower_lod_response`, `upper_lod_response`,
#'   `lower_lod_log10_conc`, `upper_lod_log10_conc`, `inflect_x`. Missing entries
#'   count as unknown.
#' @param draws Optional named list (or data frame) of parameter draws, one
#'   element per term (`calib_draws` rows with `param_scope = 'curve'`). A
#'   length-1 element is a fixed parameter and is recycled.
#' @param fixed_a Numeric scalar, the fixed lower asymptote, or `NULL`.
#' @param noise Optional list for Bayesian fits, per draw and aligned with
#'   `draws`: `sigma_obs`, `nu`, and for the additive + proportional model
#'   `sigma_add`, `cv` (with `noise_model = "add_prop"`; see [noise_sd()]). When given, each draw's
#'   response is perturbed by the fitted observation noise before inversion, as
#'   the engine does with `include_measurement_error = TRUE`.
#' @param is_log_independent Logical: is x (concentration) on the log10 scale in
#'   the fit?
#' @param is_log_response Logical: is the response on the log10 scale in the fit?
#'   Only used by the additive + proportional noise model (frequentist terms
#'   `sigma_add`, `cv`, or Bayesian `noise$noise_model = "add_prop"`).
#' @param nominal Optional numeric vector: the points' nominal concentration on
#'   the fit's x scale (standards). When given, `recovery_pct` and a bias-aware
#'   `pcov_rmse` are added via [standard_recovery()].
#' @param cv_x_max,pcov_threshold As in the engines (defaults 150 and 20).
#' @param decision_cutoff Decision cutoff on the dilution-corrected natural
#'   scale, or `NULL`.
#' @param seed Integer seed for the measurement-noise draws (only used with
#'   `noise`). The caller's random-number state is left unchanged.
#'
#' @return A data frame with one row per point and the estimate columns of
#'   [point_estimates_contract()] (`recovery_pct` only when `nominal` is given).
#'   Attribute `se_source` records, per point, where se came from: `"draws"`
#'   (Bayesian draws, or frequentist draws outside the profile), `"profile"`
#'   (the stored precision profile) or `"none"`.
#' @seealso [point_estimates_contract()], [standard_recovery()]
#' @export
estimate_points_from_stored <- function(response, dilution = 1, model_name, params,
                                        method = c("frequentist", "bayesian"),
                                        grid = NULL, limits = list(), draws = NULL,
                                        fixed_a = NULL, noise = NULL,
                                        is_log_independent = TRUE, nominal = NULL,
                                        is_log_response = TRUE,
                                        cv_x_max = 150, pcov_threshold = 20,
                                        decision_cutoff = NULL, seed = 1L) {
  method <- match.arg(method)
  n <- length(response)
  response <- suppressWarnings(as.numeric(response))
  if (length(dilution) == 1L) dilution <- rep(dilution, n)
  dilution <- suppressWarnings(as.numeric(dilution))
  if (length(dilution) != n) stop("estimate_points_from_stored(): dilution must have one value per point.")
  if (!is.null(nominal) && length(nominal) != n)
    stop("estimate_points_from_stored(): nominal must have one value per point.")
  lim <- function(nm) {
    v <- suppressWarnings(as.numeric(limits[[nm]] %||% NA_real_))
    if (length(v) == 1L) v else NA_real_
  }
  params <- unlist(params)
  storage.mode(params) <- "double"
  # Frequentist residual SE (calib_param term 'sigma_resid'): the response noise
  # the engine adds to every back-calculation. Not a curve parameter.
  # With the additive + proportional model (terms sigma_add, cv) the SD depends
  # on the response: one SD per point, from the observed response as the proxy
  # for its mean (as the engines do for samples).
  se_response <- rep(0, n)
  if ("sigma_resid" %in% names(params)) {
    if (method == "frequentist" && is.finite(params[["sigma_resid"]]))
      se_response <- rep(params[["sigma_resid"]], n)
    params <- params[names(params) != "sigma_resid"]
  }
  if (all(c("sigma_add", "cv") %in% names(params))) {
    if (method == "frequentist" && all(is.finite(params[c("sigma_add", "cv")])))
      se_response <- noise_sd(response, "add_prop", sigma_add = params[["sigma_add"]],
                              cv = params[["cv"]], is_log_response = is_log_response)
    params <- params[!names(params) %in% c("sigma_add", "cv")]
  }
  se_response[!is.finite(se_response)] <- 0
  use_fixed <- !is.null(fixed_a) && length(fixed_a) == 1L && is.finite(fixed_a)
  if (use_fixed) params <- params[names(params) != "a"]
  fns <- tryCatch(make_inv_and_grad_fixed(model_name, fixed_a = if (use_fixed) fixed_a else NULL),
                  error = function(e) NULL)
  if (is.null(fns)) stop("estimate_points_from_stored(): unknown model '", model_name, "'.")
  inv1 <- function(y, p) suppressWarnings(tryCatch(fns$inv(y, p), error = function(e) NA_real_))
  a_of <- function(p) if (use_fixed) fixed_a else (p[["a"]] %||% NA_real_)

  # ---- draws: a draws x terms matrix ----------------------------------------
  P <- NULL
  if (!is.null(draws) && length(draws)) {
    dl <- lapply(as.list(draws), function(v) suppressWarnings(as.numeric(unlist(v))))
    S  <- max(lengths(dl))
    if (S > 1L) {
      dl <- lapply(dl, function(v) if (length(v) == 1L) rep(v, S) else v)
      if (all(lengths(dl) == S)) {
        P <- do.call(cbind, dl)
        P <- P[, !colnames(P) %in% c("sigma_resid", "sigma_add", "cv"), drop = FALSE]
        if (use_fixed) P <- P[, colnames(P) != "a", drop = FALSE]
      }
    }
  }
  # response used by each draw: the observed one perturbed by the fitted noise
  # (frequentist: N(0, sigma_resid), as the engine's intervals; bayesian: the
  # observation model, when `noise` is given)
  Y <- NULL
  if (!is.null(P)) {
    Y <- matrix(response, nrow = n, ncol = nrow(P))
    if (method == "frequentist" && any(se_response > 0)) {
      Y <- .with_local_seed(seed, Y + matrix(stats::rnorm(n * nrow(P), 0, rep(se_response, nrow(P))), nrow = n))
    } else if (method == "bayesian" && !is.null(noise) && length(noise$sigma_obs) == nrow(P)) {
      addp <- identical(noise$noise_model, "add_prop") && !is.null(noise$sigma_add) && !is.null(noise$cv)
      Y <- .with_local_seed(seed, {
        for (s in seq_len(nrow(P))) {
          sig <- if (addp) noise_sd(response, "add_prop", sigma_add = noise$sigma_add[s],
                                    cv = noise$cv[s], is_log_response = is_log_response)
                 else noise$sigma_obs[s]
          Y[, s] <- response + sig * stats::rt(n, df = noise$nu[s])
        }
        Y
      })
    }
  }

  # ---- point estimate ----------------------------------------------------------
  x_mat <- NULL
  if (!is.null(P) && method == "bayesian") {
    x_mat <- matrix(NA_real_, nrow = n, ncol = nrow(P))
    for (s in seq_len(nrow(P))) {
      p <- P[s, ]
      for (i in seq_len(n)) if (is.finite(Y[i, s])) x_mat[i, s] <- inv1(Y[i, s], p)
    }
    x_mat[!is.finite(x_mat)] <- NA_real_
  }
  x_est <- vapply(seq_len(n), function(i) if (is.finite(response[i])) inv1(response[i], params) else NA_real_,
                  numeric(1))
  if (method == "bayesian" && !is.null(x_mat))
    x_est <- suppressWarnings(apply(x_mat, 1, stats::median, na.rm = TRUE))
  x_est[!is.finite(x_est)] <- NA_real_

  # ---- se ------------------------------------------------------------------------
  # bayesian with draws: SD of the inverted draws (the engine's definition).
  # Otherwise, in order:
  # 1. the stored precision profile, for points inside its range (frequentist:
  #    the engine's delta-method value to grid resolution, verified within 0.1%);
  # 2. the stored draws: frequentist delta method with the draws' covariance
  #    plus the response noise sigma_resid (the engine's formula, up to Monte
  #    Carlo error);
  # 3. unknown (pcov at the cap).
  se <- rep(NA_real_, n); se_source <- rep("none", n)
  if (method == "bayesian" && !is.null(x_mat)) {
    # The engine's own definition: SD of the (noise-perturbed) inverted draws.
    se <- suppressWarnings(apply(x_mat, 1, stats::sd, na.rm = TRUE))
    se_source[is.finite(se)] <- "draws"
  }
  if (any(!is.finite(se)) && is.data.frame(grid) && nrow(grid)) {
    fill <- !is.finite(se)
    gx <- suppressWarnings(as.numeric(if (isTRUE(is_log_independent)) grid$log10_concentration
                                      else grid$concentration %||% grid$log10_concentration))
    gs <- suppressWarnings(as.numeric(grid$se_concentration))
    ok <- is.finite(gx) & is.finite(gs)
    if (sum(ok) >= 2L) {
      o <- order(gx[ok])
      prof <- stats::approx(gx[ok][o], gs[ok][o], xout = x_est, rule = 1, ties = mean)$y
      se[fill] <- prof[fill]
      se_source[fill & is.finite(se)] <- "profile"
    }
  }
  todo <- which(!is.finite(se) & is.finite(x_est))
  if (length(todo) && !is.null(P) && method == "frequentist") {
    {
      terms <- intersect(colnames(P), names(params))
      Sigma <- if (length(terms) > 1L) stats::cov(P[, terms, drop = FALSE]) else NULL
      if (!is.null(Sigma)) for (i in todo) {
        g <- tryCatch(fns$grad(response[i], params), error = function(e) NULL)
        if (is.null(g) || !all(is.finite(g))) next
        cm <- intersect(names(g), terms)
        v <- as.numeric(t(g[cm]) %*% Sigma[cm, cm, drop = FALSE] %*% g[cm])
        gy <- tryCatch(fns$grad_y(response[i], params), error = function(e) NA_real_)
        if (se_response[i] > 0 && is.finite(gy)) v <- v + gy^2 * se_response[i]^2
        if (is.finite(v) && v >= 0) se[i] <- sqrt(v)
      }
    }
    se_source[todo[is.finite(se[todo])]] <- "draws"
  }
  se[!is.finite(x_est)] <- NA_real_
  se_source[!is.finite(se)] <- "none"

  pcov <- vapply(seq_len(n), function(i) {
    if (!is.finite(se[i])) return(cv_x_max)
    raw <- if (isTRUE(is_log_independent)) se[i] * log(10) * 100
           else if (abs(x_est[i]) > 1e-10) se[i] / abs(x_est[i]) * 100 else Inf
    min(raw, cv_x_max)
  }, numeric(1))

  x_log10 <- if (isTRUE(is_log_independent)) x_est
             else ifelse(is.finite(x_est) & x_est > 0, log10(pmax(x_est, 1e-300)), NA_real_)
  final <- if (isTRUE(is_log_independent)) 10^x_est * dilution else x_est * dilution

  out <- data.frame(predicted_concentration = x_est,
                    predicted_log10_concentration = x_log10,
                    final_concentration = final,
                    se_concentration = se, pcov = pcov, pcov_rmse = pcov)

  gate <- classify_pcov_gate(x_log10, pcov, pcov_threshold,
                             lloq_log10 = lim("lloq_log10"), uloq_log10 = lim("uloq_log10"),
                             inflect_x = lim("inflect_x"))
  out$pcov_pass <- gate$pcov_pass
  out$pcov_gate_class <- as.character(gate$pcov_gate_class)

  cens <- classify_censoring(response, x_log10, dilution,
                             lloq_log10 = lim("lloq_log10"), uloq_log10 = lim("uloq_log10"),
                             lod_response = lim("lower_lod_response"),
                             ulod_response = lim("upper_lod_response"),
                             lod_log10 = lim("lower_lod_log10_conc"),
                             ulod_log10 = lim("upper_lod_log10_conc"),
                             a = a_of(params), d = params[["d"]] %||% NA_real_)
  out <- cbind(out, cens)

  if (!is.null(P)) {
    X <- matrix(NA_real_, nrow = n, ncol = nrow(P))
    for (s in seq_len(nrow(P))) {
      p <- P[s, ]
      X[, s] <- invert_with_bounds(function(yy) fns$inv(yy, p), Y[, s], a_of(p), p[["d"]])
    }
    X[!is.finite(response), ] <- NA_real_
    if (!isTRUE(is_log_independent)) {
      fin <- is.finite(X)
      X[fin] <- ifelse(X[fin] > 0, log10(X[fin]), -Inf)
    }
    summ <- summarize_conc_draws(X, dilution, lim("lloq_log10"), lim("uloq_log10"), decision_cutoff)
  } else {
    summ <- summarize_conc_draws(matrix(NA_real_, nrow = n, ncol = 1L), dilution)
  }
  out <- cbind(out, summ)

  if (!is.null(nominal)) {
    rc <- standard_recovery(x_est, nominal, se, is_log_independent, cv_x_max)
    out$pcov_rmse <- rc$pcov_rmse
    out$recovery_pct <- rc$recovery_pct
  }
  attr(out, "se_source") <- se_source
  out
}


# Evaluate `expr` with a fixed seed, leaving the caller's RNG state untouched.
.with_local_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", old, envir = globalenv())
  })
  set.seed(seed)
  expr
}
