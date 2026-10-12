# =============================================================================
# censored_draws.R -- predictive concentration draws that keep their
# out-of-range information, and the per-sample summaries built from them.
#
# Both engines back-calculate a sample many times (Bayesian: once per
# posterior draw; frequentist: once per Monte Carlo parameter draw). A draw
# whose response lies beyond one of that draw's asymptotes has no finite
# concentration. The legacy summaries (se_concentration / pcov) drop those
# draws, which is deliberately left unchanged (pcov is frozen, CONTRIBUTING
# §3.3). The functions here instead keep such a draw as -Inf (beyond the
# low-concentration asymptote: concentration -> 0) or +Inf (beyond the
# high-concentration asymptote: concentration -> unbounded), so quantiles and
# probabilities computed from the full draw set treat it as censored rather
# than silently conditioning on invertibility.
# =============================================================================


#' Invert a response with out-of-range responses kept as censored draws
#'
#' Wraps a model inverse so that a response beyond the low-concentration
#' asymptote returns `-Inf` (on the log10 scale: concentration 0) and one
#' beyond the high-concentration asymptote returns `+Inf`, instead of the
#' `NA`/`NaN`/accidental `Inf` the bare inverses produce. Which of `a`/`d` is the
#' low-concentration asymptote is probed from the inverse itself, so the
#' wrapper is correct for decreasing curves and any sign of `b`.
#'
#' @param inv_fn Function of one argument `y` (vectorised) returning `x` for a
#'   single parameter set, e.g. `function(y) fns$inv(y, theta)`.
#' @param y Numeric vector of responses on the fitting scale.
#' @param a,d Numeric scalars: the two asymptotes for this parameter set.
#' @param tol Numeric scalar: a response within `tol` (as a fraction of
#'   `d - a`) of an asymptote counts as beyond it. Default `1e-9`.
#'
#' @return Numeric vector the length of `y`: finite `x`, `-Inf`, `+Inf`, or
#'   `NA` (missing `y`, or unusable asymptotes).
#' @keywords internal
#' @export
invert_with_bounds <- function(inv_fn, y, a, d, tol = 1e-9) {
  out <- rep(NA_real_, length(y))
  if (length(a) != 1L || length(d) != 1L || !is.finite(a) || !is.finite(d) || a == d)
    return(out)

  t <- (y - a) / (d - a)                       # 0 at a, 1 at d
  x_q1 <- suppressWarnings(tryCatch(inv_fn(a + 0.25 * (d - a)), error = function(e) NA_real_))
  x_q3 <- suppressWarnings(tryCatch(inv_fn(a + 0.75 * (d - a)), error = function(e) NA_real_))
  a_is_low <- if (is.finite(x_q1) && is.finite(x_q3)) x_q1 < x_q3 else TRUE
  low_val  <- if (a_is_low) -Inf else Inf

  ok       <- !is.na(t)
  beyond_a <- ok & t <= tol
  beyond_d <- ok & t >= 1 - tol
  inside   <- ok & !beyond_a & !beyond_d

  if (any(inside))
    out[inside] <- suppressWarnings(tryCatch(inv_fn(y[inside]),
                                             error = function(e) rep(NA_real_, sum(inside))))
  out[beyond_a] <- low_val
  out[beyond_d] <- -low_val

  # Numerically non-finite inside the open interval (very close to an
  # asymptote): assign to the nearer asymptote rather than dropping it.
  bad <- inside & !is.finite(out)
  if (any(bad)) out[bad] <- ifelse(t[bad] < 0.5, low_val, -low_val)
  out
}


#' Summarise predictive concentration draws for test samples
#'
#' Turns a samples x draws matrix of back-calculated log10 concentrations,
#' where draws beyond an asymptote are `-Inf`/`+Inf` (see
#' [invert_with_bounds()]), into the per-sample interval, threshold-probability
#' and censored-draw columns shared by both engines.
#'
#' Quantiles use `type = 1` (inverse empirical CDF), which is well defined when
#' some draws are infinite. Concentrations are reported on the natural scale and
#' multiplied by each sample's dilution, the same scale as `final_concentration`.
#' `p_below_lloq` / `p_above_uloq` compare pre-dilution log10 concentration with
#' the curve's LLOQ/ULOQ (the same footing as [classify_pcov_gate()]);
#' `p_below_cutoff` compares the dilution-corrected concentration with
#' `decision_cutoff`.
#'
#' @param x_log10 Numeric matrix, one row per sample, one column per draw:
#'   pre-dilution log10 concentration, `-Inf`/`+Inf` for censored draws, `NA`
#'   for draws that could not be evaluated at all (ignored).
#' @param dilution Numeric vector, one per row (`NA` gives `NA` natural-scale
#'   columns and `p_below_cutoff`).
#' @param lloq_log10,uloq_log10 Numeric scalars (log10, pre-dilution) or `NA`.
#' @param decision_cutoff Numeric scalar on the dilution-corrected natural scale,
#'   or `NULL`/`NA` for none.
#' @param probs Numeric length-3 vector of quantile levels. Default
#'   `c(0.025, 0.5, 0.975)`.
#'
#' @return A data.frame with one row per sample and columns `conc_q_lo`,
#'   `conc_q_med`, `conc_q_hi` (an unbounded upper quantile is `NA`; a lower
#'   quantile at the low asymptote is `0`), `p_below_lloq`, `p_above_uloq`,
#'   `p_below_cutoff`, `frac_draws_below_a` (share of draws beyond the
#'   low-concentration asymptote) and `frac_draws_above_d` (share beyond the
#'   high-concentration asymptote).
#'
#' @seealso [classify_censoring()]
#' @export
summarize_conc_draws <- function(x_log10, dilution, lloq_log10 = NA_real_,
                                 uloq_log10 = NA_real_, decision_cutoff = NULL,
                                 probs = c(0.025, 0.5, 0.975)) {
  if (!is.matrix(x_log10)) x_log10 <- matrix(x_log10, nrow = 1L)
  n <- nrow(x_log10)
  if (length(dilution) == 1L) dilution <- rep(dilution, n)
  if (length(dilution) != n)
    stop("summarize_conc_draws(): dilution must have one value per row of x_log10.")
  if (length(probs) != 3L) stop("summarize_conc_draws(): probs must have length 3.")

  dil <- suppressWarnings(as.numeric(dilution))
  dil[!is.finite(dil) | dil <= 0] <- NA_real_
  have_lloq   <- length(lloq_log10) == 1L && is.finite(lloq_log10)
  have_uloq   <- length(uloq_log10) == 1L && is.finite(uloq_log10)
  have_cutoff <- length(decision_cutoff) == 1L && is.finite(decision_cutoff) &&
                 decision_cutoff > 0

  out <- data.frame(
    conc_q_lo = rep(NA_real_, n), conc_q_med = NA_real_, conc_q_hi = NA_real_,
    p_below_lloq = NA_real_, p_above_uloq = NA_real_, p_below_cutoff = NA_real_,
    frac_draws_below_a = NA_real_, frac_draws_above_d = NA_real_)

  for (i in seq_len(n)) {
    x <- x_log10[i, ]
    x <- x[!is.na(x)]                         # NaN is.na too
    if (!length(x)) next

    out$frac_draws_below_a[i] <- mean(x == -Inf)
    out$frac_draws_above_d[i] <- mean(x ==  Inf)
    if (have_lloq) out$p_below_lloq[i] <- mean(x < lloq_log10)
    if (have_uloq) out$p_above_uloq[i] <- mean(x > uloq_log10)

    if (!is.na(dil[i])) {
      q  <- stats::quantile(x, probs = probs, type = 1, names = FALSE)
      qn <- 10^q * dil[i]                      # -Inf -> 0, +Inf -> Inf
      qn[!is.finite(qn)] <- NA_real_
      out$conc_q_lo[i]  <- qn[1]
      out$conc_q_med[i] <- qn[2]
      out$conc_q_hi[i]  <- qn[3]
      if (have_cutoff)
        out$p_below_cutoff[i] <- mean(x < log10(decision_cutoff / dil[i]))
    }
  }
  out
}
