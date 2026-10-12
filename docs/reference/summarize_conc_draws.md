# Summarise predictive concentration draws for test samples

Turns a samples x draws matrix of back-calculated log10 concentrations,
where draws beyond an asymptote are `-Inf`/`+Inf` (see
[`invert_with_bounds()`](https://immunoplex.github.io/curveRcore/reference/invert_with_bounds.md)),
into the per-sample interval, threshold-probability and censored-draw
columns shared by both engines.

## Usage

``` r
summarize_conc_draws(
  x_log10,
  dilution,
  lloq_log10 = NA_real_,
  uloq_log10 = NA_real_,
  decision_cutoff = NULL,
  probs = c(0.025, 0.5, 0.975)
)
```

## Arguments

- x_log10:

  Numeric matrix, one row per sample, one column per draw: pre-dilution
  log10 concentration, `-Inf`/`+Inf` for censored draws, `NA` for draws
  that could not be evaluated at all (ignored).

- dilution:

  Numeric vector, one per row (`NA` gives `NA` natural-scale columns and
  `p_below_cutoff`).

- lloq_log10, uloq_log10:

  Numeric scalars (log10, pre-dilution) or `NA`.

- decision_cutoff:

  Numeric scalar on the dilution-corrected natural scale, or `NULL`/`NA`
  for none.

- probs:

  Numeric length-3 vector of quantile levels. Default
  `c(0.025, 0.5, 0.975)`.

## Value

A data.frame with one row per sample and columns `conc_q_lo`,
`conc_q_med`, `conc_q_hi` (an unbounded upper quantile is `NA`; a lower
quantile at the low asymptote is `0`), `p_below_lloq`, `p_above_uloq`,
`p_below_cutoff`, `frac_draws_below_a` (share of draws beyond the
low-concentration asymptote) and `frac_draws_above_d` (share beyond the
high-concentration asymptote).

## Details

Quantiles use `type = 1` (inverse empirical CDF), which is well defined
when some draws are infinite. Concentrations are reported on the natural
scale and multiplied by each sample's dilution, the same scale as
`final_concentration`. `p_below_lloq` / `p_above_uloq` compare
pre-dilution log10 concentration with the curve's LLOQ/ULOQ (the same
footing as
[`classify_pcov_gate()`](https://immunoplex.github.io/curveRcore/reference/classify_pcov_gate.md));
`p_below_cutoff` compares the dilution-corrected concentration with
`decision_cutoff`.

## See also

[`classify_censoring()`](https://immunoplex.github.io/curveRcore/reference/classify_censoring.md)
