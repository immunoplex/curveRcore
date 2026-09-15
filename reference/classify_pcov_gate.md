# Classify pcov gate membership for a set of points on one curve

Classify pcov gate membership for a set of points on one curve

## Usage

``` r
classify_pcov_gate(
  x_log10,
  pcov,
  pcov_threshold,
  lloq_log10 = NA_real_,
  uloq_log10 = NA_real_,
  inflect_x = NA_real_
)
```

## Arguments

- x_log10:

  Numeric vector. log10-scale position: `predicted_log10_ concentration`
  for test samples, `log10_concentration` for grid points.

- pcov:

  Numeric vector. Percent CV at each `x_log10`. Only used in the
  fallback branch (when LLOQ and/or ULOQ are undefined).

- pcov_threshold:

  Numeric scalar. Only used in the fallback branch.

- lloq_log10, uloq_log10:

  Numeric scalars (log10 scale), or `NA` if undefined for this
  curve/model.

- inflect_x:

  Numeric scalar (log10 scale). Only used in the fallback branch, to
  decide which side of the curve a failing point falls on. If `NA` and a
  point fails, that point's `pcov_gate_class` stays `NA` (we know it
  fails, we just can't say which side) while `pcov_pass` is still
  correctly `FALSE`.

## Value

A data.frame with columns:

- pcov_gate_class:

  Factor: `belowLLOQ` / `inDynamicRange` / `aboveULOQ` / `NA`. In the
  fallback branch, a failing point with no `inflect_x` available is `NA`
  here (we know it fails, not which side).

- pcov_pass:

  Logical, computed independently of `pcov_gate_class` (NOT derived from
  it): `TRUE`/`FALSE` whenever the point could be evaluated at all, `NA`
  only when `x_log10` (or, in the fallback branch, `pcov`) itself is
  unusable. A failing point with unknown side is `pcov_gate_class = NA`
  but `pcov_pass = FALSE`, never `NA`.

Row count and order match the input vectors exactly.
