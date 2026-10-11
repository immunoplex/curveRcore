# Invert a response with out-of-range responses kept as censored draws

Wraps a model inverse so that a response beyond the low-concentration
asymptote returns `-Inf` (on the log10 scale: concentration 0) and one
beyond the high-concentration asymptote returns `+Inf`, instead of the
`NA`/`NaN`/accidental `Inf` the bare inverses produce. Which of `a`/`d`
is the low-concentration asymptote is probed from the inverse itself, so
the wrapper is correct for decreasing curves and any sign of `b`.

## Usage

``` r
invert_with_bounds(inv_fn, y, a, d, tol = 1e-09)
```

## Arguments

- inv_fn:

  Function of one argument `y` (vectorised) returning `x` for a single
  parameter set, e.g. `function(y) fns$inv(y, theta)`.

- y:

  Numeric vector of responses on the fitting scale.

- a, d:

  Numeric scalars: the two asymptotes for this parameter set.

- tol:

  Numeric scalar: a response within `tol` (as a fraction of `d - a`) of
  an asymptote counts as beyond it. Default `1e-9`.

## Value

Numeric vector the length of `y`: finite `x`, `-Inf`, `+Inf`, or `NA`
(missing `y`, or unusable asymptotes).
