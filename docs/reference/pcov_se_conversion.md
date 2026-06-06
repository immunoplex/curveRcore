# Convert between posterior CV (pcov) and the log10-scale concentration SD

In the curveR ecosystem the back-calculated concentration is reported on
the log10 scale (when `is_log_independent = TRUE`). `se_concentration`
is the delta-method standard deviation of that log10 concentration, and
`pcov` is the percent coefficient of variation derived from it and then
capped:

## Usage

``` r
pcov_from_se(se)

se_from_pcov(pcov)
```

## Arguments

- se:

  Numeric vector of `se_concentration` values (log10-scale SD).

- pcov:

  Numeric vector of `pcov` values (percent).

## Value

A numeric vector of the same length.

## Details

\$\$\mathrm{pcov} = \mathrm{se\\concentration} \times \ln(10) \times
100, \quad \text{then capped at } cv\\x\\max.\$\$

`se_concentration` is therefore the *uncapped* modelling-scale SD;
`pcov` is a censored percent. Downstream variance/weight models should
consume `se_concentration`, not `pcov` (the cap destroys the precision
gradient).

These helpers are the single canonical implementation of the
relationship. Do not reimplement it elsewhere.
