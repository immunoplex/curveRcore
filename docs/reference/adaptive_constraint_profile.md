# Build an Adaptive Constraint Profile from Observed Data

Inspects the response range, dynamic range, and scale to choose
appropriate bounds for nonlinear optimisation. The returned profile is
consumed by per-model constraint builders in curveRfreq and curveRbayes.

## Usage

``` r
adaptive_constraint_profile(
  data,
  response_variable,
  is_log_response,
  antigen_settings
)
```

## Arguments

- data:

  Data frame. Must contain response and `concentration` cols.

- response_variable:

  Character. Response column name.

- is_log_response:

  Logical. Is the response already log10-transformed?

- antigen_settings:

  List with `l_asy_min_constraint` and `l_asy_max_constraint`.

## Value

A named list: `y_min`, `y_max`, `dynamic_range`, `conc_range`,
`scale_class`, `slope_min`, `slope_max`, `g_min`, `g_max`,
`conc_pad_frac`, `d_margin_frac`.

## Details

Three scale classes are recognised:

- high:

  MFI-like (log-max \> 2.5 or raw max \> 1000)

- medium:

  Intermediate signals

- low:

  OD/absorbance-like (narrow dynamic range)

Narrower dynamic ranges receive wider slope and asymmetry bounds to
avoid near-singular Jacobians.
