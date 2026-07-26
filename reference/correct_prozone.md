# Correct the Prozone (Hook) Effect

At very high concentrations the measured signal can decrease (hook
effect). This function compresses the post-peak delta toward the peak
value.

## Usage

``` r
correct_prozone(
  stdframe,
  prop_diff = 0.1,
  dil_scale = 2,
  response_variable = "mfi",
  independent_variable = "concentration",
  include_col = "included",
  verbose = FALSE
)
```

## Arguments

- stdframe:

  Data frame of standard curve data.

- prop_diff:

  Numeric. Dampening factor (e.g. 0.1).

- dil_scale:

  Numeric. Dilution scale factor (e.g. 2).

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

- include_col:

  Character. Logical column marking fitted rows. Absent = all rows
  included (backward compatible).

- verbose:

  Logical.

## Value

`stdframe` (all rows) with post-peak response values adjusted. The peak
reference is attached as `attr(., "prozone_peak_response")` and
`attr(., "prozone_logc_at_peak")`.

## Details

The peak (`max_response` and `logc_at_max`) is a **set-level statistic**
computed from the *included* points only. The post-peak reflection is
then applied to **all** rows relative to that peak, so masked points
beyond the hook are dampened onto the same reference as the fitted
points. Rows are never dropped (grain is preserved); rows with a missing
response or concentration are passed through untouched.
