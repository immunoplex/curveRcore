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

- verbose:

  Logical.

## Value

`stdframe` with post-peak response values adjusted.
