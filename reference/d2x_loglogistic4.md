# Second Derivative of the loglogistic4 Model (Numerical)

Second Derivative of the loglogistic4 Model (Numerical)

## Usage

``` r
d2x_loglogistic4(x, a, b, c, d, h = 1e-05)
```

## Arguments

- x:

  Numeric vector. Independent variable (typically log10-concentration).

- a:

  Numeric scalar. Lower asymptote (baseline response).

- b:

  Numeric scalar. Scale parameter (\\b \> 0\\); controls steepness.

- c:

  Numeric scalar. Inflection-point location on the x-axis.

- d:

  Numeric scalar. Upper asymptote (saturation response).

- h:

  Step size. Default `1e-5`.

## Value

Numeric vector of d²y/dx² values.
