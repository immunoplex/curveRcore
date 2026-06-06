# Second Derivative of the 4PL Model

\$\$\frac{d^2y}{dx^2} = -\frac{(d - a)\\u\\(1 - u)}{b^2\\(1 + u)^3}\$\$

## Usage

``` r
d2x_logistic4(x, a, b, c, d)
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

## Value

Numeric vector of d²y/dx² values.
