# Second Derivative of the Gompertz Model

\$\$\frac{d^2y}{dx^2} = b^2\\(d - a)\\v\\(v - 1)\\\exp(-v)
\quad\text{where } v = \exp(-b(x - c))\$\$

## Usage

``` r
d2x_gompertz4(x, a, b, c, d)
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
