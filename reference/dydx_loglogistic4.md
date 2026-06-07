# First Derivative of the loglogistic4 Model

\$\$\frac{dy}{dx} = \frac{b\\(d - a)\\r}{x\\(1 + r)^2} \quad\text{where
} r = (c/x)^b\$\$

## Usage

``` r
dydx_loglogistic4(x, a, b, c, d)
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

Numeric vector of dy/dx values.

## See also

Other derivatives:
[`dydx_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/dydx_gompertz4.md),
[`dydx_logistic4()`](https://immunoplex.github.io/curveRcore/reference/dydx_logistic4.md),
[`dydx_logistic5()`](https://immunoplex.github.io/curveRcore/reference/dydx_logistic5.md),
[`dydx_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/dydx_loglogistic5.md)
