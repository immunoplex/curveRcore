# First Derivative of the 4PL Model

\$\$\frac{dy}{dx} = \frac{(d - a)\\u}{b\\(1 + u)^2} \quad\text{where } u
= \exp\\\left(-\frac{x - c}{b}\right)\$\$

## Usage

``` r
dydx_logistic4(x, a, b, c, d)
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
[`dydx_logistic5()`](https://immunoplex.github.io/curveRcore/reference/dydx_logistic5.md),
[`dydx_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/dydx_loglogistic4.md),
[`dydx_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/dydx_loglogistic5.md)
