# First Derivative of the 5PL Model

\$\$\frac{dy}{dx} = \frac{g\\(d - a)\\u}{b\\(1 + u)^{g+1}}\$\$

## Usage

``` r
dydx_logistic5(x, a, b, c, d, g)
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

- g:

  Numeric scalar. Asymmetry parameter (\\g \> 0\\).

## Value

Numeric vector of dy/dx values.

## See also

Other derivatives:
[`dydx_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/dydx_gompertz4.md),
[`dydx_logistic4()`](https://immunoplex.github.io/curveRcore/reference/dydx_logistic4.md),
[`dydx_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/dydx_loglogistic4.md),
[`dydx_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/dydx_loglogistic5.md)
