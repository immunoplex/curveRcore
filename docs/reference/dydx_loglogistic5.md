# First Derivative of the loglogistic5 Model

\$\$\frac{dy}{dx} = b\\(d - a)\\\exp(-b(x-c))\\ u^{-1/g - 1}
\quad\text{where } u = 1 + g\\\exp(-b(x-c))\$\$

## Usage

``` r
dydx_loglogistic5(x, a, b, c, d, g)
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

  Numeric scalar. Asymmetry (Richards) parameter (\\g \> 0\\).

## Value

Numeric vector of dy/dx values.

## See also

Other derivatives:
[`dydx_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/dydx_gompertz4.md),
[`dydx_logistic4()`](https://immunoplex.github.io/curveRcore/reference/dydx_logistic4.md),
[`dydx_logistic5()`](https://immunoplex.github.io/curveRcore/reference/dydx_logistic5.md),
[`dydx_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/dydx_loglogistic4.md)
