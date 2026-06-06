# Four-Parameter Log-Logistic (Dose-Response) Forward Function

Classical Hill-equation parameterisation for positive-valued x: \$\$y =
a + \frac{d - a}{1 + (c / x)^b}\$\$

## Usage

``` r
loglogistic4(x, a, b, c, d)
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

Numeric vector of predicted response values.

## Details

Here \\c\\ is the EC50 (inflection on the concentration scale) and \\b
\> 0\\ is the Hill coefficient. Requires \\x \> 0\\ and \\c \> 0\\.

## Domain restriction

This model requires positive `x`. When the independent variable is
log10-transformed (and can be negative), use
[`logistic4()`](https://immunoplex.github.io/curveRcore/reference/logistic4.md)
instead — on the log scale the two models are algebraically equivalent.

## See also

Other forward-models:
[`gompertz4()`](https://immunoplex.github.io/curveRcore/reference/gompertz4.md),
[`logistic4()`](https://immunoplex.github.io/curveRcore/reference/logistic4.md),
[`logistic5()`](https://immunoplex.github.io/curveRcore/reference/logistic5.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/loglogistic5.md)
