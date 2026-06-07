# Four-Parameter Logistic (4PL) Forward Function

Computes the response for a four-parameter logistic curve: \$\$y = a +
\frac{d - a}{1 + \exp\\\left(-\frac{x - c}{b}\right)}\$\$

## Usage

``` r
logistic4(x, a, b, c, d)
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

Numeric vector of predicted response values, same length as `x`.

## Details

Symmetric about its inflection point at \\(c, (a+d)/2)\\. Always
monotonically increasing when \\b \> 0\\ and \\a \< d\\.

## See also

Other forward-models:
[`gompertz4()`](https://immunoplex.github.io/curveRcore/reference/gompertz4.md),
[`logistic5()`](https://immunoplex.github.io/curveRcore/reference/logistic5.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/loglogistic4.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/loglogistic5.md)
