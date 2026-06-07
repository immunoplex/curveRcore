# Four-Parameter Gompertz Forward Function

Computes the response for a Gompertz saturation curve: \$\$y = a + (d -
a)\\\exp\\\bigl(-\exp(-b\\(x - c))\bigr)\$\$

## Usage

``` r
gompertz4(x, a, b, c, d)
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

The Gompertz is intrinsically asymmetric (skewed toward the upper
asymptote) without requiring a fifth parameter.

## See also

Other forward-models:
[`logistic4()`](https://immunoplex.github.io/curveRcore/reference/logistic4.md),
[`logistic5()`](https://immunoplex.github.io/curveRcore/reference/logistic5.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/loglogistic4.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/loglogistic5.md)
