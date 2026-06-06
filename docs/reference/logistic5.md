# Five-Parameter Logistic (5PL) Forward Function

Extends
[`logistic4()`](https://immunoplex.github.io/curveRcore/reference/logistic4.md)
with an asymmetry parameter \\g\\: \$\$y = a + \frac{d - a}{\bigl(1 +
\exp\\\bigl(-\frac{x - c}{b}\bigr)\bigr)^g}\$\$

## Usage

``` r
logistic5(x, a, b, c, d, g)
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

Numeric vector of predicted response values.

## Details

When \\g = 1\\ this reduces to
[`logistic4()`](https://immunoplex.github.io/curveRcore/reference/logistic4.md).
\\g \> 1\\ skews toward the upper asymptote; \\0 \< g \< 1\\ skews
toward the lower asymptote.

## See also

Other forward-models:
[`gompertz4()`](https://immunoplex.github.io/curveRcore/reference/gompertz4.md),
[`logistic4()`](https://immunoplex.github.io/curveRcore/reference/logistic4.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/loglogistic4.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/loglogistic5.md)
