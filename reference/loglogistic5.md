# Five-Parameter Generalised Logistic (Richards) Forward Function

A five-parameter generalised logistic curve with asymmetry: \$\$y = a +
(d - a)\\\bigl(1 + g\\\exp(-b\\(x - c))\bigr)^{-1/g}\$\$

## Usage

``` r
loglogistic5(x, a, b, c, d, g)
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

Numeric vector of predicted response values.

## Details

When \\g = 1\\ this reduces to
[`logistic4()`](https://immunoplex.github.io/curveRcore/reference/logistic4.md)
(after rescaling \\b\\). Works with any real \\x\\ (including
log10-transformed).

## See also

Other forward-models:
[`gompertz4()`](https://immunoplex.github.io/curveRcore/reference/gompertz4.md),
[`logistic4()`](https://immunoplex.github.io/curveRcore/reference/logistic4.md),
[`logistic5()`](https://immunoplex.github.io/curveRcore/reference/logistic5.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/loglogistic4.md)
