# Create Antigen-Level Constraint Settings

Constructs and validates the antigen-specific constraint metadata that
controls lower-asymptote handling, concentration scaling, and pcov
thresholds.

## Usage

``` r
new_antigen_constraints(
  antigen,
  l_asy_min = 0,
  l_asy_max = 0,
  l_asy_method = "default",
  std_curve_conc = 10000,
  pcov_threshold = 15,
  std_error_blank = NULL
)
```

## Arguments

- antigen:

  Character. Antigen identifier.

- l_asy_min:

  Numeric. Lower bound for the `a` parameter.

- l_asy_max:

  Numeric. Upper bound for the `a` parameter.

- l_asy_method:

  Character. One of `"default"`, `"user_defined"`, `"range_of_blanks"`,
  `"geometric_mean_of_blanks"`.

- std_curve_conc:

  Numeric. Undiluted standard concentration (e.g. 10000).

- pcov_threshold:

  Numeric. Percent CV acceptance threshold (e.g. 15 or 20).

- std_error_blank:

  Numeric or NULL. SE of blank wells.

## Value

A named list of class `antigen_constraints`.
