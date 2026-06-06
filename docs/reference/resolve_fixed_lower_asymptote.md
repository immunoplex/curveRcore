# Resolve the Fixed Lower Asymptote Value

Determines whether the lower asymptote `a` should be fixed based on the
constraint method. Returns the fixed value (on the **raw** scale, before
any log-transform) or NULL if `a` should be estimated freely.

## Usage

``` r
resolve_fixed_lower_asymptote(l_asy_constraints)
```

## Arguments

- l_asy_constraints:

  Named list. Must contain `l_asy_constraint_method` and, for
  non-default methods, `l_asy_min_constraint` and/or
  `l_asy_max_constraint`. Typically an `antigen_constraints` object or a
  plain list.

## Value

Numeric scalar (the fixed value on the raw scale) or NULL (parameter is
free).

## Details

The constraint method controls behaviour:

- `"default"`:

  `a` is always free (returns NULL).

- `"user_defined"`:

  `a` is fixed at `l_asy_min_constraint` (which should equal
  `l_asy_max_constraint`).

- `"range_of_blanks"`:

  `a` is fixed at `l_asy_min_constraint`. The caller is responsible for
  computing this from blank data upstream.

- `"geometric_mean_of_blanks"`:

  `a` is fixed at `l_asy_min_constraint`. The caller is responsible for
  computing this from blank data upstream.
