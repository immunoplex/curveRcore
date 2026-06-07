# Build Inverse, Gradient, and grad_y Closures for a Model

Returns three closures (`inv`, `grad`, `grad_y`) that evaluate the
inverse function, parameter gradient, and response-derivative for a
given model. Used by error propagation via the delta method.

## Usage

``` r
make_inv_and_grad_fixed(model, fixed_a = NULL)
```

## Arguments

- model:

  Character. One of `"logistic4"`, `"logistic5"`, `"loglogistic4"`,
  `"loglogistic5"`, `"gompertz4"`.

- fixed_a:

  Numeric scalar or NULL. If non-NULL, `a` is treated as a known
  constant and excluded from the parameter gradient.

## Value

A list with closures: `inv(y, p)`, `grad(y, p)`, `grad_y(y, p)`. Each
accepts a response value `y` and named parameter vector `p` (from
`coef(fit)`).
