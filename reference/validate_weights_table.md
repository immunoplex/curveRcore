# Check rows against the weights contract before writing them

Verifies that a data frame destined for `calib_weights` or
`calib_weights_fit` has a valid `model`, contains no column the contract
does not know, has every key and non-nullable column filled, and leaves
empty the columns that belong to the other model. Stops with a message
listing the problems, so a bad batch is caught before it reaches the
database.

## Usage

``` r
validate_weights_table(
  x,
  table = c("sample", "fit"),
  model = c("multiplicative", "additive")
)
```

## Arguments

- x:

  Data frame of rows to be written.

- table:

  `"sample"` or `"fit"`.

- model:

  `"multiplicative"` or `"additive"`: the model that produced the rows.
  If `x` has a `model` column it must agree.

## Value

`x`, invisibly, if it conforms.

## See also

[`weights_contract()`](https://immunoplex.github.io/curveRcore/reference/weights_contract.md)
