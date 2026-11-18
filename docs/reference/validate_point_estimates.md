# Check control or standard rows against the contract before writing them

Verifies that a data frame destined for `calib_controls` or
`calib_standards` contains no column the contract does not know, has
every key and non-nullable column filled, and that the controlled
vocabularies (`exclusion_reason`, `pcov_gate_class`, `censor_class`,
`censor_type`) hold only allowed values. Stops with a message listing
the problems, so a bad batch is caught before it reaches the database.
`created_at` is filled by the database and may be absent.

## Usage

``` r
validate_point_estimates(
  x,
  table = c("controls", "standards", "control_reference")
)
```

## Arguments

- x:

  Data frame of rows to be written. For `control_reference`, the
  generated `reference_id` and the defaulted `assigned_at`,
  `valid_from`, `study_accession`, `control_lot` and `feature` may be
  absent.

- table:

  `"controls"`, `"standards"` or `"control_reference"`.

## Value

`x`, invisibly, if it conforms.

## See also

[`point_estimates_contract()`](https://immunoplex.github.io/curveRcore/reference/point_estimates_contract.md)
