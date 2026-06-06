# Filter a Dataset List by Curve ID

Subsets all data-frame elements of a dataset list (standards, blanks,
samples, curve_id_lookup) to rows matching a specific `curve_id`.
Non-data-frame elements and data frames without a `curve_id` column are
passed through unchanged.

## Usage

``` r
filter_by_curve_id(
  loaded_data,
  curve_id,
  target_names = c("standards", "blanks", "samples", "curve_id_lookup"),
  verbose = FALSE
)
```

## Arguments

- loaded_data:

  A named list as returned by `data(bead_assay_example)` or
  `data(elisa_assay_example)`. Must contain data frames with a
  `curve_id` column.

- curve_id:

  Scalar. The `curve_id` value to filter on.

- target_names:

  Character vector. Which list elements to filter. Default:
  `c("standards", "blanks", "samples", "curve_id_lookup")`.

- verbose:

  Logical. Emit messages about skipped elements.

## Value

A copy of `loaded_data` with the target data frames filtered to the
specified `curve_id`, plus `$whole_standards` and
`$curve_id_whole_lookup` preserving the full unfiltered versions.

## Details

Also attaches `$whole_standards` and `$curve_id_whole_lookup` as
unfiltered copies for downstream cross-plate operations.
