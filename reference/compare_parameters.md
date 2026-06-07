# Compare Parameters Between Two Calibration Results

Produces a side-by-side parameter comparison table for the best model
from each result.

## Usage

``` r
compare_parameters(result_a, result_b, plate_id_b = 1L)
```

## Arguments

- result_a, result_b:

  `calibration_result` objects.

- plate_id_b:

  Integer. For Bayesian results with multiple plates, which plate's
  parameters to extract. Default 1.

## Value

Data frame with columns: `term`, `a_estimate`, `b_estimate`, `a_se`,
`b_se`, `diff`, `rel_diff`.
