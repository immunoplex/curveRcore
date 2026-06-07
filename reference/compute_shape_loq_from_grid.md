# Compute curvature-based (shape) LOQs from an enriched grid

Identifies the lower shape-LOQ as the concentration at the global
maximum of `d2y_dx2` and the upper shape-LOQ as the concentration at the
global minimum. Vertex positions are refined by 3-point parabolic
interpolation.

## Usage

``` r
compute_shape_loq_from_grid(grid)
```

## Arguments

- grid:

  Data frame with columns `log10_concentration`, `predicted_response`,
  and `d2y_dx2`.

## Value

A named list with elements:

- `shape_lloq_log10`:

  log10 concentration at the global maximum of d2y_dx2 (lower
  shape-LOQ).

- `shape_uloq_log10`:

  log10 concentration at the global minimum of d2y_dx2 (upper
  shape-LOQ).

- `shape_lloq_conc`:

  Natural-scale lower shape-LOQ.

- `shape_uloq_conc`:

  Natural-scale upper shape-LOQ.

- `shape_lloq_response`:

  Predicted response at shape-LLOQ.

- `shape_uloq_response`:

  Predicted response at shape-ULOQ.

## Details

Designed to be called inside
[`assess_model_eligibility()`](https://immunoplex.github.io/curveRcore/reference/assess_model_eligibility.md)
at Step 3, returning a named list that is merged into the `$eligibility`
slot.
