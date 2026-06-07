# Construct a Multi-Plate Calibration Result

Wraps multiple single-plate `calibration_result` objects into a
multi-plate container.

## Usage

``` r
new_calibration_result_multiplate(meta, plates)
```

## Arguments

- meta:

  Named list. Multi-plate metadata (must include `plates` character
  vector).

- plates:

  Named list of `calibration_result` objects, one per plate.

## Value

An object of class `calibration_result_multiplate`.
