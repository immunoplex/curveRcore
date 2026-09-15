# Reclassify pcov_pass (grid + samples) and add pcov_gate_class (samples), for every curve in a fitted result, using each curve's own LLOQ/ULOQ

Call this AFTER fit_calibration_freq()/fit_calibration_bayes() and AFTER
compute_detection_limits_multiplate() – typically the very next line in
the worker pipeline. Overwrites `cr$grid$pcov_pass` and
`cr$samples$pcov_pass` for every curve, and adds
`cr$samples$ pcov_gate_class`. `calib_grid` intentionally gets NO
gate_class column (grid pcov_pass is recomputed for display symmetry
only; it plays no role in deriving LLOQ/ULOQ, which happens upstream
from the raw pcov profile before this function ever runs, so there is no
circularity).

## Usage

``` r
classify_pcov_gates_multiplate(mp)
```

## Arguments

- mp:

  A `calibration_result` or `calibration_result_multiplate`.

## Value

The same class of object, with `pcov_pass`/`pcov_gate_class`
(re)computed in place.

## Details

Also computes each curve's inflection point exactly once and caches it
at `cr$detection_limits$inflection` (list(x=, y=)), so
`flatten_and_save.R`'s diagnostics-row construction can read the SAME
cached value instead of recomputing it independently – guarantees
`calib_diagnostics.inflect_x/y` can never drift from whatever this
function used for gating.
