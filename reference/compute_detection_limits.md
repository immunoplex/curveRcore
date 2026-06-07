# Compute and attach detection limits to a calibration_result

Computes LODs, MDC, and RDL for the best eligible model (or a specified
model) and returns the calibration_result with a new `$detection_limits`
element.

## Usage

``` r
compute_detection_limits(cr, model_name = NULL, alpha = 0.05, verbose = FALSE)
```

## Arguments

- cr:

  A `calibration_result` object.

- model_name:

  Character; model to use. If `NULL` (default), uses
  `cr$selection$best_model_name`.

- alpha:

  Numeric; significance level for CI-based LODs (default 0.05).

- verbose:

  Logical; emit diagnostic messages (default `FALSE`).

## Value

The input `cr` with `$detection_limits` populated.

## Details

This can be called:

- Automatically inside the fitting pipeline (after Step 4), or

- By the user post-hoc on an existing `calibration_result`.
