# Package index

## Forward Models

Five canonical sigmoid/saturation model functions.

- [`available_models()`](https://immunoplex.github.io/curveRcore/reference/available_models.md)
  : List all available model names
- [`model_params()`](https://immunoplex.github.io/curveRcore/reference/model_params.md)
  : Model registry: parameter names for each model family
- [`build_nls_formulas()`](https://immunoplex.github.io/curveRcore/reference/build_nls_formulas.md)
  : Build NLS Formulas for Candidate Models
- [`gompertz4()`](https://immunoplex.github.io/curveRcore/reference/gompertz4.md)
  : Four-Parameter Gompertz Forward Function
- [`logistic4()`](https://immunoplex.github.io/curveRcore/reference/logistic4.md)
  : Four-Parameter Logistic (4PL) Forward Function
- [`logistic5()`](https://immunoplex.github.io/curveRcore/reference/logistic5.md)
  : Five-Parameter Logistic (5PL) Forward Function
- [`loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/loglogistic4.md)
  : Four-Parameter Log-Logistic (Dose-Response) Forward Function
- [`loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/loglogistic5.md)
  : Five-Parameter Generalised Logistic (Richards) Forward Function
- [`adaptive_constraint_profile()`](https://immunoplex.github.io/curveRcore/reference/adaptive_constraint_profile.md)
  : Build an Adaptive Constraint Profile from Observed Data
- [`validate_fixed_lower_asymptote()`](https://immunoplex.github.io/curveRcore/reference/validate_fixed_lower_asymptote.md)
  : Validate a Fixed Lower Asymptote Before Log Transformation

## Settings & Configuration

Constructor functions for fit settings objects.

- [`new_antigen_constraints()`](https://immunoplex.github.io/curveRcore/reference/new_antigen_constraints.md)
  : Create Antigen-Level Constraint Settings
- [`new_calibration_result()`](https://immunoplex.github.io/curveRcore/reference/new_calibration_result.md)
  : Construct a Calibration Result Object
- [`new_calibration_result_multiplate()`](https://immunoplex.github.io/curveRcore/reference/new_calibration_result_multiplate.md)
  : Construct a Multi-Plate Calibration Result
- [`new_fit_options()`](https://immunoplex.github.io/curveRcore/reference/new_fit_options.md)
  : Create Model Fitting and Grid Options
- [`new_study_params()`](https://immunoplex.github.io/curveRcore/reference/new_study_params.md)
  : Create Study-Level Fitting Parameters
- [`resolve_effective_models()`](https://immunoplex.github.io/curveRcore/reference/resolve_effective_models.md)
  : Resolve Effective Models for a Given Configuration
- [`resolve_fixed_lower_asymptote()`](https://immunoplex.github.io/curveRcore/reference/resolve_fixed_lower_asymptote.md)
  : Resolve the Fixed Lower Asymptote Value
- [`resolve_response_col()`](https://immunoplex.github.io/curveRcore/reference/resolve_response_col.md)
  : Resolve the Response Column Name

## Calibration Result Class

The shared calibration_result S3 class.

- [`agreement_metrics()`](https://immunoplex.github.io/curveRcore/reference/agreement_metrics.md)
  : Compute Agreement Metrics Between Paired Predictions
- [`compare_calibrations()`](https://immunoplex.github.io/curveRcore/reference/compare_calibrations.md)
  : Compare Two Calibration Results (Grid Predictions)
- [`compare_parameters()`](https://immunoplex.github.io/curveRcore/reference/compare_parameters.md)
  : Compare Parameters Between Two Calibration Results
- [`compare_samples()`](https://immunoplex.github.io/curveRcore/reference/compare_samples.md)
  : Compare Sample Predictions Between Two Calibration Results
- [`new_calibration_result()`](https://immunoplex.github.io/curveRcore/reference/new_calibration_result.md)
  : Construct a Calibration Result Object
- [`new_calibration_result_multiplate()`](https://immunoplex.github.io/curveRcore/reference/new_calibration_result_multiplate.md)
  : Construct a Multi-Plate Calibration Result

## Tidy Extractors

Broom-style helpers for extracting tabular data.

- [`pcov_from_se()`](https://immunoplex.github.io/curveRcore/reference/pcov_se_conversion.md)
  [`se_from_pcov()`](https://immunoplex.github.io/curveRcore/reference/pcov_se_conversion.md)
  : Convert between posterior CV (pcov) and the log10-scale
  concentration SD
- [`tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.md)
  : Tidy the precision grid from a calibration result
- [`tidy_samples()`](https://immunoplex.github.io/curveRcore/reference/tidy_samples.md)
  : Tidy the per-sample predictions from a calibration result

## Eligibility & Selection

Model gating and best-model selection.

- [`assess_model_eligibility()`](https://immunoplex.github.io/curveRcore/reference/assess_model_eligibility.md)
  : Assess Model Eligibility for Quantification
- [`select_best_eligible()`](https://immunoplex.github.io/curveRcore/reference/select_best_eligible.md)
  : Select the Best Eligible Model

## Detection & Quantification Limits

LOD, MDC, RDL, and shape-LOQ functions.

- [`compute_detection_limits()`](https://immunoplex.github.io/curveRcore/reference/compute_detection_limits.md)
  : Compute and attach detection limits to a calibration_result
- [`compute_detection_limits_multiplate()`](https://immunoplex.github.io/curveRcore/reference/compute_detection_limits_multiplate.md)
  : Compute detection limits for all plates in a multiplate result
- [`compute_shape_loq_from_grid()`](https://immunoplex.github.io/curveRcore/reference/compute_shape_loq_from_grid.md)
  : Compute curvature-based (shape) LOQs from an enriched grid
- [`.empty_detection_limits()`](https://immunoplex.github.io/curveRcore/reference/dot-empty_detection_limits.md)
  : Empty detection_limits list for non-converged models
- [`enrich_grid_with_d2y()`](https://immunoplex.github.io/curveRcore/reference/enrich_grid_with_d2y.md)
  : Add a d2y_dx2 column to an existing prediction grid

## Grid & Predictions

Prediction grid generation and curve confidence intervals.

- [`compute_curve_ci()`](https://immunoplex.github.io/curveRcore/reference/compute_curve_ci.md)
  : Compute Confidence Interval for Fitted Curve
- [`generate_prediction_grid()`](https://immunoplex.github.io/curveRcore/reference/generate_prediction_grid.md)
  : Generate a Prediction Grid of Concentrations
- [`predict_grid_response()`](https://immunoplex.github.io/curveRcore/reference/predict_grid_response.md)
  : Compute Predicted Response for a Grid

## Transforms

Data transforms applied before fitting.

- [`preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.md)
  : Full Preprocessing Pipeline for Standard Curve Data (mask-aware)
- [`perform_blank_operation()`](https://immunoplex.github.io/curveRcore/reference/perform_blank_operation.md)
  : Apply a Blank Operation to Standard Curve Data
- [`compute_concentration()`](https://immunoplex.github.io/curveRcore/reference/compute_concentration.md)
  : Compute Concentration from Dilution and Undiluted Standard
- [`compute_log_response()`](https://immunoplex.github.io/curveRcore/reference/compute_log_response.md)
  : Log10-Transform the Assay Response (mask-aware)
- [`correct_prozone()`](https://immunoplex.github.io/curveRcore/reference/correct_prozone.md)
  : Correct the Prozone (Hook) Effect

## Inverses

Back-calculation functions for all five models.

- [`inv_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/inv_gompertz4.md)
  [`inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRcore/reference/inv_gompertz4.md)
  : Inverse of the Gompertz Model
- [`inv_logistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic4.md)
  [`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic4.md)
  : Inverse of the 4PL Model
- [`inv_logistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic5.md)
  [`inv_logistic5_fixed()`](https://immunoplex.github.io/curveRcore/reference/inv_logistic5.md)
  : Inverse of the 5PL Model
- [`inv_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic4.md)
  [`inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic4.md)
  : Inverse of the loglogistic4 Model
- [`inv_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic5.md)
  [`inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRcore/reference/inv_loglogistic5.md)
  : Inverse of the loglogistic5 Model

## Derivatives & Gradients

Analytical derivatives and gradient vectors.

- [`dydx_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/dydx_gompertz4.md)
  : First Derivative of the Gompertz Model
- [`dydx_logistic4()`](https://immunoplex.github.io/curveRcore/reference/dydx_logistic4.md)
  : First Derivative of the 4PL Model
- [`dydx_logistic5()`](https://immunoplex.github.io/curveRcore/reference/dydx_logistic5.md)
  : First Derivative of the 5PL Model
- [`dydx_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/dydx_loglogistic4.md)
  : First Derivative of the loglogistic4 Model
- [`dydx_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/dydx_loglogistic5.md)
  : First Derivative of the loglogistic5 Model
- [`d2x_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/d2x_gompertz4.md)
  : Second Derivative of the Gompertz Model
- [`d2x_logistic4()`](https://immunoplex.github.io/curveRcore/reference/d2x_logistic4.md)
  : Second Derivative of the 4PL Model
- [`d2x_logistic5()`](https://immunoplex.github.io/curveRcore/reference/d2x_logistic5.md)
  : Second Derivative of the 5PL Model (Numerical)
- [`d2x_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/d2x_loglogistic4.md)
  : Second Derivative of the loglogistic4 Model (Numerical)
- [`d2x_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/d2x_loglogistic5.md)
  : Second Derivative of the loglogistic5 Model (Numerical)
- [`grad_gompertz4()`](https://immunoplex.github.io/curveRcore/reference/grad_gompertz4.md)
  : Analytical Gradient of the Inverse Gompertz
- [`grad_logistic4()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic4.md)
  : Analytical Gradient of the Inverse 4PL
- [`grad_logistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_logistic5.md)
  : Analytical Gradient of the Inverse 5PL
- [`grad_loglogistic4()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic4.md)
  : Analytical Gradient of the Inverse loglogistic4
- [`grad_loglogistic5()`](https://immunoplex.github.io/curveRcore/reference/grad_loglogistic5.md)
  : Analytical Gradient of the Inverse loglogistic5
- [`make_inv_and_grad_fixed()`](https://immunoplex.github.io/curveRcore/reference/make_inv_and_grad_fixed.md)
  : Build Inverse, Gradient, and grad_y Closures for a Model

## Example data sets

Data from immune assays

- [`bead_assay_example`](https://immunoplex.github.io/curveRcore/reference/bead_assay_example.md)
  : Bead-based immunoassay example dataset
- [`elisa_assay_example`](https://immunoplex.github.io/curveRcore/reference/elisa_assay_example.md)
  : ELISA example dataset

## Utilities

Internal helpers and shared utilities.

- [`filter_by_curve_id()`](https://immunoplex.github.io/curveRcore/reference/filter_by_curve_id.md)
  : Filter a Dataset List by Curve ID
- [`safe_unique()`](https://immunoplex.github.io/curveRcore/reference/safe_unique.md)
  : Safe unique with NA handling
- [`geom_mean()`](https://immunoplex.github.io/curveRcore/reference/geom_mean.md)
  : Geometric mean
- [`.extract_param_ci()`](https://immunoplex.github.io/curveRcore/reference/dot-extract_param_ci.md)
  : Extract parameter CIs from an ensemble entry
- [`.fmt_na()`](https://immunoplex.github.io/curveRcore/reference/dot-fmt_na.md)
  : Format a number for messages, NA-safe
- [`.interp_response()`](https://immunoplex.github.io/curveRcore/reference/dot-interp_response.md)
  : Linear interpolation of predicted_response at a log10 concentration
- [`.parabolic_refine()`](https://immunoplex.github.io/curveRcore/reference/dot-parabolic_refine.md)
  : 3-point parabolic interpolation for a set of extremum indices
- [`.safe_invert_model()`](https://immunoplex.github.io/curveRcore/reference/dot-safe_invert_model.md)
  : Safely invert a named model at a response value
- [`.safe_pow10()`](https://immunoplex.github.io/curveRcore/reference/dot-safe_pow10.md)
  : Safe 10^x, NA-propagating
