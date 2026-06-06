# curveRcore 0.1.0

* Initial release.
* Five canonical forward models: `logistic4`, `logistic5`, `loglogistic4`,
  `loglogistic5`, `gompertz4`.
* Shared `calibration_result` S3 class.
* Eligibility gating via `assess_model_eligibility()` and
  `select_best_eligible()`.
* Full detection/quantification limit suite (LOD, LLOQ, ULOQ).
