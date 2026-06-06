Foundation package for the curveR ecosystem. Provides the five canonical
forward model functions (logistic4, logistic5, loglogistic4,
loglogistic5, gompertz4), their analytical inverses, first and second
derivatives, gradient closures for delta-method error propagation, data
transformation utilities, settings schemas, and the shared
calibration_result S3 output class. Both curveRfreq (frequentist) and
curveRbayes (Bayesian) depend on this package for model math, ensuring
exact parameterisation parity across statistical frameworks.
