# Back-transform Fitted Responses to the Raw Scale

Inverse of
[`transform_sample_response()`](https://immunoplex.github.io/curveRcore/reference/transform_sample_response.md)
for model-scale responses (e.g. predicted curves, asymptotes):
`10^y - offset + blank_value` on the log path, `y + blank_value`
otherwise.

## Usage

``` r
untransform_response(y, is_log_response, adjustment = NULL)
```

## Arguments

- y:

  Numeric. Responses on the fitting scale.

- is_log_response:

  Logical.

- adjustment:

  List or NULL (NULL = no subtraction, no offset).

## Value

Numeric vector on the raw response scale.
