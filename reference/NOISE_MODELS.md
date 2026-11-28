# Response Noise Models

The noise models the engines accept: `"constant"` (one SD on the fitting
scale; the default) and `"add_prop"` (additive + proportional,
\\\sigma^2 = \sigma\_{add}^2 + (CV \mu)^2\\ on the raw response).

## Usage

``` r
NOISE_MODELS
```

## Format

An object of class `character` of length 2.
