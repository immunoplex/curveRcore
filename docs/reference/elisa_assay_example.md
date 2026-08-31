# ELISA example dataset

A named list of simulated multi-plate ELISA data for a single analyte
(alpha), spanning 6 plates. Plates 1-3 use a 5-parameter logistic
standard curve; plates 4-6 use Gompertz, reflecting realistic
between-plate variability in curve shape.

A named list of simulated multi-plate ELISA data for a single analyte
(alpha), spanning 6 plates. Plates 1-3 use a 5-parameter logistic
standard curve; plates 4-6 use Gompertz, reflecting realistic
between-plate variability in curve shape.

## Usage

``` r
elisa_assay_example

elisa_assay_example
```

## Format

A named list with six elements:

- standards:

  Data frame (60 rows x 8 cols). Standard curve data with columns:
  `curve_id`, `stype`, `sampleid`, `well`, `dilution`, `od`,
  `assay_response_variable`, `assay_independent_variable`.

- blanks:

  Data frame (24 rows x 7 cols). Blank well measurements.

- samples:

  Data frame (120 rows x 12 cols). Patient samples at fixed serum
  dilution 1:400.

- curve_id_lookup:

  Data frame (6 rows x 5 cols). Maps `curve_id` to `antigen`,
  `study_accession`, `experiment_accession`, `plate`.

- response_var:

  Character: `"od"`.

- indep_var:

  Character: `"concentration"`.

A named list with six elements:

- standards:

  Data frame (60 rows x 8 cols). Standard curve data with columns:
  `curve_id`, `stype`, `sampleid`, `well`, `dilution`, `od`,
  `assay_response_variable`, `assay_independent_variable`.

- blanks:

  Data frame (24 rows x 7 cols). Blank well measurements.

- samples:

  Data frame (120 rows x 12 cols). Patient samples at fixed serum
  dilution 1:400.

- curve_id_lookup:

  Data frame (6 rows x 5 cols). Maps `curve_id` to `antigen`,
  `study_accession`, `experiment_accession`, `plate`.

- response_var:

  Character: `"od"`.

- indep_var:

  Character: `"concentration"`.

## Source

Synthetic data with biologically plausible parameters and proportional
CV matching typical ELISA plate-reader behaviour.

Synthetic data with biologically plausible parameters and proportional
CV matching typical ELISA plate-reader behaviour.

## Examples

``` r
data(elisa_assay_example)
str(elisa_assay_example, max.level = 1)
#> List of 6
#>  $ standards      :'data.frame': 60 obs. of  8 variables:
#>  $ blanks         :'data.frame': 24 obs. of  7 variables:
#>  $ samples        :'data.frame': 120 obs. of  12 variables:
#>  $ curve_id_lookup:'data.frame': 6 obs. of  5 variables:
#>  $ response_var   : chr "od"
#>  $ indep_var      : chr "concentration"
head(elisa_assay_example$standards)
#>   curve_id stype sampleid well    dilution     od assay_response_variable
#> 1        1     S   STD_01   A1 1000.000000 0.0955                      od
#> 2        1     S   STD_02   B1  333.333333 0.0865                      od
#> 3        1     S   STD_03   C1  100.000000 0.0672                      od
#> 4        1     S   STD_04   D1   33.333333 0.0937                      od
#> 5        1     S   STD_05   E1   10.000000 0.1885                      od
#> 6        1     S   STD_06   F1    3.333333 0.7303                      od
#>   assay_independent_variable
#> 1              concentration
#> 2              concentration
#> 3              concentration
#> 4              concentration
#> 5              concentration
#> 6              concentration
data(elisa_assay_example)
str(elisa_assay_example, max.level = 1)
#> List of 6
#>  $ standards      :'data.frame': 60 obs. of  8 variables:
#>  $ blanks         :'data.frame': 24 obs. of  7 variables:
#>  $ samples        :'data.frame': 120 obs. of  12 variables:
#>  $ curve_id_lookup:'data.frame': 6 obs. of  5 variables:
#>  $ response_var   : chr "od"
#>  $ indep_var      : chr "concentration"
head(elisa_assay_example$standards)
#>   curve_id stype sampleid well    dilution     od assay_response_variable
#> 1        1     S   STD_01   A1 1000.000000 0.0955                      od
#> 2        1     S   STD_02   B1  333.333333 0.0865                      od
#> 3        1     S   STD_03   C1  100.000000 0.0672                      od
#> 4        1     S   STD_04   D1   33.333333 0.0937                      od
#> 5        1     S   STD_05   E1   10.000000 0.1885                      od
#> 6        1     S   STD_06   F1    3.333333 0.7303                      od
#>   assay_independent_variable
#> 1              concentration
#> 2              concentration
#> 3              concentration
#> 4              concentration
#> 5              concentration
#> 6              concentration
```
