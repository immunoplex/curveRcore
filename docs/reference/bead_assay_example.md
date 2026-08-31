# Bead-based immunoassay example dataset

A named list of simulated multi-plate bead-based immunoassay data
spanning 6 plates across 2 analytes (alpha and beta), measured on 3
replicate plates each. Contains 20 patient samples per plate across 3
timepoints (baseline, month3, month6) and 2 treatment groups.

A named list of simulated multi-plate bead-based immunoassay data
spanning 6 plates across 2 analytes (alpha and beta), measured on 3
replicate plates each. Contains 20 patient samples per plate across 3
timepoints (baseline, month3, month6) and 2 treatment groups.

## Usage

``` r
bead_assay_example

bead_assay_example
```

## Format

A named list with six elements:

- standards:

  Data frame (60 rows x 8 cols). Standard curve data with columns:
  `curve_id`, `stype`, `sampleid`, `well`, `dilution`, `mfi`,
  `assay_response_variable`, `assay_independent_variable`.

- blanks:

  Data frame (24 rows x 7 cols). Blank well measurements (4 per plate).

- samples:

  Data frame (120 rows x 13 cols). Patient samples with columns
  including `timeperiod`, `patientid`, `agroup`, `dilution`, `mfi`.

- curve_id_lookup:

  Data frame (6 rows x 5 cols). Maps `curve_id` to `antigen`,
  `study_accession`, `experiment_accession`, `plate`.

- response_var:

  Character: `"mfi"`.

- indep_var:

  Character: `"concentration"`.

A named list with six elements:

- standards:

  Data frame (60 rows x 8 cols). Standard curve data with columns:
  `curve_id`, `stype`, `sampleid`, `well`, `dilution`, `mfi`,
  `assay_response_variable`, `assay_independent_variable`.

- blanks:

  Data frame (24 rows x 7 cols). Blank well measurements (4 per plate).

- samples:

  Data frame (120 rows x 13 cols). Patient samples with columns
  including `timeperiod`, `patientid`, `agroup`, `dilution`, `mfi`.

- curve_id_lookup:

  Data frame (6 rows x 5 cols). Maps `curve_id` to `antigen`,
  `study_accession`, `experiment_accession`, `plate`.

- response_var:

  Character: `"mfi"`.

- indep_var:

  Character: `"concentration"`.

## Source

Synthetic data. Standard curves simulated using fitted parameters from
real Luminex data with realistic plate-to-plate variability and
heteroscedastic noise. See `data-raw/migrate_datasets.R`.

Synthetic data. Standard curves simulated using fitted parameters from
real Luminex data with realistic plate-to-plate variability and
heteroscedastic noise. See `data-raw/migrate_datasets.R`.

## Examples

``` r
data(bead_assay_example)
str(bead_assay_example, max.level = 1)
#> List of 6
#>  $ standards      :'data.frame': 60 obs. of  8 variables:
#>  $ blanks         :'data.frame': 24 obs. of  7 variables:
#>  $ samples        :'data.frame': 120 obs. of  13 variables:
#>  $ curve_id_lookup:'data.frame': 6 obs. of  5 variables:
#>  $ response_var   : chr "mfi"
#>  $ indep_var      : chr "concentration"
head(bead_assay_example$standards)
#>   curve_id stype sampleid well    dilution     mfi assay_response_variable
#> 1        1     S   STD_01   A1 1000.000000   109.4                     mfi
#> 2        1     S   STD_02   B1  333.333333   316.9                     mfi
#> 3        1     S   STD_03   C1  100.000000  1133.0                     mfi
#> 4        1     S   STD_04   D1   33.333333  4156.1                     mfi
#> 5        1     S   STD_05   E1   10.000000 12458.1                     mfi
#> 6        1     S   STD_06   F1    3.333333 18933.4                     mfi
#>   assay_independent_variable
#> 1              concentration
#> 2              concentration
#> 3              concentration
#> 4              concentration
#> 5              concentration
#> 6              concentration
data(bead_assay_example)
str(bead_assay_example, max.level = 1)
#> List of 6
#>  $ standards      :'data.frame': 60 obs. of  8 variables:
#>  $ blanks         :'data.frame': 24 obs. of  7 variables:
#>  $ samples        :'data.frame': 120 obs. of  13 variables:
#>  $ curve_id_lookup:'data.frame': 6 obs. of  5 variables:
#>  $ response_var   : chr "mfi"
#>  $ indep_var      : chr "concentration"
head(bead_assay_example$standards)
#>   curve_id stype sampleid well    dilution     mfi assay_response_variable
#> 1        1     S   STD_01   A1 1000.000000   109.4                     mfi
#> 2        1     S   STD_02   B1  333.333333   316.9                     mfi
#> 3        1     S   STD_03   C1  100.000000  1133.0                     mfi
#> 4        1     S   STD_04   D1   33.333333  4156.1                     mfi
#> 5        1     S   STD_05   E1   10.000000 12458.1                     mfi
#> 6        1     S   STD_06   F1    3.333333 18933.4                     mfi
#>   assay_independent_variable
#> 1              concentration
#> 2              concentration
#> 3              concentration
#> 4              concentration
#> 5              concentration
#> 6              concentration
```
