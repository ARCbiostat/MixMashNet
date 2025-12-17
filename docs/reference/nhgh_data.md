# NHGH dataset

Example dataset used in MixMashNet examples. This dataset contains 15
variables derived from the National Health and Nutrition Examination
Survey (NHANES)

## Usage

``` r
data(nhgh_data)
```

## Format

A data frame with 5621 rows and 15 variables:

- wt:

  Weight (numeric).

- ht:

  Height (numeric).

- bmi:

  Body mass index (numeric).

- leg:

  Leg length (numeric).

- arml:

  Arm length (numeric).

- armc:

  Arm circumference (numeric).

- tri:

  Triceps skinfold (numeric).

- sub:

  Subscapular skinfold (numeric).

- gh:

  Glychohemoglobin (numeric).

- albumin:

  Albumin (numeric).

- bun:

  Blood urea nitrogen (numeric).

- SCr:

  Serum creatinine (numeric).

- age:

  Age (numeric).

- sex:

  Sex with 0=woman and 1=man.

- re:

  Race with 1=Mexican American, 2=Other Hispanic, 3=Non-Hispanic White,
  4=Non-Hispanic Black, 5=Other Race Including Multi-Racial.

## Examples

``` r
data(nhgh_data)
str(nhgh_data)
#> 'data.frame':    5621 obs. of  15 variables:
#>  $ wt     : num  87.4 72.3 116.8 97.6 86.7 ...
#>  $ ht     : num  165 181 166 173 168 ...
#>  $ bmi    : num  32.2 22 42.4 32.6 30.6 ...
#>  $ leg    : num  41.5 42 35.3 41.7 37.5 42.8 43 39.8 39.2 38 ...
#>  $ arml   : num  40 39.5 39 38.7 36.1 40 41.7 38.1 33.4 34.7 ...
#>  $ armc   : num  36.4 26.6 42.2 37 33.3 30.2 33.3 33.4 23 31.5 ...
#>  $ tri    : num  16.4 10.2 29.6 19 30.3 8.6 19.4 12.4 13.8 18.6 ...
#>  $ sub    : num  24.9 10.5 35.6 23.2 28 15.2 26.2 15 7.6 7.7 ...
#>  $ gh     : num  5.2 5.7 6 5.1 5.3 5.4 6.8 5.1 5.6 5.1 ...
#>  $ albumin: num  4.8 4.6 3.9 4.2 4.3 4.3 4.3 4.7 4.3 4.4 ...
#>  $ bun    : num  6 9 10 8 13 16 16 11 10 6 ...
#>  $ SCr    : num  0.94 0.89 1.11 0.8 0.79 0.83 0.9 1 0.46 0.86 ...
#>  $ age    : num  34.2 16.8 60.2 26.1 49.7 ...
#>  $ sex    : int  0 0 1 0 1 0 0 0 0 1 ...
#>  $ re     : int  3 4 4 1 3 3 3 2 4 3 ...
```
