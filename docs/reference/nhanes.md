# NHANES dataset

Example dataset used in MixMashNet examples. This dataset contains 29
variables derived from the National Health and Nutrition Examination
Survey (NHANES)

## Usage

``` r
data(nhanes)
```

## Format

A data frame with 2759 rows and 29 variables:

- TotChol:

  Total Cholesterol (numeric).

- HDL:

  High-density lipoprotein cholesterol (numeric).

- Creatinine:

  Creatinine (numeric).

- UricAcid:

  Uric acid (numeric).

- ALT:

  Alanine aminotransferase (numeric).

- AST:

  Aspartate aminotransferase (numeric).

- GGT:

  Gamma-glutamyl transferase (numeric).

- Bilirubin:

  Bilirubin (numeric).

- Albumin:

  Albumin (numeric).

- TotProtein:

  Total protein (numeric).

- HbA1c:

  Glycated hemoglobin (numeric).

- hsCRP:

  High-sensitivity C-reactive protein (numeric).

- BMI:

  Body mass index (numeric).

- Waist:

  waist circumference (numeric).

- Height:

  Height (numeric).

- ArmCirc:

  Arm circumference (numeric).

- HipCirc:

  Hip circumference (numeric).

- LegLength:

  Leg length (numeric).

- ArmLength:

  Arm Length (numeric).

- TroubleSleep:

  Trouble Sleep with 0=no and 1=yes.

- PhysicalActivity:

  Physical Activity with 0=no and 1=yes.

- Smoke:

  Smoking with 0=no and 1=yes.

- Drug:

  Drug use with 0=no and 1=yes.

- Diet:

  Dietary quality with 1=Poor, 2=Fair, 3=Good, 4=Very good, 5=Excellent.

- Alcohol:

  Alcohol consumption with 0=no and 1=yes.

- Work:

  Employment status (1–4: working; employed but absent; seeking work;
  not working).

- MonInc:

  Monthly income category (1–12), where higher values indicate higher
  income.

- Gender:

  Gender with 0=male and 1=female.

- Age:

  Age (numeric).

## Examples

``` r
data(nhanes)
str(nhanes)
#> 'data.frame':    2759 obs. of  29 variables:
#>  $ TotChol         : num  157 238 182 184 230 225 213 122 184 202 ...
#>  $ HDL             : num  60 72 48 48 42 39 53 45 78 54 ...
#>  $ Creatinine      : num  0.92 1.13 0.77 1.13 0.45 1.13 0.81 0.67 0.84 1.03 ...
#>  $ UricAcid        : num  5.8 4.2 5.8 5.7 6.5 5.5 6.6 4.9 4.9 4.6 ...
#>  $ ALT             : num  16 20 46 18 22 27 13 53 19 17 ...
#>  $ AST             : num  20 23 35 18 16 20 17 36 28 20 ...
#>  $ GGT             : num  21 19 11 26 85 27 16 38 16 12 ...
#>  $ Bilirubin       : num  0.6 0.3 0.8 0.3 0.3 0.3 0.4 0.3 0.2 0.4 ...
#>  $ Albumin         : num  4.4 4 4.8 4.3 3.4 4.1 4.3 4.1 4.1 4.2 ...
#>  $ TotProtein      : num  7.3 7.1 8.1 6.4 7 7.5 7.6 7.6 6.6 6.7 ...
#>  $ HbA1c           : num  6.2 5.7 5.4 5.6 12.7 5.6 5.1 5.2 5.8 5.3 ...
#>  $ hsCRP           : num  2.72 0.82 0.37 1.66 5.71 1.75 0.8 3.04 0.41 1.02 ...
#>  $ BMI             : num  31.7 21.3 19.7 23.5 39.9 30.7 24.5 35.9 23.8 22.4 ...
#>  $ Waist           : num  101.8 86.6 72 99.7 118.4 ...
#>  $ Height          : num  158 171 173 179 148 ...
#>  $ ArmCirc         : num  32 30.8 28.7 30.6 34.7 34.2 29.9 37.1 26.5 30 ...
#>  $ HipCirc         : num  110 90.7 88.2 91 133.1 ...
#>  $ LegLength       : num  37 40.1 44.5 39.1 26 37.4 44 34.4 34 42.7 ...
#>  $ ArmLength       : num  36 37.2 37.2 41.4 32 32.6 41.4 35.7 33.4 37.9 ...
#>  $ TroubleSleep    : num  0 1 0 1 1 1 1 0 1 0 ...
#>  $ PhysicalActivity: num  0 1 1 0 0 1 0 0 0 1 ...
#>  $ Smoke           : num  1 0 1 1 1 1 1 0 0 0 ...
#>  $ Drug            : num  0 0 0 0 0 0 1 0 0 0 ...
#>  $ Diet            : num  3 4 2 2 2 5 4 3 5 3 ...
#>  $ Alcohol         : num  0 0 0 0 0 0 1 0 0 0 ...
#>  $ Work            : num  1 1 1 4 2 1 4 1 1 1 ...
#>  $ MonInc          : num  2 12 4 7 6 11 1 7 5 8 ...
#>  $ Gender          : num  1 0 0 0 1 0 0 1 1 0 ...
#>  $ Age             : num  66 56 18 67 54 61 22 60 60 64 ...
```
