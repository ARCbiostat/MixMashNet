# Bacteremia dataset

Example dataset used in MixMashNet examples. This dataset contains 7240
patients with clinical suspicion of bacteremia who underwent blood
culture testing at the Vienna General Hospital

## Usage

``` r
data(bacteremia)
```

## Format

A data frame with 7420 rows and 16 variables:

- AGE:

  Age (numeric).

- WBC:

  White blood cell (numeric).

- NEU:

  Neutrophil counts (numeric).

- HGB:

  Hemoglobin (numeric).

- PLT:

  Platelet count (numeric).

- CRP:

  C-reactive protein (numeric).

- APTT:

  Activated partial thromboplastin time (numeric).

- FIB:

  Fibrinogen (numeric).

- CREA:

  Creatinine (numeric).

- BUN:

  Blood urea nitrogen (numeric).

- GLU:

  Glucose (numeric).

- ALAT:

  High-sensitivity C-reactive protein (numeric).

- GBIL:

  Total bilirubin (numeric).

- ALB:

  Albumin (numeric).

- SEX:

  Sex with 0=male and 1=female.

- BloodCulture:

  Positive blood culture results with 0=no and 1=yes.

## Examples

``` r
data(bacteremia)
str(bacteremia)
#> 'data.frame':    7420 obs. of  16 variables:
#>  $ AGE         : num  62 46 84 38 68 55 52 47 29 59 ...
#>  $ WBC         : num  24.1 17.45 11.58 9.86 9.94 ...
#>  $ NEU         : num  22 14.7 9.7 8.4 6.8 1.2 3.8 8.2 3.8 0.5 ...
#>  $ HGB         : num  11.5 7.4 10.3 13.7 15.7 10.8 10.3 9.1 7.3 10 ...
#>  $ PLT         : num  307 64 309 183 144 38 105 216 188 92 ...
#>  $ CRP         : num  3.94 12.09 3.78 11.17 5.89 ...
#>  $ APTT        : num  28.8 36.3 38.2 33.1 41.8 36.3 28.1 28.5 38.3 43.6 ...
#>  $ FIB         : num  578 313 487 490 400 413 407 604 476 369 ...
#>  $ CREA        : num  0.65 1.25 2.78 0.65 0.82 1.77 0.58 1.02 6.75 1.15 ...
#>  $ BUN         : num  5.7 50.6 47.5 8.5 15.3 29.8 14 18.6 46.3 77.1 ...
#>  $ GLU         : num  107 107 105 93 89 96 104 104 102 161 ...
#>  $ ALAT        : num  14 135 72 22 11 32 156 63 23 27 ...
#>  $ GBIL        : num  0.59 8.42 0.35 0.42 2.4 0.45 2.46 3.21 0.63 0.98 ...
#>  $ ALB         : num  36.7 22.1 33.2 43.8 30.1 43.6 24.8 26.2 36.1 29.6 ...
#>  $ SEX         : int  1 0 0 1 0 0 1 0 0 1 ...
#>  $ BloodCulture: int  0 0 0 0 0 0 0 0 0 0 ...
```
