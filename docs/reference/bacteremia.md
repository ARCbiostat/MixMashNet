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

## References

Ratzinger, F., Dedeyan, M., Rammerstorfer, M., Perkmann, T., Burgmann,
H., Makristathis, A., Dorffner, G., Loetsch, F., Blacky, A., &
Ramharter, M. (2014). A risk prediction model for screening bacteremic
patients: A cross sectional study. *PLoS ONE*, 9(9), e106765.
[doi:10.1371/journal.pone.0106765](https://doi.org/10.1371/journal.pone.0106765)
