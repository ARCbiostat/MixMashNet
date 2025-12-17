# Node-by-community stability table

Returns a data frame with node-by-community stability proportions,
optionally rounded; zero entries are replaced with NA for readability.

## Usage

``` r
membershipStab_table(stab_obj, digits = 3)
```

## Arguments

- stab_obj:

  An object from [`membershipStab()`](membershipStab.md).

- digits:

  Integer; decimal places for rounding.

## Value

A data.frame with rows = nodes and columns = communities (D1..DK).
