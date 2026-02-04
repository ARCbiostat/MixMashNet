# Print method for objects of class `"bridge_profiles"`

Print method for objects of class `"bridge_profiles"`

## Usage

``` r
# S3 method for class 'bridge_profiles'
print(
  x,
  statistic = c("bridge_strength", "bridge_ei1", "bridge_ei2", "bridge_closeness",
    "bridge_betweenness"),
  digits = 3,
  ...
)
```

## Arguments

- x:

  An object of class `"bridge_profiles"`.

- statistic:

  Character string indicating which bridge profile to print. If missing,
  all available profiles are printed.

- digits:

  Number of decimal digits used when printing numeric results.

- ...:

  Further arguments passed to or from other methods (currently unused).

## Value

The input object `x`, invisibly.
