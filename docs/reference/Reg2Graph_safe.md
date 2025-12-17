# Safe graph assembly for mgm outputs (modified from mgm)

Robust post-processing to assemble the pairwise (and nodewise) weighted
adjacency matrices from an `mgm`-style object. This variant tolerates
masked designs and empty parameter blocks by using conservative
fallbacks and `na.rm = TRUE` where appropriate. It preserves the
structure of the original mgm output.

## Usage

``` r
Reg2Graph_safe(mgmobj, thresholding = TRUE)
```

## Arguments

- mgmobj:

  An object produced by a nodewise mgm-style routine.

- thresholding:

  Logical; if `TRUE`, applies the thresholding rule recorded in
  `mgmobj$call$threshold`.

## Value

The input object with populated `$pairwise` matrices and interaction
lists.
