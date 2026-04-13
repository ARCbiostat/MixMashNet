# Summarize top intralayer edges for single-layer MixMashNet fits

Returns the top 10 intralayer edges for fitted objects returned by
[`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md),
in the same long-format structure used for edge summaries.

## Usage

``` r
# S3 method for class 'mixMN_fit'
summary(object, top_n = 10L, digits = NULL, drop_na_boot = TRUE, ...)
```

## Arguments

- object:

  An object of class `"mixMN_fit"`.

- top_n:

  Number of top edges to retain. Default is `10`.

- digits:

  Optional number of digits used to round numeric columns.

- drop_na_boot:

  Logical. If `TRUE` (default), bootstrap-related columns that are
  entirely `NA` are removed.

- ...:

  Further arguments for S3 compatibility.

## Value

An object of class `"summary.mixMN_fit"` containing the top intralayer
edges.
