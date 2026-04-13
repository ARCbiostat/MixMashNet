# Extract edge-level summaries

Extracts edge-level summaries from fitted objects returned by
[`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md)
and
[`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md)
in a long-format data frame.

For single layer fits of class `"mixMN_fit"`, only intralayer edges are
available.

For multilayer fits of class `"multimixMN_fit"`, `what = "intra"`
returns intralayer edges, whereas `what = "inter"` returns interlayer
edges. If `what` is not specified for a multilayer fit, both scopes are
returned by default, unless `layer` or `pairs` imply a specific scope.

The function returns the original edge weights and, when available,
bootstrap means, standard errors, and bootstrap quantile regions.

## Usage

``` r
get_edges(object, ...)

# S3 method for class 'mixMN_fit'
get_edges(object, what = "intra", digits = NULL, drop_na_boot = TRUE, ...)

# S3 method for class 'multimixMN_fit'
get_edges(
  object,
  what = c("intra", "inter"),
  layer = NULL,
  pairs = NULL,
  digits = NULL,
  drop_na_boot = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted object of class `"mixMN_fit"` or `"multimixMN_fit"`.

- ...:

  Further arguments passed to methods.

- what:

  Character string indicating which edge-level summaries to extract:

  - `"intra"`: intralayer edges;

  - `"inter"`: interlayer edges (available only for `"multimixMN_fit"`
    objects).

  For single layer fits, only `"intra"` is allowed. For multilayer fits,
  if omitted, both scopes are returned by default unless `layer` or
  `pairs` imply a specific scope.

- digits:

  Optional number of digits used to round numeric columns.

- drop_na_boot:

  Logical. If `TRUE` (default), bootstrap-related columns that are
  entirely `NA` are removed from the output.

- layer:

  Optional character vector of layer names to subset. Relevant for
  intralayer output in multilayer fits.

- pairs:

  Optional character vector of layer-pair names to subset. Relevant for
  interlayer output in multilayer fits.

## Value

A tibble in long format with one row per edge. It contains the columns:

- `edge`

- `layer` for intralayer edges

- `pairs` for interlayer edges

- `scope`

- `estimated`

When available, the output also contains bootstrap summary columns:

- `mean.bootstrap`

- `SE.bootstrap`

- `quantile.lower.bootstrap`

- `quantile.upper.bootstrap`

The quantile level used to compute the bootstrap quantile region is
stored as the `"quantile_level"` attribute of the returned tibble.

## Details

The returned data frame is in long format, with one row per edge.

For single layer fits, only `what = "intra"` is available.

For multilayer fits, `layer` can be used to subset intralayer output,
whereas `pairs` can be used to subset interlayer output.
