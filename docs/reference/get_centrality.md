# Extract node-level centrality indices

Extracts node-level centrality indices from fitted objects returned by
[`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md)
and
[`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md)
in a long-format data frame.

For single layer fits of class `"mixMN_fit"`, only intralayer node-level
indices are available.

For multilayer fits of class `"multimixMN_fit"`, `what = "intra"`
returns intralayer node-level indices, whereas `what = "inter"` returns
node-level indices computed on the interlayer-only graph. If `what` is
not specified for a multilayer fit, both intralayer and interlayer
node-level indices are returned by default, unless `layer` or `pairs`
imply a specific scope.

The function returns the original estimates and, when available,
bootstrap means, standard errors, and bootstrap quantile regions.

## Usage

``` r
get_centrality(object, ...)

# S3 method for class 'mixMN_fit'
get_centrality(
  object,
  what = "intra",
  statistics = NULL,
  digits = NULL,
  drop_na_boot = TRUE,
  ...
)

# S3 method for class 'multimixMN_fit'
get_centrality(
  object,
  what = c("intra", "inter"),
  statistics = NULL,
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

  Character string indicating which node-level indices to extract:

  - `"intra"`: intralayer node-level indices;

  - `"inter"`: interlayer-only node-level indices (available only for
    `"multimixMN_fit"` objects).

  For single layer fits, only `"intra"` is allowed. For multilayer fits,
  if omitted, both scopes are returned by default unless `layer` or
  `pairs` imply a specific scope.

- statistics:

  Character vector specifying which node-level statistics to include.

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

A tibble in long format with one row per node-statistic combination. It
contains the columns:

- `node`

- `layer`

- `scope`

- `metric`

- `estimated`

When available, the output also contains bootstrap summary columns:

- `mean.bootstrap`

- `SE.bootstrap`

- `quantile.lower.bootstrap`

- `quantile.upper.bootstrap`

The quantile level used to compute the bootstrap quantile region is
stored as the `"quantile_level"` attribute of the returned tibble.

## Details

The returned data frame is in long format, with one row per
node-statistic combination.

For single layer fits, only `what = "intra"` is available.

For multilayer fits, `layer` can be used to subset intralayer output,
whereas `pairs` can be used to subset interlayer output.

The set of admissible statistics depends on `what` and on the class of
`object`. In particular, bridge-related indices are available only for
intralayer output.
