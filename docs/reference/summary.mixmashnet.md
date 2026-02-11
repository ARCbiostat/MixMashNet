# Summarize MixMashNet fits (single and multilayer) in long format

Summarizes fitted MixMashNet objects (single and multilayer). The
summary includes the original estimates and, when available, bootstrap
means, standard errors, and quantile regions.

## Usage

``` r
# S3 method for class 'mixmashnet'
summary(
  object,
  what = c("intra", "inter"),
  statistics = NULL,
  layer = NULL,
  pairs = NULL,
  digits = 3,
  ...
)
```

## Arguments

- object:

  An object of class `"mixmashnet"` returned by [`mixMN()`](mixMN.md) or
  [`multimixMN()`](multimixMN.md).

- what:

  Character string indicating which part of the model to summarize:

  - `"intra"`: intralayer quantities (node-level indices and/or
    intralayer edges);

  - `"inter"`: interlayer quantities (node-level indices on the
    interlayer-only graph and/or cross-layer edges; multilayer fits
    only).

- statistics:

  Character vector specifying which statistics to include. For
  `what = "intra"`, valid values are:
  `c("edges", "strength", "expected_influence", "closeness", "betweenness", "bridge_strength", "bridge_closeness", "bridge_betweenness", "bridge_ei1", "bridge_ei2", "bridge_strength_excluded", "bridge_betweenness_excluded", "bridge_closeness_excluded", "bridge_ei1_excluded", "bridge_ei2_excluded")`.

  For `what = "inter"`, valid values are:
  `c("edges", "strength", "expected_influence", "closeness", "betweenness")`.

  If `statistics = NULL`, then:

  - for `what = "intra"`, all available intralayer statistics (including
    `"edges"`) are returned;

  - for `what = "inter"`, all available interlayer statistics (including
    `"edges"`) are returned.

- layer:

  Optional character vector of layer names to subset. Used for
  `what = "intra"` in multilayer fits. Ignored for single layer fits.

- pairs:

  Optional character vector of layer-pair names (e.g. `"bio_dis"`) used
  for `what = "inter"` when summarizing interlayer edges. If `NULL`, all
  available layer pairs are included. For interlayer node indices,
  `pairs` can be used to restrict the summary to nodes belonging to one
  of the layers in the given pair.

- digits:

  Number of digits to round numeric summaries.

- ...:

  Not used (for S3 compatibility).

## Value

A list (class `"summary.mixmashnet"`) with up to four data frames
(`$index`, `$edges`, `$interlayer_index`, `$interlayer_edges`) and the
quantile region used to compute the bootstrap quantile regions
(`$quantile_level`).

Depending on `what` and `statistics`, some of these elements may be
`NULL`.
