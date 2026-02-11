# Plot method for MixMashNet objects

Unified plotting interface for objects returned by [`mixMN()`](mixMN.md)
and [`multimixMN()`](multimixMN.md). Depending on `what`, it can:

- `what = "network"`: plot the estimated network (single layer or
  multilayer);

- `what = "intra"`: plot intralayer node/edge statistics with bootstrap
  quantile regions at the level stored in the object (centrality and
  bridge metrics);

- `what = "inter"`: plot interlayer node metrics or interlayer edge
  weights with bootstrap quantile regions at the level stored in the
  object (multilayer only), and the chosen `statistics`;

- `what = "stability"`: plot node stability within communities based on
  bootstrap community assignments.

## Usage

``` r
# S3 method for class 'mixmashnet'
plot(x, what = c("network", "intra", "inter", "stability"), layer = NULL, ...)
```

## Arguments

- x:

  An object of class `mixmashnet`, as returned by [`mixMN()`](mixMN.md)
  or [`multimixMN()`](multimixMN.md).

- what:

  Type of plot to produce. One of
  `c("network","intra","inter","stability")`.

- layer:

  Optional layer name. For `what = "intra"` or `what = "stability"` on a
  `multimixMN_fit` object, this selects which layer-specific fit to use.

- ...:

  Additional arguments passed to the underlying helpers:

  - For `what = "intra"`: forwarded to `plotCentrality()`, e.g.
    `statistics`, `ordering`, `standardize`, `edges_top_n`,
    `exclude_nodes`, `color_by_community`.

  - For `what = "inter"`: forwarded to `plotInterlayer()`, e.g.
    `statistics`, `pairs`, `edges_top_n`, `ordering`, `standardize`,
    `nodes_layer`.

  - For `what = "stability"`: forwarded to `membershipStab_plot()`, e.g.
    `title`.

  - For `what = "network"`: forwarded to internal network plotting
    helpers `.plot_network_single()` or `.plot_network_multi()`, e.g.
    `color_by`, `edge_color_by`, `vertex_size`.

## Value

If `what != "network"`, the function returns a `ggplot` object. If
`what = "network"`, the network is plotted directly.
