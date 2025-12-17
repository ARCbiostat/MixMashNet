# Plot method for MixMashNet objects

Unified plotting interface for objects returned by [`mixMN()`](mixMN.md)
and [`multimixMN()`](multimixMN.md). Depending on `what`, it can:

- `what = "network"`: plot the estimated network (single-layer or
  multilayer);

- `what = "intra"`: plot intra-layer node/edge statistics with 95\\

- `what = "inter"`: plot interlayer node metrics or interlayer edge
  weights with 95\\ `plotInterlayer()` and the chosen `statistics`;

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
  `c("network","intra","inter","stability")`. If missing, a default is
  chosen based on the presence of centrality-related arguments and on
  whether `x` is single-layer or multilayer (see Details).

- layer:

  Optional layer name. For `what = "intra"` or `what = "stability"` on a
  `multimixMN_fit` object, this selects which layer-specific fit to use.
  If `NULL`, the behaviour depends on `what`:

  - `what = "network"`: plots the global multilayer network;

  - `what = "intra"` or `"stability"` on multilayer: either all layers
    are plotted in a combined layout or an error is raised if the
    context is ambiguous.

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

For `what != "network"`, a `ggplot` object is returned. For
`what = "network"`, the corresponding network plotting helper is called
for its side-effect and `x` is returned invisibly.

## Details

When `what` is missing, the function inspects `...` to decide whether
the user is requesting centrality/edge statistics (e.g., by passing
`statistics`, `ordering`, etc.). In that case:

- for single-layer fits, the default is `what = "intra"`;

- for multilayer fits, an informative error is raised if neither `what`
  nor `layer` is specified and the request is ambiguous (intra-layer vs
  interlayer statistics).
