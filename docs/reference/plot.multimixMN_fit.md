# Plot method for multilayer MixMashNet objects

Plotting interface for objects returned by
[`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md).

Depending on `what`, the method can:

- `what = "network"`: plot the estimated multilayer network, or a single
  layer if `layer` is specified;

- `what = "intra"`: plot intralayer node-level metrics or edge weights
  with bootstrap quantile regions;

- `what = "inter"`: plot interlayer node metrics or interlayer edge
  weights with bootstrap quantile regions;

- `what = "stability"`: plot node stability within communities based on
  bootstrap community assignments.

## Usage

``` r
# S3 method for class 'multimixMN_fit'
plot(x, what = c("network", "intra", "inter", "stability"), layer = NULL, ...)
```

## Arguments

- x:

  An object of class `"multimixMN_fit"`, as returned by
  [`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md).

- what:

  Type of plot to produce. One of
  `c("network","intra","inter","stability")`.

- layer:

  Optional layer name. For `what = "intra"` or `what = "stability"`,
  this selects which layer-specific fit to use. For `what = "network"`,
  if provided, the selected layer is plotted as a single layer network.

- ...:

  Additional arguments. Supported arguments depend on `what`: see the
  details below.

## Value

If `what != "network"`, the function returns a `ggplot` object. If
`what = "network"`, the network is plotted directly.

## Details

**Network plots (`what = "network"`):** Supported arguments (via `...`):

- `color_by`:

  Node coloring: `c("layer","community","none")`.

- `edge_color_by`:

  Edge coloring: `c("sign","none")`.

- `edge_scale`:

  Numeric scaling factor for edge widths (multiplied by `abs(weight)`).

- `graphics::plot.igraph` arguments:

  e.g., `vertex.size`, `vertex.label.cex`, `edge.width`,
  `vertex.label.color`, etc.

**Intralayer statistics (`what = "intra"`):** Plots node-level metrics
or edge weights with bootstrap quantile regions. If `layer` is provided,
only that layer is plotted. If `layer` is `NULL`, all layers are
plotted, one panel per layer.

Supported arguments (via `...`):

- `statistics`:

  Character vector of metrics. Options include: `"strength"`,
  `"expected_influence"`, `"closeness"`, `"betweenness"`, bridge metrics
  `"bridge_strength"`, `"bridge_ei1"`, `"bridge_ei2"`,
  `"bridge_closeness"`, `"bridge_betweenness"`, excluded bridge metrics
  `"bridge_strength_excluded"`, `"bridge_ei1_excluded"`,
  `"bridge_ei2_excluded"`, `"bridge_closeness_excluded"`,
  `"bridge_betweenness_excluded"`, and `"edges"`. Different metric
  families cannot be mixed in the same call (e.g., `"edges"` cannot be
  combined with node metrics).

- `ordering`:

  Node ordering: `c("value","alphabetical","community")`.

- `standardize`:

  Logical; if `TRUE`, z-standardize the displayed values within each
  panel.

- `exclude_nodes`:

  Optional character vector of node names to remove before plotting.

- `color_by_community`:

  Logical; if `TRUE`, color nodes by community when available.

- `edges_top_n`:

  Integer; when `statistics = "edges"`, keep the top edges by absolute
  weight.

- `title`:

  Optional plot title. If omitted and multiple layers are shown,
  layer-specific titles are added automatically.

**Interlayer summaries (`what = "inter"`):** Plots interlayer node
metrics or interlayer edge weights with bootstrap quantile regions.

Supported arguments (via `...`):

- `statistics`:

  Character vector. Node metrics:
  `c("strength","expected_influence","closeness","betweenness")`, or
  `"edges"` for interlayer edge weights. Node metrics and `"edges"`
  cannot be combined.

- `pairs`:

  Layer pairs to show. Either `"*"` (all available) or a character
  vector of pair keys like `"bio_dis"` (order-insensitive).

- `edges_top_n`:

  Integer; keep the top interlayer edges by absolute weight.

- `ordering`:

  Ordering within panels: `c("value","alphabetical")`.

- `standardize`:

  Logical; if `TRUE`, z-standardize values (node metrics by metric,
  edges by pair).

- `exclude_nodes`:

  Optional character vector; removes nodes and incident interlayer
  edges.

- `nodes_layer`:

  Optional layer name to restrict node metrics to nodes belonging to
  that layer.

- `title`:

  Optional plot title.

**Community membership stability (`what = "stability"`):** Plots node
stability by community. If `layer` is provided, only that layer is
shown. Otherwise, stability plots are shown for all layers.

Supported arguments (via `...`):

- `title`:

  Plot title. Default: `"Node Stability by Community"`.

- `cutoff`:

  Optional numeric threshold in \\\[0,1\]\\ shown as a dashed vertical
  line. Use `NULL` to hide the line. Default: 0.7.

The quantile region level is taken from the fitted object
(`x$settings$quantile_level`); if missing or invalid, a default of 0.95
is used.

## See also

[`multimixMN`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md)
