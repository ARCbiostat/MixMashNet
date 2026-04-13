# Plot method for single layer MixMashNet objects

Plotting interface for objects returned by
[`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md).

Depending on `what`, the method can:

- `what = "network"`: plot the estimated single-layer network;

- `what = "intra"`: plot node-level metrics or edge weights with
  bootstrap quantile regions at the level stored in the object;

- `what = "stability"`: plot node stability within communities based on
  bootstrap community assignments.

## Usage

``` r
# S3 method for class 'mixMN_fit'
plot(x, what = c("network", "intra", "stability"), ...)
```

## Arguments

- x:

  An object of class `"mixMN_fit"`, as returned by
  [`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md).

- what:

  Type of plot to produce. One of `c("network","intra","stability")`.

- ...:

  Additional arguments. Supported arguments depend on `what`: see the
  details below.

## Value

If `what != "network"`, the function returns a `ggplot` object. If
`what = "network"`, the network is plotted directly.

## Details

**Network plots (`what = "network"`):** Supported arguments (via `...`):

- `color_by`:

  Node coloring: `c("community","none")`.

- `edge_color_by`:

  Edge coloring: `c("sign","none")`.

- `edge_scale`:

  Numeric scaling factor for edge widths (multiplied by `abs(weight)`).

- `graphics::plot.igraph` arguments:

  e.g., `vertex.size`, `vertex.label.cex`, `edge.width`,
  `vertex.label.color`, etc.

**Within-network statistics (`what = "intra"`):** Plots node-level
metrics or edge weights with bootstrap quantile regions.

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

  Logical; if `TRUE`, z-standardize the displayed values within the
  panel.

- `exclude_nodes`:

  Optional character vector of node names to remove before plotting.

- `color_by_community`:

  Logical; if `TRUE`, color nodes by community when available.

- `edges_top_n`:

  Integer; when `statistics = "edges"`, keep the top edges by absolute
  weight.

- `title`:

  Optional plot title.

**Community membership stability (`what = "stability"`):** Plots node
stability by community.

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

[`mixMN`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md)
