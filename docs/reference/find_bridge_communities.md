# Bridge profiles of a node across communities

Identifies which communities contribute most to the bridge role of a
given node, by decomposing its bridge connectivity into
community-specific contributions, excluding its own community when
assigned. The function is designed as an interpretative companion to
[`bridge_metrics()`](https://arcbiostat.github.io/MixMashNet/reference/bridge_metrics.md)
and
[`bridge_metrics_excluded()`](https://arcbiostat.github.io/MixMashNet/reference/bridge_metrics_excluded.md),
providing the components underlying the corresponding overall bridge
indices.

Bridge connectivity is summarized using five complementary profiles:
bridge strength, bridge EI1, bridge EI2, bridge closeness, and bridge
betweenness.

For single layer fits (`mixMN_fit`), profiles are computed directly on
the supplied fitted object.

For multilayer fits (`multimixMN_fit`), profiles are computed within the
selected layer only, by applying the same single layer procedure to the
corresponding intralayer fit stored in `fit$layer_fits[[layer]]`.

## Usage

``` r
find_bridge_communities(fit, node, layer = NULL)
```

## Arguments

- fit:

  An object of class `mixMN_fit` or `multimixMN_fit`.

- node:

  Character scalar: node of interest.

- layer:

  Character scalar giving the layer of interest for `multimixMN_fit`
  objects. Ignored for `mixMN_fit` objects.

## Value

An object of class `"bridge_profiles"` (a named list) with the following
components:

- `bridge_strength`:

  Bridge strength. List with `overall`, the total value across all other
  communities, and `by_comm`, a tibble with community-specific
  contributions (`community`, `sum_abs_w`).

- `bridge_ei1`:

  Bridge expected influence (order 1). List with `overall` and `by_comm`
  (`community`, `sum_signed_w`).

- `bridge_ei2`:

  Bridge expected influence (order 2). List with `overall` and `by_comm`
  (`community`, `sum_signed_w2`).

- `bridge_closeness`:

  Bridge closeness. List with `overall` and `by_comm` (`community`,
  `inv_mean_dist`).

- `bridge_betweenness`:

  Bridge betweenness. List with `overall` and `by_pair`, a tibble with
  contributions by community pair (`Ci`, `Cj`, `hits`).

## Details

Bridge profiles are computed using only connections from the focal node
to nodes in communities different from its own. If the focal node is not
assigned to any community, i.e. excluded, connections to all assigned
nodes in communities are considered.

Bridge betweenness is computed by counting all shortest paths between
pairs of nodes in different communities that pass through the focal node
as an intermediate vertex. When multiple shortest paths exist, each path
is counted separately.
