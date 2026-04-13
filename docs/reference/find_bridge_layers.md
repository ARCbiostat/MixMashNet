# Bridge profiles of a node across layers

Identifies which layers contribute most to the interlayer bridge role of
a given node, by decomposing its interlayer connectivity into
layer-specific contributions. The function is designed as an
interpretative companion to the interlayer node-level indices returned
by
[`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md),
providing the components underlying the corresponding overall interlayer
indices.

Interlayer connectivity is summarized using four complementary profiles:
interlayer strength, interlayer expected influence (order 1), interlayer
closeness, and interlayer betweenness.

## Usage

``` r
find_bridge_layers(fit, node, layer)
```

## Arguments

- fit:

  An object of class `multimixMN_fit`.

- node:

  Character scalar: node of interest.

- layer:

  Character scalar giving the layer of the focal node.

## Value

An object of class `"bridge_layer_profiles"` (a named list) with the
following components:

- `bridge_strength`:

  List with `overall` and `by_layer`, where `by_layer` is a tibble with
  columns `target_layer` and `sum_abs_w`.

- `bridge_ei1`:

  List with `overall` and `by_layer`, where `by_layer` is a tibble with
  columns `target_layer` and `sum_signed_w`.

- `bridge_closeness`:

  List with `overall` and `by_layer`, where `by_layer` is a tibble with
  columns `target_layer` and `contribution`.

- `bridge_betweenness`:

  List with `overall` and `by_pair`, where `by_pair` is a tibble with
  columns `Li`, `Lj`, and `contribution`.

## Details

The function operates on the interlayer-only graph, i.e. on the graph
containing only edges between nodes belonging to different layers.

For a focal node in the selected `layer`, the function decomposes:

- interlayer strength into contributions toward each layer;

- interlayer expected influence (order 1) into signed contributions
  toward each layer;

- interlayer closeness into additive harmonic-distance contributions
  toward each layer;

- interlayer betweenness into additive contributions from shortest paths
  between layer pairs, using the standard fraction \\\sigma\_{st}(v) /
  \sigma\_{st}\\.

Contributions are defined so that they sum to the corresponding overall
interlayer index.
