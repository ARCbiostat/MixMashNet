# Bridge profiles of a node across communities

For a given node, compute how strongly it "bridges" to *other*
communities using five profiles:

- **strength**: sum of \\\|w\|\\ to nodes in other communities

- **ei1**: signed sum of \\w\\ to nodes in other communities

- **ei2**: signed influence \\A + A^2\\ to nodes in other communities

- **closeness**: \\1 /\\ mean shortest-path distance to nodes in other
  communities; paths computed on \\\|W\|\\ with edge length \\1/\|w\|\\

- **betweenness**: number of inter-community *shortest paths* (same
  graph as for closeness, but *unweighted* in the path search) that pass
  through the node as an *intermediate* (path multiplicity is counted)

## Usage

``` r
find_bridge_communities(fit, node)
```

## Arguments

- fit:

  An object of class `mixMN_fit` containing at least:

  - `$graph$keep_nodes_graph`: character vector of node names (graph
    vertex set);

  - `$statistics$edge$true`: data.frame with columns `edge` ("a–b") and
    `weight` (signed).

  - (optional) `$communities$groups`: named factor with communities for
    a subset of nodes.

- node:

  Character scalar: node of interest; must belong to
  `fit$graph$keep_nodes_graph`.

## Value

A list with components:

- `strength`:

  list with `overall` (numeric) and `by_comm` (tibble: `community`,
  `sum_abs_w`).

- `ei1`:

  list with `overall` and `by_comm` (tibble: `community`,
  `sum_signed_w`).

- `ei2`:

  list with `overall` and `by_comm` (tibble: `community`,
  `sum_signed_w2`).

- `closeness`:

  list with `overall` and `by_comm` (tibble: `community`,
  `inv_mean_dist`).

- `betweenness`:

  list with `overall` and `by_pair` (tibble: `Ci`, `Cj`, `hits`).

## Details

Notes:

- The signed adjacency \\W\\ is rebuilt *exactly* from
  `fit$statistics$edge$true` (symmetric, diagonal set to 0).

- "Other communities" means targets whose community differs from the
  node's community; if the node is unassigned (NA), targets are all
  assigned nodes.

- Betweenness counts *per path* (multiplicity). Endpoints are always
  *assigned* nodes in *different* communities. The focal node may be
  assigned or unassigned; it is counted only as intermediate (never as
  endpoint).
