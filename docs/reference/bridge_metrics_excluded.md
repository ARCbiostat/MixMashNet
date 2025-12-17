# Bridge metrics for nodes excluded from communities (mmn)

Computes bridge centrality metrics for nodes that are **not** assigned
to any community (treated as cluster "Z"). Uses
[`networktools::bridge`](https://rdrr.io/pkg/networktools/man/bridge.html)
on the full weighted adjacency matrix with a community vector where
unassigned nodes are labeled "Z".

## Usage

``` r
bridge_metrics_excluded(g, membership)
```

## Arguments

- g:

  An `igraph` graph. If edge attribute `weight` is missing, unweighted
  adjacency (1 for edges, 0 otherwise) is used.

- membership:

  A named vector/factor of community labels for a subset of nodes; names
  must match `V(g)$name`. Nodes not present here are treated as excluded
  ("Z").

## Value

A data.frame with rows = excluded nodes and columns: `node`,
`bridge_strength`, `bridge_closeness`, `bridge_betweenness`,
`bridge_expected_influence1`, `bridge_expected_influence2`, `cluster`.

## Details

Returns, for excluded nodes only:

- `bridge_strength`

- `bridge_closeness`

- `bridge_betweenness`

- `bridge_expected_influence1` (1-step)

- `bridge_expected_influence2` (2-step)
