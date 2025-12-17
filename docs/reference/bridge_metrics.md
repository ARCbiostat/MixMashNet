# Bridge metrics for nodes across communities (original-style)

Computes bridge centralities for nodes with an assigned community:

- bridge_strength: sum of absolute weights to nodes in other communities

- bridge_ei1: signed expected influence to nodes in other communities

- bridge_ei2: community-weighted expected influence (networktools-like)

- bridge_betweenness: \# times a node lies on shortest paths between
  different communities

- bridge_closeness: inverse mean distance to nodes in other communities

## Usage

``` r
bridge_metrics(g, membership)
```

## Arguments

- g:

  An igraph object with edge attribute `weight`.

- membership:

  Named vector/factor of community labels for a subset of nodes (names
  must match `V(g)$name`).

## Value

A data.frame with columns: node, cluster, bridge_strength, bridge_ei1,
bridge_ei2, bridge_betweenness, bridge_closeness.
