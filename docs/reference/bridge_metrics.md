# Bridge metrics for nodes across communities

Computes bridge centrality measures for nodes with an assigned
community. Specifically, the function computes bridge strength as the
sum of absolute edge weights connecting a node to nodes in other
communities; bridge expected influence of order one (EI1) as the signed
sum of direct connections to nodes in other communities; bridge expected
influence of order two (EI2) as the signed influence that propagates
indirectly to nodes in other communities via one intermediate neighbor
(i.e., through paths of length two); bridge betweenness as the number of
times a node lies on shortest paths between nodes belonging to different
communities; and bridge closeness as the inverse of the mean
shortest-path distance to nodes in other communities.

## Usage

``` r
bridge_metrics(g, membership)
```

## Arguments

- g:

  An igraph object with edge attribute \`weight\`.

- membership:

  Named vector/factor of community labels for a subset of nodes (names
  must match \`V(g)\$name\`).

## Value

A data.frame with columns: node, cluster, bridge_strength, bridge_ei1,
bridge_ei2, bridge_betweenness, bridge_closeness.

## Details

Bridge betweenness and closeness are computed on the positive-weight
subgraph only, with weights converted to distances as \\d = 1/w\\.

## References

Jones, P. J. (2025). networktools: Tools for identifying important nodes
in networks. R package version 1.6.1.
<https://github.com/paytonjjones/networktools>

Jones, P. J., Ma, R., & McNally, R. J. (2021). Bridge Centrality: A
Network Approach to Understanding Comorbidity. *Multivariate Behavioral
Research*, 56(2), 353–367.
[doi:10.1080/00273171.2019.1614898](https://doi.org/10.1080/00273171.2019.1614898)
