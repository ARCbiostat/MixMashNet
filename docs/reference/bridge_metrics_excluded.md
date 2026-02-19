# Bridge metrics for nodes excluded from communities

Computes bridge centrality measures for nodes that are not assigned to
any community. This function is used internally by
[`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md)
and
[`multimixMN()`](https://arcbiostat.github.io/MixMashNet/reference/multimixMN.md).
For these excluded nodes, the function computes bridge strength, bridge
closeness, bridge betweenness, and bridge expected influence of order
one and two (EI1 and EI2), quantifying their role in connecting nodes
across different communities.

## Usage

``` r
bridge_metrics_excluded(g, membership)
```

## Arguments

- g:

  An igraph object with edge attribute `weight`.

- membership:

  Named vector/factor of community labels for a subset of nodes (names
  must match `V(g)$name`). Nodes not present here are treated as
  excluded.

## Value

A data.frame with columns: `node`, `bridge_strength`,
`bridge_closeness`, `bridge_betweenness`, `bridge_ei1`, `bridge_ei2`.

## Details

Bridge betweenness excluded and closeness excluded are computed on the
positive-weight subgraph only, with weights converted to distances as
\\d = 1/w\\.

## References

Jones, P. J. (2025). networktools: Tools for identifying important nodes
in networks. R package version 1.6.1.
<https://github.com/paytonjjones/networktools>
