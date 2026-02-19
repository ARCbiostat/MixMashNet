# Build MixMashNet metrics from a signed weighted adjacency matrix

Low-level constructor that computes clustering, centrality, and bridge
metrics from a signed weighted adjacency matrix. It returns a
`mixMN_fit` object analogous a quello prodotto da
[`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md),
ma stimato direttamente da `wadj_signed` invece che da `mgm()`. The
arguments `reps`, `seed_boot`, and `boot_what` are stored in `$settings`
for compatibility with bootstrap-based functions, but no bootstrap is
performed inside this function.

## Usage

``` r
mixMN_from_wadj(
  wadj_signed,
  nodes,
  conf_level = 0.95,
  exclude_from_graph = NULL,
  exclude_from_cluster = NULL,
  cluster_method = c("louvain", "fast_greedy", "infomap", "walktrap", "edge_betweenness"),
  reps = 0,
  seed_boot = NULL,
  treat_singletons_as_excluded = FALSE,
  boot_what = c("general_index", "interlayer_index", "bridge_index", "excluded_index",
    "community", "community_scores", "none")
)
```

## Arguments

- wadj_signed:

  Square numeric matrix with signed edge weights (rows/columns = nodes).
  Typically `wadj * signs` from mgm.

- nodes:

  Character vector of node names (must match `rownames(wadj_signed)` and
  `colnames(wadj_signed)`).

- conf_level:

  Confidence level for percentile bootstrap CIs (default 0.95). Stored
  in `$settings` for consistency with
  [`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md).

- exclude_from_graph:

  Character vector of nodes to exclude entirely from the graph (no
  edges, no centralities).

- exclude_from_cluster:

  Character vector of nodes to exclude only from clustering (they remain
  in the graph, but have no community label).

- cluster_method:

  Community detection algorithm, one of
  `c("louvain","fast_greedy","infomap","walktrap","edge_betweenness")`.

- reps:

  Integer; number of bootstrap replications. Currently not used inside
  this function, but saved in `$settings$reps` for consistency with
  [`mixMN()`](https://arcbiostat.github.io/MixMashNet/reference/mixMN.md)
  and downstream bootstrap utilities.

- seed_boot:

  Optional integer; bootstrap seed, stored in `$settings` but not used
  here.

- treat_singletons_as_excluded:

  Logical; if `TRUE`, communities of size 1 are treated as "excluded"
  (their nodes receive `NA` community labels).

- boot_what:

  Character vector specifying which quantities are meant to be
  bootstrapped in higher-level functions: can include `"general_index"`,
  `"interlayer_index"`, `"bridge_index"`, `"excluded_index"`,
  `"community"`, `"community_scores"`, or `"none"`. The value is stored
  in `$settings$boot_what` but no bootstrap is computed here.

## Value

An object of class `c("mixmashnet","mixMN_fit")` with:

- `graph`:

  List with: `igraph` (signed network), `keep_nodes_graph`, and
  `keep_nodes_cluster`.

- `communities`:

  List with empirical community structure (`original_membership`,
  `groups`), color palette (`palette`), and empty `boot_memberships`.

- `statistics`:

  Node-level metrics in `$statistics$node$true` (strength, EI1,
  closeness, betweenness, bridge metrics, including "excluded" versions)
  and edge weights in `$statistics$edge$true`. Bootstrap slots
  (`$statistics$node$boot`, `$statistics$edge$boot`,
  `$statistics$node$ci`, `$statistics$edge$ci`) are set to `NULL`.

- `community_scores`:

  Empty containers (`obj`, `df`, `ci`, `boot`) for eventual community
  scores.

- `settings`:

  List echoing the main arguments (`reps`, `cluster_method`,
  `exclude_from_graph`, `exclude_from_cluster`,
  `treat_singletons_as_excluded`, `boot_what`).

## Details

Centrality indices are computed on the signed matrix `wadj_signed` (via
[`qgraph::centrality`](https://rdrr.io/pkg/qgraph/man/centrality.html)),
while distances for closeness and betweenness are based on \\\|w\|\\
with edge length \\1/\|w\|\\. Bridge metrics are computed separately on
the absolute and signed graphs using
[`bridge_metrics()`](https://arcbiostat.github.io/MixMashNet/reference/bridge_metrics.md)
and
[`bridge_metrics_excluded()`](https://arcbiostat.github.io/MixMashNet/reference/bridge_metrics_excluded.md).
