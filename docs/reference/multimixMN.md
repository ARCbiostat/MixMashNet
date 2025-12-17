# Multilayer MGM with bootstrap, intra-/inter-layer metrics, and CIs

Estimates a **multilayer** Mixed Graphical Model (MGM) using a mask that
enforces which cross-layer edges are allowed according to `layer_rules`.
Within each layer, the function computes community structure,
non-parametric row-bootstrap for node centralities and edges, bridge
metrics (including metrics for nodes treated as "excluded"), and
(optionally) community network scores. It also returns interlayer-only
node metrics and cross-layer edge summaries.

## Usage

``` r
multimixMN(
  data,
  type,
  level,
  layers,
  layer_rules,
  scale = TRUE,
  reps = 100,
  lambdaSel = c("CV", "EBIC"),
  lambdaFolds = 5,
  lambdaGam = 0.25,
  alphaSeq = 1,
  alphaSel = "CV",
  alphaFolds = 5,
  alphaGam = 0.25,
  k = 2,
  ruleReg = "AND",
  threshold = "LW",
  overparameterize = FALSE,
  thresholdCat = TRUE,
  conf_level = 0.95,
  exclude_from_graph = NULL,
  exclude_from_cluster = NULL,
  seed_model = NULL,
  seed_boot = NULL,
  treat_singletons_as_excluded = FALSE,
  cluster_method = c("louvain", "fast_greedy", "infomap", "walktrap", "edge_betweenness"),
  compute_community_scores = FALSE,
  boot_what = c("general_index", "interlayer_index", "bridge_index", "excluded_index",
    "community", "community_scores", "none"),
  progress = TRUE
)
```

## Arguments

- data:

  A numeric matrix or data frame (n × p) with variables in columns.

- type:

  Length-`p` vector of variable types as required by
  [`mgm::mgm`](https://rdrr.io/pkg/mgm/man/mgm.html).

- level:

  Length-`p` vector of variable levels as required by
  [`mgm::mgm`](https://rdrr.io/pkg/mgm/man/mgm.html).

- layers:

  A named vector (names = variable names) assigning each node to a layer
  (character or factor). Must cover all columns of `data`.

- layer_rules:

  A logical or numeric square matrix with row/column names equal to
  layer names. Values `TRUE` or `1` indicate that cross-layer edges are
  allowed between the corresponding layer pair. Intra-layer edges are
  always allowed; if the diagonal is `NA`, it is treated as allowed.

- scale:

  Logical; if `TRUE` (default) Gaussian variables (`type == "g"`) are
  z-standardized internally by `mgm()`. Use `scale = FALSE` if your data
  are already standardized.

- reps:

  Integer (\>= 0). Number of bootstrap replications (row resampling with
  replacement). If `reps = 0`, no bootstrap is performed.

- lambdaSel:

  Method for lambda selection in mgm: `"CV"` or `"EBIC"`.

- lambdaFolds:

  Number of folds for CV (if `lambdaSel = "CV"`).

- lambdaGam:

  EBIC gamma parameter (if `lambdaSel = "EBIC"`).

- alphaSeq:

  Alpha parameters of the elastic net penalty (values in \\(0, 1\]\\).

- alphaSel:

  Method for selecting the alpha parameter: `"CV"` or `"EBIC"`.

- alphaFolds:

  Number of folds for CV (if `alphaSel = "CV"`).

- alphaGam:

  EBIC gamma parameter (if `alphaSel = "EBIC"`).

- k:

  Integer (\>= 1). Order of modeled interactions.

- ruleReg:

  Rule to combine neighborhood estimates: `"AND"` or `"OR"`.

- threshold:

  Threshold below which edge-weights are set to zero: `"LW"`, `"HW"` or
  `"none"`.

- overparameterize:

  Logical; if `TRUE` uses the over-parameterized version of `mgm`.

- thresholdCat:

  Logical; if `FALSE` thresholds of categorical variables are set to
  zero.

- conf_level:

  Confidence level for percentile bootstrap CIs (default 0.95). Must be
  a single number between 0 and 1 (e.g., 0.90, 0.95, 0.99).

- exclude_from_graph:

  Character vector of node names. Nodes in this set are excluded from
  the global graph and from all node-level metrics.

- exclude_from_cluster:

  Character vector of node names. Nodes in this set are excluded from
  community detection in addition to `exclude_from_graph`.

- seed_model:

  Optional integer seed for reproducibility of the initial MGM fit.

- seed_boot:

  Optional integer seed passed to `future.apply` for reproducibility of
  bootstrap replications.

- treat_singletons_as_excluded:

  Logical; if `TRUE`, singleton communities (size 1) are treated as
  "excluded" when computing bridge metrics and related summaries.

- cluster_method:

  Community detection algorithm applied within each layer. One of
  `"louvain"`, `"fast_greedy"`, `"infomap"`, `"walktrap"`, or
  `"edge_betweenness"`.

- compute_community_scores:

  Logical; if `TRUE`, compute community network scores (EGAnet
  standardized `std.scores`) with bootstrap arrays and confidence
  intervals.

- boot_what:

  Character vector specifying which quantities to bootstrap. Valid
  options are: `"general_index"` (intra-layer centrality indices),
  `"interlayer_index"` (interlayer-only node metrics), `"bridge_index"`
  (bridge metrics for nodes in communities), `"excluded_index"` (bridge
  metrics for nodes treated as excluded), `"community"` (bootstrap
  community memberships), `"community_scores"` (bootstrap community
  scores, only if `compute_community_scores = TRUE`), and `"none"` (skip
  all node-level bootstrap: only edge-weight bootstrap is performed if
  `reps > 0`).

- progress:

  Logical; if `TRUE` (default), show a bootstrap progress bar when the
  progressr package is available.

## Value

An object of class `c("mixmashnet", "multimixMN_fit")`. The returned
list contains at least the following components:

- `call`: the matched function call.

- `settings`: list of main settings (e.g., `reps`, `cluster_method`,
  `exclude_from_graph`, `exclude_from_cluster`,
  `treat_singletons_as_excluded`, `boot_what`).

- `model`: list with the fitted masked `mgm` object and basic
  information on the data (`nodes`, `n`, `p`).

- `layers`: list describing the multilayer structure (assignment of
  nodes to layers and the `layer_rules` matrix used).

- `layer_fits`: named list (one element per layer) with single-layer
  fits, including community structure, node-level statistics, edge-level
  statistics, bridge metrics, and (optionally) community scores with
  bootstrap information.

- `interlayer`: list collecting interlayer-only node metrics (strength,
  expected influence, closeness, betweenness, with or without bootstrap)
  and cross-layer edge summaries for each allowed pair of layers.

- `graph`: list containing a global igraph object built on the retained
  nodes (`keep_nodes_graph`), with vertex attributes such as `name`,
  `layer`, `membership`, and edge attributes such as `weight`,
  `abs_weight`, `sign`, `type` (intra vs inter) and `layer_pair`.

## Details

This function does **not** call
[`future::plan()`](https://future.futureverse.org/reference/plan.html).
To enable parallel bootstrap, set a plan (e.g.
`future::plan(multisession)`) before calling `multimixMN()`. If `"none"`
is the only element of `boot_what` and `reps > 0`, node-level metrics
are not bootstrapped, but intra- and inter-layer edge-weight bootstrap
and the corresponding confidence intervals are still computed.

## Dependencies

Relies on internal helpers [`mgm_masked()`](mgm_masked.md),
[`mixMN_from_wadj()`](mixMN_from_wadj.md),
[`bridge_metrics()`](bridge_metrics.md), and
[`bridge_metrics_excluded()`](bridge_metrics_excluded.md), which must be
available in the same package namespace.
