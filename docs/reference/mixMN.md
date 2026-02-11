# Estimate single layer MGM network with bootstrap centrality, bridge metrics, clustering, and (optionally) community score loadings

Estimates a single layer Mixed Graphical Model (MGM) network on the
original data, using the estimation framework implemented in the mgm
package, and performs non-parametric bootstrap (row resampling) to
compute centrality indices, bridge metrics, clustering stability, and
quantile regions for node metrics and edge weights. Optionally, the
function computes community score loadings (for later prediction on new
data) and can bootstrap the corresponding loadings.

## Usage

``` r
mixMN(
  data,
  reps = 100,
  scale = TRUE,
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
  quantile_level = 0.95,
  covariates = NULL,
  exclude_from_cluster = NULL,
  treat_singletons_as_excluded = FALSE,
  seed_model = NULL,
  seed_boot = NULL,
  cluster_method = c("louvain", "fast_greedy", "infomap", "walktrap", "edge_betweenness"),
  compute_loadings = TRUE,
  boot_what = c("general_index", "bridge_index", "excluded_index", "community",
    "loadings"),
  save_data = FALSE,
  progress = TRUE
)
```

## Arguments

- data:

  A `data.frame` (n x p) with variables in columns. Variables may be
  numeric, integer, logical, or factors. Character and Date/POSIXt
  variables are not supported and must be converted prior to model
  fitting. Variable types are internally mapped to MGM types as follows:
  numeric variables are treated as Gaussian; integer variables are
  treated as Poisson unless they take only values in {0,1}, in which
  case they are treated as binary categorical; factors and logical
  variables are treated as categorical. Binary categorical variables
  (two-level factors and logical variables) are internally recoded to
  {0,1} for model fitting. The original input data are not modified.

- reps:

  Integer (\>= 0). Number of bootstrap replications.

- scale:

  Logical; if `TRUE` (default) Gaussian variables (`type == "g"`) are
  z-standardized internally by `mgm()`. Use `scale = FALSE` if your data
  are already standardized.

- lambdaSel:

  Method for lambda selection: `"CV"` or `"EBIC"`.

- lambdaFolds:

  Number of folds for CV (if `lambdaSel = "CV"`).

- lambdaGam:

  EBIC gamma parameter (if `lambdaSel = "EBIC"`).

- alphaSeq:

  Alpha parameters of the elastic net penalty (values between 0 and 1).

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

- quantile_level:

  Level of the central bootstrap quantile region (default 0.95). Must be
  a single number between 0 and 1.

- covariates:

  Character vector. Variables used as adjustment covariates in model
  estimation.

- exclude_from_cluster:

  Character vector. Nodes excluded from community detection (in addition
  to `covariates`).

- treat_singletons_as_excluded:

  Logical; if `TRUE`, singleton communities (size 1) are treated as
  excluded nodes when computing bridge metrics.

- seed_model:

  Optional integer seed for reproducibility of the initial MGM fit.

- seed_boot:

  Optional integer seed passed to `future.apply` for reproducibility of
  bootstrap replications.

- cluster_method:

  Community detection algorithm used on the network: `"louvain"`,
  `"fast_greedy"`, `"infomap"`, `"walktrap"`, or `"edge_betweenness"`.

- compute_loadings:

  Logical; if `TRUE` (default), compute network loadings (EGAnet
  net.loads) for communities.

- boot_what:

  Character vector specifying which quantities to bootstrap. Valid
  options are: `"general_index"` (centrality indices), `"bridge_index"`
  (bridge metrics for nodes in communities), `"excluded_index"` (bridge
  metrics for nodes treated as excluded), `"community"` (community
  memberships), `"loadings"` (community loadings, only if
  `compute_loadings = TRUE`), and `"none"` (skip all node-level
  bootstrap: only edge-weight bootstrap is performed if `reps > 0`).

- save_data:

  Logical; if `TRUE`, store the original data in the output object.

- progress:

  Logical; if `TRUE` (default), show a bootstrap progress bar.

## Value

An object of class `c("mixmashnet", "mixMN_fit")`, that is a list with
the following top-level components:

- `call`:

  The matched function call.

- `settings`:

  List of main settings used in the call, including `reps`,
  `cluster_method`, `covariates`, `exclude_from_cluster`,
  `treat_singletons_as_excluded`, and `boot_what`.

- `data_info`:

  List with information derived from the input data used for model
  setup: `mgm_type_level` (data frame with one row per variable,
  reporting the original R class and the inferred MGM `type` and
  `level`, as used in the call to
  [`mgm::mgm`](https://rdrr.io/pkg/mgm/man/mgm.html)), and
  `binary_recode_map` (named list describing the mapping from original
  binary labels to the internal {0,1} coding used for model fitting).

- `model`:

  List with: `mgm` (the fitted `mgm` object), `nodes` (character vector
  of all node names), `n` (number of observations), `p` (number of
  variables), and `data`.

- `graph`:

  List describing the graph: `igraph` (an igraph object built on
  `keep_nodes_graph`, with edge attributes `weight`, `abs_weight`,
  `sign` and vertex attribute `membership` for communities),
  `keep_nodes_graph` (nodes retained in the graph and all node-level
  metrics), and `keep_nodes_cluster` (nodes used for community
  detection).

- `communities`:

  List describing community structure with: `original_membership`
  (integer vector of community labels on `keep_nodes_cluster`), `groups`
  (factor of community labels actually used for bridge metrics,
  optionally with singletons treated as excluded), `palette` (named
  vector of colors per community), and `boot_memberships` (list of
  bootstrap memberships if `"community"` is requested in `boot_what`,
  otherwise an empty list).

- `statistics`:

  List with node- and edge-level summaries: `node` is a list with:
  `true` (data frame with one row per node in `keep_nodes_graph`,
  containing the node name and metrics `strength`, `ei1`, `closeness`,
  `betweenness`, `bridge_strength`, `bridge_betweenness`,
  `bridge_closeness`, `bridge_ei1`, `bridge_ei2`, and for nodes treated
  as excluded from communities also `bridge_strength_excluded`,
  `bridge_betweenness_excluded`, `bridge_closeness_excluded`,
  `bridge_ei1_excluded`, `bridge_ei2_excluded`); `boot` (list of
  bootstrap matrices for each metric, each of dimension
  `reps x length(keep_nodes_graph)`, possibly `NULL` if the metric was
  not requested or if `reps = 0`); and `quantile_region` (list of
  quantile regions for each node metric, one `p x 2` matrix per metric,
  with columns corresponding to the lower and upper quantile bounds
  implied by `quantile_level`, or `NULL` if no bootstrap was performed).

  `edge` is a list with: `true` (data frame with columns `edge` and
  `weight` for all unique undirected edges among `keep_nodes_graph`);
  `boot` (matrix of bootstrap edge weights of dimension
  `n_edges x reps`); and `quantile_region` (matrix of quantile regions
  for edge weights, `n_edges x 2`, with columns corresponding to the
  lower and upper bootstrap quantile bounds, or `NULL` if `reps = 0`).

- `community_loadings`:

  List containing community-loading information (based on
  [`EGAnet::net.loads`](https://r-ega.net/reference/net.loads.html)) for
  later community-score computation on new data: `nodes` (nodes used for
  loadings), `wc` (integer community labels aligned with `nodes`),
  `true` (matrix of standardized loadings, nodes x communities), and
  `boot` (list of bootstrap loading matrices, one per replication, or
  `NULL` if not bootstrapped).

## Details

This function does **not** call
[`future::plan()`](https://future.futureverse.org/reference/plan.html).
To enable parallel bootstrap, set a plan (e.g.
`future::plan(multisession)`) before calling `mixMN()`. If `boot_what`
is `"none"` and `reps > 0`, node-level metrics are not bootstrapped but
edge-weight bootstrap and corresponding quantile regions are still
computed.

## References

Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating
Time-Varying Mixed Graphical Models in High-Dimensional Data. *Journal
of Statistical Software*, 93(8).
[doi:10.18637/jss.v093.i08](https://doi.org/10.18637/jss.v093.i08)
