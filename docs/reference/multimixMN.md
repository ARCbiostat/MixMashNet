# Multilayer MGM with bootstrap, intra/interlayer metrics, and quantile regions

Estimates a multilayer Mixed Graphical Model (MGM) using the estimation
framework implemented in the mgm package, with a masking scheme that
enforces which cross-layer edges are allowed according to `layer_rules`.
Within each layer, the function computes community structure and
performs non-parametric row-bootstrap to obtain node centrality indices,
edge weights, and bridge metrics, including metrics for nodes treated as
excluded. Optionally, within-layer community loadings can also be
estimated and bootstrapped. The function additionally returns
interlayer-only node metrics and summaries of cross-layer edge weights.

## Usage

``` r
multimixMN(
  data,
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
  quantile_level = 0.95,
  covariates = NULL,
  exclude_from_cluster = NULL,
  seed_model = NULL,
  seed_boot = NULL,
  treat_singletons_as_excluded = FALSE,
  cluster_method = c("louvain", "fast_greedy", "infomap", "walktrap", "edge_betweenness"),
  compute_loadings = TRUE,
  boot_what = c("general_index", "interlayer_index", "bridge_index", "excluded_index",
    "community", "loadings"),
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
  continuous numeric (double) variables are treated as Gaussian; integer
  variables are treated as Poisson unless they take only values in
  {0,1}, in which case they are treated as binary categorical; factors
  and logical variables are treated as categorical. Binary categorical
  variables (two-level factors and logical variables) are internally
  recoded to {0,1} for model fitting. The original input data are not
  modified.

- layers:

  A named vector (names = variable names) assigning each node to a layer
  (character or factor). Must cover all columns of `data` except
  variables listed in `covariates` (treated as adjustment covariates).

- layer_rules:

  A square matrix (L × L), where L is the number of layers. Row and
  column names must match the layer names. Entries equal to `TRUE` or
  `1` allow cross-layer edges between the corresponding pair of layers,
  while `FALSE` or `0` disallow them. The matrix is symmetrized
  internally. Diagonal entries are ignored (intralayer edges are always
  permitted).

- scale:

  Logical; if `TRUE` (default) Gaussian variables (`type == "g"`) are
  z-standardized internally by `mgm()`. Use `scale = FALSE` if your data
  are already standardized.

- reps:

  Integer (\>= 0). Number of bootstrap replications (row resampling with
  replacement). If `reps = 0`, no bootstrap is performed.

- lambdaSel:

  Method for lambda selection in `mgm`: `"CV"` or `"EBIC"`.

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

  Threshold below which edge-weights are set to zero: Available options
  are `"LW"`, `"HW"`, or `"none"`. `"LW"` applies the threshold proposed
  by Loh & Wainwright; `"HW"` applies the threshold proposed by Haslbeck
  & Waldorp; `"none"` disables thresholding. Defaults to `"LW"`.

- overparameterize:

  Logical; controls how categorical interactions are parameterized in
  the neighborhood regressions. If `TRUE`, categorical interactions are
  represented using a fully over-parameterized design matrix (i.e., all
  category combinations are explicitly modeled). If `FALSE`, the
  standard `glmnet` parameterization is used, where one category serves
  as reference. For continuous variables, both parameterizations are
  equivalent. The default is `FALSE`. The over-parameterized option may
  be advantageous when distinguishing pairwise from higher-order
  interactions.

- thresholdCat:

  Logical; if `FALSE` thresholds of categorical variables are set to
  zero.

- quantile_level:

  Level of the central bootstrap quantile region (default `0.95`). Must
  be a single number between 0 and 1.

- covariates:

  Character vector. Variables used as adjustment covariates in model
  estimation.

- exclude_from_cluster:

  Character vector of node names. Nodes in this set are excluded from
  community detection in addition to `covariates`.

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

- compute_loadings:

  Logical; if `TRUE` (default), compute community loadings
  ([`EGAnet::net.loads`](https://r-ega.net/reference/net.loads.html)).
  Only supported for Gaussian, Poisson, and binary categorical nodes;
  otherwise loadings are skipped and the reason is stored in
  `community_loadings$reason`.

- boot_what:

  Character vector specifying which quantities to bootstrap. Valid
  options are: `"general_index"` (intralayer centrality indices),
  `"interlayer_index"` (interlayer-only node metrics), `"bridge_index"`
  (bridge metrics for nodes in communities), `"excluded_index"` (bridge
  metrics for nodes treated as excluded), `"community"` (community
  memberships), `"loadings"` (within-layer community loadings, only if
  `compute_loadings = TRUE`), and `"none"` (skip all node-level
  bootstrap: only edge-weight bootstrap is performed if `reps > 0`).

- save_data:

  Logical; if `TRUE`, store the original data in the output object.

- progress:

  Logical; if `TRUE` (default), show a bootstrap progress bar.

## Value

An object of class `c("mixmashnet", "multimixMN_fit")`. The returned
list contains at least the following components:

- `call`:

  The matched function call.

- `settings`:

  List of main settings used in the call, including `reps`,
  `cluster_method`, `covariates`, `exclude_from_cluster`,
  `treat_singletons_as_excluded`, `boot_what`).

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
  variables), and `data` (if `save_data = TRUE`))

- `layers`:

  List describing the multilayer structure (assignment of nodes to
  layers, `layer_rules` matrix used and color of each layer in
  `palette`).

- `layer_fits`:

  Named list (one element per layer) with single layer fits, including
  community structure, node-level statistics, edge-level statistics,
  bridge metrics, and (optionally) community loadings with bootstrap
  information.

- `interlayer`:

  List collecting interlayer-only node metrics (strength, expected
  influence, closeness, betweenness, with or without bootstrap) and
  cross-layer edge summaries for each allowed pair of layers.

- `graph`:

  List containing a global igraph object built on the retained nodes
  (`keep_nodes_graph`), with vertex attributes such as `name`, `layer`,
  `membership`, and edge attributes such as `weight`, `abs_weight`,
  `sign`, `type` (intra vs inter) and `layer_pair`.

## Details

This function does **not** call
[`future::plan()`](https://future.futureverse.org/reference/plan.html).
To enable parallel bootstrap, set a plan (e.g.
`future::plan(multisession)`) before calling `multimixMN()`. If `"none"`
is the only element of `boot_what` and `reps > 0`, node-level metrics
are not bootstrapped, but intra and interlayer edge-weight bootstrap and
the corresponding quantile regions are still computed.

## References

Haslbeck, J. M. B., & Waldorp, L. J. (2020). mgm: Estimating
Time-Varying Mixed Graphical Models in High-Dimensional Data. *Journal
of Statistical Software*, 93(8).
[doi:10.18637/jss.v093.i08](https://doi.org/10.18637/jss.v093.i08)
