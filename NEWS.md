# MixMashNet 1.0.0

## Major updates
* Added `get_centrality()` and `get_edges()` methods for `mixMN_fit` and `multimixMN_fit`, as well as `layer_slice()` for `multimixMN_fit`.
* Added the `find_bridge_layers()` function.
* Added support in `mixMN()` and `multimixMN()` for user-supplied functions in the `cluster_method` argument.
* Added the `cluster_args` argument to `mixMN()` and `multimixMN()`.
* The summaries of `mixMN()` and `multimixMN()` now display the top 10 intralayer and interlayer edges, respectively

# MixMashNet 0.6.0

## Major updates
* Deleted `warnings` parameter from `mgm_masked()`.

# MixMashNet 0.5.0

## Major updates
* Added examples to `mixMN()`, `multimixNM()`, `community_scores()` and `update_palette()`.

# MixMashNet 0.4.0

## Major updates
* Added `print()` and `summary()` for `community_scores()`.

# MixMashNet 0.3.0

## Major updates
* Simplified the interface of `mixMN()` and `multimixMN()` by removing the `type` and `level` arguments. Variable types are now inferred automatically.
* Community score computation is now restricted to gaussian, poisson and binary categorical variables.

# MixMashNet 0.1.0

* Development version.
