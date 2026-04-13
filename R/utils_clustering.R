#' Internal clustering utilities
#'
#' @description
#' Internal helper functions used to perform community detection across
#' MixMashNet functions. These utilities provide a unified interface for
#' resolving clustering methods, preparing method-specific arguments,
#' running clustering algorithms, and extracting membership vectors.
#'
#' @keywords internal
#' @noRd

# Default supported clustering methods
.default_cluster_methods <- function() {
  c(
    "louvain", "edge_betweenness", "fast_greedy", "infomap", "label_prop",
    "leading_eigen", "leiden", "optimal", "spinglass", "walktrap"
  )
}

# Null-coalescing helper
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Resolve a clustering method to the corresponding function
.resolve_cluster_fun <- function(cluster_method) {
  if (is.function(cluster_method)) {
    return(cluster_method)
  }

  switch(
    cluster_method,
    louvain          = igraph::cluster_louvain,
    edge_betweenness = igraph::cluster_edge_betweenness,
    fast_greedy      = igraph::cluster_fast_greedy,
    infomap          = igraph::cluster_infomap,
    label_prop       = igraph::cluster_label_prop,
    leading_eigen    = igraph::cluster_leading_eigen,
    leiden           = igraph::cluster_leiden,
    optimal          = igraph::cluster_optimal,
    spinglass        = igraph::cluster_spinglass,
    walktrap         = igraph::cluster_walktrap,
    stop("Unsupported `cluster_method`.")
  )
}

# Whether the clustering graph should be simplified before community detection
.needs_simplify_clustering_graph <- function(cluster_method) {
  if (is.function(cluster_method)) {
    return(TRUE)
  }

  cluster_method %in% .default_cluster_methods()
}

# Prepare method-specific arguments
.prepare_cluster_args <- function(graph, cluster_method, cluster_args = list()) {
  args <- cluster_args

  if (is.function(cluster_method)) {
    return(args)
  }

  if (cluster_method %in% c(
    "louvain", "fast_greedy", "walktrap", "edge_betweenness",
    "label_prop", "leading_eigen", "leiden", "optimal", "spinglass"
  )) {
    args$weights <- args$weights %||% igraph::E(graph)$weight
  }

  if (cluster_method == "infomap") {
    if (!is.null(args$weights) && is.null(args$e.weights)) {
      args$e.weights <- args$weights
    }
    args$weights <- NULL
    args$e.weights <- args$e.weights %||% igraph::E(graph)$weight
  }

  if (cluster_method == "leiden" && is.null(args$objective_function)) {
    args$objective_function <- "modularity"
  }

  args
}

# Run clustering
.run_clustering <- function(graph, cluster_method, cluster_args = list()) {
  fun <- .resolve_cluster_fun(cluster_method)
  args <- .prepare_cluster_args(graph, cluster_method, cluster_args)

  do.call(fun, c(list(graph), args))
}

# Extract membership robustly from different possible return types
.extract_membership <- function(x, nodes) {
  memb <- NULL

  if (inherits(x, "communities")) {
    memb <- x$membership
  } else if (is.list(x) && !is.null(x$membership)) {
    memb <- x$membership
  } else if (is.atomic(x)) {
    memb <- x
  } else {
    stop(
      "Custom clustering function must return either an igraph 'communities' object, ",
      "a list with component `$membership`, or a membership vector."
    )
  }

  if (length(memb) != length(nodes)) {
    stop("Returned membership has wrong length.")
  }

  stats::setNames(as.integer(memb), nodes)
}
