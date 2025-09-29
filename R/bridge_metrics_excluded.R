#' Bridge metrics for nodes excluded from communities (mmn)
#'
#' Computes bridge centrality metrics for nodes that are **not** assigned to any
#' community (treated as cluster "Z"). Uses \code{networktools::bridge} on the
#' full weighted adjacency matrix with a community vector where unassigned
#' nodes are labeled "Z".
#'
#' Returns, for excluded nodes only:
#' - `bridge_strength`
#' - `bridge_closeness`
#' - `bridge_betweenness`
#' - `bridge_expected_influence1` (1-step)
#' - `bridge_expected_influence2` (2-step)
#'
#' @param g An \code{igraph} graph. If edge attribute \code{weight} is missing,
#'          unweighted adjacency (1 for edges, 0 otherwise) is used.
#' @param membership A named vector/factor of community labels for a subset of nodes;
#'          names must match \code{V(g)$name}. Nodes not present here are treated as excluded ("Z").
#'
#' @return A data.frame with rows = excluded nodes and columns:
#'   \code{node}, \code{bridge_strength}, \code{bridge_closeness},
#'   \code{bridge_betweenness}, \code{bridge_expected_influence1},
#'   \code{bridge_expected_influence2}, \code{cluster}.
#'
#' @importFrom igraph V E vcount as_adjacency_matrix
#' @importFrom networktools bridge
#' @importFrom stats setNames
#' @export
bridge_metrics_excluded <- function(g, membership) {
  # Ensure vertex names (character)
  if (is.null(igraph::V(g)$name)) {
    warning("Vertices had no names: assigning sequential character names.")
    igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
  }
  igraph::V(g)$name <- as.character(igraph::V(g)$name)
  node_names <- igraph::V(g)$name

  # Validate membership names
  if (is.null(names(membership))) {
    stop("`membership` must be a named vector/factor with names matching V(g)$name.")
  }
  names(membership) <- as.character(names(membership))

  # Build full membership over all nodes; unassigned -> "Z"
  full_membership <- stats::setNames(rep("Z", length(node_names)), node_names)
  full_membership[names(membership)] <- as.character(membership)

  full_membership <- full_membership[node_names]

  # Weighted adjacency; fall back to unweighted if weight attr missing
  adj <- suppressWarnings(as.matrix(igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE)))
  if (!is.matrix(adj) || all(is.na(adj))) {
    adj <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
    # ensure numeric 0/1
    adj[adj != 0] <- 1
  }
  rownames(adj) <- colnames(adj) <- node_names

  adj[is.na(adj)] <- 0

  # Compute bridge metrics via networktools
  bridge_res <- networktools::bridge(network = adj, communities = full_membership)

  # Select excluded nodes ("Z")
  outside_nodes <- names(full_membership)[full_membership == "Z"]

  # If none, return empty data frame with proper columns
  if (length(outside_nodes) == 0L) {
    return(data.frame(
      node = character(0),
      bridge_strength = numeric(0),
      bridge_closeness = numeric(0),
      bridge_betweenness = numeric(0),
      bridge_expected_influence1 = numeric(0),
      bridge_expected_influence2 = numeric(0),
      cluster = character(0),
      row.names = NULL
    ))
  }

  # Assemble output (index by node names)
  out <- data.frame(
    node = outside_nodes,
    bridge_strength = bridge_res$`Bridge Strength`[outside_nodes],
    bridge_closeness = bridge_res$`Bridge Closeness`[outside_nodes],
    bridge_betweenness = bridge_res$`Bridge Betweenness`[outside_nodes],
    bridge_expected_influence1 = bridge_res$`Bridge Expected Influence (1-step)`[outside_nodes],
    bridge_expected_influence2 = bridge_res$`Bridge Expected Influence (2-step)`[outside_nodes],
    cluster = full_membership[outside_nodes],
    row.names = NULL
  )

  return(out)
}
