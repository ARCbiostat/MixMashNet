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

  # --------- Recompute betweenness for EXCLUDED nodes (Z) ONLY ----------
  directed <- igraph::is_directed(g)
  if (directed) {
    g_nt <- igraph::graph_from_adjacency_matrix(adj, mode = "directed", diag = FALSE, weighted = TRUE)
  } else {
    g_nt <- igraph::graph_from_adjacency_matrix(adj, mode = "upper",   diag = FALSE, weighted = TRUE)
  }
  igraph::E(g_nt)$weight <- 1 / igraph::E(g_nt)$weight
  g2 <- g_nt
  if (igraph::ecount(g2) > 0 && min(igraph::E(g2)$weight, na.rm = TRUE) < 0) {
    g2 <- igraph::delete_edges(g2, which(igraph::E(g2)$weight < 0))
  }

  # Endpoints = only real communities (exclude "Z")
  assigned_nodes <- names(full_membership)[full_membership != "Z"]
  assigned_comm  <- full_membership[assigned_nodes]

  # Helper: intermediate vertex ids of a path
  mid_ids <- function(path_ids) {
    ids <- as.vector(path_ids)
    if (length(ids) <= 2L) integer(0) else ids[-c(1L, length(ids))]
  }

  # Compute custom betweenness for excluded nodes, counting PER TARGET
  custom_betw <- setNames(numeric(0), character(0))
  outside_nodes <- names(full_membership)[full_membership == "Z"]

  if (length(outside_nodes) && length(assigned_nodes) && igraph::vcount(g2) > 0) {
    vnames_g2 <- igraph::V(g2)$name
    for (z in outside_nodes) {
      # If z not present (or isolated) in g2, its betweenness is 0
      if (!(z %in% vnames_g2) || igraph::degree(g2, v = z) == 0) {
        custom_betw[z] <- 0
        next
      }
      zid <- match(z, vnames_g2)
      total <- 0L

      for (src in assigned_nodes) {
        cs <- assigned_comm[[src]]
        to_nodes <- assigned_nodes[assigned_comm != cs]
        if (!length(to_nodes)) next

        for (t in to_nodes) {
          sp <- suppressWarnings(
            igraph::get.all.shortest.paths(g2, from = src, to = t, mode = if (directed) "out" else "all")
          )
          if (!length(sp$res)) next

          # +1 SE z appare in ALMENO un shortest path (conteggio per-target)
          hit <- FALSE
          for (p in sp$res) { if (any(mid_ids(p) == zid)) { hit <- TRUE; break } }
          if (hit) total <- total + 1L
        }
      }

      # Undirected correction (s–t e t–s sono la stessa coppia)
      if (!directed) total <- total / 2
      custom_betw[z] <- as.numeric(total)
    }
  }
  # ----------------------------------------------------------------------


  # Assemble output (index by node names)
  out <- data.frame(
    node = outside_nodes,
    bridge_strength = bridge_res$`Bridge Strength`[outside_nodes],
    bridge_closeness = bridge_res$`Bridge Closeness`[outside_nodes],
    bridge_betweenness = bridge_res$`Bridge Betweenness`[outside_nodes],  # placeholder
    bridge_expected_influence1 = bridge_res$`Bridge Expected Influence (1-step)`[outside_nodes],
    bridge_expected_influence2 = bridge_res$`Bridge Expected Influence (2-step)`[outside_nodes],
    cluster = full_membership[outside_nodes],
    row.names = NULL
  )

  if (length(custom_betw)) {
    out$bridge_betweenness <- ifelse(out$node %in% names(custom_betw),
                                     custom_betw[out$node],
                                     out$bridge_betweenness)
  }

  return(out)
}
