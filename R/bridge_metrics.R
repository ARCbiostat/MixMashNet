#' Bridge metrics for nodes across communities
#'
#' Computes bridge centrality measures for nodes with an assigned community.
#' This function is used internally by \code{mixMN()} and \code{multimixMN()}.
#' Specifically, the function computes bridge strength as the sum of absolute
#' edge weights connecting a node to nodes in other communities; bridge expected
#' influence of order one (EI1) as the signed sum of direct connections to nodes
#' in other communities; bridge expected influence of order two (EI2) as the
#' signed influence that propagates indirectly to nodes in other communities
#' via one intermediate neighbor (i.e., through paths of length two);
#' bridge betweenness as the number of times a node lies on shortest paths
#' between nodes belonging to different communities; and bridge closeness as the
#' inverse of the mean shortest-path distance to nodes in other communities.
#'
#' @param g An igraph object with edge attribute \code{weight}.
#' @param membership Named vector/factor of community labels for a subset of nodes (names must match \code{V(g)$name}).
#' @return A data.frame with columns: \code{node}, \code{community}, \code{bridge_strength}, \code{bridge_ei1},
#'        \code{bridge_ei2}, \code{bridge_betweenness}, \code{bridge_closeness}.
#' @details
#' Bridge betweenness and closeness are computed on the positive-weight subgraph
#' only, with weights converted to distances as \eqn{d = 1/w}.
#'
#' @references
#'
#' Jones, P. J. (2025). \pkg{networktools}: Tools for identifying important nodes
#' in networks. R package version 1.6.1.
#' \url{https://github.com/paytonjjones/networktools}
#'
#' Jones, P. J., Ma, R., & McNally, R. J. (2021).
#' Bridge Centrality: A Network Approach to Understanding Comorbidity.
#' \emph{Multivariate Behavioral Research}, 56(2), 353–367.
#' \doi{10.1080/00273171.2019.1614898}
#'
#' @importFrom igraph V E vcount as_adjacency_matrix delete_edges get.all.shortest.paths distances is_directed
#' @importFrom stats setNames
#' @export
bridge_metrics <- function(g, membership) {

  # --- Ensure vertex names ---
  if (is.null(igraph::V(g)$name)) {
    igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
  }
  nodes <- igraph::V(g)$name

  # --- Full membership over all nodes (NA if unassigned) ---
  full_membership <- stats::setNames(rep(NA, length(nodes)), nodes)
  full_membership[names(membership)] <- membership

  communities <- full_membership
  assigned_nodes <- names(communities[!is.na(communities)])
  assigned_communities <- communities[!is.na(communities)]

  # --- Weighted adjacency (all nodes) ---
  adj <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  # Fallback in case weights are missing
  if (all(is.na(adj))) {
    adj <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  }
  rownames(adj) <- colnames(adj) <- nodes

  # --- EI1 & influence helpers ---
  expectedInfBridge <- function(node_of_interest, network, nodes, communities) {
    comm_int <- communities[node_of_interest]
    other_comm <- assigned_nodes[assigned_communities != comm_int]
    included_nodes <- c(node_of_interest, other_comm)
    new_net <- network[included_nodes, included_nodes, drop = FALSE]
    new_net[node_of_interest, node_of_interest] <- 0
    sum(new_net[node_of_interest, ])
  }

  influence_on_comm_j <- function(node_of_interest, network, nodes, communities, j) {
    included_nodes <- unique(c(
      node_of_interest,
      assigned_nodes[assigned_communities == unique(assigned_communities)[j]]
    ))
    new_net <- network[included_nodes, included_nodes, drop = FALSE]
    new_net[node_of_interest, node_of_interest] <- 0
    sum(new_net[node_of_interest, ])
  }

  # --- Precompute influence for all communities (for EI2) ---
  infcomm <- vector("list", length(unique(assigned_communities)))
  for (j in seq_along(infcomm)) {
    infcomm[[j]] <- sapply(nodes, function(node) {
      influence_on_comm_j(node, adj, nodes, assigned_communities, j)
    })
    names(infcomm[[j]]) <- nodes
  }

  # --- Bridge Strength & EI1 ---
  bridge_strength <- sapply(assigned_nodes, function(node) {
    node_comm <- assigned_communities[node]
    target_nodes <- assigned_nodes[assigned_communities != node_comm]
    if (length(target_nodes) == 0) 0 else sum(abs(adj[node, target_nodes, drop = FALSE]))
  })

  bridge_ei1 <- sapply(assigned_nodes, function(node) {
    expectedInfBridge(node, adj, nodes, assigned_communities)
  })

  # --- Bridge EI2 ---
  bridge_ei2 <- sapply(assigned_nodes, function(node_of_interest) {
    comm_int <- assigned_communities[node_of_interest]
    non_comm_vec <- unique(assigned_communities)[unique(assigned_communities) != comm_int]
    ei2_vec <- numeric()
    for (i in seq_along(non_comm_vec)) {
      comm_index <- which(unique(assigned_communities) == non_comm_vec[i])
      ei2_wns <- sweep(adj, MARGIN = 2, infcomm[[comm_index]], `*`)
      if (node_of_interest %in% rownames(ei2_wns)) {
        ei2_vec[i] <- sum(ei2_wns[node_of_interest, ])
      } else {
        ei2_vec[i] <- 0
      }
    }
    ei1_node <- expectedInfBridge(node_of_interest, adj, nodes, assigned_communities)
    sum(ei2_vec) + ei1_node
  })

  # --- Graph for Betweenness/Closeness ---
  g_inv <- igraph::delete_edges(g, which(is.na(E(g)$weight) | abs(E(g)$weight) <= 1e-10))
  igraph::E(g_inv)$weight <- 1 / igraph::E(g_inv)$weight
  g_pos <- igraph::delete_edges(g_inv, which(igraph::E(g_inv)$weight < 0))

  # --- Bridge Betweenness ---
  delete.ends <- function(x) {
    nodes_path <- names(igraph::V(g_pos))
    nodes_path[as.vector(x)[-c(1, length(x))]]
  }

  short.bridge.mid.paths <- function(x) {
    b <- suppressWarnings(igraph::get.all.shortest.paths(
      g_pos,
      from = assigned_nodes[x],
      to   = assigned_nodes[assigned_communities != assigned_communities[x]],
      mode = "all"
    ))
    intermediates <- unlist(lapply(b$res, delete.ends))
    return(intermediates)
  }

  betweenness.df <- as.data.frame(table(factor(unlist(
    lapply(seq_along(assigned_nodes), short.bridge.mid.paths)
  ), levels = assigned_nodes)))
  bridge_betweenness <- betweenness.df[, 2]

  # Undirected correction
  if (!igraph::is_directed(g_pos)) {
    bridge_betweenness <- bridge_betweenness / 2
  }
  names(bridge_betweenness) <- assigned_nodes

  # --- Bridge Closeness ---
  bridge_closeness <- numeric(length(assigned_nodes))
  for (i in seq_along(assigned_nodes)) {
    node <- assigned_nodes[i]
    node_comm <- assigned_communities[node]
    target_nodes <- assigned_nodes[assigned_communities != node_comm]
    target_nodes <- target_nodes[target_nodes %in% igraph::V(g_pos)$name]

    if (length(target_nodes) == 0) {
      bridge_closeness[i] <- NA
    } else {
      dists <- igraph::distances(g_pos, v = node, to = target_nodes, mode = "all")
      if (all(is.infinite(dists))) {
        bridge_closeness[i] <- NA
      } else {
        bridge_closeness[i] <- 1 / mean(dists[is.finite(dists)])
      }
    }
  }

  # --- Output ---
  data.frame(
    node = assigned_nodes,
    community = assigned_communities,
    bridge_strength = bridge_strength,
    bridge_ei1 = bridge_ei1,
    bridge_ei2 = bridge_ei2,
    bridge_betweenness = bridge_betweenness,
    bridge_closeness = bridge_closeness,
    row.names = NULL
  )
}
