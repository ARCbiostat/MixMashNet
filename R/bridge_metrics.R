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
#' Jones, P. J., Ma, R., & McNally, R. J. (2021).
#' Bridge Centrality: A Network Approach to Understanding Comorbidity.
#' \emph{Multivariate Behavioral Research}, 56(2), 353–367.
#' \doi{10.1080/00273171.2019.1614898}
#'
#' @importFrom igraph V E vcount as_adjacency_matrix graph_from_adjacency_matrix
#' @importFrom igraph distances all_shortest_paths is_directed
#' @importFrom stats setNames
#' @export
bridge_metrics <- function(g, membership) {

  EPS <- 1e-10

  # ---------------------------------------------------------------------------
  # 1) Vertex names and membership
  # ---------------------------------------------------------------------------
  if (is.null(igraph::V(g)$name)) {
    igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
  }
  igraph::V(g)$name <- as.character(igraph::V(g)$name)
  nodes <- igraph::V(g)$name

  if (is.null(names(membership))) {
    stop("`membership` must be a named vector/factor with names matching V(g)$name.")
  }

  names(membership) <- as.character(names(membership))

  full_membership <- stats::setNames(rep(NA_character_, length(nodes)), nodes)
  keep <- intersect(names(membership), nodes)
  if (length(keep)) {
    full_membership[keep] <- as.character(membership[keep])
  }

  assigned_nodes <- names(full_membership)[!is.na(full_membership)]
  if (!length(assigned_nodes)) {
    return(data.frame(
      node = character(0),
      community = character(0),
      bridge_strength = numeric(0),
      bridge_ei1 = numeric(0),
      bridge_ei2 = numeric(0),
      bridge_betweenness = numeric(0),
      bridge_closeness = numeric(0),
      row.names = NULL
    ))
  }

  assigned_comm <- full_membership[assigned_nodes]

  # ---------------------------------------------------------------------------
  # 2) Weighted adjacency on full graph
  # ---------------------------------------------------------------------------
  W <- suppressWarnings(
    as.matrix(igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  )

  if (!is.matrix(W) || all(is.na(W))) {
    W <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
    W[W != 0] <- 1
  }

  rownames(W) <- colnames(W) <- nodes
  W[is.na(W)] <- 0
  diag(W) <- 0

  # ---------------------------------------------------------------------------
  # 3) Bridge strength and EI1 (full graph, assigned targets only)
  # ---------------------------------------------------------------------------
  bridge_strength <- numeric(length(assigned_nodes))
  bridge_ei1 <- numeric(length(assigned_nodes))
  names(bridge_strength) <- names(bridge_ei1) <- assigned_nodes

  for (i in seq_along(assigned_nodes)) {
    node_i <- assigned_nodes[i]
    comm_i <- assigned_comm[node_i]

    targets <- assigned_nodes[assigned_comm != comm_i]

    if (length(targets)) {
      w_it <- W[node_i, targets, drop = TRUE]
      bridge_strength[i] <- sum(abs(w_it))
      bridge_ei1[i] <- sum(w_it)
    } else {
      bridge_strength[i] <- 0
      bridge_ei1[i] <- 0
    }
  }

  # ---------------------------------------------------------------------------
  # 4) Bridge EI2 (direct + indirect through one intermediate neighbor)
  # ---------------------------------------------------------------------------
  bridge_ei2 <- numeric(length(assigned_nodes))
  names(bridge_ei2) <- assigned_nodes

  for (i in seq_along(assigned_nodes)) {
    node_i <- assigned_nodes[i]
    comm_i <- assigned_comm[node_i]

    targets <- assigned_nodes[assigned_comm != comm_i]

    if (!length(targets)) {
      bridge_ei2[i] <- 0
      next
    }

    # Direct bridge influence
    direct_i <- sum(W[node_i, targets])

    # Indirect bridge influence through one intermediate neighbor k:
    # sum_k W[i,k] * sum_t W[k,t], where t are assigned targets in other communities
    influence_to_targets <- rowSums(W[, targets, drop = FALSE])
    indirect_i <- sum(W[node_i, ] * influence_to_targets)

    bridge_ei2[i] <- direct_i + indirect_i
  }

  # ---------------------------------------------------------------------------
  # 5) Positive-weight graph for closeness and betweenness
  # ---------------------------------------------------------------------------
  W_pos <- W
  W_pos[W_pos <= EPS] <- 0
  diag(W_pos) <- 0

  g_pos <- igraph::graph_from_adjacency_matrix(
    W_pos,
    mode = if (igraph::is_directed(g)) "directed" else "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # distances = 1 / weight
  if (igraph::ecount(g_pos) > 0) {
    igraph::E(g_pos)$dist <- 1 / igraph::E(g_pos)$weight
  }

  # ---------------------------------------------------------------------------
  # 6) Bridge closeness
  # ---------------------------------------------------------------------------
  bridge_closeness <- numeric(length(assigned_nodes))
  names(bridge_closeness) <- assigned_nodes

  for (i in seq_along(assigned_nodes)) {
    node_i <- assigned_nodes[i]
    comm_i <- assigned_comm[node_i]

    targets <- assigned_nodes[assigned_comm != comm_i]
    targets <- intersect(targets, igraph::V(g_pos)$name)

    if (!length(targets) || igraph::ecount(g_pos) == 0) {
      bridge_closeness[i] <- NA_real_
      next
    }

    d_i <- suppressWarnings(
      igraph::distances(
        g_pos,
        v = node_i,
        to = targets,
        weights = igraph::E(g_pos)$dist,
        mode = "all"
      )
    )

    d_i <- as.numeric(d_i)
    d_i <- d_i[is.finite(d_i)]

    bridge_closeness[i] <- if (length(d_i)) 1 / mean(d_i) else NA_real_
  }

  # ---------------------------------------------------------------------------
  # 7) Bridge betweenness
  # ---------------------------------------------------------------------------
  bridge_betweenness <- stats::setNames(
    numeric(length(assigned_nodes)),
    assigned_nodes
  )

  bridge_betweenness <- stats::setNames(numeric(length(assigned_nodes)), assigned_nodes)

  # positive-weight graph for path computations
  W_pos <- W
  W_pos[W_pos <= EPS] <- 0
  diag(W_pos) <- 0

  g_pos <- igraph::graph_from_adjacency_matrix(
    W_pos,
    mode = if (igraph::is_directed(g)) "directed" else "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  if (igraph::ecount(g_pos) > 0 && length(assigned_nodes) > 1L) {

    # convert weights to distances
    igraph::E(g_pos)$dist <- 1 / igraph::E(g_pos)$weight

    vertex_names <- igraph::V(g_pos)$name

    # keep only assigned nodes that are actually in g_pos
    S <- intersect(assigned_nodes, vertex_names)
    commS <- assigned_comm[S]

    # helper: remove endpoints from a path and convert ids to names
    middle_vertices <- function(path_obj) {
      ids <- as.vector(path_obj)
      if (length(ids) <= 2L) {
        character(0)
      } else {
        vertex_names[ids[-c(1L, length(ids))]]
      }
    }

    # for undirected graphs, count each pair only once
    if (!igraph::is_directed(g_pos)) {
      ord <- order(match(S, vertex_names))
      S <- S[ord]
      commS <- commS[S]

      for (i in seq_len(length(S) - 1L)) {
        src <- S[i]
        src_comm <- commS[[src]]

        j_idx <- which(commS != src_comm)
        j_idx <- j_idx[j_idx > i]
        if (!length(j_idx)) next

        to_nodes <- S[j_idx]

        sp <- suppressWarnings(
          igraph::all_shortest_paths(
            g_pos,
            from = src,
            to = to_nodes,
            weights = igraph::E(g_pos)$dist,
            mode = "all"
          )
        )

        if (!length(sp$res)) next

        mids <- unlist(lapply(sp$res, middle_vertices), use.names = FALSE)
        if (!length(mids)) next

        # count only assigned intermediate nodes
        mids <- mids[mids %in% S]
        if (!length(mids)) next

        tb <- table(factor(mids, levels = S))
        bridge_betweenness[names(tb)] <- bridge_betweenness[names(tb)] + as.numeric(tb)
      }

    } else {
      # directed case: keep ordered source-target pairs
      for (src in S) {
        src_comm <- assigned_comm[[src]]

        to_nodes <- S[commS != src_comm]
        if (!length(to_nodes)) next

        sp <- suppressWarnings(
          igraph::all_shortest_paths(
            g_pos,
            from = src,
            to = to_nodes,
            weights = igraph::E(g_pos)$dist,
            mode = "out"
          )
        )

        if (!length(sp$res)) next

        mids <- unlist(lapply(sp$res, middle_vertices), use.names = FALSE)
        if (!length(mids)) next

        mids <- mids[mids %in% S]
        if (!length(mids)) next

        tb <- table(factor(mids, levels = S))
        bridge_betweenness[names(tb)] <- bridge_betweenness[names(tb)] + as.numeric(tb)
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 8) Output
  # ---------------------------------------------------------------------------
  out <- data.frame(
    node = assigned_nodes,
    community = unname(assigned_comm),
    bridge_strength = unname(bridge_strength[assigned_nodes]),
    bridge_ei1 = unname(bridge_ei1[assigned_nodes]),
    bridge_ei2 = unname(bridge_ei2[assigned_nodes]),
    bridge_betweenness = unname(bridge_betweenness[assigned_nodes]),
    bridge_closeness = unname(bridge_closeness[assigned_nodes]),
    row.names = NULL
  )

  out
}
