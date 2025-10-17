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
  # --- Fixed numerical tolerance for edge filtering ---
  EPS <- 1e-10

  # --- 0) Ensure vertex names are present and are character -------------------
  if (is.null(igraph::V(g)$name)) {
    warning("Vertices had no names: assigning sequential character names.")
    igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
  }
  igraph::V(g)$name <- as.character(igraph::V(g)$name)
  node_names <- igraph::V(g)$name

  # --- 1) Build full membership over all nodes; unassigned -> "Z" ------------
  if (is.null(names(membership))) {
    stop("`membership` must be a *named* vector/factor with names matching V(g)$name.")
  }
  names(membership) <- as.character(names(membership))

  full_membership <- stats::setNames(rep("Z", length(node_names)), node_names)
  in_graph <- intersect(names(membership), node_names)
  if (length(in_graph)) {
    full_membership[in_graph] <- as.character(membership[in_graph])
  }

  # --- 2) Get (weighted) adjacency; fallback to unweighted if needed ----------
  adj <- suppressWarnings(
    as.matrix(igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  )
  if (!is.matrix(adj) || all(is.na(adj))) {
    # No usable weights: use unweighted adjacency (1/0)
    adj <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
    adj[adj != 0] <- 1
  }
  rownames(adj) <- colnames(adj) <- node_names
  adj[is.na(adj)] <- 0
  diag(adj) <- 0

  # --- 3) Baseline bridge metrics from networktools ---------------------------
  bridge_res <- networktools::bridge(network = adj, communities = full_membership)

  # --- 4) Identify excluded ("Z") nodes --------------------------------------
  outside_nodes <- names(full_membership)[full_membership == "Z"]
  if (!length(outside_nodes)) {
    # Nothing to compute; return a typed empty frame
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

  # --- 5) Custom Bridge Betweenness for Z nodes -------------------------------
  # Build a path graph starting from |w|:
  # - take absolute weights
  # - drop edges with |w| <= EPS
  # - invert remaining weights (1/|w|) to represent distances
  # - BUT compute shortest paths as UNWEIGHTED counts on this edge set
  Wabs <- suppressWarnings(
    as.matrix(igraph::as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  )
  rownames(Wabs) <- colnames(Wabs) <- node_names
  Wabs[is.na(Wabs)] <- 0
  diag(Wabs) <- 0
  Wabs <- abs(Wabs)

  g_absw <- igraph::graph_from_adjacency_matrix(
    Wabs, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  drop_e <- which(igraph::E(g_absw)$weight <= EPS | is.na(igraph::E(g_absw)$weight))
  g_inv  <- if (length(drop_e)) igraph::delete_edges(g_absw, drop_e) else g_absw
  if (igraph::ecount(g_inv) > 0) igraph::E(g_inv)$weight <- 1 / igraph::E(g_inv)$weight

  # Keep names consistent
  g_pos <- g_inv
  igraph::V(g_pos)$name <- node_names

  # Assigned endpoints (S) vs excluded intermediates (Z)
  S <- intersect(names(membership), igraph::V(g_pos)$name)    # assigned endpoints
  Z <- setdiff(igraph::V(g_pos)$name, S)                      # excluded candidates

  custom_betw <- stats::setNames(numeric(length(Z)), Z)

  if (length(S) && length(Z) && igraph::ecount(g_pos) > 0) {
    # Community labels for assigned endpoints (as plain character)
    commS <- stats::setNames(as.character(membership[S]), S)

    # Order S according to the vertex order in g_pos; iterate i<j
    o <- order(match(S, igraph::V(g_pos)$name))
    S <- S[o]
    commS <- commS[S]

    nodes_path <- igraph::V(g_pos)$name

    # Helpers to extract intermediate vertices (exclude endpoints)
    mid_ids <- function(path_ids) {
      ids <- as.vector(path_ids)
      if (length(ids) <= 2L) integer(0) else ids[-c(1L, length(ids))]
    }
    mids_to_names <- function(x) nodes_path[mid_ids(x)]

    # Enumerate unweighted shortest paths between assigned endpoints
    # that belong to different communities, and count Z intermediates.
    for (i in seq_len(length(S) - 1L)) {
      src <- S[i]
      Ci  <- commS[[src]]

      # Only destinations in a different community and with j > i
      j_idx <- which(commS != Ci)
      if (!length(j_idx)) next
      j_idx <- j_idx[j_idx > i]
      if (!length(j_idx)) next

      to_nodes <- S[j_idx]

      sp <- suppressWarnings(igraph::get.all.shortest.paths(
        g_pos, from = src, to = to_nodes, mode = "all"
      ))
      if (!length(sp$res)) next

      inter <- unlist(lapply(sp$res, mids_to_names), use.names = FALSE)
      if (!length(inter)) next

      inter_Z <- inter[inter %in% Z]
      if (!length(inter_Z)) next

      tb <- table(factor(inter_Z, levels = Z))
      custom_betw[names(tb)] <- custom_betw[names(tb)] + as.numeric(tb)
    }
  }

  # Override Bridge Betweenness ONLY for excluded nodes
  bridge_res$`Bridge Betweenness`[outside_nodes] <- custom_betw[outside_nodes]

  # --- 6) Assemble output for excluded nodes ----------------------------------
  out <- data.frame(
    node = outside_nodes,
    bridge_strength              = bridge_res$`Bridge Strength`[outside_nodes],
    bridge_closeness             = bridge_res$`Bridge Closeness`[outside_nodes],
    bridge_betweenness           = bridge_res$`Bridge Betweenness`[outside_nodes],
    bridge_expected_influence1   = bridge_res$`Bridge Expected Influence (1-step)`[outside_nodes],
    bridge_expected_influence2   = bridge_res$`Bridge Expected Influence (2-step)`[outside_nodes],
    cluster = full_membership[outside_nodes],
    row.names = NULL
  )

  out
}
