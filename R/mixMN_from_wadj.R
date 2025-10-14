' Build MixMashNet metrics from a signed weighted adjacency matrix
#'
#' Computes clustering, centrality, and bridge metrics from a signed
#' weighted adjacency (wadj) selecting nodes for graph vs clustering.
#'
#' @param wadj_signed Square numeric matrix with signed weights.
#' @param nodes Character vector of node names (matching row/colnames).
#' @param exclude_from_graph Nodes to exclude from the graph metrics.
#' @param exclude_from_cluster Nodes to exclude only from clustering.
#' @param cluster_method One of "louvain","fast_greedy","infomap",
#'   "walktrap","edge_betweenness".
#' @param reps Integer, bootstrap reps (currently unused here).
#' @param seed_boot Optional integer seed.
#' @param treat_singletons_as_excluded Logical; treat size-1 groups as excluded.
#'
#' @return A list of class \code{mixMN_fit} with clustering, palettes,
#'   centrality/bridge metrics, and edge weights.
#'
#' @importFrom igraph cluster_louvain cluster_fast_greedy cluster_infomap cluster_walktrap cluster_edge_betweenness
#' @importFrom qgraph centrality
#' @importFrom grDevices hcl
#' @importFrom stats setNames
#' @importFrom utils combn
#' @export
mixMN_from_wadj2 <- function(
    wadj_signed,
    nodes,
    exclude_from_graph = NULL,
    exclude_from_cluster = NULL,
    cluster_method = c("louvain", "fast_greedy", "infomap", "walktrap", "edge_betweenness"),
    reps = 0,
    seed_boot = NULL,
    treat_singletons_as_excluded = FALSE
) {
  cluster_method <- match.arg(cluster_method)

  all_nodes <- nodes

  # --- Keep nodes for graph/cluster scopes ---
  keep_nodes_graph   <- setdiff(all_nodes, exclude_from_graph)
  keep_nodes_cluster <- setdiff(all_nodes, unique(c(exclude_from_graph, exclude_from_cluster)))

  wadj_signed_graph   <- wadj_signed[keep_nodes_graph,   keep_nodes_graph, drop = FALSE]
  wadj_signed_cluster <- wadj_signed[keep_nodes_cluster, keep_nodes_cluster, drop = FALSE]

  # --- Graphs for clustering (abs) ---
  g_cluster <- igraph::graph_from_adjacency_matrix(abs(wadj_signed_cluster),
                                                   mode = "undirected", weighted = TRUE, diag = FALSE)
  if (cluster_method %in% c("infomap", "edge_betweenness", "walktrap")) {
    g_cluster <- igraph::simplify(g_cluster, remove.multiple = TRUE, remove.loops = TRUE)
  }

  # --- Clustering ---
  cluster_fun <- function(graph) {
    switch(cluster_method,
           louvain          = igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight),
           fast_greedy      = igraph::cluster_fast_greedy(graph, weights = igraph::E(graph)$weight),
           infomap          = igraph::cluster_infomap(graph),
           walktrap         = igraph::cluster_walktrap(graph, weights = igraph::E(graph)$weight),
           edge_betweenness = igraph::cluster_edge_betweenness(graph, weights = igraph::E(graph)$weight))
  }

  .base_hue_for_nodes <- function(nodes_ref) {
    s <- paste(sort(nodes_ref), collapse = "|")
    as.numeric(sum(utf8ToInt(s)) %% 360)
  }

  .make_palette_layer <- function(unique_clusters, nodes_ref,
                                  h_spread = 120,
                                  c_range = c(50, 80),
                                  l_range = c(60, 75)) {
    k <- length(unique_clusters)
    if (k == 0) return(stats::setNames(character(0), unique_clusters))
    base_hue <- .base_hue_for_nodes(nodes_ref)
    hues <- if (k == 1) base_hue else (base_hue + seq(-h_spread/2, h_spread/2, length.out = k)) %% 360
    c_vals <- if (k > 1) seq(c_range[1], c_range[2], length.out = k) else rep(mean(c_range), k)
    l_vals <- if (k > 1) seq(l_range[1], l_range[2], length.out = k) else rep(mean(l_range), k)
    pal <- grDevices::hcl(h = hues, c = c_vals, l = l_vals)
    stats::setNames(pal, unique_clusters)
  }

  # Membership on abs graph (as in the rest of the code base)
  original_membership <- tryCatch({
    m <- cluster_fun(g_cluster)$membership
    stats::setNames(m, keep_nodes_cluster)
  }, error = function(e) rep(NA_integer_, length(keep_nodes_cluster)))

  sizes <- table(original_membership)
  if (treat_singletons_as_excluded) {
    valid_labels <- as.integer(names(sizes[sizes >= 2L]))
    in_valid <- original_membership %in% valid_labels
    unique_clusters <- sort(valid_labels)
    palette_clusters <- .make_palette_layer(unique_clusters, nodes_ref = keep_nodes_graph)
    groups <- as.factor(original_membership[in_valid])
    names(groups) <- keep_nodes_cluster[in_valid]
  } else {
    unique_clusters <- sort(unique(stats::na.omit(original_membership)))
    palette_clusters <- .make_palette_layer(unique_clusters, nodes_ref = keep_nodes_graph)
    groups <- as.factor(original_membership)
    names(groups) <- keep_nodes_cluster
  }

  # --- Classical centralities ---
  # qgraph::centrality on the SIGNED matrix -> strength (sum of abs by qgraph) + EI1 (signed)
  cent_true <- qgraph::centrality(wadj_signed_graph)

  # Distances for closeness/betweenness on ABS with dist = 1/|w|
  g_dist <- igraph::graph_from_adjacency_matrix(abs(wadj_signed_graph),
                                                mode = "undirected", weighted = TRUE, diag = FALSE)
  if (igraph::ecount(g_dist) > 0) igraph::E(g_dist)$dist <- 1 / igraph::E(g_dist)$weight

  .harmonic_closeness <- function(g) {
    if (igraph::ecount(g) == 0) return(setNames(rep(0, igraph::vcount(g)), igraph::V(g)$name))
    D <- igraph::distances(g, weights = igraph::E(g)$dist)
    diag(D) <- NA
    cl <- rowSums(1 / D, na.rm = TRUE) / (nrow(D) - 1)
    names(cl) <- igraph::V(g)$name
    cl
  }

  closeness_vec <- .harmonic_closeness(g_dist)
  betweenness_vec <- if (igraph::ecount(g_dist) > 0) {
    b <- igraph::betweenness(g_dist, weights = igraph::E(g_dist)$dist,
                             directed = FALSE, normalized = FALSE)
    names(b) <- igraph::V(g_dist)$name
    b[keep_nodes_graph]
  } else setNames(rep(0, length(keep_nodes_graph)), keep_nodes_graph)

  # --- Bridge metrics (split ABS vs SIGNED) ---
  # ABS graph: used for bridge_strength / bridge_betweenness / bridge_closeness
  g_bridge_abs <- igraph::graph_from_adjacency_matrix(abs(wadj_signed_graph),
                                                      mode = "undirected", weighted = TRUE, diag = FALSE)
  g_bridge_abs <- igraph::simplify(g_bridge_abs, remove.multiple = TRUE, remove.loops = TRUE)

  # SIGNED graph: used for bridge EI1 / EI2 (needs signs)
  g_bridge_signed <- igraph::graph_from_adjacency_matrix(wadj_signed_graph,
                                                         mode = "undirected", weighted = TRUE, diag = FALSE)
  g_bridge_signed <- igraph::simplify(g_bridge_signed, remove.multiple = TRUE, remove.loops = TRUE)

  # 1) ABS → strength/clos/BTW
  bridge_abs <- tryCatch({
    b <- bridge_metrics(g_bridge_abs, membership = groups)
    b[match(keep_nodes_graph, b$node),
      c("node","bridge_strength","bridge_betweenness","bridge_closeness")]
  }, error = function(e) {
    data.frame(node = keep_nodes_graph,
               bridge_strength = NA_real_,
               bridge_betweenness = NA_real_,
               bridge_closeness = NA_real_)
  })

  # 2) SIGNED → EI1/EI2
  bridge_signed <- tryCatch({
    b <- bridge_metrics(g_bridge_signed, membership = groups)
    # some implementations call them bridge_ei1/ei2; keep those
    b[match(keep_nodes_graph, b$node), c("node","bridge_ei1","bridge_ei2")]
  }, error = function(e) {
    data.frame(node = keep_nodes_graph,
               bridge_ei1 = NA_real_,
               bridge_ei2 = NA_real_)
  })

  # merge by node
  bridge_true <- merge(bridge_abs, bridge_signed, by = "node", all = TRUE)

  # --- Bridge (excluded) split ABS vs SIGNED as well ---
  # ABS excluded → strength / closeness / betweenness
  bridge_excl_abs <- tryCatch({
    b <- bridge_metrics_excluded(g_bridge_abs, membership = groups)
    b
  }, error = function(e) {
    data.frame(node = keep_nodes_graph,
               bridge_strength = NA_real_,
               bridge_betweenness = NA_real_,
               bridge_closeness = NA_real_)
  })

  # SIGNED excluded → EI1 / EI2 (column names can vary across implementations)
  bridge_excl_sgn <- tryCatch({
    b <- bridge_metrics_excluded(g_bridge_signed, membership = groups)
    # normalize possible column names
    if (!is.null(b$bridge_ei1) || !is.null(b$bridge_expected_influence1)) {
      ei1 <- if (!is.null(b$bridge_ei1)) b$bridge_ei1 else b$bridge_expected_influence1
    } else ei1 <- NA_real_
    if (!is.null(b$bridge_ei2) || !is.null(b$bridge_expected_influence2)) {
      ei2 <- if (!is.null(b$bridge_ei2)) b$bridge_ei2 else b$bridge_expected_influence2
    } else ei2 <- NA_real_
    data.frame(node = b$node, bridge_ei1_excluded = ei1, bridge_ei2_excluded = ei2)
  }, error = function(e) {
    data.frame(node = keep_nodes_graph,
               bridge_ei1_excluded = NA_real_,
               bridge_ei2_excluded = NA_real_)
  })

  bridge_excl_merged <- merge(
    bridge_excl_abs[, c("node","bridge_strength","bridge_betweenness","bridge_closeness")],
    bridge_excl_sgn, by = "node", all = TRUE
  )

  # --- Assemble centrality table (true) ---
  centrality_true_df <- data.frame(
    node = keep_nodes_graph,
    strength   = cent_true$OutDegree,
    ei1        = cent_true$OutExpectedInfluence,
    closeness  = closeness_vec,
    betweenness= betweenness_vec,
    # bridge (included)
    bridge_strength    = bridge_true$bridge_strength[match(keep_nodes_graph, bridge_true$node)],
    bridge_betweenness = bridge_true$bridge_betweenness[match(keep_nodes_graph, bridge_true$node)],
    bridge_closeness   = bridge_true$bridge_closeness[match(keep_nodes_graph, bridge_true$node)],
    bridge_ei1         = bridge_true$bridge_ei1[match(keep_nodes_graph, bridge_true$node)],
    bridge_ei2         = bridge_true$bridge_ei2[match(keep_nodes_graph, bridge_true$node)],
    # bridge (excluded)
    bridge_strength_excluded    = bridge_excl_merged$bridge_strength[match(keep_nodes_graph, bridge_excl_merged$node)],
    bridge_betweenness_excluded = bridge_excl_merged$bridge_betweenness[match(keep_nodes_graph, bridge_excl_merged$node)],
    bridge_closeness_excluded   = bridge_excl_merged$bridge_closeness[match(keep_nodes_graph, bridge_excl_merged$node)],
    bridge_ei1_excluded         = bridge_excl_merged$bridge_ei1_excluded[match(keep_nodes_graph, bridge_excl_merged$node)],
    bridge_ei2_excluded         = bridge_excl_merged$bridge_ei2_excluded[match(keep_nodes_graph, bridge_excl_merged$node)]
  )

  # --- Edges (signed) ---
  edge_names <- utils::combn(keep_nodes_graph, 2, FUN = function(x) paste(x[1], x[2], sep = "--"))
  edge_mat_true <- wadj_signed_graph[lower.tri(wadj_signed_graph)]
  names(edge_mat_true) <- edge_names
  edges_true_df <- data.frame(edge = edge_names, weight = edge_mat_true, row.names = NULL)

  # --- Return object skeleton (no bootstrap here) ---
  out <- list(
    original_membership = original_membership,
    groups              = groups,
    community_palette   = palette_clusters,
    boot_memberships    = list(),
    reps                = reps,
    cluster_method      = cluster_method,
    mgm_model           = NULL,

    ci_results          = NULL,

    centrality_true     = centrality_true_df,
    strength_boot       = NULL,
    ei1_boot            = NULL,
    closeness_boot      = NULL,
    betweenness_boot    = NULL,

    bridge_strength_boot        = NULL,
    bridge_betweenness_boot     = NULL,
    bridge_closeness_boot       = NULL,
    bridge_ei1_boot             = NULL,
    bridge_ei2_boot             = NULL,
    bridge_strength_excl_boot   = NULL,
    bridge_betweenness_excl_boot= NULL,
    bridge_closeness_excl_boot  = NULL,
    bridge_ei1_excl_boot        = NULL,
    bridge_ei2_excl_boot        = NULL,

    edges_true          = edges_true_df,
    edge_boot_mat       = NULL,

    keep_nodes_graph    = keep_nodes_graph,
    keep_nodes_cluster  = keep_nodes_cluster,

    community_scores_obj  = NULL,
    community_scores_df   = NULL,
    community_scores_ci   = NULL,
    community_scores_boot = NULL,

    excluded_score        = NULL,
    excluded_score_ci     = NULL,
    excluded_score_boot   = NULL
  )

  class(out) <- c("mixMN_fit")
  return(out)
}
