#' Build a single-layer mixMN_fit from a signed adjacency (wadj * signs)
#'
#' @description
#' Constructs a **single-layer** \code{mixMN_fit} object starting from a signed
#' weighted adjacency matrix (e.g., \code{wadj * signs} from an MGM). The function
#' supports excluding nodes from the graph and/or from clustering, computes
#' community structure on the clustering graph, and derives per-node centrality
#' (strength, EI1, closeness, betweenness) and bridge metrics (including an
#' "outside/excluded" variant). It also returns the intra-layer edge list with
#' true weights. Bootstrap containers are initialized but left \code{NULL} here
#' (bootstrap is handled upstream).
#'
#' @param wadj_signed Square numeric matrix of signed weights (nodes × nodes).
#'   Row/column names must correspond to \code{nodes}.
#' @param nodes Character vector with the node names (order must match the
#'   variables represented in \code{wadj_signed}).
#' @param exclude_from_graph Character vector of nodes to exclude from the
#'   **graph** (no edges/centralities for these nodes).
#' @param exclude_from_cluster Character vector of nodes excluded from
#'   **clustering** (no membership), but still allowed in the graph.
#' @param cluster_method Community detection algorithm. One of
#'   \code{c("louvain","fast_greedy","infomap","walktrap","edge_betweenness")}.
#' @param reps Integer; number of bootstrap iterations to be stored in the
#'   output object (metadata only here). Default \code{0}.
#' @param seed_boot Optional integer; placeholder for upstream bootstrap seed metadata.
#' @param treat_singletons_as_excluded Logical; if \code{TRUE}, communities of
#'   size 1 (singletons) are treated as "excluded" for community-dependent metrics
#'   and are omitted from the empirical membership used for plotting/stability.
#'
#' @details
#' Two graphs are built:
#' \enumerate{
#'   \item \strong{g_graph}: for centrality/bridge metrics (uses \code{abs(wadj_signed)}).
#'   \item \strong{g_cluster}: for community detection (same weights; simplified for
#'         methods requiring it).
#' }
#' Membership is computed on \strong{g_cluster}. If
#' \code{treat_singletons_as_excluded = TRUE}, membership labels with size 1 are
#' dropped from \code{groups} (kept in \code{original_membership} for reference).
#'
#' Centrality includes:
#' \itemize{
#'   \item Strength = \code{qgraph::centrality(wadj_signed_graph)$OutDegree}
#'   \item EI1 = \code{qgraph::centrality(wadj_signed_graph)$OutExpectedInfluence}
#'   \item Closeness = harmonic closeness on a distance graph where
#'         \code{dist = 1 / weight}
#'   \item Betweenness = igraph betweenness on the same distance graph
#' }
#'
#' Bridge metrics are computed via package-internal helpers:
#' \code{bridge_metrics()} (on included nodes) and
#' \code{bridge_metrics_excluded()} (for nodes considered outside/excluded).
#'
#' @return An object of class \code{mixMN_fit} with fields:
#' \describe{
#'   \item{\code{original_membership}}{Integer membership for clustering nodes (named).}
#'   \item{\code{groups}}{Factor membership used for downstream plots/IS; may omit singletons if requested.}
#'   \item{\code{community_palette}}{Named character vector of HEX colors for communities.}
#'   \item{\code{boot_memberships}}{Empty list (bootstrap to be populated upstream).}
#'   \item{\code{reps}, \code{cluster_method}}{Metadata echoing inputs.}
#'   \item{\code{mgm_model}}{Always \code{NULL} here (no model fit inside).}
#'   \item{\code{ci_results}}{Always \code{NULL} here (CIs computed upstream).}
#'   \item{\code{centrality_true}}{Data frame with per-node strength, EI1, closeness, betweenness,
#'     plus bridge metrics (included and excluded variants).}
#'   \item{\code{strength_boot}, \code{ei1_boot}, \code{closeness_boot}, \code{betweenness_boot}}{All \code{NULL} here.}
#'   \item{\code{bridge_*_boot}, \code{bridge_*_excl_boot}}{All \code{NULL} here.}
#'   \item{\code{edges_true}}{Data frame of lower-triangular intra-layer edges with weights.}
#'   \item{\code{edge_boot_mat}}{\code{NULL} (bootstrap stored upstream).}
#'   \item{\code{keep_nodes_graph}, \code{keep_nodes_cluster}}{Vectors with nodes retained per role.}
#' }
#'
#' @section Dependencies:
#' Uses internal helpers \code{bridge_metrics()} and
#' \code{bridge_metrics_excluded()}.
#'
#' @seealso \code{\link{mixMN_multi}}, \code{\link{membershipStab}},
#'   \code{\link{membershipStab_plot}}
#'
#' @importFrom igraph graph_from_adjacency_matrix simplify ecount E V distances betweenness
#' @importFrom igraph cluster_louvain cluster_fast_greedy cluster_infomap cluster_walktrap cluster_edge_betweenness
#' @importFrom qgraph centrality
#' @importFrom grDevices hcl
#' @importFrom stats setNames
#' @importFrom utils combn
#' @export
mixMN_from_wadj <- function(
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

  # --- Keep nodes ---
  keep_nodes_graph   <- setdiff(all_nodes, exclude_from_graph)
  keep_nodes_cluster <- setdiff(all_nodes, unique(c(exclude_from_graph, exclude_from_cluster)))

  wadj_signed_graph   <- wadj_signed[keep_nodes_graph,   keep_nodes_graph, drop = FALSE]
  wadj_signed_cluster <- wadj_signed[keep_nodes_cluster, keep_nodes_cluster, drop = FALSE]

  # --- Graphs ---
  g_graph   <- igraph::graph_from_adjacency_matrix(abs(wadj_signed_graph),
                                                   mode = "undirected", weighted = TRUE, diag = FALSE)
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

  # --- Centrality ---
  cent_true <- qgraph::centrality(wadj_signed_graph)

  # Distance graph for closeness/betweenness
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

  # --- Bridge metrics ---
  bridge_true <- bridge_metrics(g_graph, membership = groups)
  bridge_true <- bridge_true[match(keep_nodes_graph, bridge_true$node), ]
  bridge_outside_true <- bridge_metrics_excluded(g_graph, membership = groups)

  centrality_true_df <- data.frame(
    node = keep_nodes_graph,
    strength = cent_true$OutDegree,
    ei1 = cent_true$OutExpectedInfluence,
    closeness = closeness_vec,
    betweenness = betweenness_vec,
    bridge_strength = bridge_true$bridge_strength,
    bridge_betweenness = bridge_true$bridge_betweenness,
    bridge_closeness = bridge_true$bridge_closeness,
    bridge_ei1 = bridge_true$bridge_ei1,
    bridge_ei2 = bridge_true$bridge_ei2
  )

  bridge_excluded_df <- data.frame(
    node = keep_nodes_graph,
    bridge_strength_excluded = NA_real_,
    bridge_betweenness_excluded = NA_real_,
    bridge_closeness_excluded = NA_real_,
    bridge_ei1_excluded = NA_real_,
    bridge_ei2_excluded = NA_real_
  )
  idx_out <- match(bridge_outside_true$node, keep_nodes_graph)
  bridge_excluded_df[idx_out, -1] <- bridge_outside_true[, c(
    "bridge_strength", "bridge_betweenness", "bridge_closeness",
    "bridge_expected_influence1", "bridge_expected_influence2"
  )]
  centrality_true_df <- cbind(centrality_true_df, bridge_excluded_df[, -1])

  # --- Edges ---
  edge_names <- utils::combn(keep_nodes_graph, 2, FUN = function(x) paste(x[1], x[2], sep = "--"))
  edge_mat_true <- wadj_signed_graph[lower.tri(wadj_signed_graph)]
  names(edge_mat_true) <- edge_names
  edges_true_df <- data.frame(edge = edge_names, weight = edge_mat_true, row.names = NULL)

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
    bridge_strength_boot= NULL,
    bridge_betweenness_boot = NULL,
    bridge_closeness_boot   = NULL,
    bridge_ei1_boot     = NULL,
    bridge_ei2_boot     = NULL,
    bridge_strength_excl_boot = NULL,
    bridge_betweenness_excl_boot = NULL,
    bridge_closeness_excl_boot   = NULL,
    bridge_ei1_excl_boot = NULL,
    bridge_ei2_excl_boot = NULL,

    edges_true          = edges_true_df,
    edge_boot_mat       = NULL,

    keep_nodes_graph    = keep_nodes_graph,
    keep_nodes_cluster  = keep_nodes_cluster
  )

  class(out) <- c("mixMN_fit")
  return(out)
}
