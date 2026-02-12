#' Build MixMashNet metrics from a signed weighted adjacency matrix
#'
#' @description
#' Low-level constructor that builds a \code{mixMN_fit} object starting from a
#' signed weighted adjacency matrix (\code{wadj_signed}), instead of estimating
#' the network via \code{mgm::mgm}. The function computes community detection,
#' node-level centrality indices, and bridge metrics (including "excluded"
#' versions), and returns an object compatible with MixMashNet downstream
#' utilities.
#'
#' No bootstrap is performed inside this function. The arguments \code{reps},
#' \code{seed_boot}, \code{quantile_level}, and \code{boot_what} are stored in
#' \code{$settings} for interface compatibility with \code{mixMN()} and other
#' bootstrap-based functions, but the corresponding bootstrap slots in the output
#' are set to \code{NULL}.
#'
#' @param wadj_signed Square numeric matrix (\eqn{p \times p}) with signed edge
#'   weights (rows/columns = nodes). Typically \code{wadj * signs} from \pkg{mgm}.
#'   Row and column names must correspond to \code{nodes}.
#' @param nodes Character vector of node names (must match
#'   \code{rownames(wadj_signed)} and \code{colnames(wadj_signed)}).
#' @param quantile_level Quantile level for quantile regions (default 0.95).
#'   Stored in \code{$settings} for consistency with \code{mixMN()}, but no quantile regions
#'   are computed here because no bootstrap is performed.
#' @param covariates Character vector. Nodes excluded entirely from the
#'   graph (removed from adjacency matrix, edges and node-level metrics).
#' @param exclude_from_cluster Character vector. Nodes excluded from community
#'   detection (in addition to \code{covariates}). These nodes remain in
#'   the graph and in node-level metrics, but receive \code{NA} community labels.
#' @param cluster_method Community detection algorithm, one of
#'   \code{c("louvain","fast_greedy","infomap","walktrap","edge_betweenness")}.
#'   Community detection is performed on the absolute-weight graph induced by
#'   \code{keep_nodes_cluster}.
#' @param reps Integer (>= 0). Number of bootstrap replications. Not used inside
#'   this function; stored in \code{$settings$reps} for downstream bootstrap
#'   utilities.
#' @param seed_boot Optional integer. Stored in \code{$settings} for consistency
#'   with bootstrap utilities; not used here.
#' @param treat_singletons_as_excluded Logical; if \code{TRUE}, singleton
#'   communities (size 1) are treated as excluded: their membership is set to
#'   \code{NA} in \code{$communities$groups} (used for bridge metrics).
#' @param boot_what Character vector specifying which quantities are intended to
#'   be bootstrapped in higher-level functions. The value is stored in
#'   \code{$settings$boot_what} but no bootstrap is computed here. Can include
#'   \code{"general_index"}, \code{"interlayer_index"}, \code{"bridge_index"},
#'   \code{"excluded_index"}, \code{"community"}, \code{"loadings"}, or
#'   \code{"none"}.
#'
#' @return
#' An object of class \code{c("mixmashnet","mixMN_fit")} with the following
#' components:
#' \describe{
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{settings}}{
#'     List echoing the main arguments, including \code{reps},
#'     \code{cluster_method}, \code{covariates},
#'     \code{exclude_from_cluster}, \code{treat_singletons_as_excluded},
#'     \code{boot_what}, and \code{quantile_level}.
#'   }
#'   \item{\code{model}}{
#'     List with \code{mgm = NULL} and \code{nodes} (character vector of all node
#'     names). This mirrors \code{mixMN()}, but no MGM is fitted here.
#'   }
#'   \item{\code{graph}}{
#'     List describing the graph:
#'     \code{igraph} (signed \pkg{igraph} object built on \code{keep_nodes_graph},
#'     with edge attributes \code{weight}, \code{abs_weight}, \code{sign} and
#'     vertex attribute \code{membership}),
#'     \code{keep_nodes_graph} (nodes retained in the graph and all node-level
#'     metrics), and \code{keep_nodes_cluster} (nodes used for community
#'     detection).
#'   }
#'   \item{\code{communities}}{
#'     List describing community structure with:
#'     \code{original_membership} (integer vector of community labels on
#'     \code{keep_nodes_cluster}),
#'     \code{groups} (factor membership used for bridge metrics; may contain
#'     \code{NA} for excluded/singleton nodes),
#'     \code{palette} (named vector of colors per community), and
#'     \code{boot_memberships} (empty list).
#'   }
#'   \item{\code{statistics}}{
#'     List with node- and edge-level summaries:
#'     \code{node$true} is a data frame with one row per node in
#'     \code{keep_nodes_graph}, containing \code{strength}, \code{ei1},
#'     \code{closeness}, \code{betweenness}, and bridge metrics
#'     (\code{bridge_strength}, \code{bridge_betweenness}, \code{bridge_closeness},
#'     \code{bridge_ei1}, \code{bridge_ei2}, plus \code{*_excluded} versions).
#'     Bootstrap slots \code{node$boot}, \code{node$quantile_region}, \code{edge$boot},
#'     \code{edge$quantile_region} are set to \code{NULL}.
#'   }
#'   \item{\code{community_loadings}}{
#'     List container for community loadings (aligned to \code{mixMN()}):
#'     \code{nodes} (nodes used for loadings), \code{wc} (integer community labels
#'     aligned with \code{nodes}), \code{true} (loadings matrix), and \code{boot}
#'     (set to \code{NULL} here).
#'   }
#' }
#'
#' @details
#' Centrality indices are computed on the signed matrix \code{wadj_signed} via
#' \code{qgraph::centrality}. Harmonic closeness and betweenness are computed on
#' an absolute-weight distance graph with edge length \eqn{1/|w|}. Bridge metrics
#' are computed separately on the absolute and signed graphs using
#' \code{bridge_metrics()} and \code{bridge_metrics_excluded()}.
#'
#' @importFrom igraph cluster_louvain cluster_fast_greedy cluster_infomap
#' @importFrom igraph cluster_walktrap cluster_edge_betweenness
#' @importFrom qgraph centrality
#' @importFrom grDevices hcl
#' @importFrom stats setNames
#' @importFrom utils combn
#' @keywords internal
#' @noRd
mixMN_from_wadj <- function(
    wadj_signed,
    nodes,
    quantile_level = 0.95,
    covariates = NULL,
    exclude_from_cluster = NULL,
    cluster_method = c("louvain", "fast_greedy", "infomap", "walktrap", "edge_betweenness"),
    reps = 0,
    seed_boot = NULL,
    treat_singletons_as_excluded = FALSE,
    boot_what = c("general_index", "interlayer_index", "bridge_index",
                  "excluded_index", "community", "loadings", "none"),
    palette_layer = NULL
) {
  cluster_method <- match.arg(cluster_method)

  boot_what <- match.arg(
    boot_what,
    choices = c("general_index","interlayer_index","bridge_index",
                "excluded_index","community","loadings","none"),
    several.ok = TRUE
  )

  # quantile region
  if (!is.numeric(quantile_level) || length(quantile_level) != 1L ||
      is.na(quantile_level) || quantile_level <= 0 || quantile_level >= 1) {
    stop("`quantile_level` must be a single number strictly between 0 and 1 (e.g., 0.95).")
  }
  alpha <- 1 - quantile_level
  quantile_region_probs <- c(alpha/2, 1 - alpha/2)

  all_nodes <- nodes

  # --- Keep nodes for graph/cluster scopes ---
  keep_nodes_graph   <- setdiff(all_nodes, covariates)
  keep_nodes_cluster <- setdiff(all_nodes, unique(c(covariates, exclude_from_cluster)))

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

  .make_palette_layer <- function(unique_clusters, palette_layer) {
    k <- length(unique_clusters)
    if (k == 0) return(stats::setNames(character(0), unique_clusters))
    if (missing(palette_layer) || is.null(palette_layer)) {
      stop("mixMN_from_wadj(): 'palette_layer' must be provided (no fallback).", call. = FALSE)
    }
    if (is.character(palette_layer) && length(palette_layer) == 1L) {
      pal <- colorspace::qualitative_hcl(k, palette = palette_layer)
      return(stats::setNames(pal, unique_clusters))
    }
    if (is.character(palette_layer) && length(palette_layer) > 1L) {
      pal <- palette_layer
      if (length(pal) >= 2L) pal <- pal[-1]
      pal <- rep(pal, length.out = k)
      return(stats::setNames(pal, unique_clusters))
    }
    stop("mixMN_from_wadj(): 'palette_layer' must be either a colorspace palette name (length 1) or a vector of color hex codes.", call. = FALSE)
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
    palette_clusters <- .make_palette_layer(unique_clusters, palette_layer = palette_layer)
    groups <- as.factor(original_membership[in_valid])
    names(groups) <- keep_nodes_cluster[in_valid]
  } else {
    unique_clusters <- sort(unique(stats::na.omit(original_membership)))
    palette_clusters <- .make_palette_layer(unique_clusters, palette_layer)
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
    vnames <- igraph::V(g)$name
    if (is.null(vnames)) vnames <- as.character(seq_len(igraph::vcount(g)))
    if (igraph::vcount(g) <= 1) {
      out <- rep(NA_real_, igraph::vcount(g))
      return(stats::setNames(out, vnames))
    }
    if (igraph::ecount(g) == 0) {
      return(stats::setNames(rep(NA_real_, igraph::vcount(g)), vnames))
    }
    D <- igraph::distances(g, weights = igraph::E(g)$dist)
    diag(D) <- NA_real_
    cl <- rowSums(1 / D, na.rm = TRUE) / (nrow(D) - 1)
    no_reach <- apply(D, 1, function(x) !any(is.finite(x), na.rm = TRUE))
    cl[no_reach] <- NA_real_
    stats::setNames(cl, vnames)
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

  # ---- iGraph object  ----
  Wg <- wadj_signed[keep_nodes_graph, keep_nodes_graph, drop = FALSE]

  g_igraph <- igraph::graph_from_adjacency_matrix(
    Wg,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  if (igraph::ecount(g_igraph) > 0) {
    igraph::E(g_igraph)$abs_weight <- abs(igraph::E(g_igraph)$weight)
    igraph::E(g_igraph)$sign <- ifelse(igraph::E(g_igraph)$weight >= 0, 1L, -1L)
  }

  igraph::V(g_igraph)$name <- keep_nodes_graph
  if (length(groups)) {
    memb <- rep(NA_integer_, length(keep_nodes_graph))
    names(memb) <- keep_nodes_graph
    memb[names(groups)] <- as.integer(groups)
    igraph::V(g_igraph)$membership <- memb
  } else {
    igraph::V(g_igraph)$membership <- NA_integer_
  }

  nodes_comm <- character(0)
  wc_comm_int <- integer(0)
  community_loadings_true <- NULL

  # --- Return object skeleton (no bootstrap here) ---
  out <- list(
    call = match.call(),

    settings = list(
      reps                         = reps,
      cluster_method               = cluster_method,
      covariates                   = covariates,
      exclude_from_cluster         = exclude_from_cluster,
      treat_singletons_as_excluded = treat_singletons_as_excluded,
      boot_what                    = boot_what,
      quantile_level = quantile_level
    ),

    model = list(
      mgm   = NULL,
      nodes = all_nodes
    ),

    graph = list(
      igraph             = g_igraph,
      keep_nodes_graph   = keep_nodes_graph,
      keep_nodes_cluster = keep_nodes_cluster
    ),

    communities = list(
      original_membership = original_membership,
      groups              = groups,
      palette             = palette_clusters,
      boot_memberships    = list()
    ),

    statistics = list(
      node = list(
        true = centrality_true_df,
        boot = list(
          strength                    = NULL,
          ei1                         = NULL,
          closeness                   = NULL,
          betweenness                 = NULL,
          bridge_strength             = NULL,
          bridge_betweenness          = NULL,
          bridge_closeness            = NULL,
          bridge_ei1                  = NULL,
          bridge_ei2                  = NULL,
          bridge_strength_excluded    = NULL,
          bridge_betweenness_excluded = NULL,
          bridge_closeness_excluded   = NULL,
          bridge_ei1_excluded         = NULL,
          bridge_ei2_excluded         = NULL
        ),
        quantile_region = NULL
      ),
      edge = list(
        true = edges_true_df,
        boot = NULL,
        quantile_region   = NULL
      )
    ),

    community_loadings = list(
      nodes = character(0),
      wc    = integer(0),
      true  = NULL,
      boot  = NULL,
      available = FALSE,
      reason = "Loadings are computed at higher level (multimixMN) when requested.",
      non_scorable_nodes = character(0)
    )
  )

  class(out) <- c("mixmashnet", "mixMN_fit")
  return(out)
}
