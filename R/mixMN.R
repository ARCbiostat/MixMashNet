#' Estimate MGM network with bootstrap centrality, bridge metrics, clustering,
#' and (optionally) community scores with CIs.
#'
#' @description
#' Estimates a Mixed Graphical Model (MGM) network on the original data and
#' Estimate MGM network with bootstrap centrality, bridge metrics, clustering,
#' and (optionally) community scores with CIs.
#'
#' @description
#' Estimates a Mixed Graphical Model (MGM) network on the original data and
#' performs non-parametric bootstrap (row resampling) to compute centrality indices,
#' bridge metrics, clustering stability, confidence intervals for node metrics and edge weights.
#' Optionally computes community network scores with CIs and bootstrap arrays.
#'
#' @param data Matrix or data.frame (n x p) with variables in columns
#' @param type,level Vectors as required by \code{mgm::mgm}
#' @param reps Integer (>= 0), number of bootstrap replications
#' @param scale Logical; if TRUE (default) Gaussian variables (type == "g") are z-standardized internally by \code{mgm()}. Use \code{scale = FALSE} if your data are already standardized
#' @param lambdaSel Method for lambda selection: \code{"CV"} or \code{"EBIC"}
#' @param lambdaFolds Number of folds for CV (if \code{lambdaSel="CV"})
#' @param lambdaGam EBIC gamma parameter (if \code{lambdaSel="EBIC"})
#' @param alphaSeq Alpha parameters of the elastic net penalty (values from 0 to 1)
#' @param alphaSel Method for selecting alpha parameter: \code{"CV"} or \code{"EBIC"}
#' @param alphaFolds Number of folds for CV (if \code{alphaSel="CV"})
#' @param alphaGam EBIC gamma parameter (if \code{alphaSel="EBIC"})
#' @param k Integer (>= 0), order of modeled interactions
#' @param ruleReg  rule to combine neighborhood estimates: \code{"AND"} or \code{"OR"}
#' @param threshold threshold below which edge-weights are put to zero: \code{"LW"}, \code{"HW"} or \code{"none"}
#' @param overparameterize Logical; if TRUE uses the over-parameterized version
#' @param thresholdCat Logical; if FALSE thresholds of categorical variables are set to zero
#' @param exclude_from_graph Character vector: nodes excluded from graph/centrality
#' @param exclude_from_cluster Character vector: nodes excluded from clustering (in addition to \code{exclude_from_graph})
#' @param seed_model,seed_boot Seeds for reproducibility (model and bootstrap)
#' @param treat_singletons_as_excluded Logical; if TRUE, singleton communities are treated as excluded
#' @param cluster_method Community detection algorithm: \code{"louvain"}, \code{"fast_greedy"}, \code{"infomap"},
#' \code{"walktrap"}, \code{"edge_betweenness"}
#' @param compute_community_scores Logical; if TRUE, compute community network scores (EGAnet std.scores) with 95% CIs and bootstrap array.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{g_igraph}
#'   \item \code{original_membership}, \code{groups}, \code{community_palette}
#'   \item \code{boot_memberships}
#'   \item \code{reps}, \code{cluster_method}, \code{mgm_model}
#'   \item \code{ci_results} (95\% CIs for node metrics and edge weights)
#'   \item \code{centrality_true} (centrality & bridge metrics on the original network)
#'   \item bootstrap matrices for each node metric (\code{*_boot})
#'   \item \code{edges_true}, \code{edge_boot_mat}
#'   \item \code{keep_nodes_graph}, \code{keep_nodes_cluster}
#'   \item (if \code{compute_community_scores}) \code{community_scores_obj},
#'         \code{community_scores_df}, \code{community_scores_ci},
#'         \code{community_scores_boot} (list of data frames; one per bootstrap replication)
#' }
#'
#' @details
#' - This function does **not** call \code{future::plan()}. Set it beforehand to enable parallel bootstrap.
#' - It calls \code{bridge_metrics()} and \code{bridge_metrics_excluded()} internally (provide them in your package).
#' - Community scores use EGAnet **std.scores** (z-standardized per community).
#'
#' @importFrom mgm mgm
#' @importFrom EGAnet net.scores
#' @importFrom igraph graph_from_adjacency_matrix simplify E ecount V distances betweenness vcount
#' @importFrom igraph cluster_louvain cluster_fast_greedy cluster_infomap cluster_walktrap cluster_edge_betweenness
#' @importFrom qgraph centrality
#' @importFrom colorspace qualitative_hcl
#' @importFrom future.apply future_lapply
#' @importFrom stats setNames quantile
#' @importFrom utils combn capture.output
#' @export
mixMN <- function(
    data, type, level,
    scale = TRUE,
    reps = 100,
    lambdaSel = c("CV", "EBIC"),
    lambdaFolds = 5, lambdaGam = 0.25,
    alphaSeq = 1,
    alphaSel = "CV",
    alphaFolds = 5, alphaGam = 0.25,
    k = 2, ruleReg = "AND", threshold = "LW",
    overparameterize = FALSE, thresholdCat = TRUE,
    exclude_from_graph = NULL,
    exclude_from_cluster = NULL,
    seed_model = NULL, seed_boot = NULL,
    treat_singletons_as_excluded = FALSE,
    cluster_method = c("louvain", "fast_greedy", "infomap", "walktrap", "edge_betweenness"),
    compute_community_scores = FALSE
) {
  lambdaSel <- match.arg(lambdaSel)
  cluster_method <- match.arg(cluster_method)
  if (!is.null(seed_model)) set.seed(seed_model)

  all_nodes <- colnames(data)

  # assign ID
  if (is.null(rownames(data))) {
    rownames(data) <- sprintf("id_%d", seq_len(nrow(data)))
  }
  subject_ids <- rownames(data)

  # ---- helpers ----
  tiny <- 1e-10

  # --- silence EGAnet message ---
  .quiet_net_scores <- function(...) {
    withCallingHandlers(
      {
        out <- NULL
        invisible(capture.output(
          out <- EGAnet::net.scores(...)
        ))
        out
      },
      message = function(m) {
        txt <- conditionMessage(m)
        if (grepl("default 'loading.method'.*changed to \"revised\".*EGAnet", txt, ignore.case = TRUE)) {
          invokeRestart("muffleMessage")
        }
      },
      warning = function(w) {
        txt <- conditionMessage(w)
        if (grepl("default 'loading.method'.*changed to \"revised\".*EGAnet", txt, ignore.case = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  # Build distance graph from signed adjacency: abs weights, distance = 1/weight.
  .make_distance_graph <- function(W_signed) {
    W <- abs(W_signed)
    diag(W) <- 0
    W[is.na(W) | W < tiny] <- 0
    g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)
    if (igraph::ecount(g) > 0) igraph::E(g)$dist <- 1 / igraph::E(g)$weight
    g
  }

  # Harmonic closeness (robust with disconnected graphs)
  .harmonic_closeness <- function(g) {
    if (igraph::ecount(g) == 0) {
      return(stats::setNames(rep(0, igraph::vcount(g)), igraph::V(g)$name))
    }
    D <- igraph::distances(g, weights = igraph::E(g)$dist)
    diag(D) <- NA
    cl <- rowSums(1 / D, na.rm = TRUE) / (nrow(D) - 1)
    stats::setNames(cl, igraph::V(g)$name)
  }

  # ---- Fit MGM on original data ----
  mgm_args <- list(
    data = as.matrix(data),
    type = type,
    level = level,
    lambdaSel = lambdaSel,
    alphaSeq = alphaSeq,
    alphaSel = alphaSel,
    alphaFolds = alphaFolds,
    alphaGam = alphaGam,
    k = k,
    ruleReg = ruleReg,
    threshold = threshold,
    overparameterize = overparameterize,
    thresholdCat = thresholdCat,
    binarySign = TRUE,
    scale = scale,
    signInfo = FALSE
  )
  if (lambdaSel == "CV") {
    mgm_args$lambdaFolds <- lambdaFolds
  } else {
    mgm_args$lambdaGam <- lambdaGam
  }
  mgm_model <- do.call(mgm::mgm, mgm_args)

  wadj  <- mgm_model$pairwise$wadj
  signs <- mgm_model$pairwise$signs
  colnames(wadj)  <- rownames(wadj)  <- all_nodes
  colnames(signs) <- rownames(signs) <- all_nodes

  wadj_signed <- wadj
  if (!is.null(signs)) {
    idx_sign <- !is.na(signs) & (abs(signs) == 1)
    wadj_signed[idx_sign] <- wadj[idx_sign] * signs[idx_sign]
  }
  wadj_signed[is.na(wadj_signed)] <- 0

  # ---- Keep nodes for graph and clustering ----
  keep_nodes_graph   <- setdiff(all_nodes, exclude_from_graph)
  keep_nodes_cluster <- setdiff(all_nodes, unique(c(exclude_from_graph, exclude_from_cluster)))

  wadj_signed_graph   <- wadj_signed[keep_nodes_graph,   keep_nodes_graph]
  wadj_signed_cluster <- wadj_signed[keep_nodes_cluster, keep_nodes_cluster]

  # ---- Build graphs ----
  g_graph   <- igraph::graph_from_adjacency_matrix(abs(wadj_signed_graph),   mode = "undirected", weighted = TRUE, diag = FALSE)
  g_cluster <- igraph::graph_from_adjacency_matrix(abs(wadj_signed_cluster), mode = "undirected", weighted = TRUE, diag = FALSE)
  if (cluster_method %in% c("infomap", "edge_betweenness", "walktrap")) {
    g_cluster <- igraph::simplify(g_cluster, remove.multiple = TRUE, remove.loops = TRUE)
  }

  # ---- two bridge graphs ----
  g_bridge_abs    <- igraph::graph_from_adjacency_matrix(abs(wadj_signed_graph), mode = "undirected", weighted = TRUE, diag = FALSE)
  g_bridge_signed <- igraph::graph_from_adjacency_matrix(wadj_signed_graph,      mode = "undirected", weighted = TRUE, diag = FALSE)

  cluster_fun <- function(graph) {
    switch(
      cluster_method,
      louvain          = igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight),
      fast_greedy      = igraph::cluster_fast_greedy(graph, weights = igraph::E(graph)$weight),
      infomap          = igraph::cluster_infomap(graph),
      walktrap         = igraph::cluster_walktrap(graph, weights = igraph::E(graph)$weight),
      edge_betweenness = igraph::cluster_edge_betweenness(graph, weights = igraph::E(graph)$weight)
    )
  }

  if (cluster_method %in% c("louvain", "infomap") && !is.null(seed_model)) set.seed(seed_model)

  original_membership <- tryCatch({
    m <- cluster_fun(g_cluster)$membership
    stats::setNames(m, keep_nodes_cluster)
  }, error = function(e) rep(NA_integer_, length(keep_nodes_cluster)))

  sizes <- table(original_membership)

  if (treat_singletons_as_excluded) {
    valid_labels <- as.integer(names(sizes[sizes >= 2L]))
    in_valid <- original_membership %in% valid_labels

    unique_clusters <- sort(valid_labels)
    n_clusters <- length(unique_clusters)
    palette_clusters <- stats::setNames(
      if (n_clusters > 0) colorspace::qualitative_hcl(n_clusters, palette = "Dark 3") else character(0),
      unique_clusters
    )

    groups <- as.factor(original_membership[in_valid])
    names(groups) <- keep_nodes_cluster[in_valid]
  } else {
    unique_clusters <- sort(unique(stats::na.omit(original_membership)))
    n_clusters <- length(unique_clusters)
    palette_clusters <- stats::setNames(
      if (n_clusters > 0) colorspace::qualitative_hcl(n_clusters, palette = "Dark 3") else character(0),
      unique_clusters
    )

    groups <- as.factor(original_membership)
    names(groups) <- keep_nodes_cluster
  }

  # --- placeholders per community score ---
  community_scores_obj  <- NULL
  community_scores_df   <- NULL
  community_scores_true <- NULL
  community_scores_boot_arr <- NULL
  community_scores_ci   <- NULL

  # ========== Community network scores (original, conditional) ==========
  nodes_comm <- integer(0)
  wc_comm_int <- integer(0)
  if (isTRUE(compute_community_scores)) {
    nodes_comm <- intersect(names(groups)[!is.na(groups)], keep_nodes_graph)

    if (length(nodes_comm) > 0) {
      A_comm <- wadj_signed_graph[nodes_comm, nodes_comm, drop = FALSE]

      wc_comm <- groups[nodes_comm]
      wc_levels <- sort(unique(as.integer(wc_comm)))
      wc_map <- stats::setNames(seq_along(wc_levels), wc_levels)
      wc_comm_int <- unname(wc_map[as.integer(wc_comm)]) # dense 1..K

      dat_comm <- as.matrix(data[, nodes_comm, drop = FALSE])
      rownames(dat_comm) <- subject_ids

      community_scores_obj <- tryCatch(
        {
          .quiet_net_scores(
            data = dat_comm,
            A = A_comm,
            wc = wc_comm_int,
            loading.method = "revised",
            structure = "simple",
            rotation = NULL,
            scores = "components"
          )
        },
        error = function(e) NULL
      )

      if (!is.null(community_scores_obj)) {
        community_scores_true <- community_scores_obj$scores$std.scores  # std.scores
      } else {
        community_scores_true <- matrix(NA_real_, nrow = nrow(dat_comm), ncol = length(unique(wc_comm_int)))
      }

      K <- ncol(community_scores_true)
      colnames(community_scores_true) <- paste0("NS_", sort(unique(as.integer(wc_comm))))
      rownames(community_scores_true) <- subject_ids

      community_scores_df <- data.frame(
        id = subject_ids,
        as.data.frame(community_scores_true[subject_ids, , drop = FALSE]),
        row.names = NULL,
        check.names = FALSE
      )
    }
  }

  # ---- Centrality & edges on original graph ----
  cent_true <- qgraph::centrality(wadj_signed_graph)  # strength & EI1 (signed)
  g_dist <- .make_distance_graph(wadj_signed_graph)

  closeness_vec <- .harmonic_closeness(g_dist)
  betweenness_vec <- if (igraph::ecount(g_dist) > 0) {
    b <- igraph::betweenness(g_dist, weights = igraph::E(g_dist)$dist, directed = FALSE, normalized = FALSE)
    stats::setNames(b, igraph::V(g_dist)$name)[keep_nodes_graph]
  } else stats::setNames(rep(0, length(keep_nodes_graph)), keep_nodes_graph)

  edge_names <- combn(keep_nodes_graph, 2, FUN = function(x) paste(x[1], x[2], sep = "--"))
  edge_mat_true <- wadj_signed_graph[lower.tri(wadj_signed_graph)]
  names(edge_mat_true) <- edge_names

  # ---- Bridge metrics on original graph ----
  bridge_true_abs <- tryCatch({
    b <- bridge_metrics(g_bridge_abs, membership = groups)
    b[match(keep_nodes_graph, b$node), ]
  }, error = function(e) {
    data.frame(
      node = keep_nodes_graph,
      bridge_strength    = NA_real_,
      bridge_betweenness = NA_real_,
      bridge_closeness   = NA_real_,
      bridge_ei1         = NA_real_,
      bridge_ei2         = NA_real_
    )
  })

  bridge_true_signed <- tryCatch({
    b <- bridge_metrics(g_bridge_signed, membership = groups)
    b[match(keep_nodes_graph, b$node), ]
  }, error = function(e) {
    data.frame(
      node = keep_nodes_graph,
      bridge_strength    = NA_real_,
      bridge_betweenness = NA_real_,
      bridge_closeness   = NA_real_,
      bridge_ei1         = NA_real_,
      bridge_ei2         = NA_real_
    )
  })

  # Combine: take strength/betw/clo from ABS, EI1/EI2 from SIGNED
  bridge_true <- data.frame(
    node               = keep_nodes_graph,
    bridge_strength    = bridge_true_abs$bridge_strength,
    bridge_betweenness = bridge_true_abs$bridge_betweenness,
    bridge_closeness   = bridge_true_abs$bridge_closeness,
    bridge_ei1         = bridge_true_signed$bridge_ei1,
    bridge_ei2         = bridge_true_signed$bridge_ei2
  )

  # ---- "Excluded" metrics: same split ----
  bridge_outside_abs <- tryCatch({
    bridge_metrics_excluded(g_bridge_abs, membership = groups)
  }, error = function(e) {
    data.frame(
      node = keep_nodes_graph,
      bridge_strength            = NA_real_,
      bridge_betweenness         = NA_real_,
      bridge_closeness           = NA_real_,
      bridge_expected_influence1 = NA_real_,
      bridge_expected_influence2 = NA_real_
    )
  })

  bridge_outside_signed <- tryCatch({
    bridge_metrics_excluded(g_bridge_signed, membership = groups)
  }, error = function(e) {
    data.frame(
      node = keep_nodes_graph,
      bridge_strength            = NA_real_,
      bridge_betweenness         = NA_real_,
      bridge_closeness           = NA_real_,
      bridge_expected_influence1 = NA_real_,
      bridge_expected_influence2 = NA_real_
    )
  })

  bridge_outside_true <- data.frame(
    node                         = bridge_outside_abs$node,
    bridge_strength              = bridge_outside_abs$bridge_strength,
    bridge_betweenness           = bridge_outside_abs$bridge_betweenness,
    bridge_closeness             = bridge_outside_abs$bridge_closeness,
    bridge_expected_influence1   = bridge_outside_signed$bridge_expected_influence1, # SIGNED
    bridge_expected_influence2   = bridge_outside_signed$bridge_expected_influence2  # SIGNED
  )

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

  n_nodes_graph <- length(keep_nodes_graph)

  # ---- Bootstrap containers ----
  boot_memberships <- list()
  ci_results <- NULL
  strength_boot <- ei1_boot <- closeness_boot <- betweenness_boot <- NULL
  bridge_strength_boot <- bridge_ei1_boot <- bridge_ei2_boot <- NULL
  bridge_betweenness_boot <- bridge_closeness_boot <- NULL
  bridge_strength_excl_boot <- bridge_betweenness_excl_boot <- NULL
  bridge_closeness_excl_boot <- bridge_ei1_excl_boot <- bridge_ei2_excl_boot <- NULL
  edge_boot_mat <- NULL

  if (reps > 0) {
    boot_output <- future.apply::future_lapply(
      X = seq_len(reps),
      FUN = function(i) {
        index <- sample(1:nrow(data), replace = TRUE)
        boot_data <- as.matrix(data[index, , drop = FALSE])

        boot_args <- list(
          data = boot_data, type = type, level = level,
          lambdaSel = lambdaSel, alphaSeq = alphaSeq,
          alphaSel = alphaSel, alphaFolds = alphaFolds, alphaGam = alphaGam,
          ruleReg = ruleReg, threshold = threshold,
          overparameterize = overparameterize, thresholdCat = thresholdCat,
          k = k, binarySign = TRUE, scale = scale, signInfo = FALSE
        )
        if (lambdaSel == "CV") boot_args$lambdaFolds <- lambdaFolds else boot_args$lambdaGam <- lambdaGam

        boot_model <- tryCatch(do.call(mgm::mgm, boot_args), error = function(e) NULL)
        if (is.null(boot_model)) {
          return(list(
            membership   = rep(NA_integer_, length(keep_nodes_cluster)),
            centralities = rep(NA_real_, 14 * n_nodes_graph),
            edges_boot   = NULL,
            nscores_boot = NULL
          ))
        }

        boot_wadj  <- boot_model$pairwise$wadj
        boot_signs <- boot_model$pairwise$signs
        colnames(boot_wadj)  <- rownames(boot_wadj)  <- all_nodes
        colnames(boot_signs) <- rownames(boot_signs) <- all_nodes

        boot_wadj_signed <- boot_wadj
        if (!is.null(boot_signs)) {
          idxb_sign <- !is.na(boot_signs) & (abs(boot_signs) == 1)
          boot_wadj_signed[idxb_sign] <- boot_wadj[idxb_sign] * boot_signs[idxb_sign]
        }
        boot_wadj_signed[is.na(boot_wadj_signed)] <- 0

        boot_wadj_signed_graph   <- boot_wadj_signed[keep_nodes_graph,   keep_nodes_graph]
        boot_wadj_signed_cluster <- boot_wadj_signed[keep_nodes_cluster, keep_nodes_cluster]

        g_graph_boot   <- igraph::graph_from_adjacency_matrix(abs(boot_wadj_signed_graph),   mode = "undirected", weighted = TRUE, diag = FALSE)
        g_cluster_boot <- igraph::graph_from_adjacency_matrix(abs(boot_wadj_signed_cluster), mode = "undirected", weighted = TRUE, diag = FALSE)
        g_bridge_abs_boot    <- igraph::graph_from_adjacency_matrix(abs(boot_wadj_signed_graph), mode = "undirected", weighted = TRUE, diag = FALSE)
        g_bridge_signed_boot <- igraph::graph_from_adjacency_matrix(boot_wadj_signed_graph,      mode = "undirected", weighted = TRUE, diag = FALSE)
        if (cluster_method %in% c("infomap", "edge_betweenness", "walktrap")) {
          g_cluster_boot <- igraph::simplify(g_cluster_boot, remove.multiple = TRUE, remove.loops = TRUE)
        }

        boot_membership <- tryCatch({
          m <- cluster_fun(g_cluster_boot)$membership
          stats::setNames(m, keep_nodes_cluster)
        }, error = function(e) rep(NA_integer_, length(keep_nodes_cluster)))

        # centrality on bootstrap
        cent_vals <- tryCatch(
          qgraph::centrality(boot_wadj_signed_graph),
          error = function(e) list(
            OutDegree = rep(NA_real_, n_nodes_graph),
            OutExpectedInfluence = rep(NA_real_, n_nodes_graph)
          )
        )

        g_dist_boot <- .make_distance_graph(boot_wadj_signed_graph)

        igraph_closeness_boot <- tryCatch(
          .harmonic_closeness(g_dist_boot)[keep_nodes_graph],
          error = function(e) rep(NA_real_, n_nodes_graph)
        )

        igraph_betweenness_boot <- tryCatch(
          if (igraph::ecount(g_dist_boot) > 0)
            igraph::betweenness(g_dist_boot, weights = igraph::E(g_dist_boot)$dist, directed = FALSE, normalized = FALSE)[keep_nodes_graph]
          else rep(0, n_nodes_graph),
          error = function(e) rep(NA_real_, n_nodes_graph)
        )

        bridge_vals_abs <- tryCatch({
          b <- bridge_metrics(g_bridge_abs_boot, membership = groups)
          b[match(keep_nodes_graph, b$node), ]
        }, error = function(e) {
          data.frame(
            node = keep_nodes_graph,
            bridge_strength    = NA_real_,
            bridge_betweenness = NA_real_,
            bridge_closeness   = NA_real_,
            bridge_ei1         = NA_real_,
            bridge_ei2         = NA_real_
          )
        })

        bridge_vals_signed <- tryCatch({
          b <- bridge_metrics(g_bridge_signed_boot, membership = groups)
          b[match(keep_nodes_graph, b$node), ]
        }, error = function(e) {
          data.frame(
            node = keep_nodes_graph,
            bridge_strength    = NA_real_,
            bridge_betweenness = NA_real_,
            bridge_closeness   = NA_real_,
            bridge_ei1         = NA_real_,
            bridge_ei2         = NA_real_
          )
        })

        bridge_vals <- data.frame(
          node               = keep_nodes_graph,
          bridge_strength    = bridge_vals_abs$bridge_strength,
          bridge_betweenness = bridge_vals_abs$bridge_betweenness,
          bridge_closeness   = bridge_vals_abs$bridge_closeness,
          bridge_ei1         = bridge_vals_signed$bridge_ei1,
          bridge_ei2         = bridge_vals_signed$bridge_ei2
        )

        bridge_outside_boot <- tryCatch({
          abs_part <- bridge_metrics_excluded(g_bridge_abs_boot, membership = groups)
          sgn_part <- bridge_metrics_excluded(g_bridge_signed_boot, membership = groups)
          cbind(
            abs_part[match(keep_nodes_graph, abs_part$node), c("bridge_strength","bridge_betweenness","bridge_closeness")],
            sgn_part[match(keep_nodes_graph, sgn_part$node), c("bridge_expected_influence1","bridge_expected_influence2")]
          )
        }, error = function(e) {
          matrix(NA_real_, nrow = n_nodes_graph, ncol = 5)
        })

        boot_edges <- boot_wadj_signed[keep_nodes_graph, keep_nodes_graph]
        edge_values <- boot_edges[lower.tri(boot_edges)]
        names(edge_values) <- combn(keep_nodes_graph, 2, FUN = function(x) paste(x[1], x[2], sep = "--"))

        # --- community scores on bootstrap
        nscores_boot <- NULL
        if (isTRUE(compute_community_scores) && length(nodes_comm) > 0) {
          A_comm_boot <- boot_wadj_signed_graph[nodes_comm, nodes_comm, drop = FALSE]

          dat_comm_full <- as.matrix(data[, nodes_comm, drop = FALSE])

          ns_obj <- tryCatch(
            {
              .quiet_net_scores(
                data = dat_comm_full,
                A = A_comm_boot,
                wc = wc_comm_int,
                loading.method = "revised",
                structure = "simple",
                rotation = NULL,
                scores = "components"
              )
            },
            error = function(e) NULL
          )

          if (!is.null(ns_obj)) {
            nscores_boot <- ns_obj$scores$std.scores
          } else {
            nscores_boot <- matrix(NA_real_, nrow = length(subject_ids), ncol = length(unique(wc_comm_int)))
          }

          if (!is.null(community_scores_true)) {
            colnames(nscores_boot) <- colnames(community_scores_true)
          } else {
            colnames(nscores_boot) <- paste0("NS_", sort(unique(as.integer(wc_comm_int))))
          }
          rownames(nscores_boot) <- subject_ids
        }

        list(
          membership = boot_membership,
          centralities = c(
            cent_vals$OutDegree,
            cent_vals$OutExpectedInfluence,
            igraph_closeness_boot,
            igraph_betweenness_boot,
            bridge_vals$bridge_strength,
            bridge_vals$bridge_ei1,
            bridge_vals$bridge_ei2,
            bridge_vals$bridge_betweenness,
            bridge_vals$bridge_closeness,
            bridge_outside_boot[, 1],
            bridge_outside_boot[, 2],
            bridge_outside_boot[, 3],
            bridge_outside_boot[, 4],
            bridge_outside_boot[, 5]
          ),
          edges_boot = edge_values,
          nscores_boot = nscores_boot
        )
      },
      future.seed = seed_boot
    )

    # Assemble bootstrap matrices
    boot_mat <- do.call(rbind, lapply(boot_output, function(x) x$centralities))
    boot_memberships <- lapply(boot_output, function(x) x$membership)

    edge_names <- combn(keep_nodes_graph, 2, FUN = function(x) paste(x[1], x[2], sep = "--"))
    edge_boot_mat <- matrix(NA_real_, nrow = length(edge_names), ncol = reps)
    rownames(edge_boot_mat) <- edge_names
    for (i in seq_len(reps)) {
      edges_i <- boot_output[[i]]$edges_boot
      if (!is.null(edges_i)) edge_boot_mat[names(edges_i), i] <- edges_i
    }

    n <- n_nodes_graph
    strength_boot                 <- boot_mat[,  (1):(n),             drop = FALSE]
    ei1_boot                      <- boot_mat[,  (n + 1):(2 * n),     drop = FALSE]
    closeness_boot                <- boot_mat[,  (2 * n + 1):(3 * n), drop = FALSE]
    betweenness_boot              <- boot_mat[,  (3 * n + 1):(4 * n), drop = FALSE]
    bridge_strength_boot          <- boot_mat[,  (4 * n + 1):(5 * n), drop = FALSE]
    bridge_ei1_boot               <- boot_mat[,  (5 * n + 1):(6 * n), drop = FALSE]
    bridge_ei2_boot               <- boot_mat[,  (6 * n + 1):(7 * n), drop = FALSE]
    bridge_betweenness_boot       <- boot_mat[,  (7 * n + 1):(8 * n), drop = FALSE]
    bridge_closeness_boot         <- boot_mat[,  (8 * n + 1):(9 * n), drop = FALSE]
    bridge_strength_excl_boot     <- boot_mat[,  (9 * n + 1):(10 * n), drop = FALSE]
    bridge_betweenness_excl_boot  <- boot_mat[, (10 * n + 1):(11 * n), drop = FALSE]
    bridge_closeness_excl_boot    <- boot_mat[, (11 * n + 1):(12 * n), drop = FALSE]
    bridge_ei1_excl_boot          <- boot_mat[, (12 * n + 1):(13 * n), drop = FALSE]
    bridge_ei2_excl_boot          <- boot_mat[, (13 * n + 1):(14 * n), drop = FALSE]

    colnames(strength_boot)                 <- keep_nodes_graph
    colnames(ei1_boot)                      <- keep_nodes_graph
    colnames(closeness_boot)                <- keep_nodes_graph
    colnames(betweenness_boot)              <- keep_nodes_graph
    colnames(bridge_strength_boot)          <- keep_nodes_graph
    colnames(bridge_ei1_boot)               <- keep_nodes_graph
    colnames(bridge_ei2_boot)               <- keep_nodes_graph
    colnames(bridge_betweenness_boot)       <- keep_nodes_graph
    colnames(bridge_closeness_boot)         <- keep_nodes_graph
    colnames(bridge_strength_excl_boot)     <- keep_nodes_graph
    colnames(bridge_betweenness_excl_boot)  <- keep_nodes_graph
    colnames(bridge_closeness_excl_boot)    <- keep_nodes_graph
    colnames(bridge_ei1_excl_boot)          <- keep_nodes_graph
    colnames(bridge_ei2_excl_boot)          <- keep_nodes_graph

    # CIs for node metrics & edges
    calc_ci <- function(mat) {
      ci <- apply(mat, 2, function(x) {
        if (all(is.na(x))) c(`2.5%` = NA_real_, `97.5%` = NA_real_)
        else stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
      })
      t(ci)
    }

    # assemble community score bootstrap array
    if (isTRUE(compute_community_scores) && !is.null(community_scores_true) && ncol(community_scores_true) > 0) {
      n_subj <- length(subject_ids)
      comm_names <- colnames(community_scores_true)
      K <- length(comm_names)

      community_scores_boot_arr <- array(
        NA_real_, dim = c(reps, n_subj, K),
        dimnames = list(
          rep  = paste0("b", seq_len(reps)),
          id   = subject_ids,
          comm = comm_names
        )
      )

      for (i in seq_len(reps)) {
        ns_i <- boot_output[[i]]$nscores_boot
        if (is.null(ns_i)) next

        ns_i <- as.matrix(ns_i)

        miss_rows <- setdiff(subject_ids, rownames(ns_i))
        if (length(miss_rows)) {
          ns_i <- rbind(ns_i, matrix(NA_real_, nrow = length(miss_rows), ncol = ncol(ns_i),
                                     dimnames = list(miss_rows, colnames(ns_i))))
        }
        miss_cols <- setdiff(comm_names, colnames(ns_i))
        if (length(miss_cols)) {
          ns_i <- cbind(ns_i, matrix(NA_real_, nrow = nrow(ns_i), ncol = length(miss_cols),
                                     dimnames = list(rownames(ns_i), miss_cols)))
        }

        ns_i <- ns_i[subject_ids, comm_names, drop = FALSE]
        community_scores_boot_arr[i, , ] <- ns_i
      }
    }

    # ---- Convert array to list of data.frame (1 per rep)
    community_scores_boot_list <- NULL
    if (isTRUE(compute_community_scores) && !is.null(community_scores_boot_arr)) {
      reps_names <- dimnames(community_scores_boot_arr)$rep
      id_names   <- dimnames(community_scores_boot_arr)$id
      comm_names <- dimnames(community_scores_boot_arr)$comm

      community_scores_boot_list <- lapply(seq_along(reps_names), function(i) {
        # estrai la "fetta" i-esima: subjects x K
        mat <- community_scores_boot_arr[i, , , drop = FALSE][1, , ]
        mat <- as.matrix(mat)
        dimnames(mat) <- list(id_names, comm_names)
        as.data.frame(mat, check.names = FALSE)
      })
      names(community_scores_boot_list) <- reps_names
    }

    ci_results <- list(
      strength                    = calc_ci(strength_boot),
      expected_influence          = calc_ci(ei1_boot),
      closeness                   = calc_ci(closeness_boot),
      betweenness                 = calc_ci(betweenness_boot),
      bridge_strength             = calc_ci(bridge_strength_boot),
      bridge_betweenness          = calc_ci(bridge_betweenness_boot),
      bridge_closeness            = calc_ci(bridge_closeness_boot),
      bridge_ei1                  = calc_ci(bridge_ei1_boot),
      bridge_ei2                  = calc_ci(bridge_ei2_boot),
      bridge_strength_excluded    = calc_ci(bridge_strength_excl_boot),
      bridge_betweenness_excluded = calc_ci(bridge_betweenness_excl_boot),
      bridge_closeness_excluded   = calc_ci(bridge_closeness_excl_boot),
      bridge_ei1_excluded         = calc_ci(bridge_ei1_excl_boot),
      bridge_ei2_excluded         = calc_ci(bridge_ei2_excl_boot),
      edge_weights                = calc_ci(t(edge_boot_mat))
    )

    # Community score CIs (std.scores; subject × community), conditional
    if (isTRUE(compute_community_scores) && !is.null(community_scores_true)) {
      qfun <- function(x, p) if (all(is.na(x))) NA_real_ else stats::quantile(x, probs = p, na.rm = TRUE)
      cs_ci_low  <- apply(community_scores_boot_arr, c(2, 3), qfun, p = 0.025)
      cs_ci_high <- apply(community_scores_boot_arr, c(2, 3), qfun, p = 0.975)
      dimnames(cs_ci_low)  <- dimnames(community_scores_true)
      dimnames(cs_ci_high) <- dimnames(community_scores_true)
      community_scores_ci <- list(lower = cs_ci_low, upper = cs_ci_high)
    }
  } # end if reps > 0

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


  # ---- Return ----
  out <- list(
    graph_igraph        = g_igraph,
    original_membership = original_membership,
    groups              = groups,
    community_palette   = palette_clusters,
    boot_memberships    = boot_memberships,
    reps                = reps,
    cluster_method      = cluster_method,
    mgm_model           = mgm_model,
    ci_results          = ci_results,
    centrality_true     = centrality_true_df,

    strength_boot                 = strength_boot,
    ei1_boot                      = ei1_boot,
    closeness_boot                = closeness_boot,
    betweenness_boot              = betweenness_boot,
    bridge_strength_boot          = bridge_strength_boot,
    bridge_betweenness_boot       = bridge_betweenness_boot,
    bridge_closeness_boot         = bridge_closeness_boot,
    bridge_ei1_boot               = bridge_ei1_boot,
    bridge_ei2_boot               = bridge_ei2_boot,
    bridge_strength_excl_boot     = bridge_strength_excl_boot,
    bridge_betweenness_excl_boot  = bridge_betweenness_excl_boot,
    bridge_closeness_excl_boot    = bridge_closeness_excl_boot,
    bridge_ei1_excl_boot          = bridge_ei1_excl_boot,
    bridge_ei2_excl_boot          = bridge_ei2_excl_boot,

    edges_true        = edges_true_df,
    edge_boot_mat     = edge_boot_mat,
    keep_nodes_graph  = keep_nodes_graph,
    keep_nodes_cluster= keep_nodes_cluster,

    # Community scores (present only if compute_community_scores = TRUE)
    community_scores_obj  = if (isTRUE(compute_community_scores)) community_scores_obj else NULL,
    community_scores_df   = if (isTRUE(compute_community_scores)) community_scores_df  else NULL,
    community_scores_ci   = if (isTRUE(compute_community_scores) && !is.null(community_scores_ci)) community_scores_ci else NULL,
    community_scores_boot = if (isTRUE(compute_community_scores) && !is.null(community_scores_boot_list)) community_scores_boot_list else NULL
  )

  class(out) <- c("mixMN_fit")
  return(out)
}
