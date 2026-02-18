#' Estimate single layer MGM network with bootstrap centrality, bridge metrics, clustering,
#' and (optionally) community score loadings
#'
#' @description
#' Estimates a single layer Mixed Graphical Model (MGM) network on the original data,
#' using the estimation framework implemented in the \pkg{mgm} package, and performs
#' non-parametric bootstrap (row resampling) to compute centrality indices, bridge
#' metrics, clustering stability, and quantile regions for node metrics
#' and edge weights.
#' Optionally, the function computes community score loadings (for later prediction
#' on new data) and can bootstrap the corresponding loadings.
#'
#' @param data A \code{data.frame} (n x p) with variables in columns.
#'   Variables may be numeric, integer, logical, or factors.
#'   Character and Date/POSIXt variables are not supported and must be converted
#'   prior to model fitting.
#'   Variable types are internally mapped to MGM types as follows:
#'   continuous numeric (double) variables are treated as Gaussian;
#'   integer variables are treated as Poisson unless they take only values
#'   in \{0,1\}, in which case they are treated as binary categorical;
#'   factors and logical variables are treated as categorical.
#'   Binary categorical variables (two-level
#'   factors and logical variables) are internally recoded to \{0,1\} for model
#'   fitting. The original input data are not modified.
#' @param reps Integer (>= 0). Number of bootstrap replications.
#' @param scale Logical; if \code{TRUE} (default) Gaussian variables
#'   (\code{type == "g"}) are z-standardized internally by \code{mgm()}. Use
#'   \code{scale = FALSE} if your data are already standardized.
#' @param lambdaSel Method for lambda selection: \code{"CV"} or \code{"EBIC"}.
#' @param lambdaFolds Number of folds for CV (if \code{lambdaSel = "CV"}).
#' @param lambdaGam EBIC gamma parameter (if \code{lambdaSel = "EBIC"}).
#' @param alphaSeq Alpha parameters of the elastic net penalty (values between 0 and 1).
#' @param alphaSel Method for selecting the alpha parameter: \code{"CV"} or \code{"EBIC"}.
#' @param alphaFolds Number of folds for CV (if \code{alphaSel = "CV"}).
#' @param alphaGam EBIC gamma parameter (if \code{alphaSel = "EBIC"}).
#' @param k Integer (>= 1). Order of modeled interactions.
#' @param ruleReg Rule to combine neighborhood estimates: \code{"AND"} or \code{"OR"}.
#' @param threshold Threshold below which edge-weights are set to zero:
#'   Available options are \code{"LW"}, \code{"HW"}, or \code{"none"}.
#'   \code{"LW"} applies the threshold proposed by Loh & Wainwright;
#'   \code{"HW"} applies the threshold proposed by Haslbeck & Waldorp;
#'   \code{"none"} disables thresholding. Defaults to \code{"LW"}.
#' @param overparameterize Logical; controls how categorical interactions are
#'   parameterized in the neighborhood regressions. If \code{TRUE}, categorical
#'   interactions are represented using a fully over-parameterized design matrix
#'   (i.e., all category combinations are explicitly modeled). If \code{FALSE},
#'   the standard \code{glmnet} parameterization is used, where one category
#'   serves as reference. For continuous variables, both parameterizations are
#'   equivalent. The default is \code{FALSE}. The over-parameterized option may
#'   be advantageous when distinguishing pairwise from higher-order interactions.
#' @param thresholdCat Logical; if \code{FALSE} thresholds of categorical
#'   variables are set to zero.
#' @param quantile_level Level of the central bootstrap quantile region (default \code{0.95}).
#'   Must be a single number between 0 and 1.
#' @param covariates Character vector. Variables used as adjustment covariates
#'   in model estimation.
#' @param exclude_from_cluster Character vector. Nodes excluded from community
#'   detection (in addition to \code{covariates}).
#' @param treat_singletons_as_excluded Logical; if \code{TRUE}, singleton
#'   communities (size 1) are treated as excluded nodes when computing
#'   bridge metrics.
#' @param seed_model Optional integer seed for reproducibility of the initial
#'   MGM fit.
#' @param seed_boot Optional integer seed passed to \code{future.apply} for
#'   reproducibility of bootstrap replications.
#' @param cluster_method Community detection algorithm used on the network:
#'   \code{"louvain"}, \code{"fast_greedy"}, \code{"infomap"},
#'   \code{"walktrap"}, or \code{"edge_betweenness"}.
#' @param compute_loadings Logical; if \code{TRUE} (default), compute community loadings
#'   (\code{EGAnet::net.loads}). Only supported for Gaussian, Poisson, and binary
#'   categorical nodes; otherwise loadings are skipped and the reason is
#'   stored in \code{community_loadings$reason}.
#' @param boot_what Character vector specifying which quantities to bootstrap.
#'   Valid options are:
#'    \code{"general_index"} (centrality indices),
#'    \code{"bridge_index"} (bridge metrics for nodes in communities),
#'    \code{"excluded_index"} (bridge metrics for nodes treated as excluded),
#'    \code{"community"} (community memberships),
#'    \code{"loadings"} (community loadings, only if \code{compute_loadings = TRUE}),
#'    and \code{"none"} (skip all node-level bootstrap: only edge-weight
#'    bootstrap is performed if \code{reps > 0}).
#' @param save_data Logical; if \code{TRUE}, store the original data in the output object.
#' @param progress Logical; if \code{TRUE} (default), show a bootstrap progress bar.
#'
#' @return
#' An object of class \code{c("mixmashnet", "mixMN_fit")}, that is a list with
#' the following top-level components:
#' \describe{
#'   \item{\code{call}}{
#'     The matched function call.
#'   }
#'   \item{\code{settings}}{
#'     List of main settings used in the call, including
#'     \code{reps}, \code{cluster_method},
#'     \code{covariates}, \code{exclude_from_cluster},
#'     \code{treat_singletons_as_excluded}, and \code{boot_what}.
#'   }
#'  \item{\code{data_info}}{
#'     List with information derived from the input data used for model setup:
#'     \code{mgm_type_level} (data frame with one row per variable, reporting
#'     the original R class and the inferred MGM \code{type} and \code{level},
#'     as used in the call to \code{mgm::mgm}),
#'     and \code{binary_recode_map} (named list describing the mapping from
#'     original binary labels to the internal \{0,1\} coding used for model fitting).
#'   }
#'   \item{\code{model}}{
#'     List with:
#'     \code{mgm} (the fitted \code{mgm} object),
#'     \code{nodes} (character vector of all node names),
#'     \code{n} (number of observations),
#'     \code{p} (number of variables), and
#'     \code{data} (if \code{save_data = TRUE}).
#'   }
#'   \item{\code{graph}}{
#'     List describing the graph:
#'     \code{igraph} (an \pkg{igraph} object built on
#'     \code{keep_nodes_graph}, with edge attributes
#'     \code{weight}, \code{abs_weight}, \code{sign} and vertex attribute
#'     \code{membership} for communities),
#'     \code{keep_nodes_graph} (nodes retained in the graph and all node-level
#'     metrics), and \code{keep_nodes_cluster} (nodes used for community
#'     detection).
#'   }
#'   \item{\code{communities}}{
#'     List describing community structure with:
#'     \code{original_membership} (integer vector of community labels on
#'     \code{keep_nodes_cluster}),
#'     \code{groups} (factor of community labels actually used for bridge
#'     metrics, optionally with singletons treated as excluded),
#'     \code{palette} (named vector of colors per community), and
#'     \code{boot_memberships} (list of bootstrap memberships if
#'     \code{"community"} is requested in \code{boot_what}, otherwise an empty
#'     list).
#'   }
#'   \item{\code{statistics}}{
#'     List with node- and edge-level summaries:
#'     \code{node} is a list with:
#'     \code{true} (data frame with one row per node in
#'     \code{keep_nodes_graph}, containing the node name and metrics
#'     \code{strength}, \code{ei1}, \code{closeness}, \code{betweenness},
#'     \code{bridge_strength}, \code{bridge_betweenness}, \code{bridge_closeness},
#'     \code{bridge_ei1}, \code{bridge_ei2}, and for nodes treated as excluded
#'     from communities also
#'     \code{bridge_strength_excluded},
#'     \code{bridge_betweenness_excluded},
#'     \code{bridge_closeness_excluded},
#'     \code{bridge_ei1_excluded}, \code{bridge_ei2_excluded});
#'     \code{boot} (list of bootstrap matrices for each metric, each of
#'     dimension \code{reps x length(keep_nodes_graph)}, possibly \code{NULL}
#'     if the metric was not requested or if \code{reps = 0}); and
#'     \code{quantile_region} (list of quantile regions for each node metric, one
#'     \code{p x 2} matrix per metric, with columns corresponding to the lower and upper
#'     quantile bounds implied by \code{quantile_level}, or \code{NULL} if no bootstrap was performed).
#'
#'     \code{edge} is a list with:
#'     \code{true} (data frame with columns \code{edge} and \code{weight} for
#'     all unique undirected edges among \code{keep_nodes_graph});
#'     \code{boot} (matrix of bootstrap edge weights of dimension
#'     \code{n_edges x reps}); and
#'     \code{quantile_region} (matrix of quantile regions for edge weights,
#'     \code{n_edges x 2}, with columns corresponding to the lower and upper
#'     bootstrap quantile bounds, or \code{NULL} if \code{reps = 0}).
#'   }
#'   \item{\code{community_loadings}}{
#'     List containing community-loading information (based on
#'     \code{EGAnet::net.loads}) for later community-score computation on new
#'     data:
#'     \code{nodes}(nodes used for loadings),
#'     \code{wc} (integer community labels aligned with \code{nodes}),
#'     \code{true} (matrix of standardized loadings, nodes x communities,
#'     or \code{NULL} if loadings were not computed.),
#'     \code{boot} (list of bootstrap loading matrices, one per replication,
#'     or \code{NULL} if not bootstrapped),
#'     \code{available} (logical indicating whether loadings were computed),
#'     \code{reason} (character string explaining why loadings were not computed,
#'         or \code{NULL} if \code{available = TRUE}),
#'     \code{non_scorable_nodes} (character vector of nodes in the community
#'         subgraph that prevented loadings from being computed (e.g., categorical variables
#'         with >2 levels), otherwise empty).
#'     }
#'   }
#'
#' @details
#' This function does \strong{not} call \code{future::plan()}. To enable
#' parallel bootstrap, set a plan (e.g. \code{future::plan(multisession)})
#' before calling \code{mixMN()}. If \code{boot_what} is \code{"none"} and
#' \code{reps > 0}, node-level metrics are not bootstrapped but edge-weight
#' bootstrap and corresponding quantile regions are still computed.
#'
#' @references
#'
#' Haslbeck, J. M. B., & Waldorp, L. J. (2020).
#' mgm: Estimating Time-Varying Mixed Graphical Models in High-Dimensional Data.
#' \emph{Journal of Statistical Software}, 93(8).
#' \doi{10.18637/jss.v093.i08}
#'
#' Loh, P. L., & Wainwright, M. J. (2012).
#' Structure estimation for discrete graphicalmodels:
#' Generalized covariance matrices and their inverses.
#' \emph{NIPS}
#'
#' @importFrom mgm mgm
#' @importFrom EGAnet net.loads
#' @importFrom igraph graph_from_adjacency_matrix simplify E ecount V distances betweenness vcount
#' @importFrom igraph cluster_louvain cluster_fast_greedy cluster_infomap cluster_walktrap cluster_edge_betweenness
#' @importFrom qgraph centrality
#' @importFrom colorspace qualitative_hcl
#' @importFrom future.apply future_lapply
#' @importFrom stats setNames quantile
#' @importFrom utils combn capture.output
#' @importFrom progressr with_progress progressor
#' @export
mixMN <- function(
    data,
    reps = 100,
    scale = TRUE,
    lambdaSel = c("CV", "EBIC"),
    lambdaFolds = 5,
    lambdaGam = 0.25,
    alphaSeq = 1,
    alphaSel = "CV",
    alphaFolds = 5,
    alphaGam = 0.25,
    k = 2,
    ruleReg = "AND",
    threshold = "LW",
    overparameterize = FALSE,
    thresholdCat = TRUE,
    quantile_level = 0.95,
    covariates = NULL,
    exclude_from_cluster = NULL,
    treat_singletons_as_excluded = FALSE,
    seed_model = NULL,
    seed_boot = NULL,
    cluster_method = c("louvain", "fast_greedy", "infomap", "walktrap", "edge_betweenness"),
    compute_loadings = TRUE,
    boot_what = c("general_index", "bridge_index", "excluded_index",
                  "community", "loadings"),
    save_data = FALSE,
    progress = TRUE
) {
  lambdaSel <- match.arg(lambdaSel)
  cluster_method <- match.arg(cluster_method)
  if (!is.null(seed_model)) set.seed(seed_model)

  # ---- Infer mgm spec + build mgm matrix (0/1 for binary factors/logicals) ----
  spec <- infer_mgm_spec(data, recode_binary = TRUE)  # TRUE perché sotto usi binarySign=TRUE
  type  <- spec$type
  level <- spec$level
  data_mgm <- spec$data_mgm

  # for print
  var_interpretation <- spec$data_info
  binary_recode_map  <- spec$binary_recode_map

  all_nodes <- colnames(data)

  # assign ID
  if (is.null(rownames(data))) {
    rownames(data) <- sprintf("id_%d", seq_len(nrow(data)))
  }
  subject_ids <- rownames(data)

  # quantile region
  if (!is.numeric(quantile_level) || length(quantile_level) != 1L ||
      is.na(quantile_level) || quantile_level <= 0 || quantile_level >= 1) {
    stop("`quantile_level` must be a single number strictly between 0 and 1 (e.g., 0.95).")
  }
  alpha <- 1 - quantile_level
  probs <- c(alpha/2, 1 - alpha/2)


  # ---- helpers ----
  tiny <- 1e-10

  # --- silence EGAnet message ---
  .quiet_net_loads <- function(...) {
    withCallingHandlers(
      {
        out <- NULL
        invisible(capture.output(
          out <- EGAnet::net.loads(...)
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

  # ---- Fit MGM on original data ----
  mgm_args <- list(
    data = data_mgm,
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
    signInfo = FALSE,
    pbar = FALSE
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
  keep_nodes_graph   <- setdiff(all_nodes, covariates)
  keep_nodes_cluster <- setdiff(all_nodes, unique(c(covariates, exclude_from_cluster)))

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
  community_loadings_true <- NULL   # matrix p x K (std loadings)
  community_loadings_boot <- NULL   # list length reps (oppure array)
  nodes_comm <- character(0)
  wc_comm_int <- integer(0)

  # ========== Community network scores (original, conditional) ==========
  nodes_comm <- intersect(names(groups)[!is.na(groups)], keep_nodes_graph)
  A_comm <- wadj_signed_graph[nodes_comm, nodes_comm, drop = FALSE]
  wc_comm <- groups[nodes_comm]
  wc_levels <- sort(unique(as.integer(wc_comm)))
  wc_map <- stats::setNames(seq_along(wc_levels), wc_levels)
  wc_comm_int <- unname(wc_map[as.integer(wc_comm)])

  # ---- Loadings eligibility (based on inferred MGM type/level) ----
  type_vec  <- stats::setNames(type,  all_nodes)
  level_vec <- stats::setNames(level, all_nodes)

  is_scorable_node <- function(nm) {
    tt <- type_vec[[nm]]
    ll <- level_vec[[nm]]
    (tt %in% c("g","p")) || (tt == "c" && !is.na(ll) && ll == 2L)
  }

  loadings_available <- FALSE
  loadings_reason <- NULL
  non_scorable_nodes <- character(0)

  if (isTRUE(compute_loadings) && length(nodes_comm) > 1) {
    ok_nodes <- vapply(nodes_comm, is_scorable_node, logical(1))

    if (!all(ok_nodes)) {
      loadings_available <- FALSE
      non_scorable_nodes <- nodes_comm[!ok_nodes]
      loadings_reason <- paste0(
        "Loadings not computed: categorical variables with more than two levels were found. ",
        "Community scores are only supported for Gaussian, Poisson, ",
        "or binary categorical variables (i.e., categorical variables with exactly two levels)."
      )

      community_loadings_true <- NULL
      compute_loadings <- FALSE
    } else {
      loadings_available <- TRUE
    }
  }

  if (isTRUE(loadings_available) && length(nodes_comm) > 1) {
    loads_obj <- tryCatch(
      .quiet_net_loads(
        A = A_comm,
        wc = wc_comm_int,
        loading.method = "revised",
        rotation = NULL
      ),
      error = function(e) NULL
    )

    if (!is.null(loads_obj) && !is.null(loads_obj$std)) {
      community_loadings_true <- loads_obj$std
      community_loadings_true <- community_loadings_true[nodes_comm, , drop = FALSE]
    } else {
      community_loadings_true <- matrix(NA_real_, nrow = length(nodes_comm),
                                        ncol = length(unique(wc_comm_int)),
                                        dimnames = list(nodes_comm, paste0("C", seq_len(length(unique(wc_comm_int))))))
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

  pairs <- combn(keep_nodes_graph, 2)
  edge_names <- apply(pairs, 2, paste, collapse = "--")
  edge_mat_true <- apply(pairs, 2, function(x) wadj_signed_graph[x[1], x[2]])
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
      bridge_strength    = NA_real_,
      bridge_betweenness = NA_real_,
      bridge_closeness   = NA_real_,
      bridge_ei1         = NA_real_,
      bridge_ei2         = NA_real_
    )
  })

  bridge_outside_signed <- tryCatch({
    bridge_metrics_excluded(g_bridge_signed, membership = groups)
  }, error = function(e) {
    data.frame(
      node = keep_nodes_graph,
      bridge_strength       = NA_real_,
      bridge_betweenness    = NA_real_,
      bridge_closeness      = NA_real_,
      bridge_ei1            = NA_real_,
      bridge_ei2            = NA_real_
    )
  })

  bridge_outside_true <- data.frame(
    node                  = bridge_outside_abs$node,
    bridge_strength       = bridge_outside_abs$bridge_strength,
    bridge_betweenness    = bridge_outside_abs$bridge_betweenness,
    bridge_closeness      = bridge_outside_abs$bridge_closeness,
    bridge_ei1            = bridge_outside_signed$bridge_ei1, # SIGNED
    bridge_ei2            = bridge_outside_signed$bridge_ei2  # SIGNED
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
    "bridge_ei1", "bridge_ei2"
  )]
  centrality_true_df <- cbind(centrality_true_df, bridge_excluded_df[, -1])

  n_nodes_graph <- length(keep_nodes_graph)

  # --- parse 'boot_what' argument ---
  boot_what <- match.arg(
    boot_what,
    choices = c("general_index", "bridge_index", "excluded_index",
                "community", "loadings", "none"),
    several.ok = TRUE
  )

  if (length(boot_what) > 1L && "none" %in% boot_what) {
    boot_what <- setdiff(boot_what, "none")
  }

  if ("none" %in% boot_what && length(boot_what) == 1L) {
    do_general_boot     <- FALSE
    do_bridge_boot      <- FALSE
    do_excluded_boot    <- FALSE
    do_community_boot   <- FALSE
    do_loadings_boot    <- FALSE
  } else {
    do_general_boot   <- "general_index"  %in% boot_what
    do_bridge_boot    <- "bridge_index"   %in% boot_what
    do_excluded_boot  <- "excluded_index" %in% boot_what
    do_community_boot <- "community"      %in% boot_what
    do_loadings_boot  <- isTRUE(compute_loadings) && ("loadings" %in% boot_what)
  }

  # ---- Bootstrap containers ----
  boot_memberships <- list()
  quantile_region_results <- NULL
  strength_boot <- ei1_boot <- closeness_boot <- betweenness_boot <- NULL
  bridge_strength_boot <- bridge_ei1_boot <- bridge_ei2_boot <- NULL
  bridge_betweenness_boot <- bridge_closeness_boot <- NULL
  bridge_strength_excl_boot <- bridge_betweenness_excl_boot <- NULL
  bridge_closeness_excl_boot <- bridge_ei1_excl_boot <- bridge_ei2_excl_boot <- NULL
  edge_boot_mat <- NULL
  community_loadings_boot <- NULL

  use_progress <- isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)
  seq_reps <- seq_len(reps)

  .boot_rep <- function(i, p = NULL) {
    if (!is.null(p)) {
      p(sprintf("Bootstrap %d/%d", i, reps))
    }

    index <- sample(seq_len(nrow(data_mgm)), replace = TRUE)
    boot_data <- data_mgm[index, , drop = FALSE]

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
      n_blocks <- 0L
      if (isTRUE(do_general_boot))  n_blocks <- n_blocks + 4L  # strength, EI1, closeness, betweenness
      if (isTRUE(do_bridge_boot))   n_blocks <- n_blocks + 5L  # bridge strength, EI1, EI2, betw, closeness
      if (isTRUE(do_excluded_boot)) n_blocks <- n_blocks + 5L  # excluded bridge metrics

      central_len <- n_blocks * n_nodes_graph

      return(list(
        membership   = if (isTRUE(do_community_boot)) rep(NA_integer_, length(keep_nodes_cluster)) else NULL,
        centralities = if (central_len > 0) rep(NA_real_, central_len) else numeric(0),
        edges_boot   = NULL,
        loadings_boot = NULL
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

    ## membership bootstrap only if requested
    boot_membership <- NULL
    if (isTRUE(do_community_boot)) {
      boot_membership <- tryCatch({
        m <- cluster_fun(g_cluster_boot)$membership
        stats::setNames(m, keep_nodes_cluster)
      }, error = function(e) rep(NA_integer_, length(keep_nodes_cluster)))
    }

    ## centrality vector
    centralities_vec <- numeric(0)

    ## --- GENERAL INDEX ---
    if (isTRUE(do_general_boot)) {
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
          igraph::betweenness(g_dist_boot, weights = igraph::E(g_dist_boot)$dist,
                              directed = FALSE, normalized = FALSE)[keep_nodes_graph]
        else rep(0, n_nodes_graph),
        error = function(e) rep(NA_real_, n_nodes_graph)
      )

      centralities_vec <- c(
        centralities_vec,
        cent_vals$OutDegree,
        cent_vals$OutExpectedInfluence,
        igraph_closeness_boot,
        igraph_betweenness_boot
      )
    }

    ## --- BRIDGE INDEX (nodes in community) ---
    if (isTRUE(do_bridge_boot)) {
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

      centralities_vec <- c(
        centralities_vec,
        bridge_vals$bridge_strength,
        bridge_vals$bridge_ei1,
        bridge_vals$bridge_ei2,
        bridge_vals$bridge_betweenness,
        bridge_vals$bridge_closeness
      )
    }

    ## --- EXCLUDED INDEX ---
    if (isTRUE(do_excluded_boot)) {
      bridge_outside_boot <- tryCatch({
        abs_part <- bridge_metrics_excluded(g_bridge_abs_boot, membership = groups)
        sgn_part <- bridge_metrics_excluded(g_bridge_signed_boot, membership = groups)
        cbind(
          abs_part[match(keep_nodes_graph, abs_part$node),
                   c("bridge_strength","bridge_betweenness","bridge_closeness")],
          sgn_part[match(keep_nodes_graph, sgn_part$node),
                   c("bridge_ei1","bridge_ei2")]
        )
      }, error = function(e) {
        matrix(NA_real_, nrow = n_nodes_graph, ncol = 5)
      })

      centralities_vec <- c(
        centralities_vec,
        bridge_outside_boot[, 1],
        bridge_outside_boot[, 2],
        bridge_outside_boot[, 3],
        bridge_outside_boot[, 4],
        bridge_outside_boot[, 5]
      )
    }

    ## edges
    edge_values <- apply(pairs, 2, function(x) boot_wadj_signed_graph[x[1], x[2]])
    names(edge_values) <- edge_names

    # --- community loadings bootstrap (only if requested) ---
    loadings_boot <- NULL
    if (isTRUE(do_loadings_boot) && length(nodes_comm) > 1) {

      A_comm_boot <- boot_wadj_signed_graph[nodes_comm, nodes_comm, drop = FALSE]

      loads_obj_b <- tryCatch(
        .quiet_net_loads(
          A = A_comm_boot,
          wc = wc_comm_int,
          loading.method = "revised",
          rotation = NULL
        ),
        error = function(e) NULL
      )

      if (!is.null(loads_obj_b) && !is.null(loads_obj_b$std)) {
        loadings_boot <- loads_obj_b$std
        loadings_boot <- loadings_boot[nodes_comm, , drop = FALSE]
      } else {
        K <- if (!is.null(community_loadings_true)) ncol(community_loadings_true) else length(unique(wc_comm_int))
        loadings_boot <- matrix(NA_real_, nrow = length(nodes_comm), ncol = K,
                                dimnames = list(nodes_comm,
                                                if (!is.null(community_loadings_true)) colnames(community_loadings_true) else paste0("C", seq_len(K))))
      }
    }

    list(
      membership   = boot_membership,
      centralities = centralities_vec,
      edges_boot   = edge_values,
      loadings_boot = loadings_boot
    )
  }

  if (reps > 0) {
    t0 <- Sys.time()

    if (use_progress) {
      boot_output <- progressr::with_progress({
        p <- progressr::progressor(steps = reps)

        future.apply::future_lapply(
          X = seq_reps,
          FUN = .boot_rep,
          p = p,
          future.seed = seed_boot
        )
      })
    } else {
      boot_output <- future.apply::future_lapply(
        X = seq_reps,
        FUN = .boot_rep,
        p = NULL,
        future.seed = seed_boot
      )
    }

    elapsed <- difftime(Sys.time(), t0, units = "secs")
    message(sprintf(
      "Total computation time: %.1f seconds (~ %.2f minutes).",
      as.numeric(elapsed),
      as.numeric(elapsed) / 60
    ))

    ## bootstrap matrix centrality
    boot_mat <- do.call(rbind, lapply(boot_output, function(x) x$centralities))

    ## memberships (only if requested)
    if (isTRUE(do_community_boot)) {
      boot_memberships <- lapply(boot_output, function(x) x$membership)
    } else {
      boot_memberships <- list()
    }

    ## edges
    edge_boot_mat <- matrix(NA_real_, nrow = length(edge_names), ncol = reps)
    rownames(edge_boot_mat) <- edge_names
    for (i in seq_len(reps)) {
      edges_i <- boot_output[[i]]$edges_boot
      if (!is.null(edges_i)) edge_boot_mat[names(edges_i), i] <- edges_i
    }

    ## slicing boot_mat using boot_what
    n <- n_nodes_graph
    offset <- 0

    if (isTRUE(do_general_boot) && !is.null(boot_mat) && ncol(boot_mat) > 0) {
      strength_boot    <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      ei1_boot         <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      closeness_boot   <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      betweenness_boot <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n

      colnames(strength_boot)    <- keep_nodes_graph
      colnames(ei1_boot)         <- keep_nodes_graph
      colnames(closeness_boot)   <- keep_nodes_graph
      colnames(betweenness_boot) <- keep_nodes_graph
    }

    if (isTRUE(do_bridge_boot) && !is.null(boot_mat) && ncol(boot_mat) >= offset + 5 * n) {
      bridge_strength_boot    <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      bridge_ei1_boot         <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      bridge_ei2_boot         <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      bridge_betweenness_boot <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      bridge_closeness_boot   <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n

      colnames(bridge_strength_boot)    <- keep_nodes_graph
      colnames(bridge_ei1_boot)         <- keep_nodes_graph
      colnames(bridge_ei2_boot)         <- keep_nodes_graph
      colnames(bridge_betweenness_boot) <- keep_nodes_graph
      colnames(bridge_closeness_boot)   <- keep_nodes_graph
    }

    if (isTRUE(do_excluded_boot) && !is.null(boot_mat) && ncol(boot_mat) >= offset + 5 * n) {
      bridge_strength_excl_boot    <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      bridge_betweenness_excl_boot <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      bridge_closeness_excl_boot   <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      bridge_ei1_excl_boot         <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n
      bridge_ei2_excl_boot         <- boot_mat[, (offset + 1):(offset + n), drop = FALSE]; offset <- offset + n

      colnames(bridge_strength_excl_boot)    <- keep_nodes_graph
      colnames(bridge_betweenness_excl_boot) <- keep_nodes_graph
      colnames(bridge_closeness_excl_boot)   <- keep_nodes_graph
      colnames(bridge_ei1_excl_boot)         <- keep_nodes_graph
      colnames(bridge_ei2_excl_boot)         <- keep_nodes_graph
    }

    ## quantile region function
    calc_quantile_region <- function(mat, probs) {
      quantile_region <- apply(mat, 2, function(x) {
        if (all(is.na(x))) {
          setNames(c(NA_real_, NA_real_), paste0(100*probs, "%"))
        } else {
          stats::quantile(x, probs = probs, na.rm = TRUE, names = TRUE)
        }
      })
      t(quantile_region)
    }

    quantile_region_results <- list()

    if (isTRUE(do_general_boot) && !is.null(strength_boot)) {
      quantile_region_results$strength           <- calc_quantile_region(strength_boot, probs)
      quantile_region_results$expected_influence <- calc_quantile_region(ei1_boot, probs)
      quantile_region_results$closeness          <- calc_quantile_region(closeness_boot, probs)
      quantile_region_results$betweenness        <- calc_quantile_region(betweenness_boot, probs)
    }

    if (isTRUE(do_bridge_boot) && !is.null(bridge_strength_boot)) {
      quantile_region_results$bridge_strength    <- calc_quantile_region(bridge_strength_boot, probs)
      quantile_region_results$bridge_betweenness <- calc_quantile_region(bridge_betweenness_boot, probs)
      quantile_region_results$bridge_closeness   <- calc_quantile_region(bridge_closeness_boot, probs)
      quantile_region_results$bridge_ei1         <- calc_quantile_region(bridge_ei1_boot, probs)
      quantile_region_results$bridge_ei2         <- calc_quantile_region(bridge_ei2_boot, probs)
    }

    if (isTRUE(do_excluded_boot) && !is.null(bridge_strength_excl_boot)) {
      quantile_region_results$bridge_strength_excluded    <- calc_quantile_region(bridge_strength_excl_boot, probs)
      quantile_region_results$bridge_betweenness_excluded <- calc_quantile_region(bridge_betweenness_excl_boot, probs)
      quantile_region_results$bridge_closeness_excluded   <- calc_quantile_region(bridge_closeness_excl_boot, probs)
      quantile_region_results$bridge_ei1_excluded         <- calc_quantile_region(bridge_ei1_excl_boot, probs)
      quantile_region_results$bridge_ei2_excluded         <- calc_quantile_region(bridge_ei2_excl_boot, probs)
    }

    ## edges
    quantile_region_results$edge_weights <- calc_quantile_region(t(edge_boot_mat), probs)

    # --- collect bootstrap loadings (one matrix per replication) ---
    if (isTRUE(do_loadings_boot) && !is.null(community_loadings_true)) {
      community_loadings_boot <- lapply(boot_output, function(x) x$loadings_boot)
      community_loadings_boot <- lapply(community_loadings_boot, function(L) {
        if (is.null(L)) return(community_loadings_true * NA_real_)
        L
      })
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

  # ---- quantile regions organized for node/edge (NULL if no bootstrap) ----
  node_quantile_region <- if (!is.null(quantile_region_results)) {
    list(
      strength                    = quantile_region_results$strength,
      expected_influence          = quantile_region_results$expected_influence,
      closeness                   = quantile_region_results$closeness,
      betweenness                 = quantile_region_results$betweenness,
      bridge_strength             = quantile_region_results$bridge_strength,
      bridge_betweenness          = quantile_region_results$bridge_betweenness,
      bridge_closeness            = quantile_region_results$bridge_closeness,
      bridge_ei1                  = quantile_region_results$bridge_ei1,
      bridge_ei2                  = quantile_region_results$bridge_ei2,
      bridge_strength_excluded    = quantile_region_results$bridge_strength_excluded,
      bridge_betweenness_excluded = quantile_region_results$bridge_betweenness_excluded,
      bridge_closeness_excluded   = quantile_region_results$bridge_closeness_excluded,
      bridge_ei1_excluded         = quantile_region_results$bridge_ei1_excluded,
      bridge_ei2_excluded         = quantile_region_results$bridge_ei2_excluded
    )
  } else {
    NULL
  }

  edge_quantile_region <- if (!is.null(quantile_region_results)) quantile_region_results$edge_weights else NULL

  # ---- Return ----
  fit <- list(
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

    data_info = list(
      mgm_type_level     = var_interpretation,
      binary_recode_map  = binary_recode_map
    ),

    model = list(
      mgm   = mgm_model,
      nodes = all_nodes,
      n     = nrow(data),
      p     = ncol(data),
      data  = if (isTRUE(save_data)) data else NULL
    ),

    graph = list(
      igraph            = g_igraph,
      keep_nodes_graph   = keep_nodes_graph,
      keep_nodes_cluster = keep_nodes_cluster
    ),

    communities = list(
      original_membership = original_membership,   # on keep_nodes_cluster
      groups              = groups,               # factor (without singletons)
      palette             = palette_clusters,
      boot_memberships    = boot_memberships
    ),

    statistics = list(
      node = list(
        true = centrality_true_df,
        boot = list(
          strength                    = strength_boot,
          ei1                         = ei1_boot,
          closeness                   = closeness_boot,
          betweenness                 = betweenness_boot,
          bridge_strength             = bridge_strength_boot,
          bridge_betweenness          = bridge_betweenness_boot,
          bridge_closeness            = bridge_closeness_boot,
          bridge_ei1                  = bridge_ei1_boot,
          bridge_ei2                  = bridge_ei2_boot,
          bridge_strength_excluded    = bridge_strength_excl_boot,
          bridge_betweenness_excluded = bridge_betweenness_excl_boot,
          bridge_closeness_excluded   = bridge_closeness_excl_boot,
          bridge_ei1_excluded         = bridge_ei1_excl_boot,
          bridge_ei2_excluded         = bridge_ei2_excl_boot
        ),
        quantile_region = node_quantile_region
      ),
      edge = list(
        true = edges_true_df,
        boot = edge_boot_mat,
        quantile_region   = edge_quantile_region
      )
    ),

    community_loadings = list(
      nodes = nodes_comm,
      wc    = wc_comm_int,
      true  = community_loadings_true,
      boot  = community_loadings_boot,
      available = isTRUE(loadings_available),
      reason = loadings_reason,
      non_scorable_nodes = non_scorable_nodes
    )
  )

  class(fit) <- c("mixmashnet", "mixMN_fit")
  return(fit)
}
