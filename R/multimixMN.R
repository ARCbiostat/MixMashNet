#' Multilayer MGM with bootstrap, intra/interlayer metrics, and CIs
#'
#' @description
#' Estimates a multilayer Mixed Graphical Model (MGM) using the estimation
#' framework implemented in the \pkg{mgm} package, with a masking scheme that
#' enforces which cross-layer edges are allowed according to \code{layer_rules}.
#' Within each layer, the function computes community structure and performs
#' non-parametric row-bootstrap to obtain node centrality indices, edge weights,
#' and bridge metrics, including metrics for nodes treated as excluded. Optionally,
#' within-layer community loadings can also be estimated and bootstrapped.
#' The function additionally returns interlayer-only node metrics and summaries
#' of cross-layer edge weights.
#'
#' @param data A data.frame (n x p) with variables in columns.
#'   Character variables are not allowed and must be converted to \code{factor}
#'   or \code{numeric} types before fitting.
#' @param layers A named vector (names = variable names) assigning each node to a
#'   layer (character or factor). Must cover all columns of \code{data}
#'   except variables listed in \code{covariates} (treated as adjustment covariates).
#' @param layer_rules A logical or numeric square matrix with row/column names
#'   equal to layer names. Values \code{TRUE} or \code{1} indicate that
#'   cross-layer edges are allowed between the corresponding layer pair.
#'   Intralayer edges are always allowed.
#' @param scale Logical; if \code{TRUE} (default) Gaussian variables
#'   (\code{type == "g"}) are z-standardized internally by \code{mgm()}. Use
#'   \code{scale = FALSE} if your data are already standardized.
#' @param reps Integer (>= 0). Number of bootstrap replications (row resampling
#'   with replacement). If \code{reps = 0}, no bootstrap is performed.
#' @param lambdaSel Method for lambda selection in \code{mgm}:
#'   \code{"CV"} or \code{"EBIC"}.
#' @param lambdaFolds Number of folds for CV (if \code{lambdaSel = "CV"}).
#' @param lambdaGam EBIC gamma parameter (if \code{lambdaSel = "EBIC"}).
#' @param alphaSeq Alpha parameters of the elastic net penalty (values between 0 and 1).
#' @param alphaSel Method for selecting the alpha parameter:
#'   \code{"CV"} or \code{"EBIC"}.
#' @param alphaFolds Number of folds for CV (if \code{alphaSel = "CV"}).
#' @param alphaGam EBIC gamma parameter (if \code{alphaSel = "EBIC"}).
#' @param k Integer (>= 1). Order of modeled interactions.
#' @param ruleReg Rule to combine neighborhood estimates:
#'   \code{"AND"} or \code{"OR"}.
#' @param threshold Threshold below which edge-weights are set to zero:
#'   \code{"LW"}, \code{"HW"} or \code{"none"}.
#' @param overparameterize Logical; if \code{TRUE} uses the over-parameterized
#'   version of \code{mgm}.
#' @param thresholdCat Logical; if \code{FALSE} thresholds of categorical
#'   variables are set to zero.
#' @param conf_level Confidence level for percentile bootstrap CIs (default 0.95).
#'   Must be a single number between 0 and 1 (e.g., 0.90, 0.95, 0.99).
#' @param covariates Character vector. Variables used as adjustment covariates
#'   in model estimation.
#' @param exclude_from_cluster Character vector of node names. Nodes in this set
#'   are excluded from community detection in addition to \code{covariates}.
#' @param seed_model Optional integer seed for reproducibility of the initial
#'   MGM fit.
#' @param seed_boot Optional integer seed passed to \code{future.apply} for
#'   reproducibility of bootstrap replications.
#' @param treat_singletons_as_excluded Logical; if \code{TRUE}, singleton
#'   communities (size 1) are treated as "excluded" when computing bridge
#'   metrics and related summaries.
#' @param cluster_method Community detection algorithm applied within each
#'   layer. One of \code{"louvain"}, \code{"fast_greedy"}, \code{"infomap"},
#'   \code{"walktrap"}, or \code{"edge_betweenness"}.
#' @param compute_loadings Logical; if \code{TRUE} (default),
#'   compute network loadings (EGAnet net.loads) for communities.
#' @param boot_what Character vector specifying which quantities to bootstrap.
#'   Valid options are:
#'   \code{"general_index"} (intralayer centrality indices),
#'   \code{"interlayer_index"} (interlayer-only node metrics),
#'   \code{"bridge_index"} (bridge metrics for nodes in communities),
#'   \code{"excluded_index"} (bridge metrics for nodes treated as excluded),
#'   \code{"community"} (community memberships),
#'   \code{"loadings"} (within-layer community loadings, only if
#'   \code{compute_loadings = TRUE}),
#'   and \code{"none"} (skip all node-level bootstrap: only edge-weight
#'   bootstrap is performed if \code{reps > 0}).
#' @param save_data Logical; if \code{TRUE}, store the original data in the output object.
#' @param progress Logical; if \code{TRUE} (default), show a bootstrap progress bar.
#'
#' @return
#' An object of class \code{c("mixmashnet", "multimixMN_fit")}. The returned
#' list contains at least the following components:
#' \describe{
#'   \item{\code{call}}{
#'   The matched function call.
#'   }
#'  \item{\code{data_info}}{
#'     List with information derived from the input data used for model setup:
#'     \code{mgm_type_level} (data frame with one row per variable, reporting
#'     the original R class and the inferred MGM \code{type} and \code{level},
#'     as used in the call to \code{mgm::mgm}),
#'     and \code{binary_recode_map} (named list describing the mapping from
#'     original binary labels to the internal \{0,1\} coding used for model fitting).
#'   }
#'   \item{\code{settings}}{
#'   List of main settings used in the call, including
#'    \code{reps}, \code{cluster_method}, \code{covariates},
#'     \code{exclude_from_cluster}, \code{treat_singletons_as_excluded},
#'     \code{boot_what}).
#'   }
#'   \item{\code{model}}{
#'     List with:
#'     \code{mgm} (the fitted \code{mgm} object),
#'     \code{nodes} (character vector of all node names),
#'     \code{n} (number of observations),
#'     \code{p} (number of variables), and
#'     \code{data}.
#'   }
#'   \item{\code{layers}}{
#'    List describing the multilayer structure
#'     (assignment of nodes to layers, \code{layer_rules} matrix used and color of each layer in \code{palette}).
#'   }
#'   \item{\code{layer_fits}}{
#'    Named list (one element per layer) with
#'     single layer fits, including community structure, node-level statistics,
#'     edge-level statistics, bridge metrics, and (optionally) community loadings
#'     with bootstrap information.
#'   }
#'   \item{\code{interlayer}}{
#'   List collecting interlayer-only node metrics
#'   (strength, expected influence, closeness, betweenness, with or without
#'   bootstrap) and cross-layer edge summaries for each allowed pair of
#'   layers.
#'   }
#'   \item{\code{graph}}{
#'   List containing a global \pkg{igraph} object built on
#'   the retained nodes (\code{keep_nodes_graph}), with vertex attributes
#'   such as \code{name}, \code{layer}, \code{membership}, and edge attributes
#'   such as \code{weight}, \code{abs_weight}, \code{sign},
#'   \code{type} (intra vs inter) and \code{layer_pair}.
#'   }
#' }
#'
#' @details
#' This function does \strong{not} call \code{future::plan()}. To enable
#' parallel bootstrap, set a plan (e.g. \code{future::plan(multisession)}) before
#' calling \code{multimixMN()}. If \code{"none"} is the only element of
#' \code{boot_what} and \code{reps > 0}, node-level metrics are not
#' bootstrapped, but intra and interlayer edge-weight bootstrap and the
#' corresponding confidence intervals are still computed.
#'
#' @references
#'
#' Haslbeck, J. M. B., & Waldorp, L. J. (2020).
#' mgm: Estimating Time-Varying Mixed Graphical Models in High-Dimensional Data.
#' \emph{Journal of Statistical Software}, 93(8).
#' \doi{10.18637/jss.v093.i08}
#'
#' @importFrom stats setNames quantile
#' @importFrom utils combn capture.output
#' @importFrom igraph graph_from_adjacency_matrix simplify ecount E V distances betweenness vcount
#' @importFrom igraph cluster_louvain cluster_fast_greedy cluster_infomap cluster_walktrap cluster_edge_betweenness
#' @importFrom future.apply future_lapply
#' @importFrom qgraph centrality
#' @importFrom progressr with_progress progressor
#' @importFrom EGAnet net.loads
#' @export
multimixMN <- function(
    data,
    layers, layer_rules,
    scale = TRUE,
    reps = 100,
    lambdaSel = c("CV", "EBIC"),
    lambdaFolds = 5, lambdaGam = 0.25,
    alphaSeq = 1,
    alphaSel = "CV",
    alphaFolds = 5, alphaGam = 0.25,
    k = 2, ruleReg = "AND", threshold = "LW",
    overparameterize = FALSE, thresholdCat = TRUE,
    conf_level = 0.95,
    covariates = NULL,
    exclude_from_cluster = NULL,
    seed_model = NULL, seed_boot = NULL,
    treat_singletons_as_excluded = FALSE,
    cluster_method = c("louvain","fast_greedy","infomap","walktrap","edge_betweenness"),
    compute_loadings = TRUE,
    boot_what = c("general_index", "interlayer_index", "bridge_index",
                  "excluded_index", "community", "loadings"),
    save_data = FALSE,
    progress = TRUE
) {
  lambdaSel <- match.arg(lambdaSel)
  cluster_method <- match.arg(cluster_method)
  if (!is.null(seed_model)) set.seed(seed_model)

  spec <- infer_mgm_spec(data, recode_binary = TRUE)

  type     <- spec$type
  level    <- spec$level
  data_mgm <- spec$data_mgm

  var_interpretation <- spec$data_info
  binary_recode_map  <- spec$binary_recode_map

  # confidence interval
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number strictly between 0 and 1 (e.g., 0.95).")
  }
  alpha <- 1 - conf_level
  probs <- c(alpha/2, 1 - alpha/2)

  # --- parse 'boot_what' argument ---
  boot_what <- match.arg(
    boot_what,
    choices = c("general_index", "interlayer_index",
                "bridge_index", "excluded_index",
                "community", "loadings", "none"),
    several.ok = TRUE
  )

  if (length(boot_what) > 1L && "none" %in% boot_what) {
    boot_what <- setdiff(boot_what, "none")
  }

  if ("none" %in% boot_what && length(boot_what) == 1L) {
    do_intra_general_boot <- FALSE
    do_interlayer_boot    <- FALSE
    do_bridge_boot        <- FALSE
    do_excluded_boot      <- FALSE
    do_community_boot     <- FALSE
    do_loadings_boot      <- FALSE
  } else {
    do_intra_general_boot <- "general_index"    %in% boot_what
    do_interlayer_boot    <- "interlayer_index" %in% boot_what
    do_bridge_boot        <- "bridge_index"     %in% boot_what
    do_excluded_boot      <- "excluded_index"   %in% boot_what
    do_community_boot     <- "community"        %in% boot_what
    do_loadings_boot      <- isTRUE(compute_loadings) && ("loadings" %in% boot_what)
  }

  all_nodes <- colnames(data); p <- ncol(data)

  if (is.null(covariates)) covariates <- character(0)
  covariates <- unique(covariates)
  covariates <- intersect(covariates, all_nodes)
  network_nodes <- setdiff(all_nodes, covariates)

  if (is.null(exclude_from_cluster)) exclude_from_cluster <- character(0)
  exclude_from_cluster <- unique(as.character(exclude_from_cluster))
  exclude_from_cluster <- setdiff(exclude_from_cluster, covariates)

  subject_ids <- rownames(data)
  if (is.null(subject_ids)) {
    subject_ids <- sprintf("id_%d", seq_len(nrow(data)))
    rownames(data) <- subject_ids
  }

  # --- Basic multilayer check
  if (length(unique(stats::na.omit(layers[network_nodes]))) < 2) {
    stop(
      paste0(
        "multimixMN is designed for MULTILAYER networks (>= 2 layers). ",
        "You provided only one layer. For single layer networks, please use the ",
        "`mixMN()` function instead."
      ),
      call. = FALSE
    )
  }

  # layers must cover all network nodes, but NOT covariates
  if (is.null(names(layers))) {
    stop("`layers` must be a named vector with names = variable names.")
  }

  if (!all(network_nodes %in% names(layers))) {
    miss <- setdiff(network_nodes, names(layers))
    stop("`layers` must cover all variables except `covariates`. Missing: ",
         paste(miss, collapse = ", "))
  }

  if (any(covariates %in% names(layers))) {
    layers <- layers[setdiff(names(layers), covariates)]
  }

  # --- layer_rules validation and normalization
  stopifnot(is.matrix(layer_rules))
  layer_order <- unique(layers[intersect(network_nodes, names(layers))])
  layer_order <- layer_order[!is.na(layer_order)]

  if (is.null(rownames(layer_rules)) || is.null(colnames(layer_rules)) ||
      !all(layer_order %in% rownames(layer_rules)) ||
      !all(layer_order %in% colnames(layer_rules))) {
    stop("`layer_rules` must have row/col names matching the provided layer names.")
  }

  layer_rules <- layer_rules[layer_order, layer_order, drop = FALSE]
  lr_bin <- (layer_rules == 1)
  layer_rules_sym <- lr_bin | t(lr_bin)

  # diagonal: always allowed (treat NA as allowed)
  diag_in <- diag(layer_rules)
  d <- diag(layer_rules_sym)
  d[is.na(diag_in)] <- TRUE
  diag(layer_rules_sym) <- d

  layer_rules <- layer_rules_sym

  # ---------- helpers ----------
  tiny <- 1e-10

  .quiet_net_loads <- function(...) {
    withCallingHandlers(
      {
        out <- NULL
        invisible(capture.output(out <- EGAnet::net.loads(...)))
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

  .make_distance_graph <- function(W_signed) {
    W <- abs(W_signed); diag(W) <- 0; W[is.na(W) | W < tiny] <- 0
    g <- igraph::graph_from_adjacency_matrix(W, mode="undirected", weighted=TRUE, diag=FALSE)
    if (igraph::ecount(g) > 0) igraph::E(g)$dist <- 1/igraph::E(g)$weight
    g
  }
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
  .calc_ci <- function(mat, probs) {
    if (is.null(mat)) return(NULL)
    ci <- apply(mat, 2, function(x) {
      if (all(is.na(x))) {
        setNames(c(NA_real_, NA_real_), paste0(100*probs, "%"))
      } else {
        stats::quantile(x, probs = probs, na.rm = TRUE, names = TRUE)
      }
    })
    t(ci)
  }
  .edge_names_lt <- function(nodes) utils::combn(nodes, 2, FUN=function(x) paste(x[1],x[2],sep="--"))
  .edge_names_cross <- function(A,B) as.vector(outer(A,B,function(a,b) paste(a,b,sep="--")))
  .align_membership <- function(membership, nodes) {
    out <- stats::setNames(rep(NA_integer_, length(nodes)), nodes)
    if (is.null(membership)) return(out)
    if (is.factor(membership)) membership <- stats::setNames(as.integer(membership), names(membership))
    if (is.null(names(membership))) names(membership) <- nodes[seq_len(min(length(nodes), length(membership)))]
    common <- intersect(names(membership), nodes)
    out[common] <- as.integer(membership[common])
    out
  }
  .pick_col <- function(df, candidates, default = NA_real_) {
    for (nm in candidates) if (nm %in% colnames(df)) return(df[[nm]])
    rep(default, nrow(df))
  }
  .as_chr <- function(x) if (is.factor(x)) as.character(x) else as.character(x)
  cluster_fun <- function(graph) {
    switch(cluster_method,
           louvain          = igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight),
           fast_greedy      = igraph::cluster_fast_greedy(graph, weights = igraph::E(graph)$weight),
           infomap          = igraph::cluster_infomap(graph),
           walktrap         = igraph::cluster_walktrap(graph, weights = igraph::E(graph)$weight),
           edge_betweenness = igraph::cluster_edge_betweenness(graph, weights = igraph::E(graph)$weight)
    )
  }

  # --- Build per-node masks from layer_rules
  mask_list <- vector("list", length(all_nodes))
  names(mask_list) <- all_nodes

  for (i in seq_along(all_nodes)) {
    v <- all_nodes[i]
    if (v %in% covariates) {
      allowed_nodes <- all_nodes
    } else {
      my_layer <- layers[v]
      allowed_layers <- names(which(layer_rules[my_layer, ]))
      allowed_layers <- unique(c(my_layer, allowed_layers))
      allowed_network <- names(layers)[layers %in% allowed_layers]
      allowed_nodes <- unique(c(allowed_network, covariates))
    }
    mask_list[[i]] <- match(allowed_nodes, all_nodes)
  }

  # --- MGM (masked) fit
  mgm_model <- mgm_masked(
    data = data_mgm, type = type, level = level,
    lambdaSel = lambdaSel, lambdaFolds = lambdaFolds, lambdaGam = lambdaGam,
    alphaSeq = alphaSeq, alphaSel = alphaSel, alphaFolds = alphaFolds, alphaGam = alphaGam,
    ruleReg = ruleReg, threshold = threshold, overparameterize = overparameterize,
    thresholdCat = thresholdCat,
    k = k, binarySign = TRUE, mask_list = mask_list, scale = scale,
    signInfo = FALSE, pbar = FALSE
  )
  wadj <- mgm_model$pairwise$wadj; signs <- mgm_model$pairwise$signs
  colnames(wadj) <- rownames(wadj) <- all_nodes
  colnames(signs) <- rownames(signs) <- all_nodes
  wadj_signed <- wadj
  if (!is.null(signs)) {
    idx_sign <- !is.na(signs) & (abs(signs) == 1)
    wadj_signed[idx_sign] <- wadj[idx_sign] * signs[idx_sign]
  }
  wadj_signed[is.na(wadj_signed)] <- 0

  # --- Inclusion sets
  keep_nodes_graph_all   <- setdiff(all_nodes, covariates)
  keep_nodes_cluster_all <- setdiff(all_nodes, unique(c(covariates, exclude_from_cluster)))
  uniq_layers <- unique(layers[network_nodes])

  layer_nodes_graph <- stats::setNames(lapply(uniq_layers, function(L) {
    intersect(names(layers)[layers == L], network_nodes)
  }), uniq_layers)

  layer_nodes_cluster <- stats::setNames(lapply(uniq_layers, function(L) {
    intersect(names(layers)[layers == L], setdiff(network_nodes, exclude_from_cluster))
  }), uniq_layers)

  # --- Per-layer fits from W (true graph, no bootstrap yet)
  layer_colors <- character(length(uniq_layers))
  names(layer_colors) <- uniq_layers
  layer_fits <- list()
  for (i in seq_along(uniq_layers)) {
    L <- uniq_layers[i]
    pal_L <- .get_layer_palette(i)
    layer_colors[L] <- pal_L[1]
    nL <- layer_nodes_graph[[L]]
    sub_w <- wadj_signed[nL, nL, drop = FALSE]
    fitL <- mixMN_from_wadj(
      wadj_signed = sub_w,
      nodes = nL,
      conf_level = conf_level,
      covariates   = NULL,
      exclude_from_cluster = intersect(exclude_from_cluster, nL),
      cluster_method = cluster_method,
      reps          = 0,
      seed_boot     = NULL,
      treat_singletons_as_excluded = treat_singletons_as_excluded,
      boot_what     = boot_what,
      palette_layer = pal_L
    )

    membL <- fitL$communities$groups
    if (is.factor(membL)) {
      membL <- stats::setNames(as.integer(membL), names(membL))
    }
    fitL$communities$groups <- membL

    # ---- COMMUNITY LOADINGS (TRUE, per-layer) ----
    fitL$community_loadings <- list(
      nodes = character(0),
      wc    = integer(0),
      true  = NULL,
      boot  = NULL
    )

    if (isTRUE(compute_loadings)) {

      groups_L <- fitL$communities$groups
      nodes_comm <- names(groups_L)[!is.na(groups_L)]
      nodes_comm <- intersect(nodes_comm, fitL$graph$keep_nodes_graph)

      if (length(nodes_comm) > 1) {

        A_comm <- sub_w[nodes_comm, nodes_comm, drop = FALSE]

        wc_fac <- groups_L[nodes_comm]
        wc_levels <- sort(unique(as.integer(wc_fac)))
        wc_map <- stats::setNames(seq_along(wc_levels), wc_levels)
        wc_int <- unname(wc_map[as.integer(wc_fac)])

        loads_obj <- tryCatch(
          .quiet_net_loads(
            A = A_comm,
            wc = wc_int,
            loading.method = "revised",
            rotation = NULL
          ),
          error = function(e) NULL
        )

        if (!is.null(loads_obj) && !is.null(loads_obj$std)) {
          L_true <- loads_obj$std
          L_true <- L_true[nodes_comm, , drop = FALSE]
        } else {
          K <- length(unique(wc_int))
          L_true <- matrix(NA_real_, nrow = length(nodes_comm), ncol = K,
                           dimnames = list(nodes_comm, paste0("C", seq_len(K))))
        }

        fitL$community_loadings <- list(
          nodes = nodes_comm,
          wc    = wc_int,
          true  = L_true,
          boot  = NULL
        )
      }
    }

    layer_fits[[L]] <- fitL
  }

  layer_nodes_excluded <- stats::setNames(
    lapply(uniq_layers, function(L) {
      membL <- layer_fits[[L]]$communities$groups
      setdiff(layer_nodes_graph[[L]], names(membL))
    }),
    uniq_layers
  )

  # --- Cross-layer pairs to track edges
  pairs <- list()
  Lvec <- uniq_layers
  if (!is.null(layer_rules)) {
    for (i in seq_along(Lvec)) for (j in seq_along(Lvec)) {
      if (j <= i) next
      if (isTRUE(layer_rules[Lvec[i], Lvec[j]])) {
        A <- layer_nodes_graph[[Lvec[i]]]; B <- layer_nodes_graph[[Lvec[j]]]
        if (length(A)>0 && length(B)>0) {
          key <- paste(Lvec[i],Lvec[j],sep="_")
          pairs[[key]] <- list(A=A, B=B, edge_names=.edge_names_cross(A,B))
        }
      }
    }
  }

  # --- Bootstrap containers (per-layer and interlayer)
  interlayer_boot <- lapply(names(pairs), function(k) {
    en <- pairs[[k]]$edge_names
    if (length(en)) matrix(NA_real_, nrow=length(en), ncol=reps, dimnames=list(en, NULL)) else NULL
  })
  names(interlayer_boot) <- names(pairs)

  layer_boot <- lapply(uniq_layers, function(L) {
    nodes_g <- layer_nodes_graph[[L]]
    membL   <- layer_fits[[L]]$communities$groups
    nodes_c <- if (!is.null(membL)) names(membL) else character(0)
    nodes_ex <- layer_nodes_excluded[[L]]
    list(
      # general indices intralayer
      strength_boot     = if (do_intra_general_boot)
        matrix(NA_real_, reps, length(nodes_g), dimnames = list(NULL, nodes_g)) else NULL,
      ei1_boot          = if (do_intra_general_boot)
        matrix(NA_real_, reps, length(nodes_g), dimnames = list(NULL, nodes_g)) else NULL,
      closeness_boot    = if (do_intra_general_boot)
        matrix(NA_real_, reps, length(nodes_g), dimnames = list(NULL, nodes_g)) else NULL,
      betweenness_boot  = if (do_intra_general_boot)
        matrix(NA_real_, reps, length(nodes_g), dimnames = list(NULL, nodes_g)) else NULL,

      # edges intralayer
      edge_boot_mat     = {
        en <- .edge_names_lt(nodes_g)
        if (length(en)) matrix(NA_real_, length(en), reps,
                               dimnames = list(en, NULL)) else NULL
      },

      # membership bootstrap (only if requested)
      boot_memberships  = if (do_community_boot) vector("list", reps) else vector("list", 0),

      # loadings bootstrap (list of matrices, one per replication)
      loadings_boot = if (isTRUE(do_loadings_boot) &&
                          !is.null(layer_fits[[L]]$community_loadings$true) &&
                          nrow(layer_fits[[L]]$community_loadings$true) > 0) {
        vector("list", reps)
      } else NULL,

      # bridge intralayer
      bridge_strength_boot     = if (do_bridge_boot && length(nodes_c) > 0)
        matrix(NA_real_, reps, length(nodes_c), dimnames = list(NULL, nodes_c)) else NULL,
      bridge_ei1_boot          = if (do_bridge_boot && length(nodes_c) > 0)
        matrix(NA_real_, reps, length(nodes_c), dimnames = list(NULL, nodes_c)) else NULL,
      bridge_ei2_boot          = if (do_bridge_boot && length(nodes_c) > 0)
        matrix(NA_real_, reps, length(nodes_c), dimnames = list(NULL, nodes_c)) else NULL,
      bridge_closeness_boot    = if (do_bridge_boot && length(nodes_c) > 0)
        matrix(NA_real_, reps, length(nodes_c), dimnames = list(NULL, nodes_c)) else NULL,
      bridge_betweenness_boot  = if (do_bridge_boot && length(nodes_c) > 0)
        matrix(NA_real_, reps, length(nodes_c), dimnames = list(NULL, nodes_c)) else NULL,

      # bridge for excluded nodes
      bridge_strength_excl_boot     = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL,
      bridge_ei1_excl_boot          = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL,
      bridge_ei2_excl_boot          = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL,
      bridge_closeness_excl_boot    = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL,
      bridge_betweenness_excl_boot  = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL
    )
  })
  names(layer_boot) <- uniq_layers

  # --- Interlayer-only TRUE metrics (from full signed W)
  nodes_int <- keep_nodes_graph_all
  W_inter_true <- matrix(0, length(nodes_int), length(nodes_int), dimnames=list(nodes_int, nodes_int))
  for (i in seq_along(nodes_int)) for (j in seq_along(nodes_int)) {
    vi <- nodes_int[i]; vj <- nodes_int[j]
    if (layers[vi] != layers[vj]) W_inter_true[vi, vj] <- wadj_signed[vi, vj]
  }
  g_inter_true <- .make_distance_graph(W_inter_true)
  inter_strength_true  <- stats::setNames(rowSums(abs(W_inter_true), na.rm=TRUE), nodes_int)
  inter_ei1_true <- tryCatch({
    cm <- qgraph::centrality(W_inter_true); stats::setNames(cm$OutExpectedInfluence, nodes_int)
  }, error=function(e) stats::setNames(rep(NA_real_, length(nodes_int)), nodes_int))
  inter_closeness_true <- .harmonic_closeness(g_inter_true)[nodes_int]
  inter_betweenness_true <- if (igraph::ecount(g_inter_true) > 0) {
    b <- igraph::betweenness(g_inter_true, weights=igraph::E(g_inter_true)$dist, directed=FALSE, normalized=FALSE)
    stats::setNames(b[nodes_int], nodes_int)
  } else stats::setNames(rep(0, length(nodes_int)), nodes_int)

  inter_strength_boot    <- if (do_interlayer_boot)
    matrix(NA_real_, reps, length(nodes_int), dimnames = list(NULL, nodes_int)) else NULL
  inter_ei1_boot         <- if (do_interlayer_boot)
    matrix(NA_real_, reps, length(nodes_int), dimnames = list(NULL, nodes_int)) else NULL
  inter_closeness_boot   <- if (do_interlayer_boot)
    matrix(NA_real_, reps, length(nodes_int), dimnames = list(NULL, nodes_int)) else NULL
  inter_betweenness_boot <- if (do_interlayer_boot)
    matrix(NA_real_, reps, length(nodes_int), dimnames = list(NULL, nodes_int)) else NULL

  use_progress <- isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)
  seq_reps <- seq_len(reps)

  .boot_rep <- function(bi, p = NULL) {
    if (!is.null(p)) {
      p(sprintf("Bootstrap %d/%d", bi, reps))
    }
    idx <- sample(seq_len(nrow(data_mgm)), replace = TRUE)
    Xb  <- data_mgm[idx, , drop = FALSE]
    boot_model <- tryCatch(mgm_masked(
      data = Xb, type = type, level = level,
      lambdaSel = lambdaSel, lambdaFolds = lambdaFolds, lambdaGam = lambdaGam,
      alphaSeq = alphaSeq, alphaSel = alphaSel, alphaFolds = alphaFolds, alphaGam = alphaGam,
      ruleReg = ruleReg, threshold = threshold, overparameterize = overparameterize,
      thresholdCat = thresholdCat,
      k = k, binarySign = TRUE, mask_list = mask_list, scale = scale,
      signInfo = FALSE
    ), error=function(e) NULL)
    if (is.null(boot_model)) return(NULL)

    bw   <- boot_model$pairwise$wadj
    bsig <- boot_model$pairwise$signs
    colnames(bw)   <- rownames(bw)   <- all_nodes
    colnames(bsig) <- rownames(bsig) <- all_nodes
    W <- bw
    if (!is.null(bsig)) {
      idxb_sign <- !is.na(bsig) & (abs(bsig) == 1)
      W[idxb_sign] <- bw[idxb_sign] * bsig[idxb_sign]
    }
    W[is.na(W)] <- 0

    per_layer <- lapply(uniq_layers, function(L) {
      nodes_g <- layer_nodes_graph[[L]]
      if (length(nodes_g) < 1) return(NULL)
      Wg <- W[nodes_g, nodes_g, drop=FALSE]
      nodes_c <- layer_nodes_cluster[[L]]
      nodes_c <- intersect(nodes_c, nodes_g)

      g_bridge_abs_boot    <- igraph::graph_from_adjacency_matrix(abs(Wg), mode="undirected", weighted=TRUE, diag=FALSE)
      g_bridge_abs_boot    <- igraph::simplify(g_bridge_abs_boot, remove.multiple=TRUE, remove.loops=TRUE)

      g_bridge_signed_boot <- igraph::graph_from_adjacency_matrix(Wg,      mode="undirected", weighted=TRUE, diag=FALSE)
      g_bridge_signed_boot <- igraph::simplify(g_bridge_signed_boot, remove.multiple=TRUE, remove.loops=TRUE)

      g_dist  <- .make_distance_graph(Wg)

      memb_boot <- NULL
      if (do_community_boot && length(nodes_c) > 0) {
        Wc <- Wg[nodes_c, nodes_c, drop = FALSE]
        g_c <- igraph::graph_from_adjacency_matrix(
          abs(Wc), mode = "undirected", weighted = TRUE, diag = FALSE
        )
        g_c <- igraph::simplify(g_c, remove.multiple = TRUE, remove.loops = TRUE)

        memb_boot <- tryCatch({
          cl <- cluster_fun(g_c)
          mb <- cl$membership
          names(mb) <- nodes_c
          mb
        }, error = function(e) {
          stats::setNames(rep(NA_integer_, length(nodes_c)), nodes_c)
        })
      }

      # --- LOADINGS bootstrap (per-layer) ---
      L_boot <- NULL
      if (isTRUE(do_loadings_boot)) {

        info <- layer_fits[[L]]$community_loadings
        nodes_comm <- info$nodes
        wc_int     <- info$wc

        if (!is.null(nodes_comm) && length(nodes_comm) > 1) {

          A_comm_boot <- Wg[nodes_comm, nodes_comm, drop = FALSE]

          loads_obj_b <- tryCatch(
            .quiet_net_loads(
              A = A_comm_boot,
              wc = wc_int,
              loading.method = "revised",
              rotation = NULL
            ),
            error = function(e) NULL
          )

          if (!is.null(loads_obj_b) && !is.null(loads_obj_b$std)) {
            L_boot <- loads_obj_b$std
            L_boot <- L_boot[nodes_comm, , drop = FALSE]
          } else {
            K <- if (!is.null(info$true)) ncol(info$true) else length(unique(wc_int))
            L_boot <- matrix(
              NA_real_,
              nrow = length(nodes_comm),
              ncol = K,
              dimnames = list(
                nodes_comm,
                if (!is.null(info$true) && !is.null(colnames(info$true))) colnames(info$true) else paste0("C", seq_len(K))
              )
            )
          }
        }
      }

      # --- general indices intralayer ---
      if (do_intra_general_boot) {
        cent <- tryCatch(
          qgraph::centrality(Wg),
          error = function(e) list(
            OutDegree            = rep(NA_real_, length(nodes_g)),
            OutExpectedInfluence = rep(NA_real_, length(nodes_g))
          )
        )
        clh <- tryCatch(
          .harmonic_closeness(g_dist)[nodes_g],
          error = function(e) rep(NA_real_, length(nodes_g))
        )
        btw <- tryCatch(
          {
            if (igraph::ecount(g_dist) > 0)
              igraph::betweenness(
                g_dist,
                weights   = igraph::E(g_dist)$dist,
                directed  = FALSE,
                normalized = FALSE
              )[nodes_g]
            else rep(0, length(nodes_g))
          },
          error = function(e) rep(NA_real_, length(nodes_g))
        )
      } else {
        cent <- list(
          OutDegree            = rep(NA_real_, length(nodes_g)),
          OutExpectedInfluence = rep(NA_real_, length(nodes_g))
        )
        clh <- rep(NA_real_, length(nodes_g))
        btw <- rep(NA_real_, length(nodes_g))
      }

      memb_orig <- layer_fits[[L]]$communities$groups
      memb_all  <- .align_membership(membership = memb_orig, nodes = nodes_g)

      ## --- BRIDGE ---
      bridge_vals <- NULL
      if (do_bridge_boot) {
        bridge_vals_abs <- tryCatch({
          b <- bridge_metrics(g_bridge_abs_boot, membership = memb_orig)
          b[match(nodes_g, b$node),
            c("bridge_strength","bridge_betweenness","bridge_closeness")]
        }, error = function(e) {
          as.data.frame(
            matrix(
              NA_real_, length(nodes_g), 3,
              dimnames = list(nodes_g, c("bridge_strength","bridge_betweenness","bridge_closeness"))
            )
          )
        })

        bridge_vals_signed <- tryCatch({
          b <- bridge_metrics(g_bridge_signed_boot, membership = memb_orig)
          b[match(nodes_g, b$node), c("bridge_ei1","bridge_ei2")]
        }, error = function(e) {
          as.data.frame(
            matrix(
              NA_real_, length(nodes_g), 2,
              dimnames = list(nodes_g, c("bridge_ei1","bridge_ei2"))
            )
          )
        })

        bridge_vals <- cbind(
          bridge_strength    = bridge_vals_abs[,"bridge_strength"],
          bridge_ei1         = bridge_vals_signed[,"bridge_ei1"],
          bridge_ei2         = bridge_vals_signed[,"bridge_ei2"],
          bridge_betweenness = bridge_vals_abs[,"bridge_betweenness"],
          bridge_closeness   = bridge_vals_abs[,"bridge_closeness"]
        )
        rownames(bridge_vals) <- nodes_g
      }

      ## --- BRIDGE "excluded"  ---
      bridge_excluded_mat <- NULL
      if (do_excluded_boot) {
        bridge_excluded_abs <- tryCatch({
          bo <- bridge_metrics_excluded(g_bridge_abs_boot, membership = memb_orig)
          idx <- match(nodes_g, bo$node)
          df_abs <- data.frame(
            node = nodes_g,
            bridge_strength    = bo$bridge_strength[idx],
            bridge_betweenness = bo$bridge_betweenness[idx],
            bridge_closeness   = bo$bridge_closeness[idx],
            stringsAsFactors   = FALSE
          )
          rownames(df_abs) <- nodes_g
          df_abs
        }, error = function(e) {
          df_abs <- data.frame(
            node = nodes_g,
            bridge_strength    = NA_real_,
            bridge_betweenness = NA_real_,
            bridge_closeness   = NA_real_,
            stringsAsFactors   = FALSE
          )
          rownames(df_abs) <- nodes_g
          df_abs
        })

        bridge_excluded_signed <- tryCatch({
          bo <- bridge_metrics_excluded(g_bridge_signed_boot, membership = memb_orig)
          ei1 <- .pick_col(bo, c("bridge_ei1","bridge_expected_influence1"))
          ei2 <- .pick_col(bo, c("bridge_ei2","bridge_expected_influence2"))
          idx <- match(nodes_g, bo$node)
          df_sgn <- data.frame(
            node = nodes_g,
            ei1  = ei1[idx],
            ei2  = ei2[idx],
            stringsAsFactors = FALSE
          )
          rownames(df_sgn) <- nodes_g
          df_sgn
        }, error = function(e) {
          data.frame(node = nodes_g, ei1 = NA_real_, ei2 = NA_real_)
        })

        bridge_excluded_mat <- matrix(
          NA_real_, nrow = length(nodes_g), ncol = 5,
          dimnames = list(nodes_g, c("strength","ei1","ei2","closeness","betweenness"))
        )

        row_map_abs <- intersect(nodes_g, rownames(bridge_excluded_abs))
        if (length(row_map_abs)) {
          bridge_excluded_mat[row_map_abs, "strength"]    <- bridge_excluded_abs[row_map_abs, "bridge_strength"]
          bridge_excluded_mat[row_map_abs, "closeness"]   <- bridge_excluded_abs[row_map_abs, "bridge_closeness"]
          bridge_excluded_mat[row_map_abs, "betweenness"] <- bridge_excluded_abs[row_map_abs, "bridge_betweenness"]
        }

        row_map_sgn <- intersect(nodes_g, bridge_excluded_signed$node)
        if (length(row_map_sgn)) {
          idx2 <- match(row_map_sgn, bridge_excluded_signed$node)
          bridge_excluded_mat[row_map_sgn, "ei1"] <- bridge_excluded_signed$ei1[idx2]
          bridge_excluded_mat[row_map_sgn, "ei2"] <- bridge_excluded_signed$ei2[idx2]
        }
      }

      en <- .edge_names_lt(nodes_g)
      evec <- if (length(en)) { v <- Wg[lower.tri(Wg)]; names(v) <- en; v } else NULL

      list(
        strength = cent$OutDegree, ei1 = cent$OutExpectedInfluence,
        closeness = clh, betweenness = btw, edges = evec,
        bridge_vals = bridge_vals, bridge_excluded_mat = bridge_excluded_mat,
        membership = memb_all,
        membership_boot = if (do_community_boot) memb_boot else NULL,
        loadings_boot = L_boot
      )
    })
    names(per_layer) <- uniq_layers

    # Interlayer-only metrics for this bootstrap
    inter_strength   <- NULL
    ei1_inter        <- NULL
    inter_clo        <- NULL
    inter_btw        <- NULL

    if (do_interlayer_boot) {
      W_inter <- matrix(0, length(nodes_int), length(nodes_int),
                        dimnames = list(nodes_int, nodes_int))
      for (i in seq_along(nodes_int)) for (j in seq_along(nodes_int)) {
        vi <- nodes_int[i]; vj <- nodes_int[j]
        if (layers[vi] != layers[vj]) W_inter[vi, vj] <- W[vi, vj]
      }
      g_inter <- .make_distance_graph(W_inter)
      inter_strength <- stats::setNames(rowSums(abs(W_inter), na.rm = TRUE), nodes_int)
      ei1_inter <- tryCatch({
        cm <- qgraph::centrality(W_inter)
        stats::setNames(cm$OutExpectedInfluence, nodes_int)
      }, error = function(e) stats::setNames(rep(NA_real_, length(nodes_int)), nodes_int))
      inter_clo <- .harmonic_closeness(g_inter)[nodes_int]
      inter_btw <- if (igraph::ecount(g_inter) > 0) {
        b <- igraph::betweenness(
          g_inter,
          weights   = igraph::E(g_inter)$dist,
          directed  = FALSE,
          normalized = FALSE
        )
        stats::setNames(b[nodes_int], nodes_int)
      } else {
        stats::setNames(rep(0, length(nodes_int)), nodes_int)
      }
    }

    # Cross-layer edges per allowed pairs
    per_pairs <- lapply(names(pairs), function(key) {
      A <- pairs[[key]]$A; B <- pairs[[key]]$B
      if (length(A)==0 || length(B)==0) return(NULL)
      cross <- W[A, B, drop=FALSE]
      v <- as.vector(cross); names(v) <- .edge_names_cross(A,B); v
    })
    names(per_pairs) <- names(pairs)

    list(
      per_layer = per_layer,
      inter_strength   = inter_strength,
      inter_ei1        = ei1_inter,
      inter_closeness  = inter_clo,
      inter_betweenness = inter_btw,
      per_pairs = per_pairs
    )
  }

  # ---------- BOOTSTRAP ----------
  if (reps > 0) {
    t0 <- Sys.time()

    if (use_progress) {
      boot_res <- progressr::with_progress({
        p <- progressr::progressor(steps = reps)
        future.apply::future_lapply(
          X = seq_reps,
          FUN = .boot_rep,
          p = p,
          future.seed = seed_boot
        )
      })
    } else {
      boot_res <- future.apply::future_lapply(
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

    # Collect bootstrap results
    for (bi in seq_len(reps)) {
      if (is.null(boot_res[[bi]])) next

      for (L in uniq_layers) {
        if (is.null(layer_boot[[L]])) next
        pl <- boot_res[[bi]]$per_layer[[L]]; if (is.null(pl)) next

        if (!is.null(layer_boot[[L]]$loadings_boot)) {
          layer_boot[[L]]$loadings_boot[[bi]] <- pl$loadings_boot
        }

        nodes_g <- layer_nodes_graph[[L]]
        nodes_c <- names(layer_fits[[L]]$communities$groups)
        if (is.null(nodes_c)) nodes_c <- character(0)
        nodes_ex <- layer_nodes_excluded[[L]]

        if (do_intra_general_boot && !is.null(layer_boot[[L]]$strength_boot)) {
          layer_boot[[L]]$strength_boot   [bi, nodes_g] <- pl$strength[nodes_g]
          layer_boot[[L]]$ei1_boot        [bi, nodes_g] <- pl$ei1[nodes_g]
          layer_boot[[L]]$closeness_boot  [bi, nodes_g] <- pl$closeness[nodes_g]
          layer_boot[[L]]$betweenness_boot[bi, nodes_g] <- pl$betweenness[nodes_g]
        }

        if (!is.null(layer_boot[[L]]$edge_boot_mat) && !is.null(pl$edges)) {
          en <- rownames(layer_boot[[L]]$edge_boot_mat)
          layer_boot[[L]]$edge_boot_mat[en, bi] <- pl$edges[en]
        }

        if (!is.null(layer_boot[[L]]$bridge_strength_boot)) {
          bdf <- pl$bridge_vals
          layer_boot[[L]]$bridge_strength_boot    [bi, nodes_c] <- bdf[nodes_c, "bridge_strength"]
          layer_boot[[L]]$bridge_ei1_boot         [bi, nodes_c] <- bdf[nodes_c, "bridge_ei1"]
          layer_boot[[L]]$bridge_ei2_boot         [bi, nodes_c] <- bdf[nodes_c, "bridge_ei2"]
          layer_boot[[L]]$bridge_closeness_boot   [bi, nodes_c] <- bdf[nodes_c, "bridge_closeness"]
          layer_boot[[L]]$bridge_betweenness_boot [bi, nodes_c] <- bdf[nodes_c, "bridge_betweenness"]
        }
        if (!is.null(layer_boot[[L]]$bridge_strength_excl_boot) && !is.null(pl$bridge_excluded_mat)) {
          bem <- pl$bridge_excluded_mat
          nn <- intersect(nodes_ex, rownames(bem))
          if (length(nn)) {
            layer_boot[[L]]$bridge_strength_excl_boot    [bi, nn] <- bem[nn, "strength"]
            layer_boot[[L]]$bridge_ei1_excl_boot         [bi, nn] <- bem[nn, "ei1"]
            layer_boot[[L]]$bridge_ei2_excl_boot         [bi, nn] <- bem[nn, "ei2"]
            layer_boot[[L]]$bridge_closeness_excl_boot   [bi, nn] <- bem[nn, "closeness"]
            layer_boot[[L]]$bridge_betweenness_excl_boot [bi, nn] <- bem[nn, "betweenness"]
          }
        }

        if (do_community_boot) {
          if (length(nodes_c)) {
            mm <- pl$membership_boot
            if (is.null(mm)) {
              mm <- stats::setNames(rep(NA_integer_, length(nodes_c)), nodes_c)
            }
            if (is.factor(mm)) {
              mm <- stats::setNames(as.integer(mm), names(mm))
            }
            layer_boot[[L]]$boot_memberships[[bi]] <- mm
          } else {
            layer_boot[[L]]$boot_memberships[[bi]] <- NULL
          }
        }
      }

      for (key in names(pairs)) {
        if (is.null(interlayer_boot[[key]])) next
        vp <- boot_res[[bi]]$per_pairs[[key]]; if (is.null(vp)) next
        en <- rownames(interlayer_boot[[key]])
        interlayer_boot[[key]][en, bi] <- vp[en]
      }

      # interlayer-only metrics
      if (do_interlayer_boot && !is.null(inter_strength_boot)) {
        inter_strength_boot   [bi, ] <- boot_res[[bi]]$inter_strength[nodes_int]
        inter_ei1_boot        [bi, ] <- boot_res[[bi]]$inter_ei1[nodes_int]
        inter_closeness_boot  [bi, ] <- boot_res[[bi]]$inter_closeness[nodes_int]
        inter_betweenness_boot[bi, ] <- boot_res[[bi]]$inter_betweenness[nodes_int]
      }
    }
  } # end bootstrap

  # ---------- CIs per layer ----------
  for (L in uniq_layers) {
    LB <- layer_boot[[L]]; if (is.null(LB)) next

    ci_list <- list(
      strength           = if (do_intra_general_boot && !is.null(LB$strength_boot))
        .calc_ci(LB$strength_boot, probs) else NULL,
      expected_influence = if (do_intra_general_boot && !is.null(LB$ei1_boot))
        .calc_ci(LB$ei1_boot, probs) else NULL,
      closeness          = if (do_intra_general_boot && !is.null(LB$closeness_boot))
        .calc_ci(LB$closeness_boot, probs) else NULL,
      betweenness        = if (do_intra_general_boot && !is.null(LB$betweenness_boot))
        .calc_ci(LB$betweenness_boot, probs) else NULL,

      edge_weights       = if (!is.null(LB$edge_boot_mat))
        .calc_ci(t(LB$edge_boot_mat), probs) else NULL,

      bridge_strength    = if (do_bridge_boot && !is.null(LB$bridge_strength_boot))
        .calc_ci(LB$bridge_strength_boot, probs) else NULL,
      bridge_betweenness = if (do_bridge_boot && !is.null(LB$bridge_betweenness_boot))
        .calc_ci(LB$bridge_betweenness_boot, probs) else NULL,
      bridge_closeness   = if (do_bridge_boot && !is.null(LB$bridge_closeness_boot))
        .calc_ci(LB$bridge_closeness_boot, probs) else NULL,
      bridge_ei1         = if (do_bridge_boot && !is.null(LB$bridge_ei1_boot))
        .calc_ci(LB$bridge_ei1_boot, probs) else NULL,
      bridge_ei2         = if (do_bridge_boot && !is.null(LB$bridge_ei2_boot))
        .calc_ci(LB$bridge_ei2_boot, probs) else NULL,

      bridge_strength_excluded    = if (do_excluded_boot && !is.null(LB$bridge_strength_excl_boot))
        .calc_ci(LB$bridge_strength_excl_boot, probs) else NULL,
      bridge_betweenness_excluded = if (do_excluded_boot && !is.null(LB$bridge_betweenness_excl_boot))
        .calc_ci(LB$bridge_betweenness_excl_boot, probs) else NULL,
      bridge_closeness_excluded   = if (do_excluded_boot && !is.null(LB$bridge_closeness_excl_boot))
        .calc_ci(LB$bridge_closeness_excl_boot, probs) else NULL,
      bridge_ei1_excluded         = if (do_excluded_boot && !is.null(LB$bridge_ei1_excl_boot))
        .calc_ci(LB$bridge_ei1_excl_boot, probs) else NULL,
      bridge_ei2_excluded         = if (do_excluded_boot && !is.null(LB$bridge_ei2_excl_boot))
        .calc_ci(LB$bridge_ei2_excl_boot, probs) else NULL
    )

    layer_fits[[L]]$settings$reps <- reps

    # --- LOADINGS: boot ---
    if (!is.null(LB$loadings_boot) && length(LB$loadings_boot) > 0) {
      layer_fits[[L]]$community_loadings$boot <- LB$loadings_boot
    }

    # --- NODE METRICS: boot ---
    layer_fits[[L]]$statistics$node$boot <- list(
      strength                    = LB$strength_boot,
      ei1                         = LB$ei1_boot,
      closeness                   = LB$closeness_boot,
      betweenness                 = LB$betweenness_boot,
      bridge_strength             = LB$bridge_strength_boot,
      bridge_betweenness          = LB$bridge_betweenness_boot,
      bridge_closeness            = LB$bridge_closeness_boot,
      bridge_ei1                  = LB$bridge_ei1_boot,
      bridge_ei2                  = LB$bridge_ei2_boot,
      bridge_strength_excluded    = LB$bridge_strength_excl_boot,
      bridge_betweenness_excluded = LB$bridge_betweenness_excl_boot,
      bridge_closeness_excluded   = LB$bridge_closeness_excl_boot,
      bridge_ei1_excluded         = LB$bridge_ei1_excl_boot,
      bridge_ei2_excluded         = LB$bridge_ei2_excl_boot
    )

    # --- NODE METRICS: CI ---
    node_ci <- ci_list[c(
      "strength","expected_influence","closeness","betweenness",
      "bridge_strength","bridge_betweenness","bridge_closeness",
      "bridge_ei1","bridge_ei2",
      "bridge_strength_excluded","bridge_betweenness_excluded",
      "bridge_closeness_excluded","bridge_ei1_excluded","bridge_ei2_excluded"
    )]
    layer_fits[[L]]$statistics$node$ci <- node_ci

    # --- EDGES: boot + CI ---
    layer_fits[[L]]$statistics$edge$boot <- LB$edge_boot_mat
    layer_fits[[L]]$statistics$edge$ci   <- ci_list$edge_weights

    # --- MEMBERSHIP BOOT ---
    layer_fits[[L]]$communities$boot_memberships <- LB$boot_memberships
  }

  # ---------- interlayer_fits ----------
  interlayer_fits <- list()
  for (key in names(pairs)) {
    A <- pairs[[key]]$A
    B <- pairs[[key]]$B

    cross_true <- wadj_signed[A, B, drop = FALSE]
    en <- .edge_names_cross(A, B)
    edges_true <- as.vector(cross_true); names(edges_true) <- en

    edges_true_df <- data.frame(
      edge   = en,
      weight = edges_true,
      row.names = NULL
    )

    eb <- interlayer_boot[[key]]
    ci_edges <- if (!is.null(eb)) .calc_ci(t(eb), probs) else NULL

    interlayer_fits[[key]] <- list(
      edges   = list(
        true = edges_true_df,   # cross-layer edges
        boot = eb,              # bootstrap matrix (edge x reps)
        ci   = ci_edges         # CIs edge
      )
    )
  }

  # ---------- interlayer-only node metrics TRUE + CI ----------
  inter_node_true <- data.frame(
    node = nodes_int,
    strength_interlayer    = inter_strength_true[nodes_int],
    ei1_interlayer         = inter_ei1_true[nodes_int],
    closeness_interlayer   = inter_closeness_true[nodes_int],
    betweenness_interlayer = inter_betweenness_true[nodes_int],
    row.names = NULL
  )
  inter_node_ci <- list(
    strength_interlayer    = if (do_interlayer_boot && !is.null(inter_strength_boot))
      .calc_ci(inter_strength_boot, probs) else NULL,
    ei1_interlayer         = if (do_interlayer_boot && !is.null(inter_ei1_boot))
      .calc_ci(inter_ei1_boot, probs) else NULL,
    closeness_interlayer   = if (do_interlayer_boot && !is.null(inter_closeness_boot))
      .calc_ci(inter_closeness_boot, probs) else NULL,
    betweenness_interlayer = if (do_interlayer_boot && !is.null(inter_betweenness_boot))
      .calc_ci(inter_betweenness_boot, probs) else NULL
  )
  inter_centrality_true <- data.frame(
    node        = inter_node_true$node,
    strength    = inter_node_true$strength_interlayer,
    ei1         = inter_node_true$ei1_interlayer,
    closeness   = inter_node_true$closeness_interlayer,
    betweenness = inter_node_true$betweenness_interlayer,
    row.names = NULL
  )
  inter_ci_results <- list(
    strength           = inter_node_ci$strength_interlayer,
    ei1                = inter_node_ci$ei1_interlayer,
    closeness          = inter_node_ci$closeness_interlayer,
    betweenness        = inter_node_ci$betweenness_interlayer
  )
  interlayer_centrality <- list(
    true       = inter_centrality_true,
    ci_results = inter_ci_results,
    boot       = list(
      strength    = inter_strength_boot,
      ei1         = inter_ei1_boot,
      closeness   = inter_closeness_boot,
      betweenness = inter_betweenness_boot
    )
  )
  interlayer <- c(
    list(centrality = interlayer_centrality),
    interlayer_fits
  )

  # ---------- global igraph ----------
  W_all <- wadj_signed[keep_nodes_graph_all, keep_nodes_graph_all, drop = FALSE]

  g_multi <- igraph::graph_from_adjacency_matrix(
    W_all,
    mode    = "undirected",
    weighted = TRUE,
    diag    = FALSE
  )

  if (igraph::ecount(g_multi) > 0) {
    igraph::E(g_multi)$abs_weight <- abs(igraph::E(g_multi)$weight)
    igraph::E(g_multi)$sign       <- ifelse(igraph::E(g_multi)$weight >= 0, 1L, -1L)
  }

  igraph::V(g_multi)$name  <- keep_nodes_graph_all
  igraph::V(g_multi)$layer <- layers[keep_nodes_graph_all]

  memb_all <- rep(NA_integer_, length(keep_nodes_graph_all))
  names(memb_all) <- keep_nodes_graph_all

  for (L in uniq_layers) {
    grpL <- layer_fits[[L]]$communities$groups
    if (!is.null(grpL) && length(grpL)) {
      grpL_int <- as.integer(grpL)
      names(grpL_int) <- names(grpL)
      common <- intersect(names(grpL_int), keep_nodes_graph_all)
      memb_all[common] <- grpL_int[common]
    }
  }
  igraph::V(g_multi)$membership <- memb_all

  if (igraph::ecount(g_multi) > 0) {
    ends_mat <- igraph::ends(g_multi, igraph::E(g_multi))
    L1 <- layers[ends_mat[, 1]]
    L2 <- layers[ends_mat[, 2]]

    igraph::E(g_multi)$type       <- ifelse(L1 == L2, "intra", "inter")
    igraph::E(g_multi)$layer_pair <- paste(pmin(L1, L2), pmax(L1, L2), sep = "_")
  }

  # ---------- output ----------
  out <- list(
    call = match.call(),

    settings = list(
      reps                         = reps,
      cluster_method               = cluster_method,
      covariates                   = covariates,
      exclude_from_cluster         = exclude_from_cluster,
      treat_singletons_as_excluded = treat_singletons_as_excluded,
      boot_what                    = boot_what,
      conf_level = conf_level
    ),

    data_info = list(
      mgm_type_level    = var_interpretation,
      binary_recode_map = binary_recode_map
    ),

    model = list(
      mgm   = mgm_model,
      nodes = all_nodes,
      n     = nrow(data),
      p     = ncol(data),
      data  = if (isTRUE(save_data)) data else NULL
    ),

    layers = list(
      assignment = layers,      # named vector node -> layer
      rules      = layer_rules,
      palette    = layer_colors
    ),

    layer_fits = layer_fits,

    interlayer = interlayer,

    graph = list(
      igraph            = g_multi,
      keep_nodes_graph  = keep_nodes_graph_all,
      keep_nodes_cluster= keep_nodes_cluster_all
    )
  )

  class(out) <- c("mixmashnet", "multimixMN_fit")
  return(out)
}
