#' Multilayer MGM with bootstrap, intra-/inter-layer metrics, and CIs
#'
#' @description
#' Estimates a \strong{multilayer} Mixed Graphical Model (MGM) using a mask that
#' enforces which cross-layer edges are allowed according to \code{layer_rules}.
#' Within each layer, the function computes community structure, non-parametric
#' row-bootstrap for node centralities and edges, bridge metrics
#' (including metrics for nodes treated as "excluded"), and (optionally)
#' community network scores. It also returns interlayer-only node metrics and
#' cross-layer edge summaries.
#'
#' @param data A numeric matrix or data frame (n × p) with variables in columns.
#' @param type Length-\code{p} vector of variable types as required by
#'   \code{mgm::mgm}.
#' @param level Length-\code{p} vector of variable levels as required by
#'   \code{mgm::mgm}.
#' @param layers A named vector (names = variable names) assigning each node to a
#'   layer (character or factor). Must cover all columns of \code{data}.
#' @param layer_rules A logical or numeric square matrix with row/column names
#'   equal to layer names. Values \code{TRUE} or \code{1} indicate that
#'   cross-layer edges are allowed between the corresponding layer pair.
#'   Intra-layer edges are always allowed; if the diagonal is \code{NA}, it is
#'   treated as allowed.
#' @param scale Logical; if \code{TRUE} (default) Gaussian variables
#'   (\code{type == "g"}) are z-standardized internally by \code{mgm()}. Use
#'   \code{scale = FALSE} if your data are already standardized.
#' @param reps Integer (>= 0). Number of bootstrap replications (row resampling
#'   with replacement). If \code{reps = 0}, no bootstrap is performed.
#' @param lambdaSel Method for lambda selection in \pkg{mgm}:
#'   \code{"CV"} or \code{"EBIC"}.
#' @param lambdaFolds Number of folds for CV (if \code{lambdaSel = "CV"}).
#' @param lambdaGam EBIC gamma parameter (if \code{lambdaSel = "EBIC"}).
#' @param alphaSeq Alpha parameters of the elastic net penalty (values in
#'   \eqn{(0, 1]}).
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
#' @param exclude_from_graph Character vector of node names. Nodes in this set
#'   are excluded from the global graph and from all node-level metrics.
#' @param exclude_from_cluster Character vector of node names. Nodes in this set
#'   are excluded from community detection in addition to \code{exclude_from_graph}.
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
#' @param compute_community_scores Logical; if \code{TRUE}, compute community
#'   network scores (EGAnet standardized \code{std.scores}) with bootstrap
#'   arrays and confidence intervals.
#' @param boot_what Character vector specifying which quantities to bootstrap.
#'   Valid options are:
#'   \code{"general_index"} (intra-layer centrality indices),
#'   \code{"interlayer_index"} (interlayer-only node metrics),
#'   \code{"bridge_index"} (bridge metrics for nodes in communities),
#'   \code{"excluded_index"} (bridge metrics for nodes treated as excluded),
#'   \code{"community"} (bootstrap community memberships),
#'   \code{"community_scores"} (bootstrap community scores, only if
#'   \code{compute_community_scores = TRUE}),
#'   and \code{"none"} (skip all node-level bootstrap: only edge-weight
#'   bootstrap is performed if \code{reps > 0}).
#' @param progress Logical; if \code{TRUE} (default), show a bootstrap progress
#'   bar when the \pkg{progressr} package is available.
#'
#' @return
#' An object of class \code{c("mixmashnet", "multimixMN_fit")}. The returned
#' list contains at least the following components:
#' \itemize{
#'   \item \code{call}: the matched function call.
#'   \item \code{settings}: list of main settings (e.g., \code{reps},
#'     \code{cluster_method}, \code{exclude_from_graph},
#'     \code{exclude_from_cluster}, \code{treat_singletons_as_excluded},
#'     \code{boot_what}).
#'   \item \code{model}: list with the fitted masked \code{mgm} object and basic
#'     information on the data (\code{nodes}, \code{n}, \code{p}).
#'   \item \code{layers}: list describing the multilayer structure
#'     (assignment of nodes to layers and the \code{layer_rules} matrix used).
#'   \item \code{layer_fits}: named list (one element per layer) with
#'     single-layer fits, including community structure, node-level statistics,
#'     edge-level statistics, bridge metrics, and (optionally) community scores
#'     with bootstrap information.
#'   \item \code{interlayer}: list collecting interlayer-only node metrics
#'     (strength, expected influence, closeness, betweenness, with or without
#'     bootstrap) and cross-layer edge summaries for each allowed pair of
#'     layers.
#'   \item \code{graph}: list containing a global \pkg{igraph} object built on
#'     the retained nodes (\code{keep_nodes_graph}), with vertex attributes
#'     such as \code{name}, \code{layer}, \code{membership}, and edge attributes
#'     such as \code{weight}, \code{abs_weight}, \code{sign},
#'     \code{type} (intra vs inter) and \code{layer_pair}.
#' }
#'
#' @details
#' This function does \strong{not} call \code{future::plan()}. To enable
#' parallel bootstrap, set a plan (e.g. \code{future::plan(multisession)}) before
#' calling \code{multimixMN()}. If \code{"none"} is the only element of
#' \code{boot_what} and \code{reps > 0}, node-level metrics are not
#' bootstrapped, but intra- and inter-layer edge-weight bootstrap and the
#' corresponding confidence intervals are still computed.
#'
#' @section Dependencies:
#' Relies on internal helpers \code{mgm_masked()}, \code{mixMN_from_wadj()},
#' \code{bridge_metrics()}, and \code{bridge_metrics_excluded()}, which must be
#' available in the same package namespace.
#'
#' @importFrom stats setNames quantile
#' @importFrom utils combn capture.output
#' @importFrom igraph graph_from_adjacency_matrix simplify ecount E V distances betweenness vcount
#' @importFrom igraph cluster_louvain cluster_fast_greedy cluster_infomap cluster_walktrap cluster_edge_betweenness
#' @importFrom future.apply future_lapply
#' @importFrom qgraph centrality
#' @importFrom progressr with_progress progressor
#' @importFrom EGAnet net.scores
#' @export
multimixMN <- function(
    data, type, level,
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
    exclude_from_graph = NULL,
    exclude_from_cluster = NULL,
    seed_model = NULL, seed_boot = NULL,
    treat_singletons_as_excluded = FALSE,
    cluster_method = c("louvain","fast_greedy","infomap","walktrap","edge_betweenness"),
    compute_community_scores = FALSE,
    boot_what = c("general_index", "interlayer_index", "bridge_index",
                  "excluded_index", "community", "community_scores", "none"),
    progress = TRUE
) {
  lambdaSel <- match.arg(lambdaSel)
  cluster_method <- match.arg(cluster_method)
  if (!is.null(seed_model)) set.seed(seed_model)

  # --- parse 'boot_what' argument ---
  boot_what <- match.arg(
    boot_what,
    choices = c("general_index", "interlayer_index",
                "bridge_index", "excluded_index",
                "community", "community_scores", "none"),
    several.ok = TRUE
  )

  if ("none" %in% boot_what && length(boot_what) == 1L) {
    do_intra_general_boot <- FALSE  # intra-layer strength/EI1/close/betw
    do_interlayer_boot    <- FALSE  # interlayer-only strength/EI1/close/betw
    do_bridge_boot        <- FALSE
    do_excluded_boot      <- FALSE
    do_community_boot     <- FALSE
    do_comm_scores_boot   <- FALSE
  } else {
    do_intra_general_boot <- "general_index"    %in% boot_what
    do_interlayer_boot    <- "interlayer_index" %in% boot_what
    do_bridge_boot        <- "bridge_index"     %in% boot_what
    do_excluded_boot      <- "excluded_index"   %in% boot_what
    do_community_boot     <- "community"        %in% boot_what
    do_comm_scores_boot   <- isTRUE(compute_community_scores) &&
      ("community_scores" %in% boot_what)
  }

  all_nodes <- colnames(data); p <- ncol(data)

  subject_ids <- rownames(data)
  if (is.null(subject_ids)) {
    subject_ids <- sprintf("id_%d", seq_len(nrow(data)))
    rownames(data) <- subject_ids
  }

  # --- Basic multilayer check
  if (length(unique(layers)) < 2) {
    stop(
      paste0(
        "multimixMN is designed for MULTILAYER networks (>= 2 layers). ",
        "You provided only one layer. For single-layer networks, please use the ",
        "`mixMN()` function instead."
      ),
      call. = FALSE
    )
  }

  # --- layer_rules validation and normalization
  stopifnot(is.matrix(layer_rules))
  layer_order <- unique(layers)

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

  .quiet_net_scores <- function(...) {
    withCallingHandlers({
      out <- NULL
      invisible(capture.output(out <- EGAnet::net.scores(...)))
      out
    }, message = function(m) {
      txt <- conditionMessage(m)
      if (grepl("default 'loading.method'.*revised", txt, ignore.case = TRUE)) invokeRestart("muffleMessage")
    }, warning = function(w) {
      txt <- conditionMessage(w)
      if (grepl("default 'loading.method'.*revised", txt, ignore.case = TRUE)) invokeRestart("muffleWarning")
    })
  }

  .make_distance_graph <- function(W_signed) {
    W <- abs(W_signed); diag(W) <- 0; W[is.na(W) | W < tiny] <- 0
    g <- igraph::graph_from_adjacency_matrix(W, mode="undirected", weighted=TRUE, diag=FALSE)
    if (igraph::ecount(g) > 0) igraph::E(g)$dist <- 1/igraph::E(g)$weight
    g
  }
  .harmonic_closeness <- function(g) {
    if (igraph::ecount(g) == 0) return(stats::setNames(rep(0, igraph::vcount(g)), igraph::V(g)$name))
    D <- igraph::distances(g, weights = igraph::E(g)$dist); diag(D) <- NA
    cl <- rowSums(1/D, na.rm=TRUE) / (nrow(D)-1); names(cl) <- igraph::V(g)$name; cl
  }
  .calc_ci <- function(mat) {
    if (is.null(mat)) return(NULL)
    ci <- apply(mat, 2, function(x)
      if (all(is.na(x))) c(`2.5%`=NA_real_, `97.5%`=NA_real_)
      else stats::quantile(x, probs=c(0.025,0.975), na.rm=TRUE))
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
  mask_list <- vector("list", p); names(mask_list) <- all_nodes
  for (i in seq_along(all_nodes)) {
    v <- all_nodes[i]; my_layer <- layers[v]
    allowed_layers <- names(which(layer_rules[my_layer, ]))
    allowed_layers <- unique(c(my_layer, allowed_layers)) # intra always allowed
    allowed_nodes <- names(layers)[layers %in% allowed_layers]
    mask_list[[i]] <- match(allowed_nodes, all_nodes)
  }

  # --- MGM (masked) fit
  mgm_model <- mgm_masked(
    data = as.matrix(data), type = type, level = level,
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
  keep_nodes_graph_all   <- setdiff(all_nodes, exclude_from_graph)
  keep_nodes_cluster_all <- setdiff(all_nodes, unique(c(exclude_from_graph, exclude_from_cluster)))
  uniq_layers <- unique(layers)

  layer_nodes_graph   <- stats::setNames(lapply(uniq_layers, function(L)
    intersect(names(layers)[layers==L], keep_nodes_graph_all)), uniq_layers)

  layer_nodes_cluster <- stats::setNames(lapply(uniq_layers, function(L)
    intersect(names(layers)[layers==L], keep_nodes_cluster_all)), uniq_layers)

  # --- Per-layer fits from W (true graph, no bootstrap yet)
  layer_fits <- list()
  for (L in uniq_layers) {
    nL <- layer_nodes_graph[[L]]
    sub_w <- wadj_signed[nL, nL, drop=FALSE]
    fitL <- mixMN_from_wadj(
      wadj_signed = sub_w, nodes = nL,
      exclude_from_graph   = NULL,
      exclude_from_cluster = intersect(exclude_from_cluster, nL),
      cluster_method = cluster_method,
      reps          = 0,
      seed_boot     = NULL,
      treat_singletons_as_excluded = treat_singletons_as_excluded,
      boot_what     = boot_what
    )

    membL <- fitL$communities$groups
    if (is.factor(membL)) {
      membL <- stats::setNames(as.integer(membL), names(membL))
    }
    fitL$communities$groups <- membL

    ## --- COMMUNITY SCORES (TRUE) ---
    if (isTRUE(compute_community_scores)) {
      membL <- fitL$communities$groups
      nodes_comm <- names(membL)[!is.na(membL)]
      nodes_comm <- intersect(nodes_comm, nL)

      if (length(nodes_comm) > 0) {
        A_comm <- sub_w[nodes_comm, nodes_comm, drop = FALSE]

        wc_fac <- membL[nodes_comm]
        wc_levels <- sort(unique(as.integer(wc_fac)))
        wc_map <- stats::setNames(seq_along(wc_levels), wc_levels)
        wc_int <- unname(wc_map[as.integer(wc_fac)])

        dat_comm <- as.matrix(data[, nodes_comm, drop = FALSE])
        subject_ids <- rownames(dat_comm)
        if (is.null(subject_ids)) subject_ids <- sprintf("id_%d", seq_len(nrow(dat_comm)))
        rownames(dat_comm) <- subject_ids

        ns_obj <- tryCatch(
          .quiet_net_scores(
            data = dat_comm,
            A    = A_comm,
            wc   = wc_int,
            loading.method = "revised",
            structure = "simple",
            rotation  = NULL,
            scores    = "components"
          ),
          error = function(e) NULL
        )

        if (!is.null(ns_obj) && !is.null(ns_obj$scores$std.scores)) {
          cs_true <- ns_obj$scores$std.scores
        } else {
          cs_true <- matrix(
            NA_real_,
            nrow = nrow(dat_comm),
            ncol = length(unique(wc_int))
          )
        }
        colnames(cs_true) <- paste0("NS_", sort(unique(as.integer(wc_fac))))
        rownames(cs_true) <- subject_ids

        fitL$community_scores$obj <- ns_obj
        fitL$community_scores$df  <- data.frame(
          id = subject_ids,
          as.data.frame(cs_true),
          check.names = FALSE
        )
        # info per il bootstrap
        fitL$community_scores$.nodes_comm <- nodes_comm
        fitL$community_scores$.wc_int     <- wc_int
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
      # indici generali intra-layer
      strength_boot     = if (do_intra_general_boot)
        matrix(NA_real_, reps, length(nodes_g), dimnames = list(NULL, nodes_g)) else NULL,
      ei1_boot          = if (do_intra_general_boot)
        matrix(NA_real_, reps, length(nodes_g), dimnames = list(NULL, nodes_g)) else NULL,
      closeness_boot    = if (do_intra_general_boot)
        matrix(NA_real_, reps, length(nodes_g), dimnames = list(NULL, nodes_g)) else NULL,
      betweenness_boot  = if (do_intra_general_boot)
        matrix(NA_real_, reps, length(nodes_g), dimnames = list(NULL, nodes_g)) else NULL,

      # edges intra-layer: li teniamo sempre
      edge_boot_mat     = {
        en <- .edge_names_lt(nodes_g)
        if (length(en)) matrix(NA_real_, length(en), reps,
                               dimnames = list(en, NULL)) else NULL
      },

      # membership bootstrap solo se richiesto
      boot_memberships  = if (do_community_boot) vector("list", reps) else vector("list", 0),

      # bridge intra-layer
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

      # bridge per nodi esclusi
      bridge_strength_excl_boot     = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL,
      bridge_ei1_excl_boot          = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL,
      bridge_ei2_excl_boot          = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL,
      bridge_closeness_excl_boot    = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL,
      bridge_betweenness_excl_boot  = if (do_excluded_boot && length(nodes_ex) > 0)
        matrix(NA_real_, reps, length(nodes_ex), dimnames = list(NULL, nodes_ex)) else NULL,

      # community scores bootstrap per layer (se richiesti)
      community_scores_boot = {
        if (reps > 0 &&
            isTRUE(do_comm_scores_boot) &&
            isTRUE(compute_community_scores) &&
            !is.null(layer_fits[[L]]$community_scores$.nodes_comm)) {

          subj <- rownames(data)
          if (is.null(subj)) subj <- sprintf("id_%d", seq_len(nrow(data)))
          K <- length(unique(layer_fits[[L]]$community_scores$.wc_int))

          array(
            NA_real_,
            dim = c(reps, length(subj), K),
            dimnames = list(
              rep  = paste0("b", seq_len(reps)),
              id   = subj,
              comm = paste0("NS_", sort(unique(layer_fits[[L]]$community_scores$.wc_int)))
            )
          )
        } else NULL
      }
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
    idx <- sample(seq_len(nrow(data)), replace=TRUE)
    Xb  <- as.matrix(data[idx, , drop=FALSE])
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

      # --- indici generali intra-layer ---
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

      ## --- ADD: COMMUNITY SCORES on bootstrap (per layer) ---
      cs_boot <- NULL
      if (isTRUE(do_comm_scores_boot)) {
        nodes_comm  <- layer_fits[[L]]$community_scores$.nodes_comm
        wc_comm_int <- layer_fits[[L]]$community_scores$.wc_int
        if (!is.null(nodes_comm) && length(nodes_comm) > 0) {
          A_comm_boot <- Wg[nodes_comm, nodes_comm, drop = FALSE]
          dat_comm_full <- as.matrix(data[, nodes_comm, drop = FALSE])
          rownames(dat_comm_full) <- subject_ids
          ns_obj <- tryCatch(
            .quiet_net_scores(
              data = dat_comm_full,
              A    = A_comm_boot,
              wc   = wc_comm_int,
              loading.method = "revised",
              structure = "simple",
              rotation  = NULL,
              scores    = "components"
            ),
            error = function(e) NULL
          )
          if (!is.null(ns_obj) && !is.null(ns_obj$scores$std.scores)) {
            cs_boot <- ns_obj$scores$std.scores
            colnames(cs_boot) <- paste0("NS_", sort(unique(wc_comm_int)))
            rownames(cs_boot) <- subject_ids
          } else {
            cs_boot <- matrix(NA_real_, nrow = length(subject_ids), ncol = length(unique(wc_comm_int)),
                              dimnames = list(subject_ids, paste0("NS_", sort(unique(wc_comm_int)))))
          }
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
        community_scores_boot = cs_boot
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

    # Cross-layer edges per allowed pairs (questi li puoi lasciare sempre)
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

    # Collect bootstrap results
    for (bi in seq_len(reps)) {
      if (is.null(boot_res[[bi]])) next

      for (L in uniq_layers) {
        if (is.null(layer_boot[[L]])) next
        pl <- boot_res[[bi]]$per_layer[[L]]; if (is.null(pl)) next

        if (!is.null(layer_boot[[L]]$community_scores_boot) && !is.null(pl$community_scores_boot)) {
          # array [rep, id, comm]
          ids   <- dimnames(layer_boot[[L]]$community_scores_boot)$id
          comms <- dimnames(layer_boot[[L]]$community_scores_boot)$comm
          mat <- pl$community_scores_boot
          # allinea
          miss_rows <- setdiff(ids, rownames(mat))
          if (length(miss_rows)) {
            mat <- rbind(mat, matrix(NA_real_, nrow = length(miss_rows), ncol = ncol(mat),
                                     dimnames = list(miss_rows, colnames(mat))))
          }
          miss_cols <- setdiff(comms, colnames(mat))
          if (length(miss_cols)) {
            mat <- cbind(mat, matrix(NA_real_, nrow = nrow(mat), ncol = length(miss_cols),
                                     dimnames = list(rownames(mat), miss_cols)))
          }
          mat <- mat[ids, comms, drop = FALSE]
          layer_boot[[L]]$community_scores_boot[bi, , ] <- mat
        }

        nodes_g <- layer_nodes_graph[[L]]
        nodes_c <- names(layer_fits[[L]]$communities$groups)        # fixed
        if (is.null(nodes_c)) nodes_c <- character(0)
        nodes_ex <- layer_nodes_excluded[[L]]           # fixed

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
        .calc_ci(LB$strength_boot) else NULL,
      expected_influence = if (do_intra_general_boot && !is.null(LB$ei1_boot))
        .calc_ci(LB$ei1_boot) else NULL,
      closeness          = if (do_intra_general_boot && !is.null(LB$closeness_boot))
        .calc_ci(LB$closeness_boot) else NULL,
      betweenness        = if (do_intra_general_boot && !is.null(LB$betweenness_boot))
        .calc_ci(LB$betweenness_boot) else NULL,

      edge_weights       = if (!is.null(LB$edge_boot_mat))
        .calc_ci(t(LB$edge_boot_mat)) else NULL,

      bridge_strength    = if (do_bridge_boot && !is.null(LB$bridge_strength_boot))
        .calc_ci(LB$bridge_strength_boot) else NULL,
      bridge_betweenness = if (do_bridge_boot && !is.null(LB$bridge_betweenness_boot))
        .calc_ci(LB$bridge_betweenness_boot) else NULL,
      bridge_closeness   = if (do_bridge_boot && !is.null(LB$bridge_closeness_boot))
        .calc_ci(LB$bridge_closeness_boot) else NULL,
      bridge_ei1         = if (do_bridge_boot && !is.null(LB$bridge_ei1_boot))
        .calc_ci(LB$bridge_ei1_boot) else NULL,
      bridge_ei2         = if (do_bridge_boot && !is.null(LB$bridge_ei2_boot))
        .calc_ci(LB$bridge_ei2_boot) else NULL,

      bridge_strength_excluded    = if (do_excluded_boot && !is.null(LB$bridge_strength_excl_boot))
        .calc_ci(LB$bridge_strength_excl_boot) else NULL,
      bridge_betweenness_excluded = if (do_excluded_boot && !is.null(LB$bridge_betweenness_excl_boot))
        .calc_ci(LB$bridge_betweenness_excl_boot) else NULL,
      bridge_closeness_excluded   = if (do_excluded_boot && !is.null(LB$bridge_closeness_excl_boot))
        .calc_ci(LB$bridge_closeness_excl_boot) else NULL,
      bridge_ei1_excluded         = if (do_excluded_boot && !is.null(LB$bridge_ei1_excl_boot))
        .calc_ci(LB$bridge_ei1_excl_boot) else NULL,
      bridge_ei2_excluded         = if (do_excluded_boot && !is.null(LB$bridge_ei2_excl_boot))
        .calc_ci(LB$bridge_ei2_excl_boot) else NULL,

      community_scores = {
        CSB <- LB$community_scores_boot
        if (isTRUE(do_comm_scores_boot) && !is.null(CSB)) {
          qfun <- function(x, p) if (all(is.na(x))) NA_real_ else
            stats::quantile(x, probs = p, na.rm = TRUE)
          low  <- apply(CSB, c(2, 3), qfun, p = 0.025)
          up   <- apply(CSB, c(2, 3), qfun, p = 0.975)
          list(lower = low, upper = up)
        } else NULL
      }
    )

    # attach standard stuff
    layer_fits[[L]]$settings$reps <- reps

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

    # --- COMMUNITY SCORES: boot + CI ---
    if (!is.null(LB$community_scores_boot) && !is.null(layer_fits[[L]]$community_scores$df)) {
      reps_names <- dimnames(LB$community_scores_boot)$rep
      id_names   <- dimnames(LB$community_scores_boot)$id
      comm_names <- dimnames(LB$community_scores_boot)$comm

      cs_list <- lapply(seq_along(reps_names), function(i) {
        mat <- LB$community_scores_boot[i, , , drop = FALSE][1, , ]
        mat <- as.matrix(mat)
        dimnames(mat) <- list(id_names, comm_names)
        as.data.frame(mat, check.names = FALSE)
      })
      names(cs_list) <- reps_names

      layer_fits[[L]]$community_scores$boot <- cs_list
      layer_fits[[L]]$community_scores$ci   <- ci_list$community_scores
    }

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
    ci_edges <- if (!is.null(eb)) .calc_ci(t(eb)) else NULL

    interlayer_fits[[key]] <- list(
      edges   = list(
        true = edges_true_df,   # pesi “osservati” cross-layer
        boot = eb,              # matrice bootstrap (edge × reps)
        ci   = ci_edges         # CIs per edge
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
      .calc_ci(inter_strength_boot) else NULL,
    ei1_interlayer         = if (do_interlayer_boot && !is.null(inter_ei1_boot))
      .calc_ci(inter_ei1_boot) else NULL,
    closeness_interlayer   = if (do_interlayer_boot && !is.null(inter_closeness_boot))
      .calc_ci(inter_closeness_boot) else NULL,
    betweenness_interlayer = if (do_interlayer_boot && !is.null(inter_betweenness_boot))
      .calc_ci(inter_betweenness_boot) else NULL
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
      exclude_from_graph           = exclude_from_graph,
      exclude_from_cluster         = exclude_from_cluster,
      treat_singletons_as_excluded = treat_singletons_as_excluded,
      boot_what                    = boot_what
    ),

    model = list(
      mgm   = mgm_model,
      nodes = all_nodes,
      n     = nrow(data),
      p     = ncol(data)
    ),

    layers = list(
      assignment = layers,      # named vector node -> layer
      rules      = layer_rules
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
