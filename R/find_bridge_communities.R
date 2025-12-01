#' Bridge profiles of a node across communities
#'
#' For a given node, compute how strongly it "bridges" to \emph{other}
#' communities using five profiles:
#' \itemize{
#'   \item \strong{strength}: sum of \eqn{|w|} to nodes in other communities
#'   \item \strong{ei1}: signed sum of \eqn{w} to nodes in other communities
#'   \item \strong{ei2}: signed influence \eqn{A + A^2} to nodes in other communities
#'   \item \strong{closeness}: \eqn{1 /} mean shortest-path distance to nodes in other communities;
#'     paths computed on \eqn{|W|} with edge length \eqn{1/|w|}
#'   \item \strong{betweenness}: number of inter-community \emph{shortest paths}
#'     (same graph as for closeness, but \emph{unweighted} in the path search)
#'     that pass through the node as an \emph{intermediate} (path multiplicity is counted)
#' }
#'
#' Notes:
#' \itemize{
#'   \item The signed adjacency \eqn{W} is rebuilt \emph{exactly} from \code{fit$statistics$edge$true}
#'         (symmetric, diagonal set to 0).
#'   \item "Other communities" means targets whose community differs from the
#'         node's community; if the node is unassigned (NA), targets are all assigned nodes.
#'   \item Betweenness counts \emph{per path} (multiplicity). Endpoints are always
#'         \emph{assigned} nodes in \emph{different} communities. The focal node may be
#'         assigned or unassigned; it is counted only as intermediate (never as endpoint).
#' }
#'
#' @param fit An object of class \code{mixMN_fit} containing at least:
#'   \itemize{
#'     \item \code{$graph$keep_nodes_graph}: character vector of node names (graph vertex set);
#'     \item \code{$statistics$edge$true}: data.frame with columns \code{edge} ("a--b") and \code{weight} (signed).
#'     \item (optional) \code{$communities$groups}: named factor with communities for a subset of nodes.
#'   }
#' @param node Character scalar: node of interest; must belong to
#'   \code{fit$graph$keep_nodes_graph}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{strength}}{list with \code{overall} (numeric) and \code{by_comm} (tibble: \code{community}, \code{sum_abs_w}).}
#'   \item{\code{ei1}}{list with \code{overall} and \code{by_comm} (tibble: \code{community}, \code{sum_signed_w}).}
#'   \item{\code{ei2}}{list with \code{overall} and \code{by_comm} (tibble: \code{community}, \code{sum_signed_w2}).}
#'   \item{\code{closeness}}{list with \code{overall} and \code{by_comm} (tibble: \code{community}, \code{inv_mean_dist}).}
#'   \item{\code{betweenness}}{list with \code{overall} and \code{by_pair} (tibble: \code{Ci}, \code{Cj}, \code{hits}).}
#' }
#'
#' @importFrom igraph graph_from_adjacency_matrix delete_edges distances is_directed degree get.all.shortest.paths
#' @importFrom dplyr bind_rows arrange desc
#' @importFrom tibble tibble
#' @importFrom tidyr separate
#' @export
find_bridge_communities <- function(fit, node) {
  # ---- guardrails ----
  stopifnot(inherits(fit, "mixMN_fit"))
  if (is.null(fit$graph$keep_nodes_graph) || is.null(fit$statistics$edge$true)) {
    stop("`fit` must contain `graph$keep_nodes_graph` and `statistics$edge$true`.")
  }

  nodes <- fit$graph$keep_nodes_graph
  if (!length(nodes)) stop("`graph$keep_nodes_graph` is empty.")
  if (!(node %in% nodes)) stop("Requested 'node' is not in `graph$keep_nodes_graph`.")

  Edf <- fit$statistics$edge$true
  if (!all(c("edge", "weight") %in% names(Edf))) {
    stop("`statistics$edge$true` must have columns 'edge' and 'weight'.")
  }

  # ---- rebuild signed W from edges_true ----
  split_edge <- function(s, sep = "--") strsplit(s, sep, fixed = TRUE)[[1]]

  W <- matrix(0, nrow = length(nodes), ncol = length(nodes),
              dimnames = list(nodes, nodes))
  for (k in seq_len(nrow(Edf))) {
    ab <- split_edge(Edf$edge[k]); if (length(ab) != 2L) next
    a <- ab[1]; b <- ab[2]
    if (!(a %in% nodes && b %in% nodes)) next
    w <- Edf$weight[k]
    W[a, b] <- w
    W[b, a] <- w
  }
  diag(W) <- 0
  W <- (W + t(W)) / 2

  # ---- membership aligned on `nodes` ----
  if (!is.null(fit$communities) && !is.null(fit$communities$groups)) {
    grp <- fit$communities$groups  # factor, named
    if (is.null(names(grp))) stop("`communities$groups` must be a named factor.")
    # inizialmente tutti NA
    comm_full <- setNames(rep(NA_integer_, length(nodes)), nodes)
    in_idx <- names(grp)[names(grp) %in% nodes]
    if (length(in_idx)) comm_full[in_idx] <- as.integer(grp[in_idx])
  } else {
    comm_full <- setNames(rep(NA_integer_, length(nodes)), nodes)
  }

  assigned_idx <- which(!is.na(comm_full))
  if (!length(assigned_idx)) {
    empty <- tibble::tibble()
    return(list(
      strength    = list(overall = NA_real_, by_comm = empty),
      ei1         = list(overall = NA_real_, by_comm = empty),
      ei2         = list(overall = NA_real_, by_comm = empty),
      closeness   = list(overall = NA_real_, by_comm = empty),
      betweenness = list(overall = 0,       by_comm = empty, by_pair = empty)
    ))
  }

  v_comm <- comm_full[[node]]

  other_comm_targets <- function() {
    if (is.na(v_comm)) setdiff(assigned_idx, match(node, nodes))
    else setdiff(assigned_idx[comm_full[assigned_idx] != v_comm], match(node, nodes))
  }
  targets <- other_comm_targets()

  # ================== PROFILES WITHOUT PATHS (on signed W) ==================
  # strength
  strength_list <- if (!length(targets)) {
    list(overall = NA_real_, by_comm = tibble::tibble())
  } else {
    tgt_comm <- comm_full[targets]
    idx_by_comm <- split(targets, tgt_comm)
    by_comm <- lapply(idx_by_comm, function(idx) {
      tibble::tibble(
        community = as.integer(unique(comm_full[idx])),
        sum_abs_w = sum(abs(W[node, idx]))
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::arrange(dplyr::desc(sum_abs_w))
    list(overall = sum(abs(W[node, targets])), by_comm = by_comm)
  }

  # ei1
  ei1_list <- if (!length(targets)) {
    list(overall = NA_real_, by_comm = tibble::tibble())
  } else {
    tgt_comm <- comm_full[targets]
    idx_by_comm <- split(targets, tgt_comm)
    by_comm <- lapply(idx_by_comm, function(idx) {
      tibble::tibble(
        community    = as.integer(unique(comm_full[idx])),
        sum_signed_w = sum(W[node, idx])
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::arrange(dplyr::desc(sum_signed_w))
    list(overall = sum(W[node, targets]), by_comm = by_comm)
  }

  # ei2 = (A + A^2) signed influence
  ei2_list <- if (!length(targets)) {
    list(overall = NA_real_, by_comm = tibble::tibble())
  } else {
    A  <- W; diag(A) <- 0
    A2 <- A %*% A
    infl2 <- A[node, ] + A2[node, ]
    tgt_comm <- comm_full[targets]
    idx_by_comm <- split(targets, tgt_comm)
    by_comm <- lapply(idx_by_comm, function(idx) {
      tibble::tibble(
        community     = as.integer(unique(comm_full[idx])),
        sum_signed_w2 = sum(infl2[idx])
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::arrange(dplyr::desc(sum_signed_w2))
    list(overall = sum(infl2[targets]), by_comm = by_comm)
  }

  # ================== PATHS ON |W| (edge length = 1/|w|) ==================
  EPS  <- 1e-10
  Wabs <- abs(W)
  g_absw <- igraph::graph_from_adjacency_matrix(Wabs, mode = "undirected", weighted = TRUE, diag = FALSE)
  if (is.null(igraph::E(g_absw)$weight)) igraph::E(g_absw)$weight <- 1
  igraph::E(g_absw)$weight <- abs(igraph::E(g_absw)$weight)

  # remove edges with NA or |w| <= EPS
  g_inv <- igraph::delete_edges(g_absw, which(is.na(igraph::E(g_absw)$weight) | igraph::E(g_absw)$weight <= EPS))

  # invert: distance = 1/|w|
  if (igraph::ecount(g_inv) > 0) igraph::E(g_inv)$weight <- 1 / igraph::E(g_inv)$weight

  # remove negative inverted edges (coherent with original-style)
  g_pos <- if (igraph::ecount(g_inv) > 0) {
    igraph::delete_edges(g_inv, which(igraph::E(g_inv)$weight < 0))
  } else g_inv

  # ---- CLOSENESS (weighted on 1/|w|) ----
  closeness_list <- if (!length(targets) || igraph::ecount(g_pos) == 0) {
    list(overall = NA_real_, by_comm = tibble::tibble())
  } else {
    v_id <- match(node, igraph::V(g_pos)$name)
    d_all <- suppressWarnings(igraph::distances(
      g_pos, v = v_id, to = igraph::V(g_pos)$name[targets],
      weights = igraph::E(g_pos)$weight, mode = "all"
    ))
    d_all <- as.numeric(d_all)
    overall <- if (all(is.infinite(d_all))) NA_real_ else 1 / mean(d_all[is.finite(d_all)])

    tgt_comm <- comm_full[targets]
    by_comm <- lapply(split(seq_along(targets), tgt_comm), function(ix) {
      di <- d_all[ix]
      tibble::tibble(
        community     = as.integer(unique(tgt_comm[ix])),
        inv_mean_dist = if (any(is.finite(di))) 1 / mean(di[is.finite(di)]) else NA_real_
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::arrange(dplyr::desc(inv_mean_dist))
    list(overall = overall, by_comm = by_comm)
  }

  # sanity checks on vertex sets (avoid name mismatches)
  stopifnot(identical(sort(nodes), sort(igraph::V(g_pos)$name)))
  stopifnot(all(names(comm_full) == nodes))

  # ---- BETWEENNESS (UNWEIGHTED on g_pos; counts path multiplicity) ----
  assigned_nodes       <- names(comm_full[!is.na(comm_full)])
  assigned_communities <- comm_full[!is.na(comm_full)]

  v_id <- match(node, igraph::V(g_pos)$name)
  if (is.na(v_id) || igraph::degree(g_pos, v = v_id) == 0) {
    betweenness_list <- list(
      overall = 0,
      by_comm = tibble::tibble(),
      by_pair = tibble::tibble(Ci = integer(), Cj = integer(), hits = integer())
    )
  } else {
    total_hits <- 0L
    by_pair_tab <- list()

    mid_ids <- function(path_ids) {
      ids <- as.vector(path_ids)
      if (length(ids) <= 2L) integer(0) else ids[-c(1L, length(ids))]
    }

    if (igraph::ecount(g_pos) > 0 && length(assigned_nodes)) {
      for (src in assigned_nodes) {
        src_comm <- assigned_communities[[src]]
        to_nodes <- assigned_nodes[assigned_communities != src_comm]
        if (!length(to_nodes)) next

        for (t in to_nodes) {
          sp <- suppressWarnings(igraph::get.all.shortest.paths(g_pos, from = src, to = t, mode = "all"))
          if (!length(sp$res)) next

          # count multiplicity: number of shortest paths that include v_id as intermediate
          k_hit <- 0L
          for (p in sp$res) if (any(mid_ids(p) == v_id)) k_hit <- k_hit + 1L

          if (k_hit > 0L) {
            total_hits <- total_hits + k_hit
            key <- paste(sort(c(as.integer(src_comm), as.integer(assigned_communities[[t]]))), collapse = "-")
            by_pair_tab[[key]] <- (if (is.null(by_pair_tab[[key]])) 0L else by_pair_tab[[key]]) + k_hit
          }
        }
      }
    }

    # undirected correction
    if (!igraph::is_directed(g_pos)) {
      total_hits <- total_hits / 2
      if (length(by_pair_tab)) by_pair_tab <- lapply(by_pair_tab, function(x) x / 2)
    }

    betweenness_by_pair <- if (length(by_pair_tab)) {
      tibble::tibble(
        pair = names(by_pair_tab),
        hits = as.numeric(unlist(by_pair_tab))
      ) |>
        tidyr::separate(pair, into = c("Ci", "Cj"), sep = "-", convert = TRUE) |>
        dplyr::arrange(dplyr::desc(hits))
    } else {
      tibble::tibble(Ci = integer(), Cj = integer(), hits = integer())
    }

    betweenness_list <- list(
      overall = as.numeric(total_hits),
      by_pair = betweenness_by_pair
    )
  }

  # ---- output ----
  list(
    strength    = strength_list,   # |W|
    ei1         = ei1_list,        # W (signed)
    ei2         = ei2_list,        # A + A^2 (signed)
    closeness   = closeness_list,  # weighted distances on 1/|w|
    betweenness = betweenness_list # unweighted shortest paths on same edge set
  )
}
