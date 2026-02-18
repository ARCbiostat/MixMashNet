#' Bridge profiles of a node across communities
#'
#' Identifies which communities contribute most to the bridge role of a
#' given node, by decomposing its bridge connectivity into community-specific
#' contributions, excluding its own community when assigned. The function is
#' designed as an interpretative companion to \code{bridge_metrics()} and
#' \code{bridge_metrics_excluded()}, providing the components underlying
#' the corresponding overall bridge indices.
#'
#' Bridge connectivity is summarized using five complementary profiles: bridge
#' strength, bridge EI1, bridge EI2, bridge closeness, and bridge betweenness.
#'
#' Notes:
#' \itemize{
#'   \item Bridge profiles are computed using only connections from the focal node to
#'         nodes in communities different from its own. If the focal node is not
#'         assigned to any community, i.e. excluded, connections to all assigned nodes in
#'         communities are considered.
#'   \item Bridge betweenness is computed by counting all shortest paths between
#'         pairs of nodes in different communities that pass through the focal
#'         node as an intermediate vertex. When multiple shortest paths exist,
#'         each path is counted separately.
#' }
#'
#' @param fit An object of class \code{mixMN_fit}.
#' @param node Character scalar: node of interest; must belong to
#'   \code{fit$graph$keep_nodes_graph}.
#'
#' @return An object of class \code{"bridge_profiles"} (a named list) with the
#'   following components:
#' \describe{
#'   \item{\code{bridge_strength}}{Bridge strength. List with \code{overall}, the total
#'     value across all other communities, and \code{by_comm}, a tibble with
#'     community-specific contributions (\code{community}, \code{sum_abs_w}).}
#'   \item{\code{bridge_ei1}}{Bridge expected influence (order 1). List with
#'     \code{overall} and \code{by_comm} (\code{community}, \code{sum_signed_w}).}
#'   \item{\code{bridge_ei2}}{Bridge expected influence (order 2). List with
#'     \code{overall} and \code{by_comm} (\code{community}, \code{sum_signed_w2}).}
#'   \item{\code{bridge_closeness}}{Bridge closeness. List with \code{overall} and
#'     \code{by_comm} (\code{community}, \code{inv_mean_dist}).}
#'   \item{\code{bridge_betweenness}}{Bridge betweenness. List with \code{overall} and
#'     \code{by_pair}, a tibble with contributions by community pair
#'     (\code{Ci}, \code{Cj}, \code{hits}).}
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
    comm_full <- setNames(rep(NA_integer_, length(nodes)), nodes)
    in_idx <- names(grp)[names(grp) %in% nodes]
    if (length(in_idx)) comm_full[in_idx] <- as.integer(grp[in_idx])
  } else {
    comm_full <- setNames(rep(NA_integer_, length(nodes)), nodes)
  }

  assigned_idx <- which(!is.na(comm_full))
  if (!length(assigned_idx)) {
    empty <- tibble::tibble()
    out <- list(
      bridge_strength    = list(overall = 0,       by_comm = empty),
      bridge_ei1         = list(overall = 0,       by_comm = empty),
      bridge_ei2         = list(overall = 0,       by_comm = empty),
      bridge_closeness   = list(overall = NA_real_, by_comm = empty),
      bridge_betweenness = list(overall = 0,       by_pair = empty)
    )
    class(out) <- c("bridge_profiles", class(out))
    return(out)
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
  out <- list(
    bridge_strength    = strength_list,
    bridge_ei1         = ei1_list,
    bridge_ei2         = ei2_list,
    bridge_closeness   = closeness_list,
    bridge_betweenness = betweenness_list
  )

  class(out) <- c("bridge_profiles", class(out))
  out
}

#' @export
print.bridge_profiles <- function(x,
                                  statistic = c("bridge_strength", "bridge_ei1", "bridge_ei2",
                                                "bridge_closeness", "bridge_betweenness"),
                                  digits = 3,
                                  ...) {
  if (missing(statistic)) {
    for (s in c("bridge_strength", "bridge_ei1", "bridge_ei2",
                "bridge_closeness", "bridge_betweenness")) {
      print.bridge_profiles(x, statistic = s, digits = digits)
    }
    return(invisible(x))
  }

  statistic <- match.arg(statistic)

  cat("\nBridge profile:", statistic, "\n")
  cat(strrep("=", 60), "\n")

  obj <- x[[statistic]]

  if (is.null(obj) || length(obj) == 0) {
    cat("No results available.\n")
    return(invisible(x))
  }

  # ---- OVERALL ----
  cat("\nOverall:\n")
  if (is.na(obj$overall)) {
    cat("  NA\n")
  } else {
    cat("  ", round(obj$overall, digits), "\n")
  }

  # labels used ONLY for printing
  pretty_labels <- c(
    community      = "Community",
    sum_abs_w      = "Contribution (|w|)",
    sum_signed_w   = "Contribution (signed)",
    sum_signed_w2  = "Contribution (EI2)",
    inv_mean_dist  = "Inverse mean distance",
    Ci             = "Community i",
    Cj             = "Community j",
    hits           = "Shortest-path hits"
  )

  prettify <- function(df) {
    df |>
      dplyr::mutate(
        dplyr::across(dplyr::where(is.numeric), \(v) round(v, digits))
      ) |>
      dplyr::rename_with(
        .fn   = \(nm) unname(pretty_labels[nm]),
        .cols = dplyr::any_of(names(pretty_labels))
      )
  }

  # ---- BY COMMUNITY ----
  if (!is.null(obj$by_comm) && nrow(obj$by_comm) > 0) {
    cat("\nBy community:\n")
    print(prettify(obj$by_comm), row.names = FALSE)
  }

  # ---- BY COMMUNITY PAIR ----
  if (!is.null(obj$by_pair) && nrow(obj$by_pair) > 0) {
    cat("\nBy community pair:\n")
    print(prettify(obj$by_pair), row.names = FALSE)
  }

  cat("\n")
  invisible(x)
}
