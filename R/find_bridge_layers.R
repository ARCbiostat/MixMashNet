#' Bridge profiles of a node across layers
#'
#' @description
#' Identifies which layers contribute most to the interlayer bridge role of a
#' given node, by decomposing its interlayer connectivity into layer-specific
#' contributions. The function is designed as an interpretative companion to the
#' interlayer node-level indices returned by \code{multimixMN()}, providing the
#' components underlying the corresponding overall interlayer indices.
#'
#' Interlayer connectivity is summarized using four complementary profiles:
#' interlayer strength, interlayer expected influence (order 1), interlayer
#' closeness, and interlayer betweenness.
#'
#' @details
#' The function operates on the interlayer-only graph, i.e. on the graph
#' containing only edges between nodes belonging to different layers.
#'
#' For a focal node in the selected \code{layer}, the function decomposes:
#' \itemize{
#'   \item interlayer strength into contributions toward each layer;
#'   \item interlayer expected influence (order 1) into signed contributions
#'         toward each layer;
#'   \item interlayer closeness into additive harmonic-distance contributions
#'         toward each layer;
#'   \item interlayer betweenness into additive contributions from shortest
#'         paths between layer pairs, using the standard fraction
#'         \eqn{\sigma_{st}(v) / \sigma_{st}}.
#' }
#'
#' Contributions are defined so that they sum to the corresponding overall
#' interlayer index.
#'
#' The returned object has class \code{"bridge_layer_profiles"} and provides a
#' dedicated \code{print()} method. By default, all interlayer profiles are
#' displayed; a specific profile can be selected through the \code{statistic}
#' argument.
#'
#' @param fit An object of class \code{multimixMN_fit}.
#' @param node Character scalar: node of interest.
#' @param layer Character scalar giving the layer of the focal node.
#'
#' @return An object of class \code{"bridge_layer_profiles"} (a named list) with
#'   the following components:
#' \describe{
#'   \item{\code{bridge_strength}}{List with \code{overall} and \code{by_layer},
#'   where \code{by_layer} is a tibble with columns
#'   \code{target_layer} and \code{sum_abs_w}.}
#'   \item{\code{bridge_ei1}}{List with \code{overall} and \code{by_layer},
#'   where \code{by_layer} is a tibble with columns
#'   \code{target_layer} and \code{sum_signed_w}.}
#'   \item{\code{bridge_closeness}}{List with \code{overall} and
#'   \code{by_layer}, where \code{by_layer} is a tibble with columns
#'   \code{target_layer} and \code{contribution}.}
#'   \item{\code{bridge_betweenness}}{List with \code{overall} and
#'   \code{by_pair}, where \code{by_pair} is a tibble with columns
#'   \code{Li}, \code{Lj}, and \code{contribution}.}
#' }
#'
#' @importFrom igraph graph_from_adjacency_matrix delete_edges distances is_directed degree get.all.shortest.paths E V ecount
#' @importFrom dplyr bind_rows arrange desc mutate across where rename_with any_of
#' @importFrom tibble tibble
#' @importFrom tidyr separate
#' @importFrom utils combn
#' @importFrom stats setNames
#' @export
find_bridge_layers <- function(fit, node, layer) {
  # ---- guardrails ----
  stopifnot(inherits(fit, "multimixMN_fit"))

  if (missing(layer) || is.null(layer) || length(layer) != 1L || !is.character(layer)) {
    stop("`layer` must be a single character string.")
  }

  if (missing(node) || is.null(node) || length(node) != 1L || !is.character(node)) {
    stop("`node` must be a single character string.")
  }

  if (is.null(fit$graph$keep_nodes_graph) ||
      is.null(fit$layers$assignment) ||
      is.null(fit$model$nodes) ||
      is.null(fit$model$mgm$pairwise$wadj)) {
    stop(
      "`fit` must contain `graph$keep_nodes_graph`, `layers$assignment`, ",
      "`model$nodes`, and `model$mgm$pairwise$wadj`."
    )
  }

  nodes <- fit$graph$keep_nodes_graph
  if (!length(nodes)) {
    stop("`graph$keep_nodes_graph` is empty.")
  }

  if (!(node %in% nodes)) {
    stop("Requested `node` is not in `graph$keep_nodes_graph`.")
  }

  layer_assign <- fit$layers$assignment
  if (is.null(names(layer_assign))) {
    stop("`fit$layers$assignment` must be a named vector.")
  }

  if (!(node %in% names(layer_assign))) {
    stop("Requested `node` is not present in `fit$layers$assignment`.")
  }

  node_layer <- as.character(layer_assign[[node]])
  if (is.na(node_layer) || !nzchar(node_layer)) {
    stop("The layer of the requested `node` is missing.")
  }

  if (!identical(node_layer, layer)) {
    stop("Requested `node` does not belong to the specified `layer`.")
  }

  # ---- rebuild signed global W on kept graph nodes ----
  all_nodes <- fit$model$nodes
  wadj <- fit$model$mgm$pairwise$wadj
  signs <- fit$model$mgm$pairwise$signs

  colnames(wadj) <- rownames(wadj) <- all_nodes
  colnames(signs) <- rownames(signs) <- all_nodes

  W_signed <- wadj
  if (!is.null(signs)) {
    idx_sign <- !is.na(signs) & (abs(signs) == 1)
    W_signed[idx_sign] <- wadj[idx_sign] * signs[idx_sign]
  }
  W_signed[is.na(W_signed)] <- 0

  W_signed <- W_signed[nodes, nodes, drop = FALSE]

  # ---- keep only interlayer edges ----
  node_layers <- as.character(layer_assign[nodes])
  node_layers_named <- stats::setNames(node_layers, nodes)

  W_inter <- matrix(
    0,
    nrow = length(nodes),
    ncol = length(nodes),
    dimnames = list(nodes, nodes)
  )

  for (i in seq_along(nodes)) {
    for (j in seq_along(nodes)) {
      if (node_layers[i] != node_layers[j]) {
        W_inter[i, j] <- W_signed[i, j]
      }
    }
  }

  diag(W_inter) <- 0
  W_inter <- (W_inter + t(W_inter)) / 2

  node_pos <- match(node, nodes)
  target_idx <- setdiff(seq_along(nodes), node_pos)

  if (!length(target_idx)) {
    empty <- tibble::tibble()
    out <- list(
      bridge_strength    = list(overall = 0,        by_layer = empty),
      bridge_ei1         = list(overall = 0,        by_layer = empty),
      bridge_closeness   = list(overall = NA_real_, by_layer = empty),
      bridge_betweenness = list(overall = 0,        by_pair = empty)
    )
    class(out) <- c("bridge_layer_profiles", class(out))
    return(out)
  }

  target_layers <- node_layers[target_idx]

  # ================== PROFILES WITHOUT PATHS ==================

  # strength
  strength_list <- {
    idx_by_layer <- split(target_idx, target_layers)
    by_layer <- lapply(idx_by_layer, function(idx) {
      tibble::tibble(
        target_layer = unique(node_layers[idx]),
        sum_abs_w    = sum(abs(W_inter[node, idx]))
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::arrange(dplyr::desc(sum_abs_w))

    list(
      overall = sum(abs(W_inter[node, target_idx])),
      by_layer = by_layer
    )
  }

  # ei1
  ei1_list <- {
    idx_by_layer <- split(target_idx, target_layers)
    by_layer <- lapply(idx_by_layer, function(idx) {
      tibble::tibble(
        target_layer = unique(node_layers[idx]),
        sum_signed_w = sum(W_inter[node, idx])
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::arrange(dplyr::desc(sum_signed_w))

    list(
      overall = sum(W_inter[node, target_idx]),
      by_layer = by_layer
    )
  }

  # ================== PATHS ON |W_inter| ==================
  EPS  <- 1e-10
  Wabs <- abs(W_inter)

  g_absw <- igraph::graph_from_adjacency_matrix(
    Wabs,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  if (is.null(igraph::E(g_absw)$weight)) {
    igraph::E(g_absw)$weight <- 1
  }
  igraph::E(g_absw)$weight <- abs(igraph::E(g_absw)$weight)

  # remove edges with NA or |w| <= EPS
  g_inv <- igraph::delete_edges(
    g_absw,
    which(is.na(igraph::E(g_absw)$weight) | igraph::E(g_absw)$weight <= EPS)
  )

  # invert to get distances
  if (igraph::ecount(g_inv) > 0) {
    igraph::E(g_inv)$weight <- 1 / igraph::E(g_inv)$weight
  }

  g_pos <- if (igraph::ecount(g_inv) > 0) {
    igraph::delete_edges(g_inv, which(igraph::E(g_inv)$weight < 0))
  } else {
    g_inv
  }

  # ---- CLOSENESS (additive decomposition of harmonic closeness) ----
  closeness_list <- if (igraph::ecount(g_pos) == 0) {
    list(overall = NA_real_, by_layer = tibble::tibble())
  } else {
    v_id <- match(node, igraph::V(g_pos)$name)

    d_all <- suppressWarnings(
      igraph::distances(
        g_pos,
        v = v_id,
        to = igraph::V(g_pos),
        weights = igraph::E(g_pos)$weight,
        mode = "all"
      )
    )
    d_all <- as.numeric(d_all)
    names(d_all) <- igraph::V(g_pos)$name

    d_all[node] <- NA_real_

    inv_d <- rep(0, length(d_all))
    names(inv_d) <- names(d_all)
    ok <- is.finite(d_all) & !is.na(d_all) & d_all > 0
    inv_d[ok] <- 1 / d_all[ok]

    denom <- length(nodes) - 1L

    overall <- if (denom <= 0) {
      NA_real_
    } else {
      sum(inv_d, na.rm = TRUE) / denom
    }

    other_nodes <- setdiff(nodes, node)

    by_layer <- lapply(split(other_nodes, node_layers_named[other_nodes]), function(nn) {
      tibble::tibble(
        target_layer = unique(node_layers_named[nn]),
        contribution = if (denom <= 0) NA_real_ else sum(inv_d[nn], na.rm = TRUE) / denom
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::arrange(dplyr::desc(contribution))

    list(
      overall = overall,
      by_layer = by_layer
    )
  }

  # ---- BETWEENNESS (additive decomposition of weighted betweenness) ----
  v_id <- match(node, igraph::V(g_pos)$name)

  if (is.na(v_id) || igraph::degree(g_pos, v = v_id) == 0) {
    betweenness_list <- list(
      overall = 0,
      by_pair = tibble::tibble(Li = character(), Lj = character(), contribution = numeric())
    )
  } else {
    v_names <- igraph::V(g_pos)$name
    other_vertices <- setdiff(v_names, node)

    if (length(other_vertices) < 2L) {
      pair_mat <- matrix(character(0), nrow = 2L)
    } else {
      pair_mat <- utils::combn(other_vertices, 2, simplify = TRUE)
    }

    total_contrib <- 0
    by_pair_tab <- list()

    mid_ids <- function(path_ids) {
      ids <- as.vector(path_ids)
      if (length(ids) <= 2L) integer(0) else ids[-c(1L, length(ids))]
    }

    if (ncol(pair_mat) > 0) {
      for (k in seq_len(ncol(pair_mat))) {
        src <- pair_mat[1, k]
        tgt <- pair_mat[2, k]

        sp <- suppressWarnings(
          igraph::get.all.shortest.paths(
            g_pos,
            from = src,
            to = tgt,
            weights = igraph::E(g_pos)$weight,
            mode = "all"
          )
        )

        if (!length(sp$res)) next

        k_tot <- length(sp$res)
        k_hit <- 0L

        for (p in sp$res) {
          if (any(mid_ids(p) == v_id)) {
            k_hit <- k_hit + 1L
          }
        }

        if (k_hit > 0L) {
          contrib <- k_hit / k_tot
          total_contrib <- total_contrib + contrib

          lp <- sort(c(node_layers_named[[src]], node_layers_named[[tgt]]))
          key <- paste(lp, collapse = "--")

          by_pair_tab[[key]] <- (if (is.null(by_pair_tab[[key]])) 0 else by_pair_tab[[key]]) + contrib
        }
      }
    }

    betweenness_by_pair <- if (length(by_pair_tab)) {
      tibble::tibble(
        pair = names(by_pair_tab),
        contribution = as.numeric(unlist(by_pair_tab))
      ) |>
        tidyr::separate(pair, into = c("Li", "Lj"), sep = "--", convert = FALSE) |>
        dplyr::arrange(dplyr::desc(contribution))
    } else {
      tibble::tibble(Li = character(), Lj = character(), contribution = numeric())
    }

    betweenness_list <- list(
      overall = as.numeric(total_contrib),
      by_pair = betweenness_by_pair
    )
  }

  # ---- output ----
  out <- list(
    bridge_strength    = strength_list,
    bridge_ei1         = ei1_list,
    bridge_closeness   = closeness_list,
    bridge_betweenness = betweenness_list
  )

  class(out) <- c("bridge_layer_profiles", class(out))
  out
}

#' @export
print.bridge_layer_profiles <- function(
    x,
    statistic = c("bridge_strength", "bridge_ei1",
                  "bridge_closeness", "bridge_betweenness"),
    digits = 3,
    ...
) {
  stat_labels <- c(
    bridge_strength    = "Interlayer Strength",
    bridge_ei1         = "Interlayer Expected Influence (EI1)",
    bridge_closeness   = "Interlayer Closeness",
    bridge_betweenness = "Interlayer Betweenness"
  )

  if (missing(statistic)) {
    for (s in names(stat_labels)) {
      print.bridge_layer_profiles(x, statistic = s, digits = digits)
    }
    return(invisible(x))
  }

  statistic <- match.arg(statistic)

  cat("\n", stat_labels[[statistic]], "\n", sep = "")
  cat(strrep("=", 60), "\n")

  obj <- x[[statistic]]

  if (is.null(obj) || length(obj) == 0) {
    cat("No results available.\n")
    return(invisible(x))
  }

  cat("\nOverall:\n")
  if (is.na(obj$overall)) {
    cat("  NA\n")
  } else {
    cat("  ", round(obj$overall, digits), "\n")
  }

  pretty_labels <- c(
    target_layer = "Target layer",
    sum_abs_w    = "Contribution (|w|)",
    sum_signed_w = "Contribution (signed)",
    contribution = "Contribution",
    Li           = "Layer i",
    Lj           = "Layer j"
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

  if (!is.null(obj$by_layer) && nrow(obj$by_layer) > 0) {

    df <- obj$by_layer
    num_col <- names(df)[sapply(df, is.numeric)]

    if (length(num_col) == 1L) {
      df_nonzero <- df[!is.na(df[[num_col]]) & abs(df[[num_col]]) > .Machine$double.eps, , drop = FALSE]

      if (nrow(df_nonzero) > 0) {
        df <- df_nonzero
      }
    }

    cat("\nBy layer:\n")
    print(prettify(df), row.names = FALSE)
  }

  if (!is.null(obj$by_pair) && nrow(obj$by_pair) > 0) {
    cat("\nBy layer pair:\n")
    print(prettify(obj$by_pair), row.names = FALSE)
  }

  cat("\n")
  invisible(x)
}
