#' List available interlayer pair keys
#'
#' @description
#' Returns the available interlayer pair keys (layer pairs) from a
#' \code{multimixMN_fit} object, normalized as \code{"A_B"} with alphabetical order.
#'
#' @param object An object returned by \code{multimixMN()} (\code{multimixMN_fit})
#'   or its \code{$interlayer} sublist.
#'
#' @return A character vector of pair keys like \code{"bio_dis"}.
#' @export
interlayerPairs <- function(object) {
  # Accept either full multimixMN_fit or directly the interlayer sublist
  if (!is.null(object$interlayer) && !is.null(object$layers)) {
    object <- object$interlayer
  }
  if (is.null(object$edges)) return(character(0))

  raw <- names(object$edges)
  if (is.null(raw) || !length(raw)) return(character(0))

  # Normalize pair keys as "A_B" with A < B
  norm_pair <- function(pk) {
    sp <- strsplit(pk, "_", fixed = TRUE)[[1]]
    if (length(sp) != 2) return(pk)
    paste(sort(sp), collapse = "_")
  }
  unique(vapply(raw, norm_pair, character(1)))
}


#' Plot interlayer node metrics or interlayer edge weights with 95% CIs
#'
#' @description
#' Works on objects returned by \code{multimixMN()} (class \code{multimixMN_fit}).
#' You can plot either:
#' \itemize{
#'   \item \strong{mode="nodes"}: interlayer-only node metrics with percentile CIs
#'         (\code{strength}, \code{ei1}, \code{closeness}, \code{betweenness});
#'   \item \strong{mode="edges"}: interlayer edge weights by layer-pair with percentile CIs.
#' }
#'
#' @param fit_multi A \code{multimixMN_fit} object (the full output of \code{multimixMN()}).
#' @param mode One of \code{"nodes"} or \code{"edges"}.
#' @param metrics For \code{mode="nodes"}: any of
#'   \code{c("strength","ei1","closeness","betweenness")}.
#' @param pairs For \code{mode="edges"}: character vector of pair keys \code{"A_B"}
#'   (order-insensitive; \code{"bio_dis" == "dis_bio"}), or \code{"*"} to plot all.
#' @param edges_top_n Keep top-N edges by \code{|weight|} after filtering (set \code{NULL} to keep all).
#' @param ordering \code{"value"} (descending observed) or \code{"alphabetical"}.
#' @param standardize If \code{TRUE}, z-standardizes observed and CIs per panel.
#' @param plot_nonzero_edges_only If \code{TRUE}, drops edges with observed weight == 0.
#' @param title Plot title (auto-chosen if \code{NULL}).
#' @param exclude_nodes Optional character vector of node names to exclude (both modes).
#' @param nodes_layer (nodes-mode) Optional single layer name to keep only nodes from that layer.
#'        The node→layer map is taken from \code{fit_multi$layers}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_errorbar scale_color_manual
#'   geom_point coord_flip labs facet_wrap theme_minimal theme element_text
#'   element_rect
#' @importFrom dplyr mutate arrange desc slice n select filter bind_rows
#'   distinct group_by ungroup
#' @importFrom rlang .data
#' @importFrom stats sd
#' @export
plotInterlayer <- function(
    fit_multi,
    mode = c("nodes", "edges"),
    metrics = c("strength","ei1","closeness","betweenness"),
    pairs = "*",
    edges_top_n = 60,
    ordering = c("value","alphabetical"),
    standardize = FALSE,
    plot_nonzero_edges_only = TRUE,
    title = NULL,
    exclude_nodes = NULL,
    nodes_layer = NULL
) {

  mode <- match.arg(mode)
  ordering <- match.arg(ordering)

  if (!is.list(fit_multi) || !"multimixMN_fit" %in% class(fit_multi)) {
    stop("`fit_multi` must be a `multimixMN_fit` object returned by multimixMN().")
  }

  interlayer <- fit_multi$interlayer
  if (is.null(interlayer)) stop("`fit_multi$interlayer` is missing.")

  if (is.null(title)) {
    title <- if (mode == "nodes") "Interlayer node metrics (95% CI)"
    else "Interlayer edge weights (95% CI)"
  }

  # ---------- helpers ----------
  norm_pair_vec <- function(v) {
    vapply(strsplit(v, "_", fixed = TRUE), function(sp) {
      if (length(sp) != 2) return(NA_character_)
      paste(sort(sp), collapse = "_")
    }, character(1))
  }

  # ------------------- NODES -------------------
  if (mode == "nodes") {
    ct <- interlayer$centrality_true
    if (is.null(ct) || !nrow(ct)) stop("No interlayer$centrality_true found.")

    # Exclude nodes if requested
    if (!is.null(exclude_nodes) && length(exclude_nodes)) {
      ct <- ct[!(ct$node %in% exclude_nodes), , drop = FALSE]
      if (!nrow(ct)) stop("No nodes left after exclude_nodes filter.")
    }

    # Filter by layer if requested
    if (!is.null(nodes_layer)) {
      if (is.null(fit_multi$layers)) stop("`fit_multi$layers` is missing; cannot filter by `nodes_layer`.")
      keep <- fit_multi$layers[ct$node] == nodes_layer
      ct <- ct[keep %in% TRUE, , drop = FALSE]
      if (!nrow(ct)) stop(sprintf("No nodes found in layer '%s' after filters.", nodes_layer))
    }

    ci_all <- interlayer$ci_results
    ci_name_map <- c(
      strength   = "strength",
      ei1        = "expected_influence",
      closeness  = "closeness",
      betweenness= "betweenness"
    )
    metrics <- metrics[metrics %in% names(ci_name_map)]
    if (!length(metrics)) stop("No valid metrics selected for nodes.")

    # Build long df for chosen metrics
    dfs <- lapply(metrics, function(m) {
      ci <- ci_all[[ ci_name_map[[m]] ]]
      if (is.null(ci)) return(NULL)
      ci <- ci[match(ct$node, rownames(ci)), , drop = FALSE]
      data.frame(
        node     = ct$node,
        metric   = m,
        observed = ct[[m]],
        lower    = ci[, "2.5%"],
        upper    = ci[, "97.5%"],
        stringsAsFactors = FALSE
      )
    })
    df <- dplyr::bind_rows(dfs)
    if (is.null(df) || !nrow(df)) stop("No data for selected node metrics after filtering.")
    df$metric <- factor(df$metric, levels = metrics)

    # Optional standardization per metric panel
    if (isTRUE(standardize)) {
      df <- df |>
        dplyr::group_by(.data$metric) |>
        dplyr::mutate(
          m = mean(.data$observed, na.rm = TRUE),
          s = stats::sd(.data$observed, na.rm = TRUE),
          observed = (.data$observed - .data$m)/.data$s,
          lower    = (.data$lower    - .data$m)/.data$s,
          upper    = (.data$upper    - .data$m)/.data$s
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-.data$m, -.data$s)
    }

    # Ordering
    if (ordering == "value") {
      ref <- df[df$metric == metrics[1], , drop = FALSE]
      node_order <- ref$node[order(-ref$observed)]
    } else {
      node_order <- sort(unique(df$node))
    }
    df$includes_zero <- ifelse(is.na(df$lower) | is.na(df$upper), NA,
                               df$lower <= 0 & df$upper >= 0)
    df$node_f <- factor(df$node, levels = rev(node_order))

    # Plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$node_f, y = .data$observed)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "gray50", linewidth = 0.4) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper, color = .data$includes_zero),
        width = 0.2, na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "gray60"),
                                  na.value = "gray80", guide = "none") +
      ggplot2::geom_point(shape = 21, fill = "gray60", color = "black",
                          size = 2.5, na.rm = TRUE) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = title,
        x = "Nodes",
        y = if (standardize) "Z-score (95% CI)" else "Observed (95% CI)"
      ) +
      ggplot2::facet_wrap(~ metric, ncol = length(metrics), scales = "free_x") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        strip.text   = ggplot2::element_text(size = 12, face = "bold"),
        axis.text.x  = ggplot2::element_text(size = 9,  face = "bold"),
        axis.text.y  = ggplot2::element_text(size = 8,  face = "bold"),
        plot.title   = ggplot2::element_text(size = 14, face = "bold"),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    return(p)
  }

  # ------------------- EDGES -------------------
  .empty_plot <- function(msg = "No interlayer edges to plot") {
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = msg) + ggplot2::theme_void()
  }

  if (is.null(interlayer$edges)) stop("No interlayer$edges found.")
  pairs_avail_raw <- names(interlayer$edges)
  if (is.null(pairs_avail_raw) || !length(pairs_avail_raw)) {
    stop("No interlayer edge containers found.")
  }
  pairs_avail <- unique(norm_pair_vec(pairs_avail_raw))

  # choose pairs (order-insensitive)
  if (length(pairs) == 1L && identical(pairs, "*")) {
    sel <- pairs_avail
  } else {
    sel <- intersect(unique(norm_pair_vec(pairs)), pairs_avail)
  }
  if (!length(sel)) {
    warning(sprintf(
      "Requested pairs not available: %s. Available pairs: %s",
      paste(pairs, collapse = ", "),
      if (length(pairs_avail)) paste(pairs_avail, collapse = ", ") else "none"
    ), call. = FALSE)
    sel <- pairs_avail
    if (!length(sel)) {
      return(.empty_plot("No interlayer pair containers available"))
    }
  }

  dl <- lapply(sel, function(pk_norm) {
    cand <- which(norm_pair_vec(names(interlayer$edges)) == pk_norm)
    if (!length(cand)) return(NULL)
    ed <- interlayer$edges[[ cand[1] ]]
    et <- ed$edges_true
    ci <- ed$ci_edges
    if (is.null(et) || is.null(ci) || !nrow(et)) return(NULL)

    # Exclude nodes (drop edges if either endpoint excluded)
    if (!is.null(exclude_nodes) && length(exclude_nodes)) {
      if (!("from" %in% names(et)) || !("to" %in% names(et))) {
        if ("edge" %in% names(et)) {
          parts <- strsplit(et$edge, "--",  fixed = TRUE)
          et$from <- vapply(parts, `[`, 1L, FUN.VALUE = character(1))
          et$to   <- vapply(parts, `[`, 2L, FUN.VALUE = character(1))
        }
      }
      if (all(c("from","to") %in% names(et))) {
        keep_edge <- !(et$from %in% exclude_nodes | et$to %in% exclude_nodes)
        et <- et[keep_edge, , drop = FALSE]
      }
      if (!nrow(et)) return(NULL)
    }

    if (isTRUE(plot_nonzero_edges_only)) {
      et <- et[et$weight != 0, , drop = FALSE]
      if (!nrow(et)) return(NULL)
    }

    ci <- ci[match(et$edge, rownames(ci)), , drop = FALSE]

    data.frame(
      pair     = pk_norm,
      edge     = et$edge,
      observed = et$weight,
      lower    = ci[, "2.5%"],
      upper    = ci[, "97.5%"],
      stringsAsFactors = FALSE
    )
  })
  df <- dplyr::bind_rows(dl)
  if (is.null(df) || !nrow(df)) {
    return(.empty_plot(sprintf(
      "No interlayer edges after filtering.\nPairs used: %s",
      if (length(sel)) paste(sel, collapse = ", ") else "none"
    )))
  }

  # Keep top-N by |observed|
  if (!is.null(edges_top_n) && is.finite(edges_top_n) && edges_top_n > 0) {
    df <- df |>
      dplyr::mutate(abs_obs = abs(.data$observed)) |>
      dplyr::arrange(dplyr::desc(.data$abs_obs)) |>
      dplyr::slice(1:min(edges_top_n, dplyr::n())) |>
      dplyr::select(-.data$abs_obs)
  }

  # Optional standardization per pair panel
  if (isTRUE(standardize)) {
    df <- df |>
      dplyr::group_by(.data$pair) |>
      dplyr::mutate(
        m = mean(.data$observed, na.rm = TRUE),
        s = stats::sd(.data$observed, na.rm = TRUE),
        observed = (.data$observed - .data$m)/.data$s,
        lower    = (.data$lower    - .data$m)/.data$s,
        upper    = (.data$upper    - .data$m)/.data$s
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-.data$m, -.data$s)
  }

  # Ordering within each pair panel
  if (ordering == "value") {
    df <- df |>
      dplyr::group_by(.data$pair) |>
      dplyr::arrange(dplyr::desc(.data$observed), .by_group = TRUE) |>
      dplyr::ungroup()
  } else {
    df <- df |>
      dplyr::group_by(.data$pair) |>
      dplyr::arrange(.data$edge, .by_group = TRUE) |>
      dplyr::ungroup()
  }

  df <- df |>
    dplyr::group_by(.data$pair) |>
    dplyr::mutate(edge_f = factor(.data$edge, levels = rev(.data$edge))) |>
    dplyr::ungroup()

  df$includes_zero <- ifelse(is.na(df$lower) | is.na(df$upper), NA,
                             df$lower <= 0 & df$upper >= 0)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$edge_f, y = .data$observed)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "gray50", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper, color = .data$includes_zero),
      width = 0.2, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "gray60"),
                                na.value = "gray80", guide = "none") +
    ggplot2::geom_point(shape = 21, fill = "gray60", color = "black",
                        size = 2.5, na.rm = TRUE) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      x = "Edges",
      y = if (standardize) "Z-score (95% CI)" else "Observed (95% CI)"
    ) +
    ggplot2::facet_wrap(~ pair, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text   = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x  = ggplot2::element_text(size = 9,  face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 7.5),
      plot.title   = ggplot2::element_text(size = 14, face = "bold"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  return(p)
}

