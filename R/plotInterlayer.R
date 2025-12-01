#' Interlayer layer-pair keys from a multimixMN_fit object
#'
#' @description
#' Internal helper to extract valid interlayer layer-pair identifiers
#' (e.g., "bio_dis") from a \code{multimixMN_fit} object or from its
#' \code{$interlayer} sublist.
#'
#' @keywords internal
#' @noRd
interlayerPairs <- function(object) {
  # Accept either full multimixMN_fit or directly the interlayer sublist
  if (!is.null(object$interlayer) && !is.null(object$layers)) {
    object <- object$interlayer
  }

  if (!is.list(object)) return(character(0))

  nms <- names(object)
  if (is.null(nms)) return(character(0))

  # prendiamo solo gli elementi che hanno $edges (cioè le coppie di layer)
  has_edges <- vapply(
    nms,
    function(nm) {
      x <- object[[nm]]
      is.list(x) && !is.null(x$edges)
    },
    logical(1)
  )

  raw <- nms[has_edges]
  if (!length(raw)) return(character(0))

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
#' Internal helper used by \code{plot.mixmashnet()} to visualize interlayer
#' node metrics or interlayer edge weights with bootstrap confidence intervals.
#'
#' @keywords internal
#' @noRd
#' @importFrom ggplot2 ggplot aes geom_hline geom_errorbar scale_color_manual
#'   geom_point coord_flip labs facet_wrap theme_minimal theme element_text
#'   element_rect
#' @importFrom dplyr mutate arrange desc slice n select filter bind_rows
#'   distinct group_by ungroup
#' @importFrom rlang .data
#' @importFrom stats sd
plotInterlayer <- function(
    fit_multi,
    statistics = c("strength","expected_influence","closeness","betweenness","edges"),
    pairs = "*",
    edges_top_n = 60,
    ordering = c("value","alphabetical"),
    standardize = FALSE,
    title = NULL,
    exclude_nodes = NULL,
    nodes_layer = NULL
) {
  # ---- checks & setup ----
  ordering   <- match.arg(ordering)
  statistics <- match.arg(statistics, several.ok = TRUE)

  if (!is.list(fit_multi) || !"multimixMN_fit" %in% class(fit_multi)) {
    stop("`fit_multi` must be a `multimixMN_fit` object returned by multimixMN().")
  }

  interlayer <- fit_multi$interlayer
  if (is.null(interlayer)) stop("`fit_multi$interlayer` is missing.")

  node_metrics <- c("strength","expected_influence","closeness","betweenness")
  edge_metric  <- "edges"

  types_selected <- c(
    node = any(statistics %in% node_metrics),
    edge = any(statistics %in% edge_metric)
  )
  if (sum(types_selected) > 1L) {
    stop(
      "Cannot combine node statistics and 'edges' in the same call.\n",
      "Please call `plotInterlayer2()` separately for node metrics and interlayer edges."
    )
  }

  # Helper for normalizing pair labels "A_B" (A < B)
  norm_pair_vec <- function(v) {
    vapply(strsplit(v, "_", fixed = TRUE), function(sp) {
      if (length(sp) != 2) return(NA_character_)
      paste(sort(sp), collapse = "_")
    }, character(1))
  }

  # ---------- NODE METRICS BRANCH ----------
  if (types_selected["node"]) {
    # Title default for nodes
    if (is.null(title)) {
      title <- "Interlayer node metrics (95% CI)"
    }

    ct <- interlayer$centrality$true
    if (is.null(ct) || !nrow(ct)) stop("No interlayer$centrality$true found.")

    # Exclude nodes if requested
    if (!is.null(exclude_nodes) && length(exclude_nodes)) {
      ct <- ct[!(ct$node %in% exclude_nodes), , drop = FALSE]
      if (!nrow(ct)) stop("No nodes left after exclude_nodes filter.")
    }

    # Filter by layer if requested
    if (!is.null(nodes_layer)) {
      if (is.null(fit_multi$layers)) {
        stop("`fit_multi$layers` is missing; cannot filter by `nodes_layer`.")
      }
      keep <- fit_multi$layers$assignment[ct$node] == nodes_layer
      ct <- ct[keep %in% TRUE, , drop = FALSE]
      if (!nrow(ct)) {
        stop(sprintf("No nodes found in layer '%s' after filters.", nodes_layer))
      }
    }

    # Map from requested statistics to names in ci_results
    ci_all <- interlayer$centrality$ci_results
    name_map <- c(
      strength           = "strength",
      expected_influence = "ei1",
      closeness          = "closeness",
      betweenness        = "betweenness"
    )
    stats_node <- statistics[statistics %in% names(name_map)]
    if (!length(stats_node)) stop("No valid node statistics selected.")

    dfs <- lapply(stats_node, function(m) {
      internal <- name_map[[m]]

      ci <- ci_all[[internal]]
      if (is.null(ci)) return(NULL)
      ci <- ci[match(ct$node, rownames(ci)), , drop = FALSE]

      data.frame(
        node     = ct$node,
        metric   = m,
        observed = ct[[internal]],
        lower    = ci[, "2.5%"],
        upper    = ci[, "97.5%"],
        stringsAsFactors = FALSE
      )
    })

    df <- dplyr::bind_rows(dfs)
    if (is.null(df) || !nrow(df)) {
      stop("No data for selected interlayer node statistics after filtering.")
    }
    df$metric <- factor(df$metric, levels = stats_node)

    # Optional standardization per metric
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
      ref <- df[df$metric == stats_node[1], , drop = FALSE]
      node_order <- ref$node[order(-ref$observed)]
    } else {
      node_order <- sort(unique(df$node))
    }

    df$includes_zero <- ifelse(is.na(df$lower) | is.na(df$upper), NA,
                               df$lower <= 0 & df$upper >= 0)
    df$node_f <- factor(df$node, levels = rev(node_order))

    # Plot nodes
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$node_f, y = .data$observed)) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = "dashed",
        color = "gray50", linewidth = 0.4
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper, color = .data$includes_zero),
        width = 0.2, na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = c("FALSE" = "black", "TRUE" = "gray60"),
        na.value = "gray80", guide = "none"
      ) +
      ggplot2::geom_point(
        shape = 21, fill = "gray60", color = "black",
        size = 2.5, na.rm = TRUE
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = title,
        x = "Nodes",
        y = if (standardize) "Z-score (95% CI)" else "Estimated (95% CI)"
      ) +
      ggplot2::facet_wrap(~ metric, ncol = length(stats_node), scales = "free_x") +
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

  # ---------- EDGE WEIGHTS BRANCH ----------
  # Helper empty plot
  .empty_plot <- function(msg = "No interlayer edges to plot") {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg) +
      ggplot2::theme_void()
  }

  # Title default per edges
  if (is.null(title)) {
    title <- "Interlayer edge weights (95% CI)"
  }

  nms <- names(interlayer)
  if (is.null(nms)) stop("No interlayer components found.")
  has_edges <- vapply(
    nms,
    function(nm) {
      x <- interlayer[[nm]]
      is.list(x) && !is.null(x$edges)
    },
    logical(1)
  )
  pairs_avail_raw <- nms[has_edges]
  if (!length(pairs_avail_raw)) {
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
    cand <- which(norm_pair_vec(pairs_avail_raw) == pk_norm)
    if (!length(cand)) return(NULL)
    nm <- pairs_avail_raw[cand[1]]
    ed <- interlayer[[nm]]$edges
    if (is.null(ed)) return(NULL)
    et <- ed$true
    ci <- ed$ci
    if (is.null(et) || is.null(ci) || !nrow(et)) return(NULL)

    # Escludi nodi se richiesto (by from/to)
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

    # Solo pesi non nulli
    et <- et[et$weight != 0, , drop = FALSE]
    if (!nrow(et)) return(NULL)

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
      dplyr::select(-"abs_obs")
  }

  # Optional standardization per pair
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

  # Plot edges
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$edge_f, y = .data$observed)) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "gray50", linewidth = 0.4
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper, color = .data$includes_zero),
      width = 0.2, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "gray60"),
      na.value = "gray80", guide = "none"
    ) +
    ggplot2::geom_point(
      shape = 21, fill = "gray60", color = "black",
      size = 2.5, na.rm = TRUE
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      x = "Edges",
      y = if (standardize) "Z-score (95% CI)" else "Estimated (95% CI)"
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
