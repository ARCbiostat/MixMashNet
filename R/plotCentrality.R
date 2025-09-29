#' Plot bootstrap centrality with 95% CIs (nodes or edges)
#'
#' @description
#' Creates lollipop/errorbar plots of bootstrap estimates with 95% CIs for:
#' - node metrics (general): strength, expected_influence, closeness, betweenness
#' - node metrics (bridge): bridge_strength, bridge_closeness, bridge_betweenness, bridge_ei1, bridge_ei2
#' - node metrics (excluded): bridge_*_excluded
#' - edge metrics: edge_weights (top-N by absolute weight, optional)
#'
#' Only one *type* of metrics can be plotted at a time (general OR bridge OR excluded OR edges).
#'
#' @param fit A \code{mixMN_fit} object from \code{mixMN()}.
#' @param metrics Character vector of metrics to plot. See details.
#' @param title Plot title.
#' @param plot_nonzero_edges_only Logical; if TRUE, keeps only non-zero observed edges.
#' @param ordering One of \code{"value"}, \code{"alphabetical"}, \code{"community"}.
#' @param exclude_nodes Optional character vector of node names to exclude.
#' @param standardize Logical; if TRUE, z-standardizes observed and CI values (per panel).
#' @param color_by_community Logical; color points by community (nodes only).
#' @param edges_top_n Integer; keep only the top-N edges by |weight| (when plotting edges).
#'
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_hline geom_errorbar geom_point
#' @importFrom ggplot2 scale_color_manual scale_fill_manual scale_x_continuous coord_flip
#' @importFrom ggplot2 labs facet_wrap theme_minimal theme element_text element_rect
#' @importFrom dplyr mutate arrange desc slice n select filter pull bind_rows distinct row_number
#' @importFrom magrittr %>%
#' @importFrom stats sd
#' @export
plotCentrality <- function(
    fit,
    metrics = c(
      "strength", "expected_influence", "closeness", "betweenness",
      "bridge_strength", "bridge_closeness", "bridge_betweenness",
      "bridge_ei1", "bridge_ei2",
      "bridge_strength_excluded", "bridge_betweenness_excluded",
      "bridge_closeness_excluded", "bridge_ei1_excluded", "bridge_ei2_excluded",
      "edge_weights"
    ),
    title = "Bootstrap Centrality with 95% CI",
    plot_nonzero_edges_only = TRUE,
    ordering = c("value", "alphabetical", "community"),
    exclude_nodes = NULL,
    standardize = FALSE,
    color_by_community = TRUE,
    edges_top_n = 60
) {
  if (!is.list(fit) || is.null(fit$ci_results))
    stop("Input 'fit' must be a 'mixMN_fit'-like list produced by mixMN().")

  metrics  <- match.arg(metrics, several.ok = TRUE)
  ordering <- match.arg(ordering, choices = c("value","alphabetical","community"))

  general_metrics  <- c("strength","expected_influence","closeness","betweenness")
  bridge_metrics   <- c("bridge_strength","bridge_closeness","bridge_betweenness","bridge_ei1","bridge_ei2")
  excluded_metrics <- c("bridge_strength_excluded","bridge_closeness_excluded",
                        "bridge_betweenness_excluded","bridge_ei1_excluded","bridge_ei2_excluded")
  edge_metric      <- "edge_weights"

  types_selected <- c(
    general  = any(metrics %in% general_metrics),
    bridge   = any(metrics %in% bridge_metrics),
    excluded = any(metrics %in% excluded_metrics),
    edge     = any(metrics %in% edge_metric)
  )
  if (sum(types_selected) > 1) {
    stop(
      "Cannot combine metrics of different types.\n",
      "- General metrics can only be plotted together.\n",
      "- Bridge metrics can only be plotted together.\n",
      "- Excluded bridge metrics can only be plotted together.\n",
      "- Edge weights can only be plotted alone."
    )
  }

  metric_map <- c(
    strength = "strength",
    expected_influence = "ei1",
    closeness = "closeness",
    betweenness = "betweenness",
    bridge_strength = "bridge_strength",
    bridge_closeness = "bridge_closeness",
    bridge_betweenness = "bridge_betweenness",
    bridge_ei1 = "bridge_ei1",
    bridge_ei2 = "bridge_ei2",
    bridge_strength_excluded = "bridge_strength_excluded",
    bridge_betweenness_excluded = "bridge_betweenness_excluded",
    bridge_closeness_excluded = "bridge_closeness_excluded",
    bridge_ei1_excluded = "bridge_ei1_excluded",
    bridge_ei2_excluded = "bridge_ei2_excluded",
    edge_weights = "edge_weights"
  )

  df_all <- list()

  for (metric in metrics) {
    df_list <- list()

    if (metric == "edge_weights") {
      edges_true <- fit$edges_true
      ci <- fit$ci_results$edge_weights
      if (is.null(edges_true) || nrow(edges_true) == 0 || is.null(ci)) next

      if (isTRUE(plot_nonzero_edges_only)) {
        edges_true <- edges_true[edges_true$weight != 0, , drop = FALSE]
      }
      ci <- ci[match(edges_true$edge, rownames(ci)), , drop = FALSE]

      df <- data.frame(
        node = edges_true$edge,
        observed = edges_true$weight,
        lower = ci[, "2.5%"],
        upper = ci[, "97.5%"],
        metric = "edge_weights",
        community = NA
      )
      df$includes_zero <- ifelse(is.na(df$lower) | is.na(df$upper), NA, df$lower <= 0 & df$upper >= 0)

      if (!is.null(edges_top_n) && is.finite(edges_top_n) && edges_top_n > 0) {
        df <- dplyr::mutate(df, abs_obs = abs(observed))
        df <- dplyr::arrange(df, dplyr::desc(abs_obs))
        df <- dplyr::slice(df, 1:min(edges_top_n, dplyr::n()))
        df <- dplyr::select(df, -abs_obs)
      }

      if (isTRUE(standardize)) {
        mean_val <- mean(df$observed, na.rm = TRUE)
        sd_val   <- stats::sd(df$observed, na.rm = TRUE)
        df$observed <- (df$observed - mean_val) / sd_val
        df$lower    <- (df$lower    - mean_val) / sd_val
        df$upper    <- (df$upper    - mean_val) / sd_val
      }

      # Complete display fields for edges (as in the original)
      df <- df %>%
        dplyr::mutate(
          lower = as.numeric(lower),
          upper = as.numeric(upper),
          community_clean  = NA_character_,
          community_factor = factor(NA_character_)
        )

      # Ordering for edges (community ordering falls back to value)
      if (ordering == "value") {
        df <- df %>% dplyr::arrange(dplyr::desc(observed)) %>% dplyr::mutate(order = dplyr::row_number())
      } else if (ordering == "alphabetical") {
        df <- df %>% dplyr::arrange(node) %>% dplyr::mutate(order = dplyr::row_number())
      } else { # "community" not meaningful for edges
        df <- df %>% dplyr::arrange(dplyr::desc(observed)) %>% dplyr::mutate(order = dplyr::row_number())
      }
      df$order_reversed <- max(df$order) - df$order + 1

      df$label_colored <- df$node

      df_list[[metric]] <- df
      df_all[[metric]]  <- do.call(rbind, df_list)
      next
    }

    # ------- Node metrics -------
    mat_boot <- switch(metric,
                       strength                    = fit$strength_boot,
                       expected_influence          = fit$ei1_boot,
                       closeness                   = fit$closeness_boot,
                       betweenness                 = fit$betweenness_boot,
                       bridge_strength             = fit$bridge_strength_boot,
                       bridge_closeness            = fit$bridge_closeness_boot,
                       bridge_betweenness          = fit$bridge_betweenness_boot,
                       bridge_ei1                  = fit$bridge_ei1_boot,
                       bridge_ei2                  = fit$bridge_ei2_boot,
                       bridge_strength_excluded    = fit$bridge_strength_excl_boot,
                       bridge_betweenness_excluded = fit$bridge_betweenness_excl_boot,
                       bridge_closeness_excluded   = fit$bridge_closeness_excl_boot,
                       bridge_ei1_excluded         = fit$bridge_ei1_excl_boot,
                       bridge_ei2_excluded         = fit$bridge_ei2_excl_boot)

    ci <- fit$ci_results[[metric]]
    if (is.null(mat_boot) || ncol(mat_boot) == 0) {
      stop(sprintf("Bootstrap matrix for metric '%s' is empty or NULL.", metric))
    }

    node_set <- colnames(mat_boot)
    ct_filtered <- fit$centrality_true %>% dplyr::filter(node %in% node_set)
    order_index <- match(node_set, ct_filtered$node)
    observed_values <- ct_filtered[[metric_map[metric]]][order_index]

    community <- if (!is.null(fit$groups)) fit$groups[node_set] else rep(NA, length(node_set))

    df <- data.frame(
      node = node_set,
      observed = observed_values,
      lower = ci[, "2.5%"],
      upper = ci[, "97.5%"],
      metric = metric,
      community = community
    )
    df$includes_zero <- ifelse(is.na(df$lower) | is.na(df$upper), NA, df$lower <= 0 & df$upper >= 0)

    if (isTRUE(standardize)) {
      mean_val <- mean(df$observed, na.rm = TRUE)
      sd_val   <- stats::sd(df$observed, na.rm = TRUE)
      df$observed <- (df$observed - mean_val) / sd_val
      df$lower    <- (df$lower    - mean_val) / sd_val
      df$upper    <- (df$upper    - mean_val) / sd_val
    }

    # keep/exclude coherent with metric type
    if (grepl("_excluded$", metric)) {
      df <- df %>% dplyr::filter(is.na(community))
    } else if (grepl("^bridge_", metric)) {
      df <- df %>% dplyr::filter(!is.na(community))
    }

    # ---- Build df_plot per original logic ----
    df_plot <- df
    if (!is.null(exclude_nodes)) {
      df_plot <- df_plot %>% dplyr::filter(!(node %in% exclude_nodes))
    }

    df_plot <- df_plot %>%
      dplyr::mutate(
        lower = as.numeric(lower),
        upper = as.numeric(upper),
        community_clean = ifelse(metric == "edge_weights",
                                 NA_character_,
                                 ifelse(is.na(community), "Excluded", as.character(community)))
      )

    community_levels_present <- unique(df_plot$community_clean)
    community_numeric <- sort(community_levels_present[community_levels_present != "Excluded"])
    community_ordered <- c(community_numeric, "Excluded")

    cmap <- c(fit$community_palette, Excluded = "gray70")
    community_colors <- cmap[community_ordered]

    df_plot$community_factor <- factor(df_plot$community_clean, levels = community_ordered)

    df_plot <- df_plot %>%
      dplyr::mutate(
        node_factor      = factor(node),
        node_order_value = observed,
        node_order_alpha = as.character(node),
        node_order_comm  = as.numeric(factor(community_factor))
      )

    if (ordering == "value") {
      df_plot <- df_plot %>% dplyr::arrange(dplyr::desc(node_order_value)) %>% dplyr::mutate(order = dplyr::row_number())
    } else if (ordering == "alphabetical") {
      df_plot <- df_plot %>% dplyr::arrange(node_order_alpha) %>% dplyr::mutate(order = dplyr::row_number())
    } else if (ordering == "community") {
      df_plot <- df_plot %>% dplyr::arrange(node_order_comm, dplyr::desc(node_order_value)) %>% dplyr::mutate(order = dplyr::row_number())
    }

    df_plot$order_reversed <- max(df_plot$order) - df_plot$order + 1

    df_plot$label_colored <- as.character(df_plot$node)

    df_list[[metric]] <- df_plot
    df_all[[metric]]  <- do.call(rbind, df_list)
  }

  # Reference order from the first metric
  reference_metric <- metrics[1]
  reference_df <- df_all[[reference_metric]]

  if (ordering == "value") {
    node_order <- reference_df %>% dplyr::arrange(dplyr::desc(observed)) %>% dplyr::pull(node)
  } else if (ordering == "alphabetical") {
    node_order <- reference_df %>% dplyr::arrange(node) %>% dplyr::pull(node)
  } else if (ordering == "community") {
    node_order <- reference_df %>% dplyr::arrange(as.numeric(factor(community_factor)), dplyr::desc(observed)) %>% dplyr::pull(node)
  }

  # Apply order to each panel and combine
  df_all <- lapply(df_all, function(df) {
    df %>%
      dplyr::mutate(
        node_factor = factor(node, levels = node_order),
        order_reversed = as.numeric(factor(node, levels = rev(node_order)))
      )
  })
  df_all_combined <- dplyr::bind_rows(df_all)
  df_all_combined$metric <- factor(df_all_combined$metric, levels = metrics)

  # Ensure label_colored exists (belt & suspenders)
  if (!"label_colored" %in% names(df_all_combined)) {
    df_all_combined$label_colored <- paste0("<b>", df_all_combined$node, "</b>")
  }

  # Palette including 'Excluded'
  cmap <- c(fit$community_palette, Excluded = "gray70")

  # Label map
  lab_map <- df_all_combined %>%
    dplyr::select(order_reversed, label_colored) %>%
    dplyr::distinct()

  # Build plot
  p <- ggplot2::ggplot(df_all_combined, ggplot2::aes(x = order_reversed, y = observed)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper, color = includes_zero),
                           width = 0.2, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "gray60"),
                                na.value = "gray80", guide = "none")

  plotting_edges <- all(df_all_combined$metric == "edge_weights")

  if (isTRUE(color_by_community) && !plotting_edges) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(fill = community_factor),
        shape = 21, color = "black", size = 2.5, na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(values = cmap, name = "Community")
  } else {
    p <- p + ggplot2::geom_point(shape = 21, fill = "gray60", color = "black", size = 2.5, na.rm = TRUE)
  }

  x_lab <- if (plotting_edges) "Edges" else "Nodes"

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = function(x) unique(df_all_combined$order_reversed),
      labels = function(x) {
        key <- setNames(lab_map$label_colored, lab_map$order_reversed)
        key[as.character(x)]
      },
      expand = c(0.01, 0.01)
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      subtitle = if (standardize)
        "Gray error bars indicate confidence intervals that include zero in the non-standardized scale."
      else NULL,
      x = x_lab,
      y = if (standardize) "Z-score of observed value with 95% CI" else "Observed value with 95% CI"
    ) +
    ggplot2::facet_wrap(~ metric, ncol = length(metrics), scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text   = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x  = ggplot2::element_text(size = 9, face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 8, face = "bold"),
      plot.title   = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = if (isTRUE(color_by_community) && !plotting_edges) "right" else "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  return(p)

}
