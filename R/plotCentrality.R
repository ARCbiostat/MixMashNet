#' Plot bootstrap centrality with CIs (nodes or edges)
#'
#' @description
#' Internal helper used by \code{plot.mixmashnet()} to visualize node and edge
#' centrality with bootstrap confidence intervals.
#'
#' @keywords internal
#' @noRd
#' @importFrom ggplot2 ggplot aes geom_errorbar geom_point
#' @importFrom ggplot2 scale_fill_manual scale_x_continuous coord_flip
#' @importFrom ggplot2 labs facet_wrap theme_minimal theme element_text element_rect
#' @importFrom dplyr mutate arrange desc slice n select filter pull bind_rows distinct row_number
#' @importFrom magrittr %>%
#' @importFrom stats sd
plotCentrality <- function(
    fit,
    statistics = c(
      "strength", "expected_influence", "closeness", "betweenness",
      "bridge_strength", "bridge_closeness", "bridge_betweenness",
      "bridge_ei1", "bridge_ei2",
      "bridge_strength_excluded", "bridge_betweenness_excluded",
      "bridge_closeness_excluded", "bridge_ei1_excluded", "bridge_ei2_excluded",
      "edges"
    ),
    title = NULL,
    ordering = c("value", "alphabetical", "community"),
    exclude_nodes = NULL,
    standardize = FALSE,
    color_by_community = TRUE,
    edges_top_n = 60
) {
  if (!inherits(fit, "mixMN_fit"))
    stop("`fit` must be an object of class 'mixMN_fit' returned by mixMN().")
  if (is.null(fit$statistics) || is.null(fit$statistics$node))
    stop("`fit$statistics$node` is missing. Did you run mixMN() with bootstrap?")

  # ---- inherit CI level from fit ----
  conf_level <- fit$settings$conf_level
  if (is.null(conf_level) || !is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    conf_level <- 0.95
  }
  alpha <- (1 - conf_level) / 2
  p_lo  <- alpha
  p_hi  <- 1 - alpha
  lab_lo <- paste0(formatC(100 * p_lo, format = "f", digits = 1), "%")
  lab_hi <- paste0(formatC(100 * p_hi, format = "f", digits = 1), "%")

  ci_txt <- paste0(round(100 * conf_level), "% CI")

  if (is.null(title)) {
    title <- paste0("Bootstrap Centrality with ", ci_txt)
  }

  statistics  <- match.arg(statistics, several.ok = TRUE)
  ordering <- match.arg(ordering, choices = c("value","alphabetical","community"))

  general_metrics  <- c("strength","expected_influence","closeness","betweenness")
  bridge_metrics   <- c("bridge_strength","bridge_closeness","bridge_betweenness","bridge_ei1","bridge_ei2")
  excluded_metrics <- c("bridge_strength_excluded","bridge_closeness_excluded",
                        "bridge_betweenness_excluded","bridge_ei1_excluded","bridge_ei2_excluded")
  edge_metric      <- "edges"

  types_selected <- c(
    general  = any(statistics %in% general_metrics),
    bridge   = any(statistics %in% bridge_metrics),
    excluded = any(statistics %in% excluded_metrics),
    edge     = any(statistics %in% edge_metric)
  )
  if (sum(types_selected) > 1) {
    stop(
      "Cannot combine statistics of different types.\n",
      "- General metrics can only be plotted together.\n",
      "- Bridge metrics can only be plotted together.\n",
      "- Excluded bridge metrics can only be plotted together.\n",
      "- Edge weights can only be plotted alone."
    )
  }

  statistic_map <- c(
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
    edges = "edges"
  )

  .pick_ci_cols <- function(ci_mat, ids, lab_lo, lab_hi) {
    if (is.null(ci_mat)) return(list(lower = rep(NA_real_, length(ids)),
                                     upper = rep(NA_real_, length(ids))))
    ci_mat <- as.matrix(ci_mat)
    if (!is.null(rownames(ci_mat))) {
      ci_mat <- ci_mat[match(ids, rownames(ci_mat)), , drop = FALSE]
    }
    rownames(ci_mat) <- NULL
    cn <- colnames(ci_mat)
    # 1) standard "x%" columns
    if (!is.null(cn) && all(c(lab_lo, lab_hi) %in% cn)) {
      return(list(lower = as.numeric(ci_mat[, lab_lo]),
                  upper = as.numeric(ci_mat[, lab_hi])))
    }
    # 2) allow "lower"/"upper"
    if (!is.null(cn) && all(c("lower", "upper") %in% cn)) {
      return(list(lower = as.numeric(ci_mat[, "lower"]),
                  upper = as.numeric(ci_mat[, "upper"])))
    }

    if (ncol(ci_mat) >= 2) {
      return(list(lower = as.numeric(ci_mat[, 1]),
                  upper = as.numeric(ci_mat[, 2])))
    }

    list(lower = rep(NA_real_, length(ids)),
         upper = rep(NA_real_, length(ids)))
  }

  .warn_outside_ci <- function(df, statistic) {
    if (is.null(df) || !nrow(df)) return(invisible(NULL))
    if (!all(c("observed", "lower", "upper") %in% names(df))) return(invisible(NULL))

    # decide id column based on statistic + columns present
    id_col <- NULL

    # if we're plotting edges, prefer "edge", otherwise "node"
    if (identical(statistic, "edges")) {
      if ("edge" %in% names(df)) {
        id_col <- "edge"
      } else if ("node" %in% names(df)) {
        id_col <- "node"   # fallback (your current edges df uses node=edge labels)
      } else {
        return(invisible(NULL))
      }
    } else {
      if ("node" %in% names(df)) {
        id_col <- "node"
      } else if ("edge" %in% names(df)) {
        id_col <- "edge"   # fallback if someone passes an edge-like df
      } else {
        return(invisible(NULL))
      }
    }

    ok_obs <- is.finite(df$observed)
    ok_ci  <- is.finite(df$lower) & is.finite(df$upper)
    outside_ci <- ok_obs & ok_ci & (df$observed < df$lower | df$observed > df$upper)

    if (!any(outside_ci, na.rm = TRUE)) return(invisible(NULL))

    bad_ids <- df[[id_col]][which(outside_ci)]
    n_bad <- length(bad_ids)

    show_max <- 15
    shown <- if (n_bad > show_max) bad_ids[1:show_max] else bad_ids
    tail_txt <- if (n_bad > show_max) paste0(" ... (+", n_bad - show_max, " more)") else ""

    obj_txt <- if (identical(statistic, "edges")) "edges" else "nodes"
    lab_txt <- if (identical(statistic, "edges")) "Edges" else "Nodes"

    warning(
      sprintf(
        "plot(): '%s' : %d %s have observed value outside the bootstrap CI. %s: %s%s",
        statistic, n_bad, obj_txt, lab_txt,
        paste(shown, collapse = ", "), tail_txt
      ),
      call. = FALSE
    )

    invisible(NULL)
  }

  df_all <- list()

  for (statistic in statistics) {
    df_list <- list()

    if (statistic == "edges") {
      edges_true <- fit$statistics$edge$true
      ci <- fit$statistics$edge$ci
      if (is.null(edges_true) || nrow(edges_true) == 0) next

      edges_true <- edges_true[edges_true$weight != 0, , drop = FALSE]

      edge_set <- edges_true$edge
      ci_pair <- if (is.null(ci)) {
        list(lower = rep(NA_real_, length(edge_set)),
             upper = rep(NA_real_, length(edge_set)))
      } else {
        .pick_ci_cols(ci, edge_set, lab_lo, lab_hi)
      }

      df <- data.frame(
        node = edge_set,
        observed = edges_true$weight,
        lower = ci_pair$lower,
        upper = ci_pair$upper,
        statistic = "edges",
        community = NA
      )

      .warn_outside_ci(df, statistic = "edges")

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

      df_list[[statistic]] <- df
      df_all[[statistic]]  <- do.call(rbind, df_list)
      next
    }

    # ------- Node metrics -------
    ct <- fit$statistics$node$true
    internal <- statistic_map[[statistic]]

    node_set <- ct$node
    observed_values <- ct[[internal]]

    ci <- fit$statistics$node$ci[[statistic]]
    ci_pair <- .pick_ci_cols(ci, node_set, lab_lo, lab_hi)
    lower_vec <- ci_pair$lower
    upper_vec <- ci_pair$upper

    ## --- observed and community ---
    ct_filtered <- fit$statistics$node$true %>%
      dplyr::filter(node %in% node_set)
    order_index <- match(node_set, ct_filtered$node)
    observed_values <- ct_filtered[[statistic_map[statistic]]][order_index]

    community <- if (!is.null(fit$communities$groups)) {
      fit$communities$groups[node_set]
    } else {
      rep(NA_integer_, length(node_set))
    }

    names(node_set)        <- NULL
    names(observed_values) <- NULL
    names(lower_vec)       <- NULL
    names(upper_vec)       <- NULL
    names(community)       <- NULL

    df <- data.frame(
      node      = node_set,
      observed  = observed_values,
      lower     = lower_vec,
      upper     = upper_vec,
      statistic    = statistic,
      community = community,
      row.names = seq_along(node_set),
      check.names = FALSE
    )

    .warn_outside_ci(df, statistic = statistic)

    # ---- WARN: CI available but observed missing ----
    miss_obs_with_ci <- is.na(df$observed) & (is.finite(df$lower) | is.finite(df$upper))
    if (any(miss_obs_with_ci, na.rm = TRUE)) {
      bad_nodes <- df$node[which(miss_obs_with_ci)]
      n_bad <- length(bad_nodes)

      # abbreviate long lists
      show_max <- 15
      shown <- if (n_bad > show_max) bad_nodes[1:show_max] else bad_nodes
      tail_txt <- if (n_bad > show_max) paste0(" ... (+", n_bad - show_max, " more)") else ""

      warning(
        sprintf(
          "plot(): '%s' : %d nodes have observed closeness = NA in the original network (isolated nodes). Bootstrap CIs may still exist because the node is connected in some resamples. Nodes: %s%s",
          statistic, n_bad, paste(shown, collapse = ", "), tail_txt
        ),
        call. = FALSE
      )
    }

    if (isTRUE(standardize)) {
      mean_val <- mean(df$observed, na.rm = TRUE)
      sd_val   <- stats::sd(df$observed, na.rm = TRUE)
      df$observed <- (df$observed - mean_val) / sd_val
      df$lower    <- (df$lower    - mean_val) / sd_val
      df$upper    <- (df$upper    - mean_val) / sd_val
    }

    # keep/exclude coherent with statistic type
    if (grepl("_excluded$", statistic)) {
      df <- df %>% dplyr::filter(is.na(community))
    } else if (grepl("^bridge_", statistic)) {
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
        community_clean = ifelse(statistic == "edges",
                                 NA_character_,
                                 ifelse(is.na(community), "Excluded", as.character(community)))
      )

    community_levels_present <- unique(df_plot$community_clean)
    community_numeric <- sort(community_levels_present[community_levels_present != "Excluded"])
    community_ordered <- c(community_numeric, "Excluded")

    cmap <- c(fit$communities$palette, Excluded = "gray70")
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

    df_list[[statistic]] <- df_plot
    df_all[[statistic]]  <- do.call(rbind, df_list)
  }

  # Reference order from the first statistic
  reference_statistic <- statistics[1]
  reference_df <- df_all[[reference_statistic]]

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
  df_all_combined$statistic <- factor(df_all_combined$statistic, levels = statistics)

  # Ensure label_colored exists (belt & suspenders)
  if (!"label_colored" %in% names(df_all_combined)) {
    df_all_combined$label_colored <- paste0("<b>", df_all_combined$node, "</b>")
  }

  # Palette including 'Excluded'
  cmap <- c(fit$communities$palette, Excluded = "gray70")

  # Label map
  lab_map <- df_all_combined %>%
    dplyr::select(order_reversed, label_colored) %>%
    dplyr::distinct()

  # Build plot
  p <- ggplot2::ggplot(df_all_combined, ggplot2::aes(x = order_reversed, y = observed)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.2,
      color = "black",
      na.rm = TRUE
    )

  plotting_edges <- all(df_all_combined$statistic == "edges")

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
      subtitle = NULL,
      x = x_lab,
      y = if (standardize)
        paste0("Z-score of estimated value with ", ci_txt)
      else
        paste0("Estimated value with ", ci_txt)
    ) +
    ggplot2::facet_wrap(~ statistic, ncol = length(statistics), scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text   = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x  = ggplot2::element_text(size = 9, face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 12, face = "bold"),
      plot.title   = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = if (isTRUE(color_by_community) && !plotting_edges) "right" else "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  return(p)

}
