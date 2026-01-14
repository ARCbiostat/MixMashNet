#' Summarize MixMashNet fits (single- and multilayer) in long format
#'
#' @param object An object of class \code{"mixmashnet"} returned by
#'   \code{mixMN()} or \code{multimixMN()}.
#' @param what Character string indicating which part of the model to summarize:
#'   \itemize{
#'     \item \code{"intra"}: intra-layer quantities (node-level indices and/or
#'       intra-layer edges);
#'     \item \code{"inter"}: interlayer quantities (node-level indices on the
#'       interlayer-only graph and/or cross-layer edges; multilayer fits only).
#'   }
#' @param statistics Character vector specifying which statistics to include.
#'   For \code{what = "intra"}, valid values are:
#'   \code{c("edges",
#'          "strength", "expected_influence", "closeness", "betweenness",
#'          "bridge_strength", "bridge_closeness", "bridge_betweenness",
#'          "bridge_ei1", "bridge_ei2",
#'          "bridge_strength_excluded", "bridge_betweenness_excluded",
#'          "bridge_closeness_excluded", "bridge_ei1_excluded",
#'          "bridge_ei2_excluded")}.
#'
#'   For \code{what = "inter"}, valid values are:
#'   \code{c("edges", "strength", "expected_influence", "closeness",
#'          "betweenness")}.
#'
#'   If \code{statistics = NULL}, then:
#'   \itemize{
#'     \item for \code{what = "intra"}, all available intra-layer statistics
#'       (including \code{"edges"}) are returned;
#'     \item for \code{what = "inter"}, all available interlayer statistics
#'       (including \code{"edges"}) are returned.
#'   }
#' @param layer Optional character vector of layer names to subset. Used for
#'   \code{what = "intra"} in multilayer fits. Ignored for single-layer fits.
#' @param pairs Optional character vector of layer-pair names (e.g.
#'   \code{"bio_dis"}) used for \code{what = "inter"} when summarizing
#'   interlayer edges. If \code{NULL}, all available layer pairs are included.
#'   For interlayer node indices, \code{pairs} can be used to restrict the
#'   summary to nodes belonging to one of the layers in the given pair.
#' @param digits Number of digits to round numeric summaries.
#' @param ... Not used (for S3 compatibility).
#'
#' @return A list (class \code{"summary.mixmashnet"}) with up to four data
#'   frames:
#'   \itemize{
#'     \item \code{$index}: intra-layer node-level indices (one row per
#'       node-metric);
#'     \item \code{$edges}: intra-layer edges (one row per edge);
#'     \item \code{$interlayer_index}: interlayer-only node indices;
#'     \item \code{$interlayer_edges}: cross-layer edges.
#'   }
#'
#'   Depending on \code{what} and \code{statistics}, some of these elements may
#'   be \code{NULL}.
#'
#' @method summary mixmashnet
#' @export
summary.mixmashnet <- function(object,
                               what       = c("intra", "inter"),
                               statistics = NULL,
                               layer      = NULL,
                               pairs      = NULL,
                               digits     = 3,
                               ...) {

  if (!inherits(object, "mixmashnet")) {
    stop("`summary.mixmashnet()` expects an object of class 'mixmashnet'.")
  }

  is_multi <- inherits(object, "multimixMN_fit")

  .get_conf_level <- function(fit, default = 0.95) {
    cl <- fit$settings$conf_level
    if (is.null(cl) || !is.numeric(cl) || length(cl) != 1L ||
        is.na(cl) || cl <= 0 || cl >= 1) cl <- default
    cl
  }

  # quantile probs from conf_level
  .ci_probs <- function(conf_level) {
    alpha <- (1 - conf_level) / 2
    c(alpha, 1 - alpha)
  }

  if (missing(what)) {
    if (!is_multi) {
      what <- "intra"

    } else {
      # multilayer
      if (!is.null(statistics)) {
        if (!is.null(layer)) {
          what <- "intra"
        } else {
          stop(
            "You are summarizing statistics on a multilayer object.\n",
            "Please specify one of:\n",
            "  - layer = \"...\" and what = \"intra\"   # intra-layer statistics for a specific layer\n",
            "  - what  = \"intra\"                     # intra-layer statistics for ALL layers\n",
            "  - what  = \"inter\"                     # interlayer statistics\n"
          )
        }
      } else {
        what <- "intra"
      }
    }
  } else {
    what <- match.arg(what)
  }

  if (what == "inter" && !is_multi) {
    stop("`what = \"inter\"` is only available for 'multimixMN_fit' objects.")
  }

  # ----- allowed statistics -----
  allowed_intra <- c(
    "edges",
    "strength", "expected_influence", "closeness", "betweenness",
    "bridge_strength", "bridge_closeness", "bridge_betweenness",
    "bridge_ei1", "bridge_ei2",
    "bridge_strength_excluded", "bridge_betweenness_excluded",
    "bridge_closeness_excluded", "bridge_ei1_excluded",
    "bridge_ei2_excluded"
  )

  allowed_inter <- c(
    "edges",
    "strength", "expected_influence", "closeness", "betweenness"
  )

  if (is.null(statistics)) {
    if (what == "intra") {
      statistics <- allowed_intra        # tutte le intra + edges
    } else {
      statistics <- allowed_inter        # tutte le inter + edges
    }
  }

  if (what == "intra") {
    invalid <- setdiff(statistics, allowed_intra)
    if (length(invalid)) {
      stop(
        "Invalid `statistics` for what = \"intra\": ",
        paste(invalid, collapse = ", "), "\n",
        "Valid values are: ",
        paste(allowed_intra, collapse = ", ")
      )
    }
    want_edges     <- "edges" %in% statistics
    stats_intra_nd <- setdiff(statistics, "edges")

  } else { # what == "inter"
    invalid <- setdiff(statistics, allowed_inter)
    if (length(invalid)) {
      stop(
        "Invalid `statistics` for what = \"inter\": ",
        paste(invalid, collapse = ", "), "\n",
        "Valid values are: ",
        paste(allowed_inter, collapse = ", ")
      )
    }
    want_edges     <- "edges" %in% statistics
    stats_inter_nd <- setdiff(statistics, "edges")
  }

  .normalize_pairs <- function(pairs_vec) {
    if (is.null(pairs_vec) || !length(pairs_vec)) return(character(0L))
    vapply(strsplit(pairs_vec, "_", fixed = TRUE), function(xx) {
      xx <- xx[xx != ""]
      if (!length(xx)) return("")
      paste(sort(xx), collapse = "_")
    }, FUN.VALUE = character(1L))
  }

  # helper: round numeric columns
  .round_numcols <- function(df, digits) {
    if (is.null(df) || !nrow(df)) return(df)
    num_cols <- vapply(df, is.numeric, logical(1L))
    df[, num_cols] <- lapply(df[, num_cols, drop = FALSE], function(x) {
      x[is.finite(x)] <- round(x[is.finite(x)], digits)
      x
    })
    df
  }

  conf_level <- .get_conf_level(object, default = 0.95)
  pr_global  <- .ci_probs(conf_level)

  # ------------------------------------------------------------------
  # 1) INTRA-LAYER NODE-LEVEL INDICES (index)
  # ------------------------------------------------------------------
  index_tab <- NULL

  if (what == "intra" && length(stats_intra_nd)) {

    build_index_long <- function(true_df, boot_list, ci_list,
                                 layer_label, stats_names, digits) {
      if (is.null(true_df) || !nrow(true_df)) return(NULL)
      if (!"node" %in% colnames(true_df)) {
        stop("`statistics$node$true` must contain a column 'node'.")
      }

      nodes <- as.character(true_df$node)

      # mapping statistic name -> column name in true_df / boot_list / ci_list
      map_intra <- c(
        expected_influence = "ei1"
      )

      metrics_block <- stats_names[
        stats_names %in% setdiff(colnames(true_df), "node") |
          stats_names %in% names(map_intra)
      ]
      if (!length(metrics_block)) return(NULL)

      rows <- vector("list", length(metrics_block))

      for (k in seq_along(metrics_block)) {
        met <- metrics_block[k]
        met_col <- if (met %in% names(map_intra)) map_intra[[met]] else met

        estimated <- true_df[[met_col]]

        mean_boot <- sd_boot <- rep(NA_real_, length(nodes))
        lower     <- upper   <- rep(NA_real_, length(nodes))

        # bootstrap (if available)
        if (!is.null(boot_list) &&
            !is.null(boot_list[[met_col]]) &&
            is.matrix(boot_list[[met_col]])) {

          boot_mat <- boot_list[[met_col]]
          # ensure columns are named by nodes
          if (is.null(colnames(boot_mat))) {
            colnames(boot_mat) <- nodes[seq_len(min(ncol(boot_mat), length(nodes)))]
          }

          boot_full <- matrix(
            NA_real_, nrow = nrow(boot_mat), ncol = length(nodes),
            dimnames = list(rownames(boot_mat), nodes)
          )
          common <- intersect(nodes, colnames(boot_mat))
          if (length(common)) {
            boot_full[, common] <- boot_mat[, common, drop = FALSE]
          }

          mean_boot <- colMeans(boot_full, na.rm = TRUE)
          sd_boot   <- apply(boot_full, 2, stats::sd, na.rm = TRUE)
        }

        # CI (if available)
        if (!is.null(ci_list) &&
            !is.null(ci_list[[met]]) &&
            is.matrix(ci_list[[met]]) &&
            ncol(ci_list[[met]]) >= 2) {
          ci_mat <- ci_list[[met]]
          if (is.null(rownames(ci_mat))) {
            k2 <- min(nrow(ci_mat), length(nodes))
            lower[seq_len(k2)] <- ci_mat[seq_len(k2), 1]
            upper[seq_len(k2)] <- ci_mat[seq_len(k2), 2]
          } else {
            common <- intersect(nodes, rownames(ci_mat))
            if (length(common)) {
              ir   <- match(common, rownames(ci_mat))
              inod <- match(common, nodes)
              lower[inod] <- ci_mat[ir, 1]
              upper[inod] <- ci_mat[ir, 2]
            }
          }
        }

        rows[[k]] <- data.frame(
          node               = nodes,
          layer              = layer_label,
          metric             = met,
          estimated          = estimated,
          mean.bootstrap     = mean_boot,
          SE.bootstrap       = sd_boot,
          CI.lower.bootstrap = lower,
          CI.upper.bootstrap = upper,
          stringsAsFactors   = FALSE
        )
      }

      out <- do.call(rbind, rows)
      rownames(out) <- NULL
      .round_numcols(out, digits)
    }

    if (is_multi) {
      all_layers <- names(object$layer_fits)
      if (is.null(layer)) {
        layers_sel <- all_layers
      } else {
        layers_sel <- intersect(all_layers, layer)
      }

      blocks <- vector("list", length(layers_sel))
      for (j in seq_along(layers_sel)) {
        L    <- layers_sel[j]
        fitL <- object$layer_fits[[L]]
        if (is.null(fitL) || is.null(fitL$statistics$node)) next

        true_df   <- fitL$statistics$node$true
        boot_list <- fitL$statistics$node$boot
        ci_list   <- fitL$statistics$node$ci

        blocks[[j]] <- build_index_long(
          true_df     = true_df,
          boot_list   = boot_list,
          ci_list     = ci_list,
          layer_label = L,
          stats_names = stats_intra_nd,
          digits      = digits
        )
      }
      blocks <- Filter(Negate(is.null), blocks)
      if (length(blocks)) index_tab <- do.call(rbind, blocks)

    } else {
      # single-layer: no concept of layer; use "1"
      true_df   <- object$statistics$node$true
      boot_list <- object$statistics$node$boot
      ci_list   <- object$statistics$node$ci

      index_tab <- build_index_long(
        true_df     = true_df,
        boot_list   = boot_list,
        ci_list     = ci_list,
        layer_label = "1",
        stats_names = stats_intra_nd,
        digits      = digits
      )
    }
  }

  # ------------------------------------------------------------------
  # 2) INTRA-LAYER EDGES (edges)
  # ------------------------------------------------------------------
  edges_tab <- NULL

  if (what == "intra" && isTRUE(want_edges)) {

    build_edge_long <- function(true_edges, boot_mat, layer_label, digits) {
      if (is.null(true_edges) || !nrow(true_edges)) return(NULL)
      if (!all(c("edge","weight") %in% colnames(true_edges))) {
        stop("`statistics$edge$true` must contain columns 'edge' and 'weight'.")
      }

      idx_keep <- !is.na(true_edges$weight) & (true_edges$weight != 0)
      if (!any(idx_keep)) return(NULL)

      true_edges <- true_edges[idx_keep, , drop = FALSE]
      edges      <- as.character(true_edges$edge)
      weight_obs <- true_edges$weight

      mean_boot <- sd_boot <- rep(NA_real_, length(edges))
      lower     <- upper   <- rep(NA_real_, length(edges))

      if (!is.null(boot_mat) && is.matrix(boot_mat)) {
        if (is.null(rownames(boot_mat))) {
          rownames(boot_mat) <- edges
        }
        boot_full <- matrix(
          NA_real_, nrow = length(edges), ncol = ncol(boot_mat),
          dimnames = list(edges, colnames(boot_mat))
        )
        common <- intersect(edges, rownames(boot_mat))
        if (length(common)) {
          boot_full[common, ] <- boot_mat[common, , drop = FALSE]
        }

        mean_boot <- rowMeans(boot_full, na.rm = TRUE)
        sd_boot   <- apply(boot_full, 1, stats::sd, na.rm = TRUE)

        pr <- .ci_probs(conf_level)
        lower <- apply(boot_full, 1, stats::quantile, probs = pr[1], na.rm = TRUE, type = 6)
        upper <- apply(boot_full, 1, stats::quantile, probs = pr[2], na.rm = TRUE, type = 6)

      }

      out <- data.frame(
        edge               = edges,
        layer              = layer_label,
        estimated          = weight_obs,
        mean.bootstrap     = mean_boot,
        SE.bootstrap       = sd_boot,
        CI.lower.bootstrap = lower,
        CI.upper.bootstrap = upper,
        stringsAsFactors   = FALSE
      )
      rownames(out) <- NULL
      .round_numcols(out, digits)
    }

    if (is_multi) {
      all_layers <- names(object$layer_fits)
      if (is.null(layer)) {
        layers_sel <- all_layers
      } else {
        layers_sel <- intersect(all_layers, layer)
      }

      blocks <- list()
      for (L in layers_sel) {
        fitL <- object$layer_fits[[L]]
        if (is.null(fitL) || is.null(fitL$statistics$edge)) next
        true_edges <- fitL$statistics$edge$true
        boot_mat   <- fitL$statistics$edge$boot
        blocks[[length(blocks) + 1L]] <- build_edge_long(
          true_edges  = true_edges,
          boot_mat    = boot_mat,
          layer_label = L,
          digits      = digits
        )
      }
      blocks <- Filter(Negate(is.null), blocks)
      if (length(blocks)) edges_tab <- do.call(rbind, blocks)

    } else {
      true_edges <- object$statistics$edge$true
      boot_mat   <- object$statistics$edge$boot
      edges_tab  <- build_edge_long(
        true_edges  = true_edges,
        boot_mat    = boot_mat,
        layer_label = "1",
        digits      = digits
      )
    }
  }

  # ------------------------------------------------------------------
  # 3) INTERLAYER NODE INDEX (interlayer_index)
  # ------------------------------------------------------------------
  inter_idx_tab <- NULL

  if (what == "inter" && is_multi && length(stats_inter_nd) &&
      !is.null(object$interlayer)) {

    inter_cent <- object$interlayer$centrality
    if (!is.null(inter_cent) && !is.null(inter_cent$true)) {

      true_df   <- inter_cent$true
      boot_list <- inter_cent$boot
      ci_list   <- inter_cent$ci_results

      if (!"node" %in% colnames(true_df)) {
        stop("`object$interlayer$centrality$true` must contain a column 'node'.")
      }
      nodes_all <- as.character(true_df$node)

      # mapping: user-level -> internal column names / CI keys
      stat_to_true <- c(
        strength           = "strength",
        expected_influence = "ei1",
        closeness          = "closeness",
        betweenness        = "betweenness"
      )
      stat_to_ci <- c(
        strength           = "strength",
        expected_influence = "ei1",
        closeness          = "closeness",
        betweenness        = "betweenness"
      )

      # subset nodes by pairs (if given)
      nodes_use <- nodes_all

      if (!is.null(pairs) && length(pairs) == 1L &&
          !is.null(object$layers) && !is.null(object$layers$assignment)) {

        parts <- strsplit(pairs, "_", fixed = TRUE)[[1]]
        parts <- parts[parts != ""]
        if (length(parts) == 2L) {
          lay_assign <- object$layers$assignment
          lay_assign <- lay_assign[nodes_all]
          nodes_use  <- nodes_all[lay_assign %in% parts]
        }

      } else if (!is.null(layer) &&
                 !is.null(object$layers) &&
                 !is.null(object$layers$assignment)) {

        lay_assign <- object$layers$assignment
        lay_assign <- lay_assign[nodes_all]
        nodes_use  <- nodes_all[lay_assign %in% layer]
      }

      if (length(nodes_use)) {

        rows <- vector("list", length(stats_inter_nd))

        for (k in seq_along(stats_inter_nd)) {
          stat_name <- stats_inter_nd[k]
          true_col  <- stat_to_true[[stat_name]]
          ci_key    <- stat_to_ci[[stat_name]]

          if (is.na(true_col) || !true_col %in% colnames(true_df)) {
            stop("Interlayer true table does not contain column '", true_col,
                 "' for statistic '", stat_name, "'.")
          }

          estimated <- true_df[match(nodes_use, nodes_all), true_col]

          mean_boot <- sd_boot <- rep(NA_real_, length(nodes_use))
          lower     <- upper   <- rep(NA_real_, length(nodes_use))

          # bootstrap (indexed by true_col)
          if (!is.null(boot_list) &&
              !is.null(boot_list[[true_col]]) &&
              is.matrix(boot_list[[true_col]])) {

            boot_mat <- boot_list[[true_col]]
            if (is.null(colnames(boot_mat))) {
              colnames(boot_mat) <- nodes_all[seq_len(min(ncol(boot_mat),
                                                          length(nodes_all)))]
            }
            boot_full <- matrix(
              NA_real_, nrow = nrow(boot_mat), ncol = length(nodes_use),
              dimnames = list(rownames(boot_mat), nodes_use)
            )
            common <- intersect(nodes_use, colnames(boot_mat))
            if (length(common)) {
              boot_full[, common] <- boot_mat[, common, drop = FALSE]
            }
            mean_boot <- colMeans(boot_full, na.rm = TRUE)
            sd_boot   <- apply(boot_full, 2, stats::sd, na.rm = TRUE)
          }

          # CI (indexed by ci_key)
          if (!is.null(ci_key) &&
              !is.null(ci_list) &&
              !is.null(ci_list[[ci_key]]) &&
              is.matrix(ci_list[[ci_key]]) &&
              ncol(ci_list[[ci_key]]) >= 2) {

            ci_mat <- ci_list[[ci_key]]
            if (is.null(rownames(ci_mat))) {
              k2 <- min(nrow(ci_mat), length(nodes_use))
              lower[seq_len(k2)] <- ci_mat[seq_len(k2), 1]
              upper[seq_len(k2)] <- ci_mat[seq_len(k2), 2]
            } else {
              common <- intersect(nodes_use, rownames(ci_mat))
              if (length(common)) {
                ir   <- match(common, rownames(ci_mat))
                inod <- match(common, nodes_use)
                lower[inod] <- ci_mat[ir, 1]
                upper[inod] <- ci_mat[ir, 2]
              }
            }
          }

          rows[[k]] <- data.frame(
            node               = nodes_use,
            metric             = stat_name,
            estimated          = estimated,
            mean.bootstrap     = mean_boot,
            SE.bootstrap       = sd_boot,
            CI.lower.bootstrap = lower,
            CI.upper.bootstrap = upper,
            stringsAsFactors   = FALSE
          )
        }

        out <- do.call(rbind, rows)
        rownames(out) <- NULL
        inter_idx_tab <- .round_numcols(out, digits)
      }
    }
  }

  # ------------------------------------------------------------------
  # 4) INTERLAYER EDGES (interlayer_edges)
  # ------------------------------------------------------------------
  inter_edges_tab <- NULL

  if (what == "inter" && is_multi &&
      isTRUE(want_edges) &&
      !is.null(object$interlayer)) {

    all_pairs <- setdiff(names(object$interlayer), "centrality")

    if (is.null(pairs)) {
      pairs_sel <- all_pairs
    } else {
      norm_all   <- .normalize_pairs(all_pairs)
      norm_pairs <- .normalize_pairs(pairs)
      pairs_sel  <- all_pairs[norm_all %in% norm_pairs]
    }

    build_inter_edges <- function(edges_obj, pair_label, digits) {
      if (is.null(edges_obj) || is.null(edges_obj$true)) return(NULL)
      true_edges <- edges_obj$true
      boot_mat   <- edges_obj$boot

      if (!all(c("edge","weight") %in% colnames(true_edges))) {
        stop("`interlayer[[pair]]$edges$true` must contain 'edge' and 'weight'.")
      }

      idx_keep <- !is.na(true_edges$weight) & (true_edges$weight != 0)
      if (!any(idx_keep)) return(NULL)

      true_edges <- true_edges[idx_keep, , drop = FALSE]
      edges      <- as.character(true_edges$edge)
      weight_obs <- true_edges$weight

      mean_boot <- sd_boot <- rep(NA_real_, length(edges))
      lower     <- upper   <- rep(NA_real_, length(edges))

      if (!is.null(boot_mat) && is.matrix(boot_mat)) {
        if (is.null(rownames(boot_mat))) {
          rownames(boot_mat) <- edges
        }
        boot_full <- matrix(
          NA_real_, nrow = length(edges), ncol = ncol(boot_mat),
          dimnames = list(edges, colnames(boot_mat))
        )
        common <- intersect(edges, rownames(boot_mat))
        if (length(common)) {
          boot_full[common, ] <- boot_mat[common, , drop = FALSE]
        }

        mean_boot <- rowMeans(boot_full, na.rm = TRUE)
        sd_boot   <- apply(boot_full, 1, stats::sd, na.rm = TRUE)

        lower <- apply(boot_full, 1, stats::quantile, probs = pr_global[1], na.rm = TRUE, type = 6)
        upper <- apply(boot_full, 1, stats::quantile, probs = pr_global[2], na.rm = TRUE, type = 6)
      }

      out <- data.frame(
        edge               = edges,
        pairs              = pair_label,
        estimated          = weight_obs,
        mean.bootstrap     = mean_boot,
        SE.bootstrap       = sd_boot,
        CI.lower.bootstrap = lower,
        CI.upper.bootstrap = upper,
        stringsAsFactors   = FALSE
      )
      rownames(out) <- NULL
      .round_numcols(out, digits)
    }

    blocks <- list()
    for (nm in pairs_sel) {
      edges_obj <- object$interlayer[[nm]]$edges
      blocks[[length(blocks) + 1L]] <- build_inter_edges(
        edges_obj  = edges_obj,
        pair_label = nm,
        digits     = digits
      )
    }
    blocks <- Filter(Negate(is.null), blocks)
    if (length(blocks)) inter_edges_tab <- do.call(rbind, blocks)
  }

  out <- list(
    index            = index_tab,
    edges            = edges_tab,
    interlayer_index = inter_idx_tab,
    interlayer_edges = inter_edges_tab,
    conf_level       = if (is_multi) .get_conf_level(object) else .get_conf_level(object)
  )
  class(out) <- c("summary.mixmashnet", "list")
  out
}


#' Print method for "summary.mixmashnet" objects
#'
#' @param x A "summary.mixmashnet" object
#' @param digits Number of digits to print
#' @param top_n Show only the top `top_n` rows per block (ranked by |estimated|).
#'   Use `Inf` to show all rows (default).
#' @param ... Unused
#'
#' @method print summary.mixmashnet
#' @export
print.summary.mixmashnet <- function(x, digits = 3, top_n = Inf, ...) {

  # rounding helper
  round_df <- function(df, digits) {
    if (is.null(df) || !nrow(df)) return(df)
    num_cols <- vapply(df, is.numeric, logical(1L))
    if (any(num_cols)) {
      df[, num_cols] <- lapply(df[, num_cols, drop = FALSE], function(z) {
        idx <- is.finite(z)
        z[idx] <- round(z[idx], digits)
        z
      })
    }
    df
  }

  # helper for robust sorting by one or more columns
  order_by <- function(df, cols) {
    cols <- intersect(cols, colnames(df))
    if (!length(cols) || !nrow(df)) return(df)
    o <- do.call(order, df[cols])
    df[o, , drop = FALSE]
  }

  # helper function to make column names more readable for printing
  prettify_colnames <- function(df, conf_level = 0.95) {
    if (is.null(df) || !nrow(df)) return(df)

    ci_pct <- paste0(round(100 * conf_level), "%")

    map <- c(
      "mean.bootstrap"      = "mean (bootstrap)",
      "SE.bootstrap"        = "SE (bootstrap)",
      "CI.lower.bootstrap"  = paste0(ci_pct, " CI lower (bootstrap)"),
      "CI.upper.bootstrap"  = paste0(ci_pct, " CI upper (bootstrap)")
    )

    cn <- colnames(df)
    cn <- ifelse(cn %in% names(map), map[cn], cn)
    cn <- gsub("\\.", " ", cn)
    colnames(df) <- cn
    df
  }

  # ---------- 1) INDEX (intra-layer, node-level) ----------
  if (!is.null(x$index) && nrow(x$index)) {
    cat("\nNode-level indices (intra-layer):\n")

    idx <- round_df(x$index, digits)
    metrics_idx <- unique(idx$metric)
    original_n  <- length(unique(idx$node))

    for (met in metrics_idx) {
      sub_met <- idx[idx$metric == met, , drop = FALSE]
      cat("\n  Metric:", met, "\n")

      if (grepl("^bridge_", met) && "estimated" %in% names(sub_met)) {
        sub_met <- sub_met[!is.na(sub_met$estimated), , drop = FALSE]
      }

      original_n_met <- length(unique(sub_met$node))

      if (!nrow(sub_met)) {
        cat("    (no non-NA values to display)\n")
        next
      }

      sub_met <- order_by(sub_met, c("node"))
      sub_met <- sub_met[, setdiff(colnames(sub_met), "metric"), drop = FALSE]
      sub_met <- prettify_colnames(sub_met, conf_level = x$conf_level %||% 0.95)

      # --- TOP_N FILTER ---
      if (is.finite(top_n) && "estimated" %in% colnames(sub_met)) {
        o <- order(abs(sub_met$estimated), decreasing = TRUE)
        sub_met <- sub_met[o, , drop = FALSE]
        sub_met <- utils::head(sub_met, top_n)
        message(
          "    Showing top ", nrow(sub_met), " of ",
          original_n_met, " nodes (ranked by |estimated|)"
        )
      }
      print(sub_met, row.names = FALSE)
    }
  }

  # ---------- 2) EDGES (intra-layer) ----------
  if (!is.null(x$edges) && nrow(x$edges)) {
    cat("\nIntra-layer edges:\n")

    ed <- round_df(x$edges, digits)

    ed <- order_by(ed, c("layer", "edge"))

    ed <- prettify_colnames(ed, conf_level = x$conf_level %||% 0.95)

    cat("\n")
    # --- TOP_N FILTER for intra edges ---
    if (is.finite(top_n) && "estimated" %in% colnames(ed)) {
      o <- order(abs(ed$estimated), decreasing = TRUE)
      ed <- ed[o, , drop = FALSE]
      ed <- utils::head(ed, top_n)
      message("    Showing top ", nrow(ed), " edges (ranked by |estimated|)")
    }

    print(ed, row.names = FALSE)
  }

  # ---------- 3) INTERLAYER INDEX ----------
  if (!is.null(x$interlayer_index) && nrow(x$interlayer_index)) {
    cat("\nInterlayer-only node indices:\n")

    idx2 <- round_df(x$interlayer_index, digits)
    metrics_int <- unique(idx2$metric)
    original_n2 <- length(unique(idx2$node))

    for (met in metrics_int) {
      sub_met <- idx2[idx2$metric == met, , drop = FALSE]
      cat("\n  Metric:", met, "\n")

      sub_met <- order_by(sub_met, c("layer", "node"))

      sub_met <- sub_met[, setdiff(colnames(sub_met), "metric"), drop = FALSE]

      sub_met <- prettify_colnames(sub_met, conf_level = x$conf_level %||% 0.95)

      # --- TOP_N FILTER for interlayer node metrics ---
      if (is.finite(top_n) && "estimated" %in% colnames(sub_met)) {
        o <- order(abs(sub_met$estimated), decreasing = TRUE)
        sub_met <- sub_met[o, , drop = FALSE]
        sub_met <- utils::head(sub_met, top_n)
        message("    Showing top ", nrow(sub_met), " of ", original_n2,
                " nodes (ranked by |estimated|)")
      }

      print(sub_met, row.names = FALSE)
    }
  }

  # ---------- 4) INTERLAYER EDGES ----------
  if (!is.null(x$interlayer_edges) && nrow(x$interlayer_edges)) {
    cat("\nInterlayer edges:\n")

    ed2 <- round_df(x$interlayer_edges, digits)

    if ("pairs" %in% colnames(ed2)) {
      pairs_vals <- unique(ed2$pairs)
      for (pp in pairs_vals) {
        sub_pp <- ed2[ed2$pairs == pp, , drop = FALSE]
        sub_pp <- order_by(sub_pp, c("pairs", "edge"))
        sub_pp <- prettify_colnames(sub_pp, conf_level = x$conf_level %||% 0.95)
        cat("\n  Pair:", pp, "\n")

        # --- TOP_N FILTER for interlayer edges (for pair) ---
        if (is.finite(top_n) && "estimated" %in% colnames(sub_pp)) {
          o <- order(abs(sub_pp$estimated), decreasing = TRUE)
          sub_pp <- sub_pp[o, , drop = FALSE]
          sub_pp <- utils::head(sub_pp, top_n)
          message("    Showing top ", nrow(sub_pp),
                  " edges (ranked by |estimated|)")
        }

        print(sub_pp, row.names = FALSE)
      }
    } else {
      ed2 <- order_by(ed2, c("edge"))
      ed2 <- prettify_colnames(ed2, conf_level = x$conf_level %||% 0.95)
      cat("\n")

      # --- TOP_N FILTER for interlayer edges (together) ---
      if (is.finite(top_n) && "estimated" %in% colnames(ed2)) {
        o <- order(abs(ed2$estimated), decreasing = TRUE)
        ed2 <- ed2[o, , drop = FALSE]
        ed2 <- utils::head(ed2, top_n)
        message("    Showing top ", nrow(ed2),
                " edges (ranked by |estimated|)")
      }

      print(ed2, row.names = FALSE)
    }
  }

  invisible(x)
}
