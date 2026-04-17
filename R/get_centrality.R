# -------------------------------------------------------------------------
# internal helpers
# -------------------------------------------------------------------------

#' Internal helper
#' @keywords internal
#' @noRd
.get_quantile_level <- function(fit, default = 0.95) {
  cl <- fit$settings$quantile_level
  if (is.null(cl) || !is.numeric(cl) || length(cl) != 1L ||
      is.na(cl) || cl <= 0 || cl >= 1) {
    cl <- default
  }
  cl
}

#' Internal helper
#' @keywords internal
#' @noRd
.round_numcols <- function(df, digits) {
  if (is.null(df) || !nrow(df) || is.null(digits)) return(df)
  num_cols <- vapply(df, is.numeric, logical(1L))
  df[, num_cols] <- lapply(df[, num_cols, drop = FALSE], function(x) {
    x[is.finite(x)] <- round(x[is.finite(x)], digits)
    x
  })
  df
}

#' Internal helper
#' @keywords internal
#' @noRd
.normalize_pairs <- function(pairs_vec) {
  if (is.null(pairs_vec) || !length(pairs_vec)) return(character(0L))
  vapply(strsplit(pairs_vec, "_", fixed = TRUE), function(xx) {
    xx <- xx[xx != ""]
    if (!length(xx)) return("")
    paste(sort(xx), collapse = "_")
  }, FUN.VALUE = character(1L))
}

#' Internal helper
#' @keywords internal
#' @noRd
.drop_all_na_boot_cols <- function(df) {
  boot_cols <- c(
    "mean.bootstrap",
    "SE.bootstrap",
    "quantile.lower.bootstrap",
    "quantile.upper.bootstrap"
  )
  keep <- setdiff(names(df), boot_cols)
  for (nm in boot_cols) {
    if (nm %in% names(df) && !all(is.na(df[[nm]]))) {
      keep <- c(keep, nm)
    }
  }
  df[, keep, drop = FALSE]
}

#' Internal helper
#' @keywords internal
#' @noRd
.build_intra_index_long <- function(true_df, boot_list, quantile_region_list,
                                    layer_label, stats_names) {
  if (is.null(true_df) || !nrow(true_df)) return(NULL)
  if (!"node" %in% colnames(true_df)) {
    stop("`statistics$node$true` must contain a column 'node'.")
  }

  nodes <- as.character(true_df$node)

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

    if (!met_col %in% colnames(true_df)) next

    estimated <- true_df[[met_col]]

    if (grepl("_excluded$", met) && all(is.na(estimated))) {
      next
    }

    mean_boot <- SE_boot <- rep(NA_real_, length(nodes))
    q_lower   <- q_upper <- rep(NA_real_, length(nodes))

    if (!is.null(boot_list) &&
        !is.null(boot_list[[met_col]]) &&
        is.matrix(boot_list[[met_col]])) {

      boot_mat <- boot_list[[met_col]]

      if (is.null(colnames(boot_mat))) {
        colnames(boot_mat) <- nodes[seq_len(min(ncol(boot_mat), length(nodes)))]
      }

      boot_full <- matrix(
        NA_real_,
        nrow = nrow(boot_mat),
        ncol = length(nodes),
        dimnames = list(rownames(boot_mat), nodes)
      )

      common <- intersect(nodes, colnames(boot_mat))
      if (length(common)) {
        boot_full[, common] <- boot_mat[, common, drop = FALSE]
      }

      mean_boot <- colMeans(boot_full, na.rm = TRUE)
      SE_boot   <- apply(boot_full, 2, stats::sd, na.rm = TRUE)
    }

    qr_key <- if (met %in% names(map_intra)) met else met

    if (!is.null(quantile_region_list) &&
        !is.null(quantile_region_list[[qr_key]]) &&
        is.matrix(quantile_region_list[[qr_key]]) &&
        ncol(quantile_region_list[[qr_key]]) >= 2) {

      qr_mat <- quantile_region_list[[qr_key]]

      if (is.null(rownames(qr_mat))) {
        k2 <- min(nrow(qr_mat), length(nodes))
        q_lower[seq_len(k2)] <- qr_mat[seq_len(k2), 1]
        q_upper[seq_len(k2)] <- qr_mat[seq_len(k2), 2]
      } else {
        common <- intersect(nodes, rownames(qr_mat))
        if (length(common)) {
          ir   <- match(common, rownames(qr_mat))
          inod <- match(common, nodes)
          q_lower[inod] <- qr_mat[ir, 1]
          q_upper[inod] <- qr_mat[ir, 2]
        }
      }
    }

    keep_rows <- rep(TRUE, length(nodes))
    if (grepl("^bridge_", met)) {
      keep_rows <- !is.na(estimated)
    }

    if (!any(keep_rows)) next

    rows[[k]] <- data.frame(
      node                     = nodes[keep_rows],
      layer                    = rep(layer_label, sum(keep_rows)),
      scope                    = rep("intra", sum(keep_rows)),
      metric                   = rep(met, sum(keep_rows)),
      estimated                = estimated[keep_rows],
      mean.bootstrap           = mean_boot[keep_rows],
      SE.bootstrap             = SE_boot[keep_rows],
      quantile.lower.bootstrap = q_lower[keep_rows],
      quantile.upper.bootstrap = q_upper[keep_rows],
      stringsAsFactors = FALSE
    )
  }

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Internal helper
#' @keywords internal
#' @noRd
.finalize_centrality <- function(out, quantile_level, what, digits,
                                         drop_na_boot, metric_order) {
  if (is.null(out) || !nrow(out)) {
    out <- tibble::tibble(
      node = character(),
      layer = character(),
      scope = character(),
      metric = character(),
      estimated = numeric()
    )
    attr(out, "quantile_level") <- quantile_level
    attr(out, "what") <- what
    class(out) <- c("get_centrality", class(out))
    return(out)
  }

  out <- .round_numcols(out, digits)

  if (isTRUE(drop_na_boot)) {
    out <- .drop_all_na_boot_cols(out)
  }

  out$metric <- factor(out$metric, levels = metric_order)
  out <- out[order(out$scope, out$layer, out$metric, out$node), , drop = FALSE]
  out <- tibble::as_tibble(out)
  out$metric <- as.character(out$metric)

  attr(out, "quantile_level") <- quantile_level
  attr(out, "what") <- what
  class(out) <- c("get_centrality", class(out))

  out
}

#' Extract node-level centrality indices
#'
#' @description
#' Extracts node-level centrality indices from fitted objects returned by
#' \code{mixMN()} and \code{multimixMN()} in a long-format data frame.
#'
#' For single layer fits of class \code{"mixMN_fit"}, only intralayer
#' node-level indices are available.
#'
#' For multilayer fits of class \code{"multimixMN_fit"}, \code{what = "intra"}
#' returns intralayer node-level indices, whereas \code{what = "inter"}
#' returns node-level indices computed on the interlayer-only graph.
#' If \code{what} is not specified for a multilayer fit, both intralayer and
#' interlayer node-level indices are returned by default, unless
#' \code{layer} or \code{pairs} imply a specific scope.
#'
#' The function returns the original estimates and, when available, bootstrap
#' means, standard errors, and bootstrap quantile regions.
#'
#' @aliases get_centrality get_centrality.mixMN_fit get_centrality.multimixMN_fit
#' @param object A fitted object of class \code{"mixMN_fit"} or
#'   \code{"multimixMN_fit"}.
#' @param what Character string indicating which node-level indices to extract:
#'   \itemize{
#'     \item \code{"intra"}: intralayer node-level indices;
#'     \item \code{"inter"}: interlayer-only node-level indices
#'       (available only for \code{"multimixMN_fit"} objects).
#'   }
#'   For single layer fits, only \code{"intra"} is allowed.
#'   For multilayer fits, if omitted, both scopes are returned by default
#'   unless \code{layer} or \code{pairs} imply a specific scope.
#' @param statistics Character vector specifying which node-level statistics
#'   to include.
#' @param layer Optional character vector of layer names to subset.
#'   Relevant for intralayer output in multilayer fits.
#' @param pairs Optional character vector of layer-pair names to subset.
#'   Relevant for interlayer output in multilayer fits.
#' @param digits Optional number of digits used to round numeric columns.
#' @param drop_na_boot Logical. If \code{TRUE} (default), bootstrap-related
#'   columns that are entirely \code{NA} are removed from the output.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' The returned data frame is in long format, with one row per
#' node-statistic combination.
#'
#' For single layer fits, only \code{what = "intra"} is available.
#'
#' For multilayer fits, \code{layer} can be used to subset intralayer output,
#' whereas \code{pairs} can be used to subset interlayer output.
#'
#' The set of admissible statistics depends on \code{what} and on the class of
#' \code{object}. In particular, bridge-related indices are available only for
#' intralayer output.
#'
#' @return
#' A tibble in long format with one row per node-statistic combination.
#' It contains the columns:
#' \itemize{
#'   \item \code{node}
#'   \item \code{layer}
#'   \item \code{scope}
#'   \item \code{metric}
#'   \item \code{estimated}
#' }
#'
#' When available, the output also contains bootstrap summary columns:
#' \itemize{
#'   \item \code{mean.bootstrap}
#'   \item \code{SE.bootstrap}
#'   \item \code{quantile.lower.bootstrap}
#'   \item \code{quantile.upper.bootstrap}
#' }
#'
#' The quantile level used to compute the bootstrap quantile region is stored
#' as the \code{"quantile_level"} attribute of the returned tibble.
#'
#' @export
get_centrality <- function(object, ...) {
  UseMethod("get_centrality")
}

#' @rdname get_centrality
#' @export
get_centrality.mixMN_fit <- function(object,
                                         what = "intra",
                                         statistics = NULL,
                                         digits = NULL,
                                         drop_na_boot = TRUE,
                                         ...) {
  allowed_intra <- c(
    "strength", "expected_influence", "closeness", "betweenness",
    "bridge_strength", "bridge_closeness", "bridge_betweenness",
    "bridge_ei1", "bridge_ei2",
    "bridge_strength_excluded", "bridge_betweenness_excluded",
    "bridge_closeness_excluded", "bridge_ei1_excluded",
    "bridge_ei2_excluded"
  )

  what <- match.arg(what, choices = "intra")

  if (is.null(statistics)) {
    statistics <- allowed_intra
  }

  invalid <- setdiff(statistics, allowed_intra)
  if (length(invalid)) {
    stop(
      "Invalid `statistics` for `mixMN_fit`: ",
      paste(invalid, collapse = ", "), "\n",
      "Valid values are: ",
      paste(allowed_intra, collapse = ", ")
    )
  }

  quantile_level <- .get_quantile_level(object)

  out <- .build_intra_index_long(
    true_df = object$statistics$node$true,
    boot_list = object$statistics$node$boot,
    quantile_region_list = object$statistics$node$quantile_region,
    layer_label = "1",
    stats_names = statistics
  )

  .finalize_centrality(
    out = out,
    quantile_level = quantile_level,
    what = what,
    digits = digits,
    drop_na_boot = drop_na_boot,
    metric_order = allowed_intra
  )
}

#' @rdname get_centrality
#' @export
get_centrality.multimixMN_fit <- function(object,
                                              what = c("intra", "inter"),
                                              statistics = NULL,
                                              layer = NULL,
                                              pairs = NULL,
                                              digits = NULL,
                                              drop_na_boot = TRUE,
                                              ...) {
  allowed_intra <- c(
    "strength", "expected_influence", "closeness", "betweenness",
    "bridge_strength", "bridge_closeness", "bridge_betweenness",
    "bridge_ei1", "bridge_ei2",
    "bridge_strength_excluded", "bridge_betweenness_excluded",
    "bridge_closeness_excluded", "bridge_ei1_excluded",
    "bridge_ei2_excluded"
  )

  allowed_inter <- c(
    "strength", "expected_influence", "closeness", "betweenness"
  )

  if (missing(what)) {
    if (!is.null(layer) && is.null(pairs)) {
      what <- "intra"
    } else if (is.null(layer) && !is.null(pairs)) {
      what <- "inter"
    } else if (is.null(layer) && is.null(pairs)) {
      what <- c("intra", "inter")
    } else {
      stop(
        "Ambiguous request for a multilayer object.\n",
        "Please specify `what = \"intra\"` or `what = \"inter\"`.\n",
        "Without `what`, `layer` suggests intralayer output and `pairs` suggests interlayer output."
      )
    }
  } else {
    what <- match.arg(what, choices = c("intra", "inter"), several.ok = TRUE)
  }

  if (is.null(statistics)) {
    if (identical(what, "intra")) {
      statistics <- allowed_intra
    } else if (identical(what, "inter")) {
      statistics <- allowed_inter
    } else {
      statistics <- union(allowed_intra, allowed_inter)
    }
  }

  if (identical(what, "intra")) {
    invalid <- setdiff(statistics, allowed_intra)
    if (length(invalid)) {
      stop(
        "Invalid `statistics` for `what = \"intra\"`: ",
        paste(invalid, collapse = ", "), "\n",
        "Valid values are: ",
        paste(allowed_intra, collapse = ", ")
      )
    }
  } else if (identical(what, "inter")) {
    invalid <- setdiff(statistics, allowed_inter)
    if (length(invalid)) {
      stop(
        "Invalid `statistics` for `what = \"inter\"`: ",
        paste(invalid, collapse = ", "), "\n",
        "Valid values are: ",
        paste(allowed_inter, collapse = ", ")
      )
    }
  } else {
    invalid <- setdiff(statistics, union(allowed_intra, allowed_inter))
    if (length(invalid)) {
      stop(
        "Invalid `statistics`: ",
        paste(invalid, collapse = ", "), "\n",
        "Valid values are: ",
        paste(union(allowed_intra, allowed_inter), collapse = ", ")
      )
    }
  }

  stats_intra <- intersect(statistics, allowed_intra)
  stats_inter <- intersect(statistics, allowed_inter)

  quantile_level <- .get_quantile_level(object)

  out_intra <- NULL
  out_inter <- NULL

  # -----------------------------------------------------------------------
  # intralayer node-level indices
  # -----------------------------------------------------------------------
  if ("intra" %in% what && length(stats_intra)) {
    all_layers <- names(object$layer_fits)
    layers_sel <- if (is.null(layer)) all_layers else intersect(all_layers, layer)

    blocks <- vector("list", length(layers_sel))

    for (j in seq_along(layers_sel)) {
      L <- layers_sel[j]
      fitL <- object$layer_fits[[L]]
      if (is.null(fitL) || is.null(fitL$statistics$node)) next

      blocks[[j]] <- .build_intra_index_long(
        true_df = fitL$statistics$node$true,
        boot_list = fitL$statistics$node$boot,
        quantile_region_list = fitL$statistics$node$quantile_region,
        layer_label = L,
        stats_names = stats_intra
      )
    }

    blocks <- Filter(Negate(is.null), blocks)
    out_intra <- if (length(blocks)) do.call(rbind, blocks) else NULL
  }

  # -----------------------------------------------------------------------
  # interlayer node-level indices
  # -----------------------------------------------------------------------
  if ("inter" %in% what && length(stats_inter)) {
    inter_cent <- object$interlayer$centrality

    if (!is.null(inter_cent) && !is.null(inter_cent$true)) {
      true_df   <- inter_cent$true
      boot_list <- inter_cent$boot
      quantile_region_list <- inter_cent$quantile_region_results

      if (!"node" %in% colnames(true_df)) {
        stop("`object$interlayer$centrality$true` must contain a column 'node'.")
      }

      nodes_all <- as.character(true_df$node)

      stat_to_true <- c(
        strength           = "strength",
        expected_influence = "ei1",
        closeness          = "closeness",
        betweenness        = "betweenness"
      )

      stat_to_qr <- c(
        strength           = "strength",
        expected_influence = "ei1",
        closeness          = "closeness",
        betweenness        = "betweenness"
      )

      node_layer <- rep(NA_character_, length(nodes_all))
      if (!is.null(object$layers) && !is.null(object$layers$assignment)) {
        lay_assign <- object$layers$assignment
        node_layer <- unname(lay_assign[nodes_all])
      }

      nodes_use <- nodes_all

      if (!is.null(pairs) && length(pairs)) {
        norm_pairs <- .normalize_pairs(pairs)
        allowed_layers <- unique(unlist(strsplit(norm_pairs, "_", fixed = TRUE)))
        nodes_use <- nodes_all[node_layer %in% allowed_layers]
      }

      if (!is.null(layer)) {
        nodes_use <- intersect(nodes_use, nodes_all[node_layer %in% layer])
      }

      if (length(nodes_use)) {
        rows <- vector("list", length(stats_inter))

        for (k in seq_along(stats_inter)) {
          stat_name <- stats_inter[k]
          true_col  <- stat_to_true[[stat_name]]
          qr_key    <- stat_to_qr[[stat_name]]

          if (is.na(true_col) || !true_col %in% colnames(true_df)) {
            stop(
              "Interlayer true table does not contain column '", true_col,
              "' for statistic '", stat_name, "'."
            )
          }

          estimated <- true_df[match(nodes_use, nodes_all), true_col]

          mean_boot <- SE_boot <- rep(NA_real_, length(nodes_use))
          q_lower   <- q_upper <- rep(NA_real_, length(nodes_use))

          if (!is.null(boot_list) &&
              !is.null(boot_list[[true_col]]) &&
              is.matrix(boot_list[[true_col]])) {

            boot_mat <- boot_list[[true_col]]

            if (is.null(colnames(boot_mat))) {
              colnames(boot_mat) <- nodes_all[seq_len(min(ncol(boot_mat), length(nodes_all)))]
            }

            boot_full <- matrix(
              NA_real_,
              nrow = nrow(boot_mat),
              ncol = length(nodes_use),
              dimnames = list(rownames(boot_mat), nodes_use)
            )

            common <- intersect(nodes_use, colnames(boot_mat))
            if (length(common)) {
              boot_full[, common] <- boot_mat[, common, drop = FALSE]
            }

            mean_boot <- colMeans(boot_full, na.rm = TRUE)
            SE_boot   <- apply(boot_full, 2, stats::sd, na.rm = TRUE)
          }

          if (!is.null(qr_key) &&
              !is.null(quantile_region_list) &&
              !is.null(quantile_region_list[[qr_key]]) &&
              is.matrix(quantile_region_list[[qr_key]]) &&
              ncol(quantile_region_list[[qr_key]]) >= 2) {

            qr_mat <- quantile_region_list[[qr_key]]

            if (is.null(rownames(qr_mat))) {
              k2 <- min(nrow(qr_mat), length(nodes_use))
              q_lower[seq_len(k2)] <- qr_mat[seq_len(k2), 1]
              q_upper[seq_len(k2)] <- qr_mat[seq_len(k2), 2]
            } else {
              common <- intersect(nodes_use, rownames(qr_mat))
              if (length(common)) {
                ir   <- match(common, rownames(qr_mat))
                inod <- match(common, nodes_use)
                q_lower[inod] <- qr_mat[ir, 1]
                q_upper[inod] <- qr_mat[ir, 2]
              }
            }
          }

          rows[[k]] <- data.frame(
            node                     = nodes_use,
            layer                    = node_layer[match(nodes_use, nodes_all)],
            scope                    = "inter",
            metric                   = stat_name,
            estimated                = estimated,
            mean.bootstrap           = mean_boot,
            SE.bootstrap             = SE_boot,
            quantile.lower.bootstrap = q_lower,
            quantile.upper.bootstrap = q_upper,
            stringsAsFactors = FALSE
          )
        }

        out_inter <- do.call(rbind, rows)
        rownames(out_inter) <- NULL
      }
    }
  }

  out <- Filter(Negate(is.null), list(out_intra, out_inter))
  out <- if (length(out)) do.call(rbind, out) else NULL

  metric_order <- if (identical(what, "inter")) {
    allowed_inter
  } else if (identical(what, "intra")) {
    allowed_intra
  } else {
    c(allowed_intra, setdiff(allowed_inter, allowed_intra))
  }

  .finalize_centrality(
    out = out,
    quantile_level = quantile_level,
    what = what,
    digits = digits,
    drop_na_boot = drop_na_boot,
    metric_order = metric_order
  )
}

# -------------------------------------------------------------------------
# print method for get_centrality output
# -------------------------------------------------------------------------

#' @export
print.get_centrality <- function(x, digits = 3, top_n = Inf, max_rows = 15, ...) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

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

  order_by <- function(df, cols) {
    cols <- intersect(cols, colnames(df))
    if (!length(cols) || !nrow(df)) return(df)
    o <- do.call(order, df[cols])
    df[o, , drop = FALSE]
  }

  drop_all_na_boot_cols <- function(df) {
    boot_cols <- c(
      "mean.bootstrap",
      "SE.bootstrap",
      "quantile.lower.bootstrap",
      "quantile.upper.bootstrap"
    )
    keep <- setdiff(names(df), boot_cols)
    for (nm in boot_cols) {
      if (nm %in% names(df) && !all(is.na(df[[nm]]))) {
        keep <- c(keep, nm)
      }
    }
    df[, keep, drop = FALSE]
  }

  prettify_colnames <- function(df, quantile_level = 0.95) {
    if (is.null(df) || !nrow(df)) return(df)

    quantile_pct <- paste0(round(100 * quantile_level), "%")

    map <- c(
      "mean.bootstrap" = "mean (bootstrap)",
      "SE.bootstrap" = "SE (bootstrap)",
      "quantile.lower.bootstrap" = paste0(quantile_pct, " quantile lower bound (bootstrap)"),
      "quantile.upper.bootstrap" = paste0(quantile_pct, " quantile upper bound (bootstrap)")
    )

    cn <- colnames(df)
    cn <- ifelse(cn %in% names(map), map[cn], cn)
    cn <- gsub("\\.", " ", cn)
    colnames(df) <- cn
    df
  }

  if (is.null(x) || !nrow(x)) {
    cat("No node-level centrality indices available.\n")
    return(invisible(x))
  }

  quantile_level <- attr(x, "quantile_level") %||% 0.95

  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x <- round_df(x, digits)

  allowed_intra <- c(
    "strength", "expected_influence", "closeness", "betweenness",
    "bridge_strength", "bridge_closeness", "bridge_betweenness",
    "bridge_ei1", "bridge_ei2",
    "bridge_strength_excluded", "bridge_betweenness_excluded",
    "bridge_closeness_excluded", "bridge_ei1_excluded",
    "bridge_ei2_excluded"
  )

  allowed_inter <- c(
    "strength", "expected_influence", "closeness", "betweenness"
  )

  print_metric_block <- function(df, metric_order, block_title,
                                 drop_cols = character(0),
                                 sort_cols = c("layer", "node")) {
    if (is.null(df) || !nrow(df)) return(invisible(NULL))

    cat("\n", block_title, ":\n", sep = "")

    metrics <- intersect(metric_order, unique(df$metric))

    for (met in metrics) {
      sub_met <- df[df$metric == met, , drop = FALSE]
      if (!nrow(sub_met)) next

      sub_met <- order_by(sub_met, sort_cols)
      n_total <- nrow(sub_met)

      if (is.finite(top_n) && "estimated" %in% colnames(sub_met)) {
        o <- order(abs(sub_met$estimated), decreasing = TRUE)
        sub_rank <- sub_met[o, , drop = FALSE]
        sub_rank <- utils::head(sub_rank, top_n)
      } else {
        sub_rank <- sub_met
      }

      n_rank <- nrow(sub_rank)

      if (is.finite(max_rows) && n_rank > max_rows) {
        sub_print <- utils::head(sub_rank, max_rows)
      } else {
        sub_print <- sub_rank
      }

      if ("layer" %in% colnames(sub_print) &&
          length(unique(sub_print$layer)) == 1L &&
          all(sub_print$layer == "1")) {
        sub_print <- sub_print[, setdiff(colnames(sub_print), "layer"), drop = FALSE]
      }

      sub_print <- sub_print[, setdiff(colnames(sub_print), c("metric", "scope", drop_cols)), drop = FALSE]
      sub_print <- drop_all_na_boot_cols(sub_print)
      sub_print <- prettify_colnames(sub_print, quantile_level = quantile_level)

      cat("\n  Metric:", met, "\n")
      print(sub_print, row.names = FALSE)

      if (n_rank > nrow(sub_print)) {
        cat("... showing ", nrow(sub_print), " of ", n_rank, " rows", sep = "")
        if (n_total > n_rank) {
          cat(" (", n_total, " total before top_n)", sep = "")
        }
        cat("\n")
      } else if (n_total > n_rank) {
        cat("... showing ", n_rank, " of ", n_total, " rows after top_n filtering\n", sep = "")
      }
    }

    invisible(NULL)
  }

  if ("intra" %in% unique(x$scope)) {
    x_intra <- x[x$scope == "intra", , drop = FALSE]
    print_metric_block(
      df = x_intra,
      metric_order = allowed_intra,
      block_title = "Node-level centrality indices (intralayer)",
      drop_cols = character(0),
      sort_cols = c("layer", "node")
    )
  }

  if ("inter" %in% unique(x$scope)) {
    x_inter <- x[x$scope == "inter", , drop = FALSE]
    print_metric_block(
      df = x_inter,
      metric_order = allowed_inter,
      block_title = "Node-level centrality indices (interlayer-only graph)",
      drop_cols = character(0),
      sort_cols = c("layer", "node")
    )
  }

  invisible(x)
}
