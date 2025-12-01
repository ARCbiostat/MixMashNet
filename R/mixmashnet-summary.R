#' Summarize MixMashNet fits (single- and multilayer) in long format
#'
#' @param object An object of class \code{"mixmashnet"} returned by
#'   \code{mixMN()} or \code{multimixMN()}.
#' @param what Character vector specifying what to summarize:
#'   \itemize{
#'     \item \code{"index"}: node-level indices on the (intra-layer) graph.
#'     \item \code{"edges"}: intra-layer edges.
#'     \item \code{"interlayer_index"}: node-level indices computed on the
#'       interlayer-only graph (from \code{object$interlayer$centrality}).
#'     \item \code{"interlayer_edges"}: cross-layer edges
#'       (from \code{object$interlayer[[pairs]]$edges}).
#'   }
#' @param metrics Optional character vector of metrics to include. If \code{NULL},
#'   all available metrics in the corresponding object are used.
#'   For \code{what = "interlayer_index"}, valid metrics are typically
#'   \code{c("strength","ei1","closeness","betweenness")}.
#' @param layer Optional character vector of layer names to subset:
#'   used for \code{what = "index"} and \code{what = "edges"} in multilayer fits.
#'   Ignored for single-layer fits.
#' @param pairs Optional character vector of layer-pair names (e.g. \code{"bio_dis"})
#'   used for \code{what = "interlayer_index"} and \code{what = "interlayer_edges"}.
#'   If \code{NULL}, all available layer pairs are included.
#'   For \code{what = "interlayer_index"}, if a single \code{pairs} is given,
#'   only nodes belonging to one of the two layers in that pair are reported.
#' @param digits Number of digits to round numeric summaries.
#' @param ... Not used (for S3 compatibility).
#'
#' @return A list (class \code{"summary.mixmashnet"}) with up to four data frames:
#'   \itemize{
#'     \item \code{$index}: node-level indices (one row per node-metric).
#'     \item \code{$edges}: intra-layer edges (one row per edge).
#'     \item \code{$interlayer_index}: interlayer-only node indices.
#'     \item \code{$interlayer_edges}: cross-layer edges.
#'   }
#' @method summary mixmashnet
#' @export
summary.mixmashnet <- function(object,
                               what       = c("index", "edges",
                                              "interlayer_index", "interlayer_edges"),
                               metrics    = NULL,
                               layer      = NULL,
                               pairs = NULL,
                               digits     = 3,
                               ...) {
  if (!inherits(object, "mixmashnet")) {
    stop("`summary.mixmashnet()` expects an object of class 'mixmashnet'.")
  }

  what     <- match.arg(what, several.ok = TRUE)
  is_multi <- inherits(object, "multimixMN_fit")

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

  # -------------------- 1) NODE-LEVEL INDEX (intra-layer) -------------------- #
  index_tab <- NULL
  if ("index" %in% what) {

    build_index_long <- function(true_df, boot_list, ci_list,
                                 layer_label, metrics, digits) {
      if (is.null(true_df) || !nrow(true_df)) return(NULL)
      if (!"node" %in% colnames(true_df)) {
        stop("`statistics$node$true` must contain a column 'node'.")
      }

      nodes <- as.character(true_df$node)

      # choose metrics
      if (is.null(metrics)) {
        metrics_block <- setdiff(colnames(true_df), "node")
      } else {
        metrics_block <- intersect(metrics, setdiff(colnames(true_df), "node"))
      }
      if (!length(metrics_block)) return(NULL)

      rows <- list()

      for (met in metrics_block) {
        estimated <- true_df[[met]]

        mean_boot <- sd_boot <- rep(NA_real_, length(nodes))
        lower     <- upper   <- rep(NA_real_, length(nodes))

        # bootstrap (if available)
        if (!is.null(boot_list) &&
            !is.null(boot_list[[met]]) &&
            is.matrix(boot_list[[met]])) {

          boot_mat <- boot_list[[met]]
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
            k <- min(nrow(ci_mat), length(nodes))
            lower[seq_len(k)] <- ci_mat[seq_len(k), 1]
            upper[seq_len(k)] <- ci_mat[seq_len(k), 2]
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

        rows[[met]] <- data.frame(
          node                 = nodes,
          layer                = layer_label,
          metric               = met,
          estimated            = estimated,
          mean.bootstrap       = mean_boot,
          SE.bootstrap         = sd_boot,
          CI.lower.bootstrap   = lower,
          CI.upper.bootstrap  = upper,
          stringsAsFactors = FALSE
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
          true_df   = true_df,
          boot_list = boot_list,
          ci_list   = ci_list,
          layer_label = L,
          metrics   = metrics,
          digits    = digits
        )
      }
      blocks <- Filter(Negate(is.null), blocks)
      if (length(blocks)) index_tab <- do.call(rbind, blocks)

    } else {
      # single-layer: no concept of layer; use NA
      true_df   <- object$statistics$node$true
      boot_list <- object$statistics$node$boot
      ci_list   <- object$statistics$node$ci

      index_tab <- build_index_long(
        true_df     = true_df,
        boot_list   = boot_list,
        ci_list     = ci_list,
        layer_label = "1",
        metrics     = metrics,
        digits      = digits
      )
    }
  }

  # ------------------------ 2) EDGES (intra-layer) ------------------------ #
  edges_tab <- NULL
  if ("edges" %in% what) {

    build_edge_long <- function(true_edges, boot_mat, layer_label, digits) {
      if (is.null(true_edges) || !nrow(true_edges)) return(NULL)
      if (!all(c("edge","weight") %in% colnames(true_edges))) {
        stop("`statistics$edge$true` must contain columns 'edge' and 'weight'.")
      }

      # solo archi con peso osservato diverso da 0
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

        lower <- apply(boot_full, 1, stats::quantile,
                       probs = 0.025, na.rm = TRUE, type = 6)
        upper <- apply(boot_full, 1, stats::quantile,
                       probs = 0.975, na.rm = TRUE, type = 6)
      }

      out <- data.frame(
        edge                = edges,
        layer               = layer_label,
        estimated           = weight_obs,
        mean.bootstrap      = mean_boot,
        SE.bootstrap        = sd_boot,
        CI.lower.bootstrap  = lower,
        CI.upper.bootstrap  = upper,
        stringsAsFactors    = FALSE
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

  # ------------------ 3) INTERLAYER NODE INDEX (interlayer_index) ------------------ #
  inter_idx_tab <- NULL
  if ("interlayer_index" %in% what && is_multi && !is.null(object$interlayer)) {

    inter_cent <- object$interlayer$centrality
    if (!is.null(inter_cent) && !is.null(inter_cent$true)) {
      true_df   <- inter_cent$true
      boot_list <- inter_cent$boot
      ci_list   <- inter_cent$ci_results

      if (!"node" %in% colnames(true_df)) {
        stop("`object$interlayer$centrality$true` must contain a column 'node'.")
      }
      nodes_all <- as.character(true_df$node)

      # mapping metric names per interlayer: colnames vs ci_list names
      # true_df columns: typically strength, ei1, closeness, betweenness
      # ci_list names:   strength, expected_influence, closeness, betweenness
      metric_map_ci <- list(
        strength  = "strength",
        ei1       = "expected_influence",
        closeness = "closeness",
        betweenness = "betweenness"
      )

      # subset nodes by pairs (if given)
      nodes_use <- nodes_all
      if (!is.null(pairs) && length(pairs) == 1L &&
          !is.null(object$layers) && !is.null(object$layers$assignment)) {

        # es. "bio_dis" -> c("bio","dis")
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

      if (!length(nodes_use)) {
        inter_idx_tab <- NULL
      } else {
        # scegli metriche
        if (is.null(metrics)) {
          metrics_int <- setdiff(colnames(true_df), "node")
        } else {
          metrics_int <- intersect(metrics, setdiff(colnames(true_df), "node"))
        }
        if (!length(metrics_int)) {
          inter_idx_tab <- NULL
        } else {
          rows <- list()
          for (met in metrics_int) {
            estimated <- true_df[match(nodes_use, nodes_all), met]

            mean_boot <- sd_boot <- rep(NA_real_, length(nodes_use))
            lower     <- upper   <- rep(NA_real_, length(nodes_use))

            # boot
            if (!is.null(boot_list) &&
                !is.null(boot_list[[met]]) &&
                is.matrix(boot_list[[met]])) {

              boot_mat <- boot_list[[met]]
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

            # CI: attenzione alla mappa dei nomi
            ci_key <- metric_map_ci[[met]]
            if (!is.null(ci_key) &&
                !is.null(ci_list) &&
                !is.null(ci_list[[ci_key]]) &&
                is.matrix(ci_list[[ci_key]]) &&
                ncol(ci_list[[ci_key]]) >= 2) {

              ci_mat <- ci_list[[ci_key]]
              if (is.null(rownames(ci_mat))) {
                k <- min(nrow(ci_mat), length(nodes_use))
                lower[seq_len(k)] <- ci_mat[seq_len(k), 1]
                upper[seq_len(k)] <- ci_mat[seq_len(k), 2]
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

            rows[[met]] <- data.frame(
              node               = nodes_use,
              metric             = met,
              estimated          = estimated,
              mean.bootstrap     = mean_boot,
              SE.bootstrap       = sd_boot,
              CI.lower.bootstrap = lower,
              CI.upper.bootstrap = upper,
              stringsAsFactors = FALSE
            )
          }
          out <- do.call(rbind, rows)
          rownames(out) <- NULL
          inter_idx_tab <- .round_numcols(out, digits)
        }
      }
    }
  }

  # ------------------ 4) INTERLAYER EDGES (interlayer_edges) ------------------ #
  inter_edges_tab <- NULL
  if ("interlayer_edges" %in% what && is_multi && !is.null(object$interlayer)) {

    # tutti i nomi pairs disponibili (escluso "centrality")
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

        lower <- apply(boot_full, 1, stats::quantile,
                       probs = 0.025, na.rm = TRUE, type = 6)
        upper <- apply(boot_full, 1, stats::quantile,
                       probs = 0.975, na.rm = TRUE, type = 6)
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
    interlayer_edges = inter_edges_tab
  )
  class(out) <- c("summary.mixmashnet", "list")
  out
}


#' Print method for \code{"summary.mixmashnet"} objects
#'
#' @param x A \code{"summary.mixmashnet"} object as returned by
#'   \code{summary.mixmashnet()}.
#' @param digits Number of digits to print.
#' @param ... Further arguments passed to or from other methods (currently unused).
#'
#' @method print summary.mixmashnet
#' @export
print.summary.mixmashnet <- function(x, digits = 3, ...) {

  # helper per arrotondare al volo (non tocca l'oggetto originale)
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

  # helper per ordinare in modo robusto per una o più colonne
  order_by <- function(df, cols) {
    cols <- intersect(cols, colnames(df))
    if (!length(cols) || !nrow(df)) return(df)
    o <- do.call(order, df[cols])
    df[o, , drop = FALSE]
  }

  # helper per rendere i nomi delle colonne più leggibili in stampa
  prettify_colnames <- function(df) {
    if (is.null(df) || !nrow(df)) return(df)

    cn <- colnames(df)

    # mapping esplicito per le colonne bootstrap
    map <- c(
      "mean.bootstrap"      = "mean (bootstrap)",
      "SE.bootstrap"        = "SE (bootstrap)",
      "CI.lower.bootstrap"  = "95% CI lower (bootstrap)",
      "CI.upper.bootstrap"  = "95% CI upper (bootstrap)"
    )

    cn <- ifelse(cn %in% names(map), map[cn], cn)

    # per il resto, solo estetica: punto -> spazio
    cn <- gsub("\\.", " ", cn)

    colnames(df) <- cn
    df
  }

  # ---------- 1) INDEX (intra-layer, node-level) ----------
  if (!is.null(x$index) && nrow(x$index)) {
    cat("\nNode-level indices (intra-layer):\n")

    idx <- round_df(x$index, digits)

    metrics_idx <- unique(idx$metric)

    for (met in metrics_idx) {
      sub_met <- idx[idx$metric == met, , drop = FALSE]
      cat("\n  Metric:", met, "\n")

      # ordina per layer, poi node (quando presenti)
      sub_met <- order_by(sub_met, c("layer", "node"))

      # non ristampo la colonna "metric" (già nel titolo)
      sub_met <- sub_met[, setdiff(colnames(sub_met), "metric"), drop = FALSE]

      # nomi colonne più leggibili
      sub_met <- prettify_colnames(sub_met)

      print(sub_met, row.names = FALSE)
    }
  }

  # ---------- 2) EDGES (intra-layer) ----------
  if (!is.null(x$edges) && nrow(x$edges)) {
    cat("\nIntra-layer edges:\n")

    ed <- round_df(x$edges, digits)

    # ordina per layer, poi edge se presenti
    ed <- order_by(ed, c("layer", "edge"))

    ed <- prettify_colnames(ed)

    cat("\n")
    print(ed, row.names = FALSE)
  }

  # ---------- 3) INTERLAYER INDEX ----------
  if (!is.null(x$interlayer_index) && nrow(x$interlayer_index)) {
    cat("\nInterlayer-only node indices:\n")

    idx2 <- round_df(x$interlayer_index, digits)

    metrics_int <- unique(idx2$metric)

    for (met in metrics_int) {
      sub_met <- idx2[idx2$metric == met, , drop = FALSE]
      cat("\n  Metric:", met, "\n")

      # se in futuro aggiungi una colonna "layer", verrà usata automaticamente
      sub_met <- order_by(sub_met, c("layer", "node"))

      sub_met <- sub_met[, setdiff(colnames(sub_met), "metric"), drop = FALSE]

      sub_met <- prettify_colnames(sub_met)

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
        sub_pp <- prettify_colnames(sub_pp)
        cat("\n  Pair:", pp, "\n")
        print(sub_pp, row.names = FALSE)
      }
    } else {
      ed2 <- order_by(ed2, c("edge"))
      ed2 <- prettify_colnames(ed2)
      cat("\n")
      print(ed2, row.names = FALSE)
    }
  }

  invisible(x)
}
