# -------------------------------------------------------------------------
# internal helpers
# -------------------------------------------------------------------------

#' Internal helper
#' @keywords internal
#' @noRd
.build_edge_long <- function(true_edges,
                             boot_mat = NULL,
                             layer_value = NA_character_,
                             pairs_value = NA_character_,
                             scope_label,
                             quantile_level) {

  if (is.null(true_edges) || !nrow(true_edges)) return(NULL)

  if (!all(c("edge", "weight") %in% colnames(true_edges))) {
    stop("Edge summary table must contain columns 'edge' and 'weight'.")
  }

  idx_keep <- !is.na(true_edges$weight) & (true_edges$weight != 0)

  if (!any(idx_keep)) return(NULL)

  true_edges <- true_edges[idx_keep, , drop = FALSE]
  edge_names <- as.character(true_edges$edge)
  estimated  <- true_edges$weight

  mean_boot <- SE_boot <- rep(NA_real_, length(edge_names))
  q_lower   <- q_upper <- rep(NA_real_, length(edge_names))

  if (!is.null(boot_mat) && is.matrix(boot_mat)) {

    if (is.null(rownames(boot_mat))) {
      rownames(boot_mat) <- edge_names[seq_len(min(nrow(boot_mat), length(edge_names)))]
    }

    boot_full <- matrix(
      NA_real_,
      nrow = length(edge_names),
      ncol = ncol(boot_mat),
      dimnames = list(edge_names, colnames(boot_mat))
    )

    common <- intersect(edge_names, rownames(boot_mat))
    if (length(common)) {
      boot_full[common, ] <- boot_mat[common, , drop = FALSE]
    }

    mean_boot <- rowMeans(boot_full, na.rm = TRUE)
    SE_boot   <- apply(boot_full, 1, stats::sd, na.rm = TRUE)

    alpha <- (1 - quantile_level) / 2
    probs <- c(alpha, 1 - alpha)

    q_lower <- apply(
      boot_full, 1, stats::quantile,
      probs = probs[1], na.rm = TRUE, type = 6
    )
    q_upper <- apply(
      boot_full, 1, stats::quantile,
      probs = probs[2], na.rm = TRUE, type = 6
    )
  }

  out <- data.frame(
    edge                     = edge_names,
    layer                    = rep(layer_value, length(edge_names)),
    pairs                    = rep(pairs_value, length(edge_names)),
    scope                    = rep(scope_label, length(edge_names)),
    estimated                = estimated,
    mean.bootstrap           = mean_boot,
    SE.bootstrap             = SE_boot,
    quantile.lower.bootstrap = q_lower,
    quantile.upper.bootstrap = q_upper,
    stringsAsFactors = FALSE
  )

  rownames(out) <- NULL
  out
}

#' Internal helper
#' @keywords internal
#' @noRd
.finalize_edges <- function(out, quantile_level, what, digits,
                                    drop_na_boot) {
  if (is.null(out) || !nrow(out)) {
    out <- tibble::tibble(
      edge = character(),
      layer = character(),
      pairs = character(),
      scope = character(),
      estimated = numeric()
    )
    attr(out, "quantile_level") <- quantile_level
    attr(out, "what") <- what
    class(out) <- c("get_edges", class(out))
    return(out)
  }

  out <- .round_numcols(out, digits)

  if (isTRUE(drop_na_boot)) {
    out <- .drop_all_na_boot_cols(out)
  }

  out <- out[order(out$scope, out$layer, out$pairs, out$edge), , drop = FALSE]

  if (identical(what, "intra")) {
    out <- out[, setdiff(colnames(out), "pairs"), drop = FALSE]
  } else if (identical(what, "inter")) {
    out <- out[, setdiff(colnames(out), "layer"), drop = FALSE]
  }

  out <- tibble::as_tibble(out)

  attr(out, "quantile_level") <- quantile_level
  attr(out, "what") <- what
  class(out) <- c("get_edges", class(out))

  out
}

#' Extract edge-level summaries
#'
#' @description
#' Extracts edge-level summaries from fitted objects returned by
#' \code{mixMN()} and \code{multimixMN()} in a long-format data frame.
#'
#' For single layer fits of class \code{"mixMN_fit"}, only intralayer
#' edges are available.
#'
#' For multilayer fits of class \code{"multimixMN_fit"}, \code{what = "intra"}
#' returns intralayer edges, whereas \code{what = "inter"} returns interlayer
#' edges. If \code{what} is not specified for a multilayer fit, both scopes are
#' returned by default, unless \code{layer} or \code{pairs} imply a specific
#' scope.
#'
#' The function returns the original edge weights and, when available,
#' bootstrap means, standard errors, and bootstrap quantile regions.
#'
#' @aliases get_edges get_edges.mixMN_fit get_edges.multimixMN_fit
#' @param object A fitted object of class \code{"mixMN_fit"} or
#'   \code{"multimixMN_fit"}.
#' @param what Character string indicating which edge-level summaries to extract:
#'   \itemize{
#'     \item \code{"intra"}: intralayer edges;
#'     \item \code{"inter"}: interlayer edges (available only for
#'       \code{"multimixMN_fit"} objects).
#'   }
#'   For single layer fits, only \code{"intra"} is allowed.
#'   For multilayer fits, if omitted, both scopes are returned by default
#'   unless \code{layer} or \code{pairs} imply a specific scope.
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
#' The returned data frame is in long format, with one row per edge.
#'
#' For single layer fits, only \code{what = "intra"} is available.
#'
#' For multilayer fits, \code{layer} can be used to subset intralayer output,
#' whereas \code{pairs} can be used to subset interlayer output.
#'
#' @return
#' A tibble in long format with one row per edge.
#' It contains the columns:
#' \itemize{
#'   \item \code{edge}
#'   \item \code{layer} for intralayer edges
#'   \item \code{pairs} for interlayer edges
#'   \item \code{scope}
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
get_edges <- function(object, ...) {
  UseMethod("get_edges")
}

#' @rdname get_edges
#' @export
get_edges.mixMN_fit <- function(object,
                                    what = "intra",
                                    digits = NULL,
                                    drop_na_boot = TRUE,
                                    ...) {
  what <- match.arg(what, choices = "intra")

  quantile_level <- .get_quantile_level(object)

  out <- .build_edge_long(
    true_edges     = object$statistics$edge$true,
    boot_mat       = object$statistics$edge$boot,
    layer_value    = "1",
    pairs_value    = NA_character_,
    scope_label    = "intra",
    quantile_level = quantile_level
  )

  .finalize_edges(
    out = out,
    quantile_level = quantile_level,
    what = what,
    digits = digits,
    drop_na_boot = drop_na_boot
  )
}

#' @rdname get_edges
#' @export
get_edges.multimixMN_fit <- function(object,
                                         what = c("intra", "inter"),
                                         layer = NULL,
                                         pairs = NULL,
                                         digits = NULL,
                                         drop_na_boot = TRUE,
                                         ...) {

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

  quantile_level <- .get_quantile_level(object)

  out_intra <- NULL
  out_inter <- NULL

  # -----------------------------------------------------------------------
  # intralayer edges
  # -----------------------------------------------------------------------
  if ("intra" %in% what) {
    all_layers <- names(object$layer_fits)
    layers_sel <- if (is.null(layer)) all_layers else intersect(all_layers, layer)

    blocks <- vector("list", length(layers_sel))

    for (j in seq_along(layers_sel)) {
      L <- layers_sel[j]
      fitL <- object$layer_fits[[L]]
      if (is.null(fitL) || is.null(fitL$statistics$edge)) next

      blocks[[j]] <- .build_edge_long(
        true_edges     = fitL$statistics$edge$true,
        boot_mat       = fitL$statistics$edge$boot,
        layer_value    = L,
        pairs_value    = NA_character_,
        scope_label    = "intra",
        quantile_level = quantile_level
      )
    }

    blocks <- Filter(Negate(is.null), blocks)
    out_intra <- if (length(blocks)) do.call(rbind, blocks) else NULL
  }

  # -----------------------------------------------------------------------
  # interlayer edges
  # -----------------------------------------------------------------------
  if ("inter" %in% what) {
    all_pairs <- setdiff(names(object$interlayer), "centrality")

    if (is.null(pairs)) {
      pairs_sel <- all_pairs
    } else {
      norm_all   <- .normalize_pairs(all_pairs)
      norm_pairs <- .normalize_pairs(pairs)
      pairs_sel  <- all_pairs[norm_all %in% norm_pairs]
    }

    blocks <- vector("list", length(pairs_sel))

    for (j in seq_along(pairs_sel)) {
      pp <- pairs_sel[j]
      edge_obj <- object$interlayer[[pp]]$edges
      if (is.null(edge_obj)) next

      blocks[[j]] <- .build_edge_long(
        true_edges     = edge_obj$true,
        boot_mat       = edge_obj$boot,
        layer_value    = NA_character_,
        pairs_value    = pp,
        scope_label    = "inter",
        quantile_level = quantile_level
      )
    }

    blocks <- Filter(Negate(is.null), blocks)
    out_inter <- if (length(blocks)) do.call(rbind, blocks) else NULL
  }

  out <- Filter(Negate(is.null), list(out_intra, out_inter))
  out <- if (length(out)) do.call(rbind, out) else NULL

  .finalize_edges(
    out = out,
    quantile_level = quantile_level,
    what = what,
    digits = digits,
    drop_na_boot = drop_na_boot
  )
}

#' @export
print.get_edges <- function(x, digits = 3, top_n = Inf, max_rows = 15, ...) {

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

  prettify_colnames <- function(df, quantile_level = 0.95) {
    if (is.null(df) || !nrow(df)) return(df)

    quantile_pct <- paste0(round(100 * quantile_level), "%")

    map <- c(
      "mean.bootstrap" = "mean (bootstrap)",
      "SE.bootstrap" = "SE (bootstrap)",
      "quantile.lower.bootstrap" =
        paste0(quantile_pct, " quantile lower bound (bootstrap)"),
      "quantile.upper.bootstrap" =
        paste0(quantile_pct, " quantile upper bound (bootstrap)")
    )

    cn <- colnames(df)
    cn <- ifelse(cn %in% names(map), map[cn], cn)
    cn <- gsub("\\.", " ", cn)
    colnames(df) <- cn
    df
  }

  if (is.null(x) || !nrow(x)) {
    cat("No edge-level summaries available.\n")
    return(invisible(x))
  }

  quantile_level <- attr(x, "quantile_level") %||% 0.95

  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x <- round_df(x, digits)

  print_edge_block <- function(df, block_title, group_var, group_label,
                               sort_cols, drop_cols = character(0)) {
    if (is.null(df) || !nrow(df)) return(invisible(NULL))

    cat("\n", block_title, ":\n", sep = "")

    groups <- unique(df[[group_var]])

    for (g in groups) {
      sub_g <- df[df[[group_var]] == g, , drop = FALSE]
      sub_g <- order_by(sub_g, sort_cols)

      n_total <- nrow(sub_g)

      if (is.finite(top_n) && "estimated" %in% colnames(sub_g)) {
        o <- order(abs(sub_g$estimated), decreasing = TRUE)
        sub_rank <- sub_g[o, , drop = FALSE]
        sub_rank <- utils::head(sub_rank, top_n)
      } else {
        sub_rank <- sub_g
      }

      n_rank <- nrow(sub_rank)

      if (is.finite(max_rows) && n_rank > max_rows) {
        sub_print <- utils::head(sub_rank, max_rows)
      } else {
        sub_print <- sub_rank
      }

      sub_print <- sub_print[, setdiff(
        colnames(sub_print),
        c("scope", group_var, drop_cols)
      ), drop = FALSE]

      sub_print <- sub_print[, setdiff(colnames(sub_print), c("scope", drop_cols)), drop = FALSE]
      sub_print <- .drop_all_na_boot_cols(sub_print)
      sub_print <- prettify_colnames(sub_print, quantile_level = quantile_level)

      cat("\n  ", group_label, ": ", g, "\n", sep = "")
      print(sub_print, row.names = FALSE)

      if (n_rank > nrow(sub_print)) {
        cat("... showing ", nrow(sub_print), " of ", n_rank, " rows", sep = "")
        if (n_total > n_rank) {
          cat(" (", n_total, " total before top_n)", sep = "")
        }
        cat("\n")
      } else if (n_total > n_rank) {
        cat("... showing ", n_rank, " of ", n_total,
            " rows after top_n filtering\n", sep = "")
      }
    }

    invisible(NULL)
  }

  if ("intra" %in% unique(x$scope)) {
    x_intra <- x[x$scope == "intra", , drop = FALSE]
    print_edge_block(
      df = x_intra,
      block_title = "Edge-level summaries (intralayer)",
      group_var = "layer",
      group_label = "Layer",
      sort_cols = c("layer", "edge"),
      drop_cols = "pairs"
    )
  }

  if ("inter" %in% unique(x$scope)) {
    x_inter <- x[x$scope == "inter", , drop = FALSE]
    print_edge_block(
      df = x_inter,
      block_title = "Edge-level summaries (interlayer)",
      group_var = "pairs",
      group_label = "Pair",
      sort_cols = c("pairs", "edge"),
      drop_cols = "layer"
    )
  }

  invisible(x)
}
