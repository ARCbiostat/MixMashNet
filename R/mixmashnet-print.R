#' Print method for MixMashNet objects
#'
#' @description
#' Compact textual summary for objects returned by \code{mixMN()} and
#' \code{multimixMN()}. The method reports:
#' \itemize{
#'   \item whether the fit is single-layer (\code{mixMN}) or multilayer
#'         (\code{multimixMN});
#'   \item number of subjects (if available) and variables;
#'   \item for multilayer fits, the number of nodes and non-zero edges per layer
#'         and, if present, per interlayer pair;
#'   \item size of the global graph (nodes and edges);
#'   \item number of communities (single-layer) or communities per layer
#'         (multilayer);
#'   \item covariates used for adjustment and nodes excluded from the graph
#'         and/or clustering;
#'   \item main settings for community detection and bootstrap.
#' }
#'
#' @param x An object of class \code{mixmashnet}, typically returned by
#'   \code{mixMN()} or \code{multimixMN()}.
#' @param ... Additional arguments.
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @seealso \code{\link{mixMN}}, \code{\link{multimixMN}}
#' @method print mixmashnet
#' @export
print.mixmashnet <- function(x, ...) {
  cls <- class(x)

  type_label <- if ("multimixMN_fit" %in% cls) {
    "Multilayer MGM (multimixMN)"
  } else if ("mixMN_fit" %in% cls) {
    "Single-layer MGM (mixMN)"
  } else {
    "Unknown MixMashNet object"
  }

  cat("MixMashNet fit\n")
  cat("  Type: ", type_label, "\n", sep = "")

  ## ---- n subjects + p variables ----
  p <- length(x$model$nodes)
  n <- x$model$n %||% NA_integer_

  if (is.finite(n)) {
    cat("  Data: ", n, " subjects x ", p, " variables\n", sep = "")
  } else {
    cat("  Variables: ", p, "\n", sep = "")
  }

  ## ---- multilayer: info on layers ----
  if ("multimixMN_fit" %in% cls && !is.null(x$layer_fits)) {
    cat("  Layers (", length(x$layer_fits), "):\n", sep = "")
    for (L in names(x$layer_fits)) {
      fitL <- x$layer_fits[[L]]

      if (!is.null(fitL$graph$keep_nodes_graph)) {
        n_nodes <- length(fitL$graph$keep_nodes_graph)
      } else if (!is.null(x$layers$assignment)) {
        n_nodes <- sum(x$layers$assignment == L)
      } else {
        n_nodes <- NA_integer_
      }

      edges_df <- NULL
      if (!is.null(fitL$statistics$edge$true)) {
        edges_df <- fitL$statistics$edge$true
      }
      n_edges_nz <- if (!is.null(edges_df) && nrow(edges_df) && "weight" %in% names(edges_df)) {
        sum(edges_df$weight != 0, na.rm = TRUE)
      } else {
        0L
      }

      cat("    - ", L, ": ", n_nodes, " nodes, ",
          n_edges_nz, " edges\n", sep = "")
    }
  }

  ## ---- interlayer info ----
  if ("multimixMN_fit" %in% cls && !is.null(x$interlayer)) {
    pair_names <- setdiff(names(x$interlayer), "centrality")
    if (length(pair_names)) {
      cat("  Interlayer edges:\n")
      for (pk in pair_names) {
        ed_obj <- x$interlayer[[pk]]$edges
        edges_df <- if (!is.null(ed_obj) && !is.null(ed_obj$true)) ed_obj$true else NULL
        n_edges_nz <- if (!is.null(edges_df) && nrow(edges_df) && "weight" %in% names(edges_df)) {
          sum(edges_df$weight != 0, na.rm = TRUE)
        } else {
          0L
        }
        cat("    - ", pk, ": ", n_edges_nz, " edges\n", sep = "")
      }
    }
  }


  ## ---- global graph ----
  if (!is.null(x$graph$igraph)) {
    g <- x$graph$igraph
    v <- igraph::vcount(g)
    e <- igraph::ecount(g)
    cat("  Graph: ", v, " nodes, ", e, " edges\n", sep = "")
  }

  ## ---- info community ----
  if ("mixMN_fit" %in% cls && !is.null(x$communities$groups)) {
    k <- length(stats::na.omit(unique(as.integer(x$communities$groups))))
    cat("  Communities: ", k, " \n", sep = "")
  } else if ("multimixMN_fit" %in% cls && !is.null(x$layer_fits)) {
    cat("  Communities per layer:\n")
    K <- vapply(
      x$layer_fits,
      function(fitL) {
        gr <- fitL$communities$groups
        if (is.null(gr)) return(0L)
        length(stats::na.omit(unique(as.integer(gr))))
      },
      integer(1)
    )
    for (nm in names(K)) {
      cat("    - ", nm, ": ", K[[nm]], "\n", sep = "")
    }
  }

  ## ---- covariate variables ----
  covars <- x$settings$covariates %||% NULL
  if (!is.null(covars) && length(covars)) {
    cat("  Covariates (adjusted for): ",
        paste(covars, collapse = ", "),
        "\n", sep = "")
  }

  ## ---- nodes excluded from graph / clustering ----
  excl_g <- x$settings$exclude_from_graph %||% NULL
  excl_c <- x$settings$exclude_from_cluster %||% NULL

  if (!is.null(excl_g) && length(excl_g)) {
    cat("  Nodes excluded from graph: ",
        paste(excl_g, collapse = ", "),
        "\n", sep = "")
  }

  if (!is.null(excl_c) && length(excl_c)) {
    excl_only_clust <- setdiff(excl_c, excl_g)
    if (length(excl_only_clust)) {
      cat("  Nodes excluded from clustering only: ",
          paste(excl_only_clust, collapse = ", "),
          "\n", sep = "")
    } else {
      cat("  Nodes excluded from clustering: ",
          paste(excl_c, collapse = ", "),
          "\n", sep = "")
    }
  }

  ## ---- settings bootstrap / cluster ----
  s <- x$settings
  if (!is.null(s$cluster_method)) {
    cat("  Community detection: ", s$cluster_method, "\n", sep = "")
  }
  if (!is.null(s$reps)) {
    cat("  Bootstrap replications: ", s$reps, "\n", sep = "")
  }
  if (!is.null(s$boot_what)) {
    cat("  Bootstrapped quantities: ",
        paste(s$boot_what, collapse = ", "),
        "\n", sep = "")
  }


  invisible(x)
}
