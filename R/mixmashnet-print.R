#' Print method for MixMashNet objects
#'
#' @description
#' Compact textual summary for objects returned by \code{mixMN()} and
#' \code{multimixMN()}. The method reports:
#' \itemize{
#'   \item whether the fit is single layer (\code{mixMN}) or multilayer
#'         (\code{multimixMN});
#'   \item number of subjects (if available) and variables;
#'   \item for multilayer fits, the number of nodes and non-zero edges per layer
#'         and, if present, per interlayer pair;
#'   \item size of the global graph (nodes and edges);
#'   \item number of communities (single layer) or communities per layer
#'         (multilayer);
#'   \item covariates used for adjustment and nodes excluded from the graph
#'         and/or clustering;
#'   \item main settings for community detection and bootstrap;
#'   \item data info.
#' }
#'
#' @param x An object of class \code{mixmashnet}, returned by
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
    "Single layer MGM (mixMN)"
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


  ## ---- nodes excluded from graph / clustering ----
  excl_g <- x$settings$covariates %||% NULL
  excl_c <- x$settings$exclude_from_cluster %||% NULL

  if (!is.null(excl_g) && length(excl_g)) {
    cat("  Covariates (adjusted for): ",
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

  reps_i <- as.integer(s$reps %||% 0L)
  if (is.finite(reps_i) && reps_i > 0L) {
    if (!is.null(s$boot_what) && length(s$boot_what)) {
      cat("  Bootstrapped quantities: ",
          paste(s$boot_what, collapse = ", "),
          "\n", sep = "")
    }
  } else {
    cat("  Bootstrapped quantities: none (reps = 0)\n", sep = "")
  }

  ## ---- data_info: inferred mgm types + ONLY variables recoded to {0,1} ----
  if (!is.null(x$data_info)) {

    tl <- x$data_info$mgm_type_level

    # --- report inferred "c" and "p" ---
    if (is.data.frame(tl) && nrow(tl) > 0 && all(c("node", "mgm_type") %in% names(tl))) {

      vars_c <- tl$node[tl$mgm_type == "c"]
      vars_p <- tl$node[tl$mgm_type == "p"]

      if (length(vars_c) || length(vars_p)) {
        cat("  Data info:\n")
        if (length(vars_c)) {
          cat("    - Inferred as 'c' (categorical): ",
              paste(vars_c, collapse = ", "), "\n", sep = "")
        }
        if (length(vars_p)) {
          cat("    - Inferred as 'p' (Poisson/count): ",
              paste(vars_p, collapse = ", "), "\n", sep = "")
        }
      }
    }

    # --- ONLY variables that were actually recoded to {0,1} ---
    brm <- x$data_info$binary_recode_map
    if (!is.null(brm) && length(brm) > 0) {

      fmt_map <- function(m) {
        # m is a named numeric vector: names = original labels, values = 0/1
        if (is.atomic(m) && !is.null(names(m))) {
          paste0(names(m), "->", as.character(m), collapse = ", ")
        } else {
          paste(as.character(m), collapse = ", ")
        }
      }

      # ensure we have the "Data info:" header even if the 'c/p' block didn't print
      if (!(is.data.frame(tl) && nrow(tl) > 0 && all(c("node", "mgm_type") %in% names(tl)) &&
            (length(tl$node[tl$mgm_type == "c"]) || length(tl$node[tl$mgm_type == "p"])))) {
        cat("  Data info:\n")
      }

      cat("    - Recoded to {0,1} for mgm fitting: ", length(brm), " variable(s)\n", sep = "")
      nm <- names(brm)
      show_nm <- nm[seq_len(min(length(nm), 10))]
      for (v in show_nm) {
        cat("      * ", v, ": ", fmt_map(brm[[v]]), "\n", sep = "")
      }
      if (length(nm) > 10) {
        cat("      * ... (", length(nm) - 10, " more)\n", sep = "")
      }
    }
  }

  invisible(x)
}
