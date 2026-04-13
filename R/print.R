# internal helper ---------------------------------------------------------

#' Internal helper for print methods
#' @keywords internal
#' @noRd
.print_mm_common_header <- function(type_label, n = NULL, p = NULL) {
  cat("MixMashNet fit\n")
  cat("  Type: ", type_label, "\n", sep = "")
  if (!is.null(n) && is.finite(n)) {
    cat("  Data: ", n, " subjects x ", p, " variables\n", sep = "")
  } else if (!is.null(p)) {
    cat("  Variables: ", p, "\n", sep = "")
  }
}

#' Internal helper for print methods
#' @keywords internal
#' @noRd
.print_mm_exclusions <- function(x) {
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
}

#' Internal helper for print methods
#' @keywords internal
#' @noRd
.print_mm_settings <- function(x) {
  s <- x$settings

  lambda_sel <- NULL
  if (!is.null(x$call$lambdaSel)) {
    lambda_sel <- as.character(x$call$lambdaSel)
  }

  if (!is.null(lambda_sel) && length(lambda_sel) == 1L) {
    cat("  Lambda selection: ", lambda_sel, "\n", sep = "")
  }

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
}

#' Internal helper for print methods
#' @keywords internal
#' @noRd
.print_mm_data_info <- function(x) {
  if (is.null(x$data_info)) return(invisible(NULL))

  tl <- x$data_info$mgm_type_level
  printed_header <- FALSE

  if (is.data.frame(tl) && nrow(tl) > 0 && all(c("node", "mgm_type") %in% names(tl))) {
    vars_c <- tl$node[tl$mgm_type == "c"]
    vars_p <- tl$node[tl$mgm_type == "p"]

    if (length(vars_c) || length(vars_p)) {
      cat("  Data info:\n")
      printed_header <- TRUE

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

  brm <- x$data_info$binary_recode_map
  if (!is.null(brm) && length(brm) > 0) {

    fmt_map <- function(m) {
      if (is.atomic(m) && !is.null(names(m))) {
        paste0(names(m), "->", as.character(m), collapse = ", ")
      } else {
        paste(as.character(m), collapse = ", ")
      }
    }

    if (!printed_header) {
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

  invisible(NULL)
}

#' Print method for single layer MixMashNet objects
#'
#' @description
#' Compact textual summary for objects returned by \code{mixMN()}.
#'
#' @param x An object of class \code{"mixMN_fit"}.
#' @param ... Additional arguments.
#'
#' @return The input object \code{x}, returned invisibly.
#' @seealso \code{\link{mixMN}}
#' @method print mixMN_fit
#' @export
print.mixMN_fit <- function(x, ...) {
  if (!inherits(x, "mixMN_fit")) {
    stop("`x` must be a 'mixMN_fit' object.", call. = FALSE)
  }

  p <- length(x$model$nodes)
  n <- x$model$n %||% NA_integer_

  .print_mm_common_header(
    type_label = "Single layer MGM (mixMN)",
    n = n,
    p = p
  )

  if (!is.null(x$graph$igraph)) {
    g <- x$graph$igraph
    cat("  Graph: ", igraph::vcount(g), " nodes, ", igraph::ecount(g), " edges\n", sep = "")
  }

  if (!is.null(x$communities$groups)) {
    k <- length(stats::na.omit(unique(as.integer(x$communities$groups))))
    cat("  Communities: ", k, "\n", sep = "")
  }

  .print_mm_exclusions(x)
  .print_mm_settings(x)
  .print_mm_data_info(x)

  invisible(x)
}

#' Print method for multilayer MixMashNet objects
#'
#' @description
#' Compact textual summary for objects returned by \code{multimixMN()}.
#'
#' @param x An object of class \code{"multimixMN_fit"}.
#' @param ... Additional arguments.
#'
#' @return The input object \code{x}, returned invisibly.
#' @seealso \code{\link{multimixMN}}
#' @method print multimixMN_fit
#' @export
print.multimixMN_fit <- function(x, ...) {
  if (!inherits(x, "multimixMN_fit")) {
    stop("`x` must be a 'multimixMN_fit' object.", call. = FALSE)
  }

  p <- length(x$model$nodes)
  n <- x$model$n %||% NA_integer_

  .print_mm_common_header(
    type_label = "Multilayer MGM (multimixMN)",
    n = n,
    p = p
  )

  if (!is.null(x$layer_fits)) {
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

      edges_df <- fitL$statistics$edge$true %||% NULL
      n_edges_nz <- if (!is.null(edges_df) && nrow(edges_df) > 0 &&
                        "weight" %in% names(edges_df)) {
        sum(edges_df$weight != 0, na.rm = TRUE)
      } else {
        0L
      }

      cat("    - ", L, ": ", n_nodes, " nodes, ", n_edges_nz, " edges\n", sep = "")
    }
  }

  if (!is.null(x$interlayer)) {
    pair_names <- setdiff(names(x$interlayer), "centrality")
    if (length(pair_names)) {
      cat("  Interlayer edges:\n")
      for (pk in pair_names) {
        ed_obj <- x$interlayer[[pk]]$edges
        edges_df <- if (!is.null(ed_obj) && !is.null(ed_obj$true)) ed_obj$true else NULL

        n_edges_nz <- if (!is.null(edges_df) && nrow(edges_df) > 0 &&
                          "weight" %in% names(edges_df)) {
          sum(edges_df$weight != 0, na.rm = TRUE)
        } else {
          0L
        }

        cat("    - ", pk, ": ", n_edges_nz, " edges\n", sep = "")
      }
    }
  }

  if (!is.null(x$graph$igraph)) {
    g <- x$graph$igraph
    cat("  Graph: ", igraph::vcount(g), " nodes, ", igraph::ecount(g), " edges\n", sep = "")
  }

  if (!is.null(x$layer_fits)) {
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

  .print_mm_exclusions(x)
  .print_mm_settings(x)
  .print_mm_data_info(x)

  invisible(x)
}
