#' Plot method for MixMashNet objects
#'
#' @description
#' Unified plotting interface for objects returned by \code{mixMN()} and
#' \code{multimixMN()}. Depending on \code{what}, it can:
#' \itemize{
#'   \item \code{what = "network"}: plot the estimated network
#'         (single-layer or multilayer);
#'   \item \code{what = "intra"}: plot intra-layer node/edge statistics with
#'         95\% bootstrap CIs (centrality and bridge metrics);
#'   \item \code{what = "inter"}: plot interlayer node metrics or interlayer
#'         edge weights with 95\% bootstrap CIs (multilayer only), via
#'         \code{plotInterlayer()} and the chosen \code{statistics};
#'   \item \code{what = "stability"}: plot node stability within communities
#'         based on bootstrap community assignments.
#' }
#'
#' @param x An object of class \code{mixmashnet}, as returned by
#'   \code{mixMN()} or \code{multimixMN()}.
#' @param what Type of plot to produce. One of
#'   \code{c("network","intra","inter","stability")}.
#'   If missing, a default is chosen based on the presence of
#'   centrality-related arguments and on whether \code{x} is single-layer
#'   or multilayer (see Details).
#' @param layer Optional layer name. For \code{what = "intra"} or
#'   \code{what = "stability"} on a \code{multimixMN_fit} object, this selects
#'   which layer-specific fit to use. If \code{NULL}, the behaviour depends
#'   on \code{what}:
#'   \itemize{
#'     \item \code{what = "network"}: plots the global multilayer network;
#'     \item \code{what = "intra"} or \code{"stability"} on multilayer:
#'           either all layers are plotted in a combined layout or an error
#'           is raised if the context is ambiguous.
#'   }
#' @param ... Additional arguments passed to the underlying helpers:
#'   \itemize{
#'     \item For \code{what = "intra"}: forwarded to \code{plotCentrality()},
#'           e.g. \code{statistics}, \code{ordering}, \code{standardize},
#'           \code{edges_top_n}, \code{exclude_nodes}, \code{color_by_community}.
#'     \item For \code{what = "inter"}: forwarded to \code{plotInterlayer()},
#'           e.g. \code{statistics}, \code{pairs}, \code{edges_top_n},
#'           \code{ordering}, \code{standardize}, \code{nodes_layer}.
#'     \item For \code{what = "stability"}: forwarded to
#'           \code{membershipStab_plot()}, e.g. \code{title}.
#'     \item For \code{what = "network"}: forwarded to internal network
#'           plotting helpers \code{.plot_network_single()} or
#'           \code{.plot_network_multi()}, e.g. \code{color_by},
#'           \code{edge_color_by}, \code{vertex_size}.
#'   }
#'
#' @details
#' When \code{what} is missing, the function inspects \code{...} to decide
#' whether the user is requesting centrality/edge statistics (e.g., by
#' passing \code{statistics}, \code{ordering}, etc.). In that case:
#' \itemize{
#'   \item for single-layer fits, the default is \code{what = "intra"};
#'   \item for multilayer fits, an informative error is raised if neither
#'         \code{what} nor \code{layer} is specified and the request is
#'         ambiguous (intra-layer vs interlayer statistics).
#' }
#'
#' @return
#' For \code{what != "network"}, a \code{ggplot} object is returned.
#' For \code{what = "network"}, the corresponding network plotting helper
#' is called for its side-effect and \code{x} is returned invisibly.
#'
#' @export
plot.mixmashnet <- function(
    x,
    what  = c("network","intra","inter","stability"),
    layer = NULL,
    ...
) {

  cls       <- class(x)
  is_multi  <- "multimixMN_fit" %in% cls
  is_single <- "mixMN_fit" %in% cls && !is_multi

  if (!is_multi && !is_single) {
    stop("`x` must be an object of class 'mixmashnet' (from mixMN() or multimixMN()).")
  }

  dots <- list(...)

  centrality_args <- c(
    "statistics", "ordering", "standardize",
    "edges_top_n", "exclude_nodes", "color_by_community"
  )
  wants_centrality <- any(names(dots) %in% centrality_args)

  # --- logica per what ---
  if (missing(what)) {
    if (is_single && wants_centrality) {
      # single-layer + statistics/order/etc → intra
      what <- "intra"

    } else if (is_multi && wants_centrality && is.null(layer)) {
      # multilayer + statistics ma SENZA layer e SENZA what esplicito → AMBIGUO
      stop(
        "You are plotting statistics on a multilayer object.\n",
        "Please specify one of:\n",
        "  - layer = \"bio\"   # intra-layer statistics for a specific layer\n",
        "  - what  = \"intra\" # intra-layer statistics for ALL layers\n",
        "  - what  = \"inter\" # interlayer statistics\n"
      )

    } else if (is_multi && !is.null(layer) && wants_centrality) {
      # multilayer + layer specificato + statistics/etc → intra
      what <- "intra"

    } else {
      # default generale: network
      what <- "network"
    }
  } else {
    # utente ha specificato what esplicitamente
    what <- match.arg(what)
  }

  ## --- helper per gestire layer in oggetti multilayer ---
  .get_layer_fit <- function(obj, layer_name) {
    if (is.null(obj$layer_fits) || !length(obj$layer_fits)) {
      stop("No 'layer_fits' found in this 'multimixMN_fit' object.")
    }
    if (is.null(layer_name)) {
      stop(
        "This is a 'multimixMN_fit' object. Please specify the ",
        "'layer' argument (e.g., layer = \"bio\")."
      )
    }
    if (!layer_name %in% names(obj$layer_fits)) {
      stop(
        "Layer '", layer_name, "' not found. Available layers: ",
        paste(names(obj$layer_fits), collapse = ", ")
      )
    }
    obj$layer_fits[[layer_name]]
  }

  ## =================== SWITCH SULLA MODALITÀ ===================
  if (what == "network") {
    # qui userai i tuoi helper per plottare la rete
    if (is_single) {
      # es: .plot_network_single(x, ...)
      .plot_network_single(x, ...)
    } else {
      # se layer specificato, puoi decidere se plottare solo quel layer:
      if (!is.null(layer)) {
        subfit <- .get_layer_fit(x, layer)
        .plot_network_single(subfit, ...)
      } else {
        # plot multilayer globale
        .plot_network_multi(x, ...)
      }
    }
    return(invisible(x))
  }

  if (what == "intra") {
    # default per intra-layer se l'utente non specifica statistics
    default_stats   <- c("strength", "expected_influence", "closeness", "betweenness")
    has_statistics  <- "statistics" %in% names(dots)

    if (!has_statistics) {
      message(
        "Using default statistics for intra-layer plots: ",
        paste(default_stats, collapse = ", "),
        ". To choose other metrics, use the `statistics` argument."
      )
    }

    if (is_single) {
      # ---- single-layer ----
      args <- dots
      if (!has_statistics) {
        args$statistics <- default_stats
      }
      args$fit <- x
      p <- do.call(plotCentrality, args)
      return(p)

    } else {
      # ---- multilayer ----
      if (!is.null(layer)) {
        # caso 1: multilayer + layer specificato → intra-layer per quel layer
        subfit <- .get_layer_fit(x, layer)
        args <- dots
        if (!has_statistics) {
          args$statistics <- default_stats
        }
        args$fit <- subfit
        p <- do.call(plotCentrality, args)
        return(p)
      }

      # caso 2: multilayer + nessun layer specificato → tutti i layer uno sotto l'altro
      if (is.null(x$layer_fits) || !length(x$layer_fits)) {
        stop("No 'layer_fits' found in this 'multimixMN_fit' object.")
      }

      layer_names <- names(x$layer_fits)

      plots <- lapply(layer_names, function(L) {
        subfit_L <- x$layer_fits[[L]]
        args_L <- dots
        if (!has_statistics) {
          args_L$statistics <- default_stats
        }
        args_L$fit <- subfit_L
        # titolo per ogni pannello, ma solo se l'utente non ha già messo un title
        if (!"title" %in% names(args_L)) {
          args_L$title <- paste0("Layer: ", L)
        }
        do.call(plotCentrality, args_L)
      })

      p <- patchwork::wrap_plots(plots, ncol = 1)
      return(p)
    }
  }

  if (what == "inter") {
    if (!is_multi) {
      stop("Interlayer node statistics are only available for 'multimixMN_fit' objects.")
    }

    default_stats  <- c("strength", "expected_influence", "closeness", "betweenness")
    has_statistics <- "statistics" %in% names(dots)

    if (!has_statistics) {
      message(
        "Using default statistics for interlayer plots: ",
        paste(default_stats, collapse = ", "),
        ". To choose other metrics, use the `statistics` argument."
      )
    }

    args <- dots
    if (!has_statistics) {
      args$statistics <- default_stats
    }

    if (!is.null(layer) && !"nodes_layer" %in% names(args)) {
      args$nodes_layer <- layer
    }

    args$fit_multi <- x

    p <- do.call(plotInterlayer, args)
    return(p)
  }

  if (what == "stability") {
    if (is_single) {
      stab <- membershipStab(x, IS.plot = FALSE)
      p <- membershipStab_plot(stab, ...)
      return(p)

    } else {
      # multilayer
      if (!is.null(layer)) {
        subfit <- .get_layer_fit(x, layer)
        stab <- membershipStab(subfit, IS.plot = FALSE)
        p <- membershipStab_plot(stab, ...)
        return(p)
      }

      # multilayer senza layer → stability per OGNI layer
      if (is.null(x$layer_fits) || !length(x$layer_fits)) {
        stop("No 'layer_fits' found in this 'multimixMN_fit' object.")
      }

      layer_names <- names(x$layer_fits)

      plots <- lapply(layer_names, function(L) {
        subfit_L <- x$layer_fits[[L]]
        stab_L   <- membershipStab(subfit_L, IS.plot = FALSE)
        membershipStab_plot(stab_L, title = paste0("Layer: ", L), ...)
      })

      # Composizione unica con patchwork (ora sempre disponibile)
      p <- patchwork::wrap_plots(plots, ncol = length(plots))
      return(p)
    }
  }

  # in teoria non ci arrivi mai:
  stop("Unsupported 'what' argument.")
}
