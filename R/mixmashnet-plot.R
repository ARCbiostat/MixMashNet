#' Plot method for MixMashNet objects
#'
#' @description
#' Unified plotting interface for objects returned by \code{mixMN()} and
#' \code{multimixMN()}. Depending on \code{what}, it can:
#' \itemize{
#'   \item \code{what = "network"}: plot the estimated network
#'         (single layer or multilayer);
#'   \item \code{what = "intra"}: plot intralayer node/edge statistics with
#'         bootstrap CIs at the level stored in the object (centrality and bridge metrics);
#'   \item \code{what = "inter"}: plot interlayer node metrics or interlayer
#'         edge weights with bootstrap CIs at the level stored in the object (multilayer only),
#'         and the chosen \code{statistics};
#'   \item \code{what = "stability"}: plot node stability within communities
#'         based on bootstrap community assignments.
#' }
#'
#' @param x An object of class \code{mixmashnet}, as returned by
#'   \code{mixMN()} or \code{multimixMN()}.
#' @param what Type of plot to produce. One of
#'   \code{c("network","intra","inter","stability")}.
#' @param layer Optional layer name. For \code{what = "intra"} or
#'   \code{what = "stability"} on a \code{multimixMN_fit} object, this selects
#'   which layer-specific fit to use.
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
#' @return
#' If \code{what != "network"}, the function returns a \code{ggplot} object.
#' If \code{what = "network"}, the network is plotted directly.
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
      # single layer + statistics/order/etc → intra
      what <- "intra"

    } else if (is_multi && wants_centrality && is.null(layer)) {
      stop(
        "You are plotting statistics on a multilayer object.\n",
        "Please specify one of:\n",
        "  - layer = \"bio\"   # intralayer statistics for a specific layer\n",
        "  - what  = \"intra\" # intralayer statistics for ALL layers\n",
        "  - what  = \"inter\" # interlayer statistics\n"
      )

    } else if (is_multi && !is.null(layer) && wants_centrality) {
      what <- "intra"

    } else {
      what <- "network"
    }
  } else {
    what <- match.arg(what)
  }

  ## --- helper for layers in multilayer objects ---
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

  ## =================== SWITCH ON MODALITY ===================
  if (what == "network") {
    if (is_single) {
      .plot_network_single(x, ...)
    } else {
      if (!is.null(layer)) {
        subfit <- .get_layer_fit(x, layer)
        .plot_network_single(subfit, ...)
      } else {
        .plot_network_multi(x, ...)
      }
    }
    return(invisible(x))
  }

  if (what == "intra") {
    default_stats   <- c("strength", "expected_influence", "closeness", "betweenness")
    has_statistics  <- "statistics" %in% names(dots)

    if (!has_statistics) {
      message(
        "Using default statistics for intralayer plots: ",
        paste(default_stats, collapse = ", "),
        ". To choose other metrics, use the `statistics` argument."
      )
    }

    if (is_single) {
      # ---- single layer ----
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
        subfit <- .get_layer_fit(x, layer)
        args <- dots
        if (!has_statistics) {
          args$statistics <- default_stats
        }
        args$fit <- subfit
        p <- do.call(plotCentrality, args)
        return(p)
      }

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

      if (is.null(x$layer_fits) || !length(x$layer_fits)) {
        stop("No 'layer_fits' found in this 'multimixMN_fit' object.")
      }

      layer_names <- names(x$layer_fits)

      plots <- lapply(layer_names, function(L) {
        subfit_L <- x$layer_fits[[L]]
        stab_L   <- membershipStab(subfit_L, IS.plot = FALSE)
        membershipStab_plot(stab_L, title = paste0("Layer: ", L), ...)
      })

      p <- patchwork::wrap_plots(plots, ncol = length(plots))
      return(p)
    }
  }

  stop("Unsupported 'what' argument.")
}
