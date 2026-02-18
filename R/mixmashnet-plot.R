#' Plot method for MixMashNet objects
#'
#' @description
#' Unified plotting interface for objects returned by \code{mixMN()} and
#' \code{multimixMN()}. Depending on \code{what}, it can:
#' \itemize{
#'   \item \code{what = "network"}: plot the estimated network
#'         (single layer or multilayer);
#'   \item \code{what = "intra"}: plot intralayer node/edge statistics with
#'         bootstrap quantile regions at the level stored in the object (centrality and bridge metrics);
#'   \item \code{what = "inter"}: plot interlayer node metrics or interlayer
#'         edge weights with bootstrap quantile regions at the level stored in the object (multilayer only),
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
#' @param ... Additional arguments. Supported arguments depend on \code{what}:
#'   see the details below.
#' @details
#' \strong{Network plots (\code{what = "network"}):}
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{color_by}}{Node coloring. Single layer: \code{c("community","none")}.
#'     Multilayer: \code{c("layer","community","none")}}.
#'   \item{\code{edge_color_by}}{Edge coloring: \code{c("sign","none")}.}
#'   \item{\code{edge_scale}}{Numeric scaling factor for edge widths (multiplied by \code{abs(weight)}).}
#'   \item{\code{graphics::plot.igraph} arguments}{e.g., \code{vertex.size},
#'     \code{vertex.label.cex}, \code{edge.width}, \code{vertex.label.color}, etc.}
#' }
#'
#' \strong{Intralayer statistics (\code{what = "intra"}):}
#' Plots node-level metrics or edge weights with bootstrap quantile regions.
#' For multilayer objects:
#' \itemize{
#'   \item if \code{layer} is provided, plots that layer only;
#'   \item if \code{layer} is \code{NULL}, plots all layers (one panel per layer).
#' }
#'
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{statistics}}{Character vector of metrics. Options include:
#'     \code{"strength"}, \code{"expected_influence"}, \code{"closeness"}, \code{"betweenness"},
#'     bridge metrics \code{"bridge_strength"}, \code{"bridge_ei1"}, \code{"bridge_ei2"},
#'     \code{"bridge_closeness"}, \code{"bridge_betweenness"},
#'     excluded bridge metrics \code{"bridge_strength_excluded"}, \code{"bridge_ei1_excluded"},
#'     \code{"bridge_ei2_excluded"}, \code{"bridge_closeness_excluded"}, \code{"bridge_betweenness_excluded"},
#'     and \code{"edges"}.
#'     Note: different metric families cannot be mixed in the same call (e.g., \code{"edges"} cannot
#'     be combined with node metrics).}
#'   \item{\code{ordering}}{Node ordering: \code{c("value","alphabetical","community")}.}
#'   \item{\code{standardize}}{Logical; if \code{TRUE}, z-standardize the displayed values (within each panel).}
#'   \item{\code{exclude_nodes}}{Optional character vector of node names to remove before plotting.}
#'   \item{\code{color_by_community}}{Logical; if \code{TRUE}, color nodes by community (when available).}
#'   \item{\code{edges_top_n}}{Integer; when \code{statistics = "edges"}, keep the top edges by
#'     absolute weight.}
#'   \item{\code{title}}{Optional plot title. In multilayer mode, if not provided a layer-specific
#'     title is added automatically.}
#' }
#'
#' \strong{Interlayer summaries (\code{what = "inter"}; multilayer only):}
#' Plots interlayer node metrics or interlayer edge weights with bootstrap quantile regions.
#'
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{statistics}}{Character vector. Node metrics:
#'     \code{c("strength","expected_influence","closeness","betweenness")}, or \code{"edges"} for
#'     interlayer edge weights. Node metrics and \code{"edges"} cannot be combined.}
#'   \item{\code{pairs}}{Layer pairs to show. Either \code{"*"} (all available) or a character vector
#'     of pair keys like \code{"bio_dis"} (order-insensitive).}
#'   \item{\code{edges_top_n}}{Integer; keep the top interlayer edges by absolute weight.}
#'   \item{\code{ordering}}{Ordering within panels: \code{c("value","alphabetical")}.}
#'   \item{\code{standardize}}{Logical; if \code{TRUE}, z-standardize values (node metrics by metric,
#'     edges by pair).}
#'   \item{\code{exclude_nodes}}{Optional character vector; removes nodes (and incident interlayer edges).}
#'   \item{\code{nodes_layer}}{Optional layer name to restrict node metrics to nodes belonging to that layer.}
#'   \item{\code{title}}{Optional plot title.}
#' }
#'
#' \strong{Community membership stability (\code{what = "stability"}):}
#' Plots node stability by community.
#' For multilayer objects, \code{layer} selects a specific layer; if \code{layer} is \code{NULL},
#' stability plots are shown for all layers.
#'
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{title}}{Plot title. Default: \code{"Node Stability by Community"}.}
#'   \item{\code{cutoff}}{Optional numeric threshold in \eqn{[0,1]} shown as a dashed vertical line.
#'     Use \code{NULL} to hide the line. Default: 0.7.}
#' }
#'
#' The quantile region level is taken from the fitted object (\code{x$settings$quantile_level});
#' if missing or invalid, a default of 0.95 is used.
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
