#' Plot method for single layer MixMashNet objects
#'
#' @description
#' Plotting interface for objects returned by \code{mixMN()}.
#'
#' Depending on \code{what}, the method can:
#' \itemize{
#'   \item \code{what = "network"}: plot the estimated single-layer network;
#'   \item \code{what = "intra"}: plot node-level metrics or edge weights with
#'         bootstrap quantile regions at the level stored in the object;
#'   \item \code{what = "stability"}: plot node stability within communities
#'         based on bootstrap community assignments.
#' }
#'
#' @param x An object of class \code{"mixMN_fit"}, as returned by
#'   \code{mixMN()}.
#' @param what Type of plot to produce. One of
#'   \code{c("network","intra","stability")}.
#' @param ... Additional arguments. Supported arguments depend on \code{what}:
#'   see the details below.
#'
#' @details
#' \strong{Network plots (\code{what = "network"}):}
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{color_by}}{Node coloring: \code{c("community","none")}.}
#'   \item{\code{edge_color_by}}{Edge coloring: \code{c("sign","none")}.}
#'   \item{\code{edge_scale}}{Numeric scaling factor for edge widths
#'     (multiplied by \code{abs(weight)}).}
#'   \item{\code{graphics::plot.igraph} arguments}{e.g., \code{vertex.size},
#'     \code{vertex.label.cex}, \code{edge.width}, \code{vertex.label.color}, etc.}
#' }
#'
#' \strong{Within-network statistics (\code{what = "intra"}):}
#' Plots node-level metrics or edge weights with bootstrap quantile regions.
#'
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{statistics}}{Character vector of metrics. Options include:
#'     \code{"strength"}, \code{"expected_influence"}, \code{"closeness"},
#'     \code{"betweenness"}, bridge metrics \code{"bridge_strength"},
#'     \code{"bridge_ei1"}, \code{"bridge_ei2"}, \code{"bridge_closeness"},
#'     \code{"bridge_betweenness"}, excluded bridge metrics
#'     \code{"bridge_strength_excluded"}, \code{"bridge_ei1_excluded"},
#'     \code{"bridge_ei2_excluded"}, \code{"bridge_closeness_excluded"},
#'     \code{"bridge_betweenness_excluded"}, and \code{"edges"}.
#'     Different metric families cannot be mixed in the same call
#'     (e.g., \code{"edges"} cannot be combined with node metrics).}
#'   \item{\code{ordering}}{Node ordering:
#'     \code{c("value","alphabetical","community")}.}
#'   \item{\code{standardize}}{Logical; if \code{TRUE}, z-standardize the
#'     displayed values within the panel.}
#'   \item{\code{exclude_nodes}}{Optional character vector of node names to
#'     remove before plotting.}
#'   \item{\code{color_by_community}}{Logical; if \code{TRUE}, color nodes by
#'     community when available.}
#'   \item{\code{edges_top_n}}{Integer; when \code{statistics = "edges"}, keep
#'     the top edges by absolute weight.}
#'   \item{\code{title}}{Optional plot title.}
#' }
#'
#' \strong{Community membership stability (\code{what = "stability"}):}
#' Plots node stability by community.
#'
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{title}}{Plot title. Default:
#'     \code{"Node Stability by Community"}.}
#'   \item{\code{cutoff}}{Optional numeric threshold in \eqn{[0,1]} shown as a
#'     dashed vertical line. Use \code{NULL} to hide the line. Default: 0.7.}
#' }
#'
#' The quantile region level is taken from the fitted object
#' (\code{x$settings$quantile_level}); if missing or invalid, a default of 0.95
#' is used.
#'
#' @return
#' If \code{what != "network"}, the function returns a \code{ggplot} object.
#' If \code{what = "network"}, the network is plotted directly.
#'
#' @seealso \code{\link{mixMN}}
#' @method plot mixMN_fit
#' @export
plot.mixMN_fit <- function(
    x,
    what = c("network", "intra", "stability"),
    ...
) {
  if (!inherits(x, "mixMN_fit")) {
    stop("`x` must be a 'mixMN_fit' object.", call. = FALSE)
  }
  dots <- list(...)

  wants_centrality <- .plot_mm_wants_centrality(dots)

  if (missing(what)) {
    what <- if (wants_centrality) "intra" else "network"
  } else {
    what <- match.arg(what)
  }

  if (what == "network") {
    .plot_network_single(x, ...)
    return(invisible(x))
  }

  if (what == "intra") {
    default_stats  <- .plot_mm_default_stats_intra()
    has_statistics <- "statistics" %in% names(dots)

    if (!has_statistics) {
      message(
        "Using default statistics for intralayer plots: ",
        paste(default_stats, collapse = ", "),
        ". To choose other metrics, use the `statistics` argument."
      )
    }

    args <- dots
    if (!has_statistics) {
      args$statistics <- default_stats
    }
    args$fit <- x

    p <- do.call(plotCentrality, args)
    return(p)
  }

  if (what == "stability") {
    stab <- membershipStab(x)
    p <- membershipStab_plot(stab, ...)
    return(p)
  }

  stop("Unsupported 'what' argument.")
}

#' Plot method for multilayer MixMashNet objects
#'
#' @description
#' Plotting interface for objects returned by \code{multimixMN()}.
#'
#' Depending on \code{what}, the method can:
#' \itemize{
#'   \item \code{what = "network"}: plot the estimated multilayer network, or a
#'         single layer if \code{layer} is specified;
#'   \item \code{what = "intra"}: plot intralayer node-level metrics or edge
#'         weights with bootstrap quantile regions;
#'   \item \code{what = "inter"}: plot interlayer node metrics or interlayer
#'         edge weights with bootstrap quantile regions;
#'   \item \code{what = "stability"}: plot node stability within communities
#'         based on bootstrap community assignments.
#' }
#'
#' @param x An object of class \code{"multimixMN_fit"}, as returned by
#'   \code{multimixMN()}.
#' @param what Type of plot to produce. One of
#'   \code{c("network","intra","inter","stability")}.
#' @param layer Optional layer name. For \code{what = "intra"} or
#'   \code{what = "stability"}, this selects which layer-specific fit to use.
#'   For \code{what = "network"}, if provided, the selected layer is plotted as
#'   a single layer network.
#' @param ... Additional arguments. Supported arguments depend on \code{what}:
#'   see the details below.
#'
#' @details
#' \strong{Network plots (\code{what = "network"}):}
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{color_by}}{Node coloring:
#'     \code{c("layer","community","none")}.}
#'   \item{\code{edge_color_by}}{Edge coloring: \code{c("sign","none")}.}
#'   \item{\code{edge_scale}}{Numeric scaling factor for edge widths
#'     (multiplied by \code{abs(weight)}).}
#'   \item{\code{graphics::plot.igraph} arguments}{e.g., \code{vertex.size},
#'     \code{vertex.label.cex}, \code{edge.width}, \code{vertex.label.color}, etc.}
#' }
#'
#' \strong{Intralayer statistics (\code{what = "intra"}):}
#' Plots node-level metrics or edge weights with bootstrap quantile regions.
#' If \code{layer} is provided, only that layer is plotted. If \code{layer} is
#' \code{NULL}, all layers are plotted, one panel per layer.
#'
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{statistics}}{Character vector of metrics. Options include:
#'     \code{"strength"}, \code{"expected_influence"}, \code{"closeness"},
#'     \code{"betweenness"}, bridge metrics \code{"bridge_strength"},
#'     \code{"bridge_ei1"}, \code{"bridge_ei2"}, \code{"bridge_closeness"},
#'     \code{"bridge_betweenness"}, excluded bridge metrics
#'     \code{"bridge_strength_excluded"}, \code{"bridge_ei1_excluded"},
#'     \code{"bridge_ei2_excluded"}, \code{"bridge_closeness_excluded"},
#'     \code{"bridge_betweenness_excluded"}, and \code{"edges"}.
#'     Different metric families cannot be mixed in the same call
#'     (e.g., \code{"edges"} cannot be combined with node metrics).}
#'   \item{\code{ordering}}{Node ordering:
#'     \code{c("value","alphabetical","community")}.}
#'   \item{\code{standardize}}{Logical; if \code{TRUE}, z-standardize the
#'     displayed values within each panel.}
#'   \item{\code{exclude_nodes}}{Optional character vector of node names to
#'     remove before plotting.}
#'   \item{\code{color_by_community}}{Logical; if \code{TRUE}, color nodes by
#'     community when available.}
#'   \item{\code{edges_top_n}}{Integer; when \code{statistics = "edges"}, keep
#'     the top edges by absolute weight.}
#'   \item{\code{title}}{Optional plot title. If omitted and multiple layers are
#'     shown, layer-specific titles are added automatically.}
#' }
#'
#' \strong{Interlayer summaries (\code{what = "inter"}):}
#' Plots interlayer node metrics or interlayer edge weights with bootstrap
#' quantile regions.
#'
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{statistics}}{Character vector. Node metrics:
#'     \code{c("strength","expected_influence","closeness","betweenness")}, or
#'     \code{"edges"} for interlayer edge weights. Node metrics and
#'     \code{"edges"} cannot be combined.}
#'   \item{\code{pairs}}{Layer pairs to show. Either \code{"*"} (all available)
#'     or a character vector of pair keys like \code{"bio_dis"}
#'     (order-insensitive).}
#'   \item{\code{edges_top_n}}{Integer; keep the top interlayer edges by
#'     absolute weight.}
#'   \item{\code{ordering}}{Ordering within panels:
#'     \code{c("value","alphabetical")}.}
#'   \item{\code{standardize}}{Logical; if \code{TRUE}, z-standardize values
#'     (node metrics by metric, edges by pair).}
#'   \item{\code{exclude_nodes}}{Optional character vector; removes nodes and
#'     incident interlayer edges.}
#'   \item{\code{nodes_layer}}{Optional layer name to restrict node metrics to
#'     nodes belonging to that layer.}
#'   \item{\code{title}}{Optional plot title.}
#' }
#'
#' \strong{Community membership stability (\code{what = "stability"}):}
#' Plots node stability by community. If \code{layer} is provided, only that
#' layer is shown. Otherwise, stability plots are shown for all layers.
#'
#' Supported arguments (via \code{...}):
#' \describe{
#'   \item{\code{title}}{Plot title. Default:
#'     \code{"Node Stability by Community"}.}
#'   \item{\code{cutoff}}{Optional numeric threshold in \eqn{[0,1]} shown as a
#'     dashed vertical line. Use \code{NULL} to hide the line. Default: 0.7.}
#' }
#'
#' The quantile region level is taken from the fitted object
#' (\code{x$settings$quantile_level}); if missing or invalid, a default of 0.95
#' is used.
#'
#' @return
#' If \code{what != "network"}, the function returns a \code{ggplot} object.
#' If \code{what = "network"}, the network is plotted directly.
#'
#' @seealso \code{\link{multimixMN}}
#' @method plot multimixMN_fit
#' @export
plot.multimixMN_fit <- function(
    x,
    what = c("network", "intra", "inter", "stability"),
    layer = NULL,
    ...
) {
  if (!inherits(x, "multimixMN_fit")) {
    stop("`x` must be a 'multimixMN_fit' object.", call. = FALSE)
  }
  dots <- list(...)

  wants_centrality <- .plot_mm_wants_centrality(dots)

  if (missing(what)) {
    if (wants_centrality && is.null(layer)) {
      stop(
        "You are plotting statistics on a multilayer object.\n",
        "Please specify one of:\n",
        "  - layer = \"bio\"   # intralayer statistics for a specific layer\n",
        "  - what  = \"intra\" # intralayer statistics for ALL layers\n",
        "  - what  = \"inter\" # interlayer statistics\n"
      )
    } else if (!is.null(layer) && wants_centrality) {
      what <- "intra"
    } else {
      what <- "network"
    }
  } else {
    what <- match.arg(what)
  }

  if (what == "network") {
    if (!is.null(layer)) {
      subfit <- .plot_mm_get_layer_fit(x, layer)
      .plot_network_single(subfit, ...)
    } else {
      .plot_network_multi(x, ...)
    }
    return(invisible(x))
  }

  if (what == "intra") {
    default_stats  <- .plot_mm_default_stats_intra()
    has_statistics <- "statistics" %in% names(dots)

    if (!has_statistics) {
      message(
        "Using default statistics for intralayer plots: ",
        paste(default_stats, collapse = ", "),
        ". To choose other metrics, use the `statistics` argument."
      )
    }

    if (!is.null(layer)) {
      subfit <- .plot_mm_get_layer_fit(x, layer)
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

  if (what == "inter") {
    default_stats  <- .plot_mm_default_stats_inter()
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
    if (!is.null(layer)) {
      subfit <- .plot_mm_get_layer_fit(x, layer)
      stab <- membershipStab(subfit)
      p <- membershipStab_plot(stab, ...)
      return(p)
    }

    if (is.null(x$layer_fits) || !length(x$layer_fits)) {
      stop("No 'layer_fits' found in this 'multimixMN_fit' object.")
    }

    layer_names <- names(x$layer_fits)

    plots <- lapply(layer_names, function(L) {
      subfit_L <- x$layer_fits[[L]]
      stab_L   <- membershipStab(subfit_L)
      membershipStab_plot(stab_L, title = paste0("Layer: ", L), ...)
    })

    p <- patchwork::wrap_plots(plots, ncol = length(plots))
    return(p)
  }

  stop("Unsupported 'what' argument.")
}
