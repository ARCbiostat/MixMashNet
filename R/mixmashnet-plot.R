#' Fruchterman–Reingold layout for multilayer graphs
#'
#' @description
#' Internal helper that computes 2D coordinates for a multilayer \pkg{igraph}
#' object by:
#' \itemize{
#'   \item running a weighted Fruchterman–Reingold layout separately within
#'         each layer (using intra-layer edges only);
#'   \item normalizing the scale of each layer;
#'   \item placing layer centroids on a larger circle so that layers are
#'         visually separated but interlayer edges remain readable.
#' }
#'
#' @param g An \code{igraph} object with a vertex attribute containing layer
#'   membership (see \code{layer_attr}).
#' @param layer_attr Character string; name of the vertex attribute that stores
#'   the layer label (default \code{"layer"}).
#'
#' @return A numeric matrix with two columns (x, y) and one row per vertex.
#'   The row order matches \code{V(g)}.
#'
#' @keywords internal
#' @importFrom igraph is_igraph V vcount induced_subgraph ecount E layout_with_fr
#' @noRd
layout_multilayer_fr <- function(g, layer_attr = "layer") {
  stopifnot(igraph::is_igraph(g))
  layer <- as.factor(igraph::V(g)$layer)

  if (is.null(layer)) stop("Vertex attribute '", layer_attr, "' not found.")

  L <- nlevels(layer)
  n         <- igraph::vcount(g)

  coords <- matrix(NA_real_, nrow = n, ncol = 2)

  # 1) layout FR separato per ciascun layer (solo archi intra-layer)
  for (i in seq_len(L)) {
    idx <- which(layer == levels(layer)[i])
    if (length(idx) == 0L) next

    subg <- igraph::induced_subgraph(g, vids = idx)

    # pesi positivi per il layout
    if (igraph::ecount(subg) > 0) {
      w <- abs(igraph::E(subg)$weight)
      w[is.na(w)] <- 0
      if (all(w == 0)) w <- rep(1, length(w))
      xy <- igraph::layout_with_fr(subg, weights = w, niter = 300)
    } else {
      # se il layer è totalmente scollegato → un piccolo cerchio
      k  <- length(idx)
      ang <- seq(0, 2*pi, length.out = k + 1L)[-1L]
      xy <- cbind(cos(ang), sin(ang))
    }

    # normalizzo dimensione del layer
    xy <- scale(xy)  # ~ raggio 1
    coords[idx, ] <- xy
  }

  # 2) posiziono i layer su un cerchio più grande
  angles_layer <- seq(0, 2*pi, length.out = L + 1L)[-1L]
  radius_layer <- 4  # distanza tra i centri dei layer

  for (i in seq_len(L)) {
    idx <- which(layer == levels(layer)[i])
    if (length(idx) == 0L) next

    center_now <- colMeans(coords[idx, , drop = FALSE])
    target     <- c(radius_layer * cos(angles_layer[i]),
                    radius_layer * sin(angles_layer[i]))
    shift      <- target - center_now
    coords[idx, ] <- sweep(coords[idx, , drop = FALSE], 2, shift, "+")
  }

  coords
}


# -------------------------------------------------------------------
# Internal network plotters (single-layer and multilayer)
# -------------------------------------------------------------------

.plot_network_single <- function(
    x,
    color_by = c("community", "none"),
    edge_color_by = c("sign", "none"),
    edge_scale = 4,
    vertex_size = 12,
    vertex_label_cex = 0.8,
    layout = c("fr"),
    ...
) {
  color_by     <- match.arg(color_by)
  edge_color_by <- match.arg(edge_color_by)
  layout       <- match.arg(layout)

  # --- recupero grafo ---
  g <- x$graph$igraph
  if (is.null(g)) stop("No igraph object found in x$graph$igraph.")

  vnames <- igraph::V(g)$name

  # --- edge weights (assicurati che abs_weight esista) ---
  if (is.null(igraph::E(g)$abs_weight)) {
    igraph::E(g)$abs_weight <- abs(igraph::E(g)$weight)
  }

  # --- layout: default fruchterman-reingold pesato ---
  lay <- igraph::layout_with_fr(g, weights = igraph::E(g)$abs_weight)

  # ======================
  #  Vertex colors
  # ======================

  if (color_by == "none") {
    vcol <- rep("skyblue", length(vnames))
  } else {
    # color_by = "community"
    memb    <- x$communities$groups   # named vector
    palette <- x$communities$palette  # named palette

    memb <- memb[vnames]  # allineo

    # default: esclusi in grigio
    vcol <- rep("grey80", length(vnames))

    idx <- !is.na(memb)
    if (any(idx)) {
      pal <- palette[as.character(memb[idx])]
      pal[is.na(pal)] <- "orange"   # fallback in caso di mismatch
      vcol[idx] <- pal
    }
  }

  # ======================
  #  Edge colors & widths
  # ======================

  ew <- edge_scale * igraph::E(g)$abs_weight

  if (edge_color_by == "sign") {
    ecol <- ifelse(igraph::E(g)$weight > 0, "darkgreen", "red")
  } else {
    ecol <- rep("grey40", igraph::ecount(g))
  }

  # ======================
  #  Drawing
  # ======================

  graphics::plot(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size  = vertex_size,
    vertex.label = vnames,
    vertex.label.cex = vertex_label_cex,
    edge.width  = ew,
    edge.color  = ecol,
    ...
  )

  invisible(x)
}


.plot_network_multi <- function(
    x,
    color_by = c("layer", "community", "none"),
    edge_color_by = c("sign", "none"),
    edge_scale = 4,
    vertex_size = 10,
    vertex_label_cex = 0.7,
    layout = c("multilayer_fr", "fr"),
    ...
) {
  color_by      <- match.arg(color_by)
  edge_color_by <- match.arg(edge_color_by)
  layout        <- match.arg(layout)

  g <- x$graph$igraph
  if (is.null(g)) {
    stop("No igraph object found in x$graph$igraph.")
  }

  vnames <- igraph::V(g)$name

  # edge weights / abs_weight
  if (is.null(igraph::E(g)$abs_weight)) {
    igraph::E(g)$abs_weight <- abs(igraph::E(g)$weight)
  }

  # layout
  lay <- switch(
    layout,
    multilayer_fr = layout_multilayer_fr(g, layer_attr = "layer"),
    fr            = igraph::layout_with_fr(g, weights = igraph::E(g)$abs_weight)
  )

  # ---------- vertex colors ----------
  vcol <- rep("grey80", igraph::vcount(g))
  vnames <- igraph::V(g)$name
  names(vcol) <- vnames

  if (color_by == "none") {
    vcol[] <- "skyblue"

  } else if (color_by == "layer") {
    layers_vec <- as.character(igraph::V(g)$layer)
    ulay       <- sort(unique(layers_vec[!is.na(layers_vec)]))
    if (length(ulay) > 0) {
      pal <- colorspace::qualitative_hcl(length(ulay), palette = "Dark 3")
      names(pal) <- ulay
      idx <- !is.na(layers_vec)
      vcol[idx] <- pal[layers_vec[idx]]
    }

  } else if (color_by == "community") {
    layers_vec <- as.character(igraph::V(g)$layer)
    memb_vec   <- igraph::V(g)$membership   # numeri di cluster

    for (i in seq_along(vnames)) {
      node  <- vnames[i]
      L     <- layers_vec[i]
      cl_id <- memb_vec[i]

      if (is.na(L) || is.na(cl_id)) next

      # palette del layer corrispondente
      fitL <- x$layer_fits[[L]]
      if (is.null(fitL)) next

      palL <- fitL$communities$palette  # named vector cluster -> colore
      if (is.null(palL)) next

      col_i <- palL[as.character(cl_id)]
      if (!is.na(col_i)) {
        vcol[i] <- col_i
      }
    }
  }

  # ---------- edge colors & width ----------
  ew <- edge_scale * igraph::E(g)$abs_weight

  if (edge_color_by == "sign") {
    ecol <- ifelse(igraph::E(g)$weight > 0, "darkgreen", "red")
  } else {
    ecol <- rep("grey40", igraph::ecount(g))
  }

  graphics::plot(
    g,
    layout = lay,
    vertex.color = vcol,
    vertex.size  = vertex_size,
    vertex.label = vnames,
    vertex.label.cex = vertex_label_cex,
    edge.width  = ew,
    edge.color  = ecol,
    ...
  )

  invisible(x)
}



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
#' @import patchwork
#' @export
plot.mixmashnet <- function(
    x,
    what  = c("network", "intra",
              "inter",
              "stability"),
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
    "edges_top_n", "exclude_nodes", "color_by_community", "title"
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
        "  - layer = \" \"                 # intra-layer statistics for a specific layer\n",
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
  get_layer_fit <- function(obj, layer_name) {
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
        subfit <- get_layer_fit(x, layer)
        .plot_network_single(subfit, layer = layer, ...)
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
        subfit <- get_layer_fit(x, layer)
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

    dots <- list(...)

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
        subfit <- get_layer_fit(x, layer)
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
