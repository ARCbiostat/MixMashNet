#' @keywords internal
#' @noRd
.get_abs_weights <- function(g) {
  w <- igraph::E(g)$weight
  if (is.null(w)) w <- rep(1, igraph::ecount(g)) else w <- abs(w)
  w[is.na(w) | w == 0] <- 1
  w
}

# ============================================================
# Network plotting helpers (igraph-based)
# ============================================================

#' @keywords internal
#' @noRd
.plot_network_single <- function(
    x,
    color_by = c("community", "none"),
    edge_color_by = c("sign", "none"),
    edge_scale = 4,
    layout_type = c("fr"),
    ...
) {
  color_by      <- match.arg(color_by)
  edge_color_by <- match.arg(edge_color_by)
  layout_type   <- match.arg(layout_type)

  g <- x$graph$igraph
  if (is.null(g)) stop("No igraph object found in x$graph$igraph.")
  vnames <- igraph::V(g)$name

  abs_w <- .get_abs_weights(g)

  # vertex colors
  if (color_by == "none") {
    vcol <- rep("skyblue", length(vnames))
  } else {
    memb <- x$communities$groups
    if (!is.null(memb) && !is.null(names(memb))) {
      memb <- memb[vnames]
    } else {
      warning("Community membership has no names; using grey80 for all nodes.")
      memb <- rep(NA_integer_, length(vnames))
    }
    palette <- x$communities$palette
    vcol <- rep("grey80", length(vnames))
    idx <- !is.na(memb)
    if (any(idx)) {
      pal <- palette[as.character(memb[idx])]
      pal[is.na(pal)] <- "orange"
      vcol[idx] <- pal
    }
  }

  # edge colors
  ecol <- if (edge_color_by == "sign") {
    w0 <- igraph::E(g)$weight
    if (is.null(w0)) w0 <- rep(1, igraph::ecount(g))
    ifelse(w0 > 0, "blue", "red")
  } else {
    rep("grey40", igraph::ecount(g))
  }

  dots <- list(...)
  if (is.null(dots$layout)) {
    lay <- switch(
      layout_type,
      fr = igraph::layout_with_fr(g, weights = abs_w)
    )
    dots$layout <- lay
  }
  if (is.null(dots$vertex.color))       dots$vertex.color <- vcol
  if (is.null(dots$vertex.label))       dots$vertex.label <- vnames
  if (is.null(dots$vertex.size))        dots$vertex.size <- 12
  if (is.null(dots$vertex.label.cex))   dots$vertex.label.cex <- 0.8
  if (is.null(dots$vertex.label.color)) dots$vertex.label.color <- "black"
  if (is.null(dots$edge.color))         dots$edge.color <- ecol
  if (is.null(dots$edge.width))         dots$edge.width <- edge_scale * abs_w

  do.call(graphics::plot, c(list(x = g), dots))
  invisible(x)
}

#' @keywords internal
#' @noRd
.plot_network_multi <- function(
    x,
    color_by = c("layer", "community", "none"),
    edge_color_by = c("sign", "none"),
    edge_scale = 4,
    layout_type = c("multilayer_fr", "fr"),
    layer_attr = "layer",
    ...
) {
  color_by      <- match.arg(color_by)
  edge_color_by <- match.arg(edge_color_by)
  layout_type   <- match.arg(layout_type)

  g <- x$graph$igraph
  if (is.null(g)) stop("No igraph object found in x$graph$igraph.")

  vnames <- igraph::V(g)$name

  # --- local abs weights (no side effects on g) ---
  abs_w <- .get_abs_weights(g)

  # ======================
  # Vertex colors
  # ======================
  vcol <- rep("grey80", igraph::vcount(g))

  if (color_by == "none") {
    vcol[] <- "skyblue"

  } else if (color_by == "layer") {
    layers_vec <- as.character(igraph::vertex_attr(g, layer_attr))
    ulay <- sort(unique(layers_vec[!is.na(layers_vec)]))
    if (length(ulay) > 0) {
      pal <- colorspace::qualitative_hcl(length(ulay), palette = "Dynamic")
      names(pal) <- ulay
      idx <- !is.na(layers_vec)
      vcol[idx] <- pal[layers_vec[idx]]
    }

  } else if (color_by == "community") {
    # assumes membership stored as vertex attribute 'membership'
    layers_vec <- as.character(igraph::vertex_attr(g, layer_attr))
    memb_vec   <- igraph::vertex_attr(g, "membership")

    # fallback if missing membership
    if (is.null(memb_vec)) {
      warning("Vertex attribute 'membership' not found; falling back to color_by = 'layer'.")
      layers_vec <- as.character(igraph::vertex_attr(g, layer_attr))
      ulay <- sort(unique(layers_vec[!is.na(layers_vec)]))
      if (length(ulay) > 0) {
        pal <- colorspace::qualitative_hcl(length(ulay), palette = "Dynamic")
        names(pal) <- ulay
        idx <- !is.na(layers_vec)
        vcol[idx] <- pal[layers_vec[idx]]
      }
    } else {
      for (i in seq_along(vnames)) {
        L     <- layers_vec[i]
        cl_id <- memb_vec[i]

        if (is.na(L) || is.na(cl_id)) next
        fitL <- x$layer_fits[[L]]
        if (is.null(fitL)) next

        palL <- fitL$communities$palette
        if (is.null(palL)) next

        col_i <- palL[as.character(cl_id)]
        if (!is.na(col_i)) vcol[i] <- col_i
      }
    }
  }

  # ======================
  # Edge colors
  # ======================
  ecol <- if (edge_color_by == "sign") {
    w0 <- igraph::E(g)$weight
    if (is.null(w0)) w0 <- rep(1, igraph::ecount(g))
    ifelse(w0 > 0, "blue", "red")
  } else {
    rep("grey40", igraph::ecount(g))
  }


  # ======================
  # Forward igraph args
  # ======================
  dots <- list(...)

  # set defaults only when absent (igraph-first behaviour)
  if (is.null(dots$layout)) {
    lay <- switch(
      layout_type,
      multilayer_fr = .layout_multilayer_fr(g, layer_attr = layer_attr),
      fr            = igraph::layout_with_fr(g, weights = abs_w)
    )
    dots$layout <- lay
  }
  if (is.null(dots$vertex.color))       dots$vertex.color <- vcol
  if (is.null(dots$vertex.label))       dots$vertex.label <- vnames
  if (is.null(dots$vertex.size))        dots$vertex.size <- 10
  if (is.null(dots$vertex.label.cex))   dots$vertex.label.cex <- 0.7
  if (is.null(dots$vertex.label.color)) dots$vertex.label.color <- "black"
  if (is.null(dots$edge.color))         dots$edge.color <- ecol
  if (is.null(dots$edge.width))         dots$edge.width <- edge_scale * abs_w

  do.call(graphics::plot, c(list(x = g), dots))
  invisible(x)
}
