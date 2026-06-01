#' @keywords internal
#' @noRd
.resolve_layout_single <- function(g, layout = "fr", weights = NULL) {

  if (is.null(layout)) layout <- "fr"

  if (is.character(layout)) {
    layout <- match.arg(layout, c("fr", "kk", "circle"))

    return(
      switch(
        layout,
        fr = igraph::layout_with_fr(g, weights = weights),
        kk = igraph::layout_with_kk(g, weights = weights),
        circle = igraph::layout_in_circle(g)
      )
    )
  }

  layout <- as.matrix(layout)

  if (nrow(layout) != igraph::vcount(g) || ncol(layout) < 2) {
    stop(
      "`layout` must be either one of c('fr', 'kk', 'circle') or a ",
      "numeric matrix with one row per vertex and at least two columns.",
      call. = FALSE
    )
  }

  layout[, 1:2, drop = FALSE]
}

#' @keywords internal
#' @noRd
.layout_multilayer <- function(
    g,
    layer_attr = "layer",
    radius_layer = 4,
    layout = "fr"
) {
  stopifnot(igraph::is_igraph(g))

  layer_raw <- igraph::vertex_attr(g, layer_attr)
  if (is.null(layer_raw)) {
    stop("Vertex attribute '", layer_attr, "' not found.", call. = FALSE)
  }

  layer <- factor(layer_raw)
  layer_levels <- levels(layer)

  n <- igraph::vcount(g)
  L <- length(layer_levels)

  coords <- matrix(NA_real_, nrow = n, ncol = 2)
  rownames(coords) <- igraph::V(g)$name

  for (i in seq_along(layer_levels)) {

    Lname <- layer_levels[i]
    idx <- which(layer == Lname)

    if (!length(idx)) next

    subg <- igraph::induced_subgraph(g, vids = idx)

    w <- igraph::E(subg)$weight
    w <- if (is.null(w)) rep(1, igraph::ecount(subg)) else abs(w)
    w[is.na(w) | w == 0] <- 1

    layout_i <- "fr"

    if (is.list(layout)) {

      if (!is.null(layout[[Lname]])) {
        layout_i <- layout[[Lname]]
      }

    } else {
      layout_i <- layout
    }

    if (igraph::vcount(subg) == 1L) {
      xy <- matrix(c(0, 0), nrow = 1)
    } else if (igraph::ecount(subg) == 0L && is.character(layout_i) && layout_i == "fr") {
      k <- igraph::vcount(subg)
      ang <- seq(0, 2 * pi, length.out = k + 1L)[-1L]
      xy <- cbind(cos(ang), sin(ang))
    } else {
      xy <- .resolve_layout_single(subg, layout = layout_i, weights = w)
    }

    if (nrow(xy) > 1L) {
      xy <- unclass(scale(xy))
      xy[is.na(xy)] <- 0
    }

    coords[idx, ] <- xy
  }

  angles_layer <- seq(0, 2 * pi, length.out = L + 1L)[-1L]

  for (i in seq_along(layer_levels)) {

    idx <- which(layer == layer_levels[i])
    if (!length(idx)) next

    target <- c(
      radius_layer * cos(angles_layer[i]),
      radius_layer * sin(angles_layer[i])
    )

    center_now <- colMeans(coords[idx, , drop = FALSE], na.rm = TRUE)

    coords[idx, ] <- sweep(
      coords[idx, , drop = FALSE],
      2,
      target - center_now,
      "+"
    )
  }

  coords
}
