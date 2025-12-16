# ============================================================
# Multilayer layouts
# ============================================================

#' Fruchterman–Reingold layout for multilayer graphs
#'
#' @description
#' Internal helper that computes 2D coordinates for a multilayer
#' \pkg{igraph} object by running a weighted Fruchterman–Reingold
#' layout separately within each layer and then arranging layer
#' centroids on a larger circle.
#'
#' @param g An \code{igraph} object with a vertex attribute
#'   containing layer membership.
#' @param layer_attr Character string giving the name of the
#'   vertex attribute that stores the layer label.
#'
#' @return A numeric matrix with two columns (x, y) and one row
#'   per vertex, ordered as \code{V(g)}.
#'
#' @keywords internal
#' @noRd
.layout_multilayer_fr <- function(g, layer_attr = "layer", radius_layer = 4) {
  stopifnot(igraph::is_igraph(g))

  layer_raw <- igraph::vertex_attr(g, layer_attr)
  if (is.null(layer_raw)) stop("Vertex attribute '", layer_attr, "' not found.")
  layer <- factor(layer_raw)

  L <- nlevels(layer)
  n <- igraph::vcount(g)
  coords <- matrix(NA_real_, nrow = n, ncol = 2)

  # 1) layout FR within each layer (intra-layer edges only)
  for (i in seq_len(L)) {
    idx <- which(layer == levels(layer)[i])
    if (length(idx) == 0L) next

    subg <- igraph::induced_subgraph(g, vids = idx)

    if (igraph::ecount(subg) > 0) {
      w <- igraph::E(subg)$weight
      w <- if (is.null(w)) rep(1, igraph::ecount(subg)) else abs(w)
      w[is.na(w) | w == 0] <- 1
      xy <- igraph::layout_with_fr(subg, weights = w, niter = 300)
    } else {
      k <- length(idx)
      ang <- seq(0, 2*pi, length.out = k + 1L)[-1L]
      xy <- cbind(cos(ang), sin(ang))
    }

    if (nrow(xy) == 1L) {
      xy <- matrix(c(0, 0), nrow = 1)
    } else {
      xy <- unclass(scale(xy))
    }

    coords[idx, ] <- xy
  }

  # 2) place layer centroids on a larger circle
  angles_layer <- seq(0, 2*pi, length.out = L + 1L)[-1L]

  for (i in seq_len(L)) {
    idx <- which(layer == levels(layer)[i])
    if (length(idx) == 0L) next
    if (anyNA(coords[idx, ])) next

    center_now <- colMeans(coords[idx, , drop = FALSE])
    target <- c(radius_layer * cos(angles_layer[i]),
                radius_layer * sin(angles_layer[i]))
    coords[idx, ] <- sweep(coords[idx, , drop = FALSE], 2, target - center_now, "+")
  }

  coords
}
