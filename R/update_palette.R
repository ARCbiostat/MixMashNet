#' Update community and layer color palettes in MixMashNet objects
#'
#' @description
#' Updates the color palettes associated with communities and/or layers in
#' \code{mixMN_fit} and \code{multimixMN_fit} objects. The function replaces only
#' the colors corresponding to the provided names, leaving all other colors
#' unchanged.
#'
#' If colors are provided for names that do not exist in the object (e.g.,
#' unknown community labels or layer names), a warning is issued and those entries
#' are ignored. If some communities or layers are not specified, their original
#' colors are preserved.
#'
#' @param fit An object of class \code{mixMN_fit} or \code{multimixMN_fit}.
#' @param community_colors Optional named character vector specifying new colors
#'   for communities. Names must correspond to existing community labels
#'   (as stored in \code{communities$palette}). Missing names are ignored.
#' @param layer_colors Optional named character vector specifying new colors for
#'   layers. Names must correspond to existing layer names
#'   (as stored in \code{layers$palette}). Only applicable to
#'   \code{multimixMN_fit} objects.
#'
#' @details
#' For \code{mixMN_fit} objects, community colors are updated in
#' \code{fit$communities$palette}.
#'
#' For \code{multimixMN_fit} objects, community colors are updated separately
#' within each layer (i.e., in \code{fit$layer_fits[[L]]$communities$palette}),
#' while layer colors are updated in \code{fit$layers$palette}.
#'
#' The function performs in-place modification of the palettes and returns the
#' updated object.
#'
#' @return
#' The input object \code{fit}, with updated community and/or layer palettes.
#'
#' @export
update_palette <- function(
    fit,
    community_colors = NULL,  # named vector: names = community labels
    layer_colors = NULL       # named vector: names = layer names
) {
  stopifnot(is.list(fit))

  is_multilayer <- inherits(fit, "multimixMN_fit") ||
    (!is.null(fit$layers) && !is.null(fit$layer_fits))
  is_single <- inherits(fit, "mixMN_fit") ||
    (!is.null(fit$communities) && is.null(fit$layers))

  if (!is_multilayer && !is_single) {
    stop("`fit` must be a mixMN_fit or multimixMN_fit object.", call. = FALSE)
  }

  .warn_extra_names <- function(given, allowed, what) {
    if (is.null(names(given))) return(invisible(NULL))
    extra <- setdiff(names(given), allowed)
    if (length(extra)) {
      warning(
        sprintf(
          "%s: ignoring %d unknown name(s): %s",
          what, length(extra), paste(extra, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    invisible(NULL)
  }

  .apply_named_update <- function(old, new_named, allowed_names, what) {
    if (is.null(new_named)) return(old)

    if (is.null(names(new_named)) || anyNA(names(new_named)) || any(names(new_named) == "")) {
      stop(sprintf(
        "%s must be a *named* vector (names = ids, values = colors).",
        what
      ), call. = FALSE)
    }

    .warn_extra_names(new_named, allowed_names, what)

    upd <- intersect(names(new_named), allowed_names)
    if (!length(upd)) return(old)

    old <- as.character(old)
    names(old) <- allowed_names
    old[upd] <- as.character(new_named[upd])
    old
  }

  # ---- communities palette ----
  if (!is.null(community_colors)) {

    # single-layer
    if (is_single && !is.null(fit$communities$palette)) {
      allowed <- names(fit$communities$palette)
      fit$communities$palette <- .apply_named_update(
        fit$communities$palette,
        community_colors,
        allowed,
        "community_colors"
      )
    }

    # multilayer: per-layer communities
    if (is_multilayer && !is.null(fit$layer_fits)) {
      for (L in names(fit$layer_fits)) {
        pal <- fit$layer_fits[[L]]$communities$palette
        if (is.null(pal)) next
        allowed <- names(pal)
        fit$layer_fits[[L]]$communities$palette <- .apply_named_update(
          pal,
          community_colors,
          allowed,
          sprintf("community_colors (layer '%s')", L)
        )
      }
    }
  }

  # ---- layers palette (multilayer only) ----
  if (is_multilayer && !is.null(layer_colors)) {
    if (is.null(fit$layers$palette)) {
      stop("This multimixMN_fit object has no `layers$palette`.", call. = FALSE)
    }
    allowed <- names(fit$layers$palette)
    fit$layers$palette <- .apply_named_update(
      fit$layers$palette,
      layer_colors,
      allowed,
      "layer_colors"
    )
  }

  fit
}
