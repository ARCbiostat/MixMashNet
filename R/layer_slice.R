#' Extract a single layer from a multilayer MixMashNet object
#'
#' @description
#' Extracts one layer from a fitted multilayer \code{multimixMN_fit} object
#' returned by \code{multimixMN()}.
#'
#' The selected layer is returned as the corresponding single layer
#' \code{mixMN_fit} object stored in \code{layer_fits}.
#'
#' @aliases layer_slice layer_slice.multimixMN_fit
#' @param object An object of class \code{"multimixMN_fit"} returned by
#'   \code{multimixMN()}.
#' @param layer Character string giving the layer to extract.
#' @param ... Further arguments passed to methods.
#'
#' @return
#' An object of class \code{"mixMN_fit"} corresponding to the
#' selected layer.
#'
#' @export
layer_slice <- function(object, ...) {
  UseMethod("layer_slice")
}

#' @rdname layer_slice
#' @export
layer_slice.multimixMN_fit <- function(object, layer, ...) {
  if (missing(layer) || length(layer) != 1L || is.na(layer)) {
    stop("`layer` must be a single non-missing character string.", call. = FALSE)
  }

  layer <- as.character(layer)

  available_layers <- names(object$layer_fits)

  if (is.null(available_layers) || !length(available_layers)) {
    stop("No layer-specific fits are available in `object$layer_fits`.", call. = FALSE)
  }

  if (!layer %in% available_layers) {
    stop(
      "`layer` not found. Available layers are: ",
      paste(available_layers, collapse = ", "),
      call. = FALSE
    )
  }

  object$layer_fits[[layer]]
}
