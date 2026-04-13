# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' Internal helper for summary methods
#' @keywords internal
#' @noRd
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

#' Internal helper for summary methods
#' @keywords internal
#' @noRd
.apply_named_update <- function(old, new_named, allowed_names, what) {
  if (is.null(new_named)) return(old)

  if (is.null(names(new_named)) ||
      anyNA(names(new_named)) ||
      any(names(new_named) == "")) {
    stop(
      sprintf(
        "%s must be a *named* vector (names = ids, values = colors).",
        what
      ),
      call. = FALSE
    )
  }

  .warn_extra_names(new_named, allowed_names, what)

  upd <- intersect(names(new_named), allowed_names)
  if (!length(upd)) return(old)

  old <- as.character(old)
  names(old) <- allowed_names
  old[upd] <- as.character(new_named[upd])
  old
}

#' Update community and layer color palettes in MixMashNet objects
#'
#' @description
#' Updates the color palettes associated with communities and/or layers in
#' fitted \code{mixMN_fit} and \code{multimixMN_fit} objects.
#'
#' For \code{mixMN_fit} objects, \code{community_colors} must be a named
#' character vector specifying colors for community labels in
#' \code{object$communities$palette}.
#'
#' For \code{multimixMN_fit} objects, \code{community_colors} must be a named
#' list whose elements correspond to layer names. Each element must be a named
#' character vector specifying colors for the community labels of that layer.
#' The list may be partial, so only the specified layers are updated.
#'
#' For \code{multimixMN_fit} objects, \code{layer_colors} updates the palette
#' stored in \code{object$layers$palette}.
#'
#' The function replaces only the colors corresponding to the provided names,
#' leaving all other colors unchanged. Unknown layer names, community labels,
#' or layer labels are ignored with a warning.
#'
#' @aliases update_palette update_palette.mixMN_fit update_palette.multimixMN_fit
#' @param object An object of class \code{mixMN_fit} or
#'   \code{multimixMN_fit}.
#' @param community_colors For \code{mixMN_fit} objects, an optional named
#'   character vector specifying new colors for communities.
#'
#'   For \code{multimixMN_fit} objects, an optional named list whose names are
#'   layer names and whose elements are named character vectors specifying new
#'   colors for communities within each layer.
#' @param layer_colors Optional named character vector specifying new colors for
#'   layers. Only applicable to \code{multimixMN_fit} objects.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' For single layer fits, only \code{community_colors} is used.
#'
#' For multilayer fits:
#' \itemize{
#'   \item \code{community_colors} updates community palettes within the
#'   specified layers;
#'   \item \code{layer_colors} updates the palette of the layers themselves.
#' }
#'
#' For multilayer fits, \code{community_colors} can be partial: layers not
#' included in the list are left unchanged.
#'
#' @return
#' The input object, with updated community and/or layer palettes.
#'
#' @examples
#' data(bacteremia)
#'
#' vars <- c("WBC", "NEU", "HGB", "PLT", "CRP")
#' df <- bacteremia[, vars]
#'
#' fit <- mixMN(
#'   data = df,
#'   lambdaSel = "EBIC",
#'   reps = 0,
#'   seed_model = 42,
#'   compute_loadings = FALSE,
#'   progress = FALSE
#' )
#'
#' fit$communities$palette
#'
#' fit2 <- update_palette(
#'   fit,
#'   community_colors = c("1" = "red", "2" = "blue")
#' )
#'
#' fit2$communities$palette
#'
#' set.seed(1)
#' plot(fit2)
#'
#' @export
update_palette <- function(object, ...) {
  UseMethod("update_palette")
}

#' @rdname update_palette
#' @export
update_palette.mixMN_fit <- function(
    object,
    community_colors = NULL,
    layer_colors = NULL,
    ...
) {
  stopifnot(is.list(object))

  if (!is.null(layer_colors)) {
    warning(
      "`layer_colors` is ignored for mixMN_fit objects.",
      call. = FALSE
    )
  }

  if (!is.null(community_colors)) {
    if (is.list(community_colors)) {
      stop(
        "`community_colors` must be a named character vector for mixMN_fit objects.",
        call. = FALSE
      )
    }

    if (is.null(object$communities$palette)) {
      stop("This mixMN_fit object has no `communities$palette`.", call. = FALSE)
    }

    allowed <- names(object$communities$palette)
    object$communities$palette <- .apply_named_update(
      object$communities$palette,
      community_colors,
      allowed,
      "community_colors"
    )
  }

  object
}

#' @rdname update_palette
#' @export
update_palette.multimixMN_fit <- function(
    object,
    community_colors = NULL,
    layer_colors = NULL,
    ...
) {
  stopifnot(is.list(object))

  # ---- communities palette ----
  if (!is.null(community_colors)) {
    if (!is.list(community_colors) || is.null(names(community_colors)) ||
        anyNA(names(community_colors)) || any(names(community_colors) == "")) {
      stop(
        paste(
          "`community_colors` must be a named list for multimixMN_fit objects.",
          "Names must correspond to layer names, and each element must be a",
          "named character vector of community colors."
        ),
        call. = FALSE
      )
    }

    if (is.null(object$layer_fits)) {
      stop("This multimixMN_fit object has no `layer_fits`.", call. = FALSE)
    }

    available_layers <- names(object$layer_fits)
    .warn_extra_names(community_colors, available_layers, "community_colors")

    layers_to_update <- intersect(names(community_colors), available_layers)

    for (L in layers_to_update) {
      pal <- object$layer_fits[[L]]$communities$palette
      if (is.null(pal)) next

      layer_update <- community_colors[[L]]

      if (is.null(layer_update)) next
      if (!is.character(layer_update)) {
        stop(
          sprintf(
            "community_colors[['%s']] must be a named character vector.",
            L
          ),
          call. = FALSE
        )
      }

      allowed <- names(pal)
      object$layer_fits[[L]]$communities$palette <- .apply_named_update(
        pal,
        layer_update,
        allowed,
        sprintf("community_colors[['%s']]", L)
      )
    }
  }

  # ---- layers palette ----
  if (!is.null(layer_colors)) {
    if (is.null(object$layers$palette)) {
      stop("This multimixMN_fit object has no `layers$palette`.", call. = FALSE)
    }
    allowed <- names(object$layers$palette)
    object$layers$palette <- .apply_named_update(
      object$layers$palette,
      layer_colors,
      allowed,
      "layer_colors"
    )
  }

  object
}
