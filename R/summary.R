# -------------------------------------------------------------------------
# internal helper
# -------------------------------------------------------------------------

#' Internal helper for summary method
#' @keywords internal
#' @noRd
.summary_top_edges <- function(x, top_n = 10L) {
  if (is.null(x) || !nrow(x)) return(x)

  if (!"estimated" %in% colnames(x)) {
    stop("Input must contain column 'estimated'.", call. = FALSE)
  }

  keep <- !is.na(x$estimated) & x$estimated != 0
  x <- x[keep, , drop = FALSE]

  if (!nrow(x)) return(x)

  x <- x[order(abs(x$estimated), decreasing = TRUE), , drop = FALSE]

  if (nrow(x) > top_n) {
    x <- utils::head(x, top_n)
  }

  rownames(x) <- NULL
  x
}

#' Summarize top intralayer edges for single-layer MixMashNet fits
#'
#' @description
#' Returns the top 10 intralayer edges for fitted objects returned by
#' \code{mixMN()}, in the same long-format structure used for edge summaries.
#'
#' @param object An object of class \code{"mixMN_fit"}.
#' @param top_n Number of top edges to retain. Default is \code{10}.
#' @param digits Optional number of digits used to round numeric columns.
#' @param drop_na_boot Logical. If \code{TRUE} (default), bootstrap-related
#'   columns that are entirely \code{NA} are removed.
#' @param ... Further arguments for S3 compatibility.
#'
#' @return
#' An object of class \code{"summary.mixMN_fit"} containing the top
#' intralayer edges.
#'
#' @method summary mixMN_fit
#' @export
summary.mixMN_fit <- function(object,
                              top_n = 10L,
                              digits = NULL,
                              drop_na_boot = TRUE,
                              ...) {
  if (!inherits(object, "mixMN_fit")) {
    stop("`object` must be a 'mixMN_fit' object.", call. = FALSE)
  }

  out <- get_edges(
    object = object,
    what = "intra",
    digits = digits,
    drop_na_boot = drop_na_boot
  )

  out <- .summary_top_edges(out, top_n = top_n)

  attr(out, "quantile_level") <- attr(out, "quantile_level")
  attr(out, "what") <- "intra"
  attr(out, "top_n") <- top_n

  class(out) <- c("summary.mixMN_fit", class(out))
  out
}

#' Summarize top interlayer edges for multilayer MixMashNet fits
#'
#' @description
#' Returns the top 10 interlayer edges for fitted objects returned by
#' \code{multimixMN()}, in the same long-format structure used for edge summaries.
#'
#' @param object An object of class \code{"multimixMN_fit"}.
#' @param top_n Number of top edges to retain. Default is \code{10}.
#' @param digits Optional number of digits used to round numeric columns.
#' @param drop_na_boot Logical. If \code{TRUE} (default), bootstrap-related
#'   columns that are entirely \code{NA} are removed.
#' @param ... Further arguments for S3 compatibility.
#'
#' @return
#' An object of class \code{"summary.multimixMN_fit"} containing the top
#' interlayer edges.
#'
#' @method summary multimixMN_fit
#' @export
summary.multimixMN_fit <- function(object,
                                   top_n = 10L,
                                   digits = NULL,
                                   drop_na_boot = TRUE,
                                   ...) {
  if (!inherits(object, "multimixMN_fit")) {
    stop("`object` must be a 'multimixMN_fit' object.", call. = FALSE)
  }

  out <- get_edges(
    object = object,
    what = "inter",
    digits = digits,
    drop_na_boot = drop_na_boot
  )

  out <- .summary_top_edges(out, top_n = top_n)

  attr(out, "quantile_level") <- attr(out, "quantile_level")
  attr(out, "what") <- "inter"
  attr(out, "top_n") <- top_n

  class(out) <- c("summary.multimixMN_fit", class(out))
  out
}

#' @method print summary.mixMN_fit
#' @export
print.summary.mixMN_fit <- function(x, digits = 3, ...) {
  if (!inherits(x, "summary.mixMN_fit")) {
    stop("`x` must be a 'summary.mixMN_fit' object.", call. = FALSE)
  }

  if (is.null(x) || !nrow(x)) {
    cat("No edge-level summaries available.\n")
    return(invisible(x))
  }

  quantile_level <- attr(x, "quantile_level")
  top_n <- attr(x, "top_n") %||% 10L

  prettify_colnames <- function(df, quantile_level = 0.95) {
    quantile_pct <- paste0(round(100 * quantile_level), "%")

    map <- c(
      "mean.bootstrap" = "mean (bootstrap)",
      "SE.bootstrap" = "SE (bootstrap)",
      "quantile.lower.bootstrap" =
        paste0(quantile_pct, " quantile lower bound (bootstrap)"),
      "quantile.upper.bootstrap" =
        paste0(quantile_pct, " quantile upper bound (bootstrap)")
    )

    cn <- colnames(df)
    cn <- ifelse(cn %in% names(map), map[cn], cn)
    cn <- gsub("\\.", " ", cn)
    colnames(df) <- cn
    df
  }

  out <- as.data.frame(x, stringsAsFactors = FALSE)

  num_cols <- vapply(out, is.numeric, logical(1L))
  if (any(num_cols)) {
    out[, num_cols] <- lapply(out[, num_cols, drop = FALSE], function(z) {
      idx <- is.finite(z)
      z[idx] <- round(z[idx], digits)
      z
    })
  }

  out <- out[, setdiff(colnames(out), c("scope", "pairs")), drop = FALSE]

  if ("layer" %in% colnames(out) &&
      length(unique(out$layer)) == 1L &&
      all(out$layer == "1")) {
    out <- out[, setdiff(colnames(out), "layer"), drop = FALSE]
  }

  out <- .drop_all_na_boot_cols(out)
  out <- prettify_colnames(out, quantile_level = quantile_level %||% 0.95)

  cat("Top ", top_n, " edge-level summaries (intralayer):\n", sep = "")
  print(out, row.names = FALSE)

  invisible(x)
}

#' @method print summary.multimixMN_fit
#' @export
print.summary.multimixMN_fit <- function(x, digits = 3, ...) {
  if (!inherits(x, "summary.multimixMN_fit")) {
    stop("`x` must be a 'summary.multimixMN_fit' object.", call. = FALSE)
  }

  if (is.null(x) || !nrow(x)) {
    cat("No edge-level summaries available.\n")
    return(invisible(x))
  }

  quantile_level <- attr(x, "quantile_level")
  top_n <- attr(x, "top_n") %||% 10L

  prettify_colnames <- function(df, quantile_level = 0.95) {
    quantile_pct <- paste0(round(100 * quantile_level), "%")

    map <- c(
      "mean.bootstrap" = "mean (bootstrap)",
      "SE.bootstrap" = "SE (bootstrap)",
      "quantile.lower.bootstrap" =
        paste0(quantile_pct, " quantile lower bound (bootstrap)"),
      "quantile.upper.bootstrap" =
        paste0(quantile_pct, " quantile upper bound (bootstrap)")
    )

    cn <- colnames(df)
    cn <- ifelse(cn %in% names(map), map[cn], cn)
    cn <- gsub("\\.", " ", cn)
    colnames(df) <- cn
    df
  }

  out <- as.data.frame(x, stringsAsFactors = FALSE)

  num_cols <- vapply(out, is.numeric, logical(1L))
  if (any(num_cols)) {
    out[, num_cols] <- lapply(out[, num_cols, drop = FALSE], function(z) {
      idx <- is.finite(z)
      z[idx] <- round(z[idx], digits)
      z
    })
  }

  out <- out[, setdiff(colnames(out), c("scope", "layer")), drop = FALSE]
  out <- .drop_all_na_boot_cols(out)
  out <- prettify_colnames(out, quantile_level = quantile_level %||% 0.95)

  cat("Top ", top_n, " edge-level summaries (interlayer):\n", sep = "")
  print(out, row.names = FALSE)

  invisible(x)
}
