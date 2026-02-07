#' Infer MGM type/level specification and prepare numeric data for \pkg{mgm}
#'
#' @description
#' Internal helper used by \code{mixMN} and \code{multimixMN} to infer \code{type} and \code{level}
#' vectors required by \code{mgm::mgm}, and to build a numeric matrix suitable for
#' model fitting. When \code{recode_binary = TRUE}, binary variables stored as
#' two-level factors (or ordered factors) and logical variables are internally
#' recoded to \{0,1\} to satisfy the requirements of \code{mgm} when
#' \code{binarySign = TRUE}. The original input \code{data} is not modified.
#'
#' @details
#' The inference rules are:
#' \itemize{
#'   \item \strong{Factor / ordered}: \code{type = "c"}, \code{level = nlevels(x)}.
#'         If \code{nlevels(x) == 2} and \code{recode_binary = TRUE}, the variable is
#'         internally recoded so that the first level maps to 0 and the second level to 1.
#'   \item \strong{Logical}: \code{type = "c"}, \code{level = 2}. If
#'         \code{recode_binary = TRUE}, \code{FALSE} is mapped to 0 and \code{TRUE} to 1.
#'   \item \strong{Integer}: if the observed non-missing values are all in \{0,1\},
#'         then \code{type = "c"}, \code{level = 2}; otherwise \code{type = "p"},
#'         \code{level = 1}.
#'   \item \strong{Numeric (double)}: \code{type = "g"}, \code{level = 1}.
#' }
#'
#' The returned \code{data_info} data frame provides a compact audit trail of the
#' inferred specification and which variables were internally recoded. The
#' \code{binary_recode_map} records the mapping from original binary labels to
#' \{0,1\} (useful for interpreting signs when \code{binarySign = TRUE}).
#'
#' @param data A \code{data.frame} (n x p) with variables in columns. Must have
#'   column names. Character and Date/POSIXt columns are not supported.
#' @param recode_binary Logical; if \code{TRUE} (default), internally recode
#'   two-level factors/ordered factors and logical variables to \{0,1\} in
#'   \code{data_mgm}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{type}}{Character vector of MGM variable types (\code{"g"}, \code{"p"}, \code{"c"}).}
#'   \item{\code{level}}{Integer vector of MGM levels (1 for \code{"g"} and \code{"p"}, number of categories for \code{"c"}).}
#'   \item{\code{data_mgm}}{Numeric matrix (n x p) to be passed to \code{mgm::mgm}.}
#'   \item{\code{data_info}}{Data frame summarizing the inferred specification and recoding flags.}
#'   \item{\code{binary_recode_map}}{Named list mapping original binary labels to \{0,1\} for recoded variables.}
#' }
#'
#' @keywords internal
#' @noRd
infer_mgm_spec <- function(data, recode_binary = TRUE) {

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (is.null(colnames(data))) stop("`data` must have column names.")

  nodes <- colnames(data)
  type  <- rep(NA_character_, length(nodes))
  level <- rep(NA_integer_,   length(nodes))

  binary_recode_map <- list()
  binary_recoded <- rep(FALSE, length(nodes))

  # this is ONLY for mgm fitting
  data_for_mgm <- data

  for (j in seq_along(nodes)) {
    x  <- data[[j]]
    nm <- nodes[j]

    if (is.character(x)) {
      stop(sprintf("Column '%s' is 'character'. Convert to factor or numeric.", nm))
    }
    if (inherits(x, "Date") || inherits(x, "POSIXt")) {
      stop(sprintf("Column '%s' is Date/POSIXt. Convert to numeric or factor.", nm))
    }

    if (is.factor(x) || is.ordered(x)) {
      type[j]  <- "c"
      level[j] <- nlevels(x)

      if (recode_binary && nlevels(x) == 2) {
        lv <- levels(x)

        if (setequal(lv, c("0","1"))) {
          data_for_mgm[[j]] <- as.integer(as.character(x))  # "0"->0, "1"->1
        } else {
          data_for_mgm[[j]] <- as.integer(x == lv[2L])  # lv[1]=0, lv[2]=1
          binary_recode_map[[nm]] <- stats::setNames(c(0, 1), lv)
          binary_recoded[j] <- TRUE
        }
      }

    } else if (is.logical(x)) {
      type[j]  <- "c"
      level[j] <- 2L

      # logical -> {0,1} for mgm matrix
      if (recode_binary) {
        data_for_mgm[[j]] <- as.integer(x)  # FALSE=0, TRUE=1
        binary_recode_map[[nm]] <- stats::setNames(c(0, 1), c("FALSE", "TRUE"))
        binary_recoded[j] <- TRUE
      }

    } else if (is.integer(x)) {
      ux <- unique(x[!is.na(x)])

      if (length(ux) > 0 && all(ux %in% c(0L, 1L))) {
        type[j]  <- "c"
        level[j] <- 2L
        # already {0,1}; do NOT mark as recoded and do NOT add to map
      } else {
        type[j]  <- "p"
        level[j] <- 1L
      }

    } else if (is.numeric(x)) {
      type[j]  <- "g"
      level[j] <- 1L

    } else {
      stop(sprintf("Column '%s' has unsupported class: %s", nm, paste(class(x), collapse = ", ")))
    }
  }

  data_info <- data.frame(
    node = nodes,
    class = vapply(data, function(x) paste(class(x), collapse = "/"), character(1)),
    mgm_type = type,
    mgm_level = level,
    n_unique = vapply(data, function(x) length(unique(x[!is.na(x)])), integer(1)),
    binary_recoded_01 = binary_recoded,
    stringsAsFactors = FALSE
  )

  list(
    type = type,
    level = level,
    data_mgm = data.matrix(data_for_mgm),   # <- ONLY for mgm()
    data_info = data_info,
    binary_recode_map = binary_recode_map
  )
}
