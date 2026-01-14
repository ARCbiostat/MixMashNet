#' Node-by-community stability table
#'
#' @description
#' Returns a data frame with node-by-community stability proportions, optionally
#' rounded; zero entries are replaced with NA for readability.
#'
#' @param stab_obj An object from \code{membershipStab()}.
#' @param digits Integer; decimal places for rounding.
#'
#' @return A data.frame with rows = nodes and columns = communities (D1..DK).
#' @export
membershipStab_table <- function(stab_obj, digits = 3) {
  prop_matrix <- stab_obj$membership.stability$all.dimensions
  if (is.null(prop_matrix)) stop("Missing membership.stability$all.dimensions.")

  prop_df <- as.data.frame(round(prop_matrix, digits))
  prop_df[prop_df == 0] <- NA

  # Order columns by numeric community index
  num_idx <- suppressWarnings(as.integer(sub("\\D+", "", colnames(prop_df))))
  ord <- order(num_idx, na.last = TRUE)
  prop_df <- prop_df[, ord, drop = FALSE]

  return(prop_df)
}
