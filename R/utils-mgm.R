#' Access internal functions from mgm safely
#'
#' @keywords internal
#' @noRd
.mgm_fun <- function(fname) {
  get(fname, envir = asNamespace("mgm"), inherits = FALSE)
}
