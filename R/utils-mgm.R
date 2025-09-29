#' Access internal functions from mgm safely
#'
#' @keywords internal
.mgm_fun <- function(fname) {
  if (!requireNamespace("mgm", quietly = TRUE)) {
    stop("Package 'mgm' is required but not installed.")
  }
  get(fname, envir = asNamespace("mgm"), inherits = FALSE)
}
