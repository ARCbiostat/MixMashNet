#' Node stability from bootstrap community assignments
#'
#' @description
#' Computes per-node stability given the empirical community structure and the
#' homogenized bootstrap memberships contained in a \code{mixMN_fit} object.
#'
#' @param fit An object returned by \code{mixMN()} (class \code{mixMN_fit}).
#' @param IS.plot Logical; if \code{TRUE}, prints the plot returned by
#'   \code{membershipStab_plot()}.
#'
#' @return An object of class \code{c("membershipStab")}.
#' @importFrom EGAnet community.homogenize
#' @export
membershipStab <- function(fit, IS.plot = FALSE) {
  # --- Extract inputs
  structure <- fit$original_membership
  boot.list <- fit$boot_memberships
  palette   <- fit$community_palette

  # --- Checks
  if (is.null(structure) || length(structure) == 0L)
    stop("original_membership is missing or empty.")
  if (any(is.na(structure)))
    stop("original_membership contains NAs. Cannot compute node stability.")

  p <- length(structure)
  if (!is.null(names(structure))) {
    var_names <- names(structure)
  } else {
    stop("original_membership must have variable names.")
  }

  # --- Bootstrap membership matrix (reps x p)
  if (!length(boot.list)) stop("boot_memberships is missing or empty.")
  boot.wc <- do.call(rbind, boot.list)
  if (ncol(boot.wc) != p)
    stop("Inconsistent bootstrap membership dimensions: expected ", p, ", got ", ncol(boot.wc), ".")
  colnames(boot.wc) <- var_names

  # --- Homogenize bootstrap labels to target structure
  homogenized <- EGAnet::community.homogenize(
    target.membership  = structure,
    convert.membership = boot.wc
  )

  # --- Max number of communities across target + bootstrap (robust to NA)
  max_communities <- max(
    max(structure, na.rm = TRUE),
    suppressWarnings(max(homogenized, na.rm = TRUE))
  )

  # --- Proportion matrix: p x K (D1..DK), per-node denominator = #valid (non-NA)
  prop_matrix <- matrix(NA_real_, nrow = p, ncol = max_communities)
  for (j in seq_len(p)) {
    colj <- homogenized[, j]
    valid <- !is.na(colj)
    if (!any(valid)) next
    counts <- tabulate(colj[valid], nbins = max_communities)
    prop_matrix[j, ] <- counts / sum(valid)
  }
  colnames(prop_matrix) <- paste0("D", seq_len(max_communities))
  rownames(prop_matrix) <- var_names

  # --- Empirical (target) community stability per node
  prop_empirical <- vapply(seq_along(structure), function(i) {
    k <- structure[i]
    if (!is.na(k) && k <= max_communities) prop_matrix[i, k] else NA_real_
  }, numeric(1))
  names(prop_empirical) <- var_names

  # --- Output object
  result <- list(
    membership = list(
      empirical  = structure,
      bootstrap  = homogenized,
      structure  = structure
    ),
    membership.stability = list(
      empirical.dimensions = prop_empirical,
      all.dimensions       = prop_matrix
    )
  )
  class(result) <- c("membershipStab")

  if (!is.null(palette)) {
    attr(result, "palette") <- palette
    result$community_palette <- palette
  }

  if (isTRUE(IS.plot)) {
    print(membershipStab_plot(result))
  }
  return(result)
}
