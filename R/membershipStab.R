#' Node stability from bootstrap community assignments
#'
#' @description
#' Computes per-node stability given the empirical community structure and the
#' homogenized bootstrap memberships contained in a \code{mixMN_fit} object.
#' This function is used internally by \code{mixMN()} and \code{multimixMN()}.
#' Stability is expressed as the proportion of bootstrap replications that
#' assign each node to its empirical (original) community.
#'
#' @param fit An object returned by \code{mixMN()} (class \code{mixMN_fit}),
#'   containing \code{$communities$original_membership} and
#'   \code{$communities$boot_memberships}. Bootstrap memberships must be
#'   available, i.e. \code{reps > 0} and \code{"community" \%in\% boot_what}.
#' @param IS.plot Logical; if \code{TRUE}, prints a stability plot via the
#'   internal helper \code{membershipStab_plot()}.
#'
#' @return An object of class \code{c("membershipStab")}, with components:
#' \describe{
#'   \item{\code{membership}}{List with:
#'     \describe{
#'       \item{\code{empirical}}{Named integer vector of empirical community labels}
#'       \item{\code{bootstrap}}{Matrix of homogenized bootstrap labels
#'         (\code{reps × p})}
#'     }
#'   }
#'   \item{\code{membership.stability}}{List with:
#'     \describe{
#'       \item{\code{empirical.dimensions}}{Named numeric vector of node-level stability
#'         (proportion assigned to empirical community)}
#'       \item{\code{all.dimensions}}{Matrix (\code{p × K}) with proportions of
#'         assignment to each community}
#'     }
#'   }
#'   \item{\code{community_palette}}{Named vector of colors for communities,
#'     if available}
#' }
#'
#' @details
#' Bootstrap community labels are first aligned to the empirical solution using
#' \code{EGAnet::community.homogenize()}. Stability is then computed node-wise as
#' the proportion of bootstrap runs in which the node's community matches its
#' empirical assignment.
#'
#' @references
#'
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the Stability of Psychological Dimensions via Bootstrap Exploratory Graph Analysis:
#' A Monte Carlo Simulation and Tutorial. \emph{Psych}, 3(3), 479–500.
#' \doi{10.3390/psych3030032}
#'
#' @importFrom EGAnet community.homogenize
#' @export
membershipStab <- function(fit, IS.plot = FALSE) {
  # --- Extract inputs
  structure <- fit$communities$original_membership
  boot.list <- fit$communities$boot_memberships
  palette   <- fit$communities$palette

  # --- Checks
  if (is.null(structure) || length(structure) == 0L)
    stop("original_membership is missing or empty.")
  if (any(is.na(structure)))
    stop("original_membership contains NAs. Cannot compute node stability.")
  if (!inherits(fit, "mixMN_fit")) {
    stop("`fit` must be an object of class 'mixMN_fit'.")
  }
  if (is.null(boot.list) || !length(boot.list)) {
    stop(
      "No bootstrap memberships found in `communities$boot_memberships`.\n",
      "Make sure `boot_what` includes \"community\" and `reps > 0`."
    )
  }

  p <- length(structure)
  if (!is.null(names(structure))) {
    var_names <- names(structure)
  } else {
    stop("original_membership must have variable names.")
  }

  # --- Bootstrap membership matrix (reps x p)
  if (!length(boot.list)) stop("boot_memberships is missing or empty.")
  boot.wc <- do.call(
    rbind,
    lapply(boot.list, function(z) {
      if (is.null(z)) {
        rep(NA_integer_, p)
      } else {
        if (is.null(names(z))) {
          stop("All bootstrap membership vectors must be named.")
        }
        z[var_names]
      }
    })
  )
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
      bootstrap  = homogenized
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
