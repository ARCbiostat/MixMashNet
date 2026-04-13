#' Node stability from bootstrap community assignments
#'
#' @description
#' Computes per-node stability given the empirical community structure and the
#' homogenized bootstrap memberships contained in a \code{mixMN_fit} object.
#' Stability is expressed as the proportion of bootstrap replications that
#' assign each node to its empirical (original) community.
#'
#' @param fit An object returned by \code{mixMN()} (class \code{mixMN_fit}),
#'   containing \code{$communities$original_membership} and
#'   \code{$communities$boot_memberships}. Bootstrap memberships must be
#'   available, i.e. \code{reps > 0} and \code{"community" \%in\% boot_what}.
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
membershipStab <- function(fit) {
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
  class(result) <- "membershipStab"

  if (!is.null(palette)) {
    attr(result, "palette") <- palette
    result$community_palette <- palette
  }

  return(result)
}

#' @export
print.membershipStab <- function(x, ...) {
  if (!inherits(x, "membershipStab")) {
    stop("`x` must be an object of class 'membershipStab'.")
  }

  empirical <- x$membership$empirical
  bootstrap <- x$membership$bootstrap
  stab <- x$membership.stability$empirical.dimensions

  n_nodes <- if (!is.null(empirical)) length(empirical) else 0L
  n_boot  <- if (!is.null(bootstrap) && is.matrix(bootstrap)) nrow(bootstrap) else 0L
  n_comm  <- if (!is.null(empirical)) length(unique(stats::na.omit(as.integer(empirical)))) else 0L

  cat("MixMashNet membership stability object\n")
  cat(strrep("=", 38), "\n", sep = "")
  cat("Nodes: ", n_nodes, "\n", sep = "")
  cat("Bootstrap replications: ", n_boot, "\n", sep = "")
  cat("Communities: ", n_comm, "\n", sep = "")

  if (!is.null(stab) && length(stab)) {
    cat(
      "Node stability range: ",
      round(min(stab, na.rm = TRUE), 3), " to ",
      round(max(stab, na.rm = TRUE), 3), "\n",
      sep = ""
    )
  }

  invisible(x)
}

#' @export
summary.membershipStab <- function(object, ...) {
  if (!inherits(object, "membershipStab")) {
    stop("`object` must be an object of class 'membershipStab'.")
  }

  empirical <- object$membership$empirical
  bootstrap <- object$membership$bootstrap
  stab <- object$membership.stability$empirical.dimensions

  n_nodes <- if (!is.null(empirical)) length(empirical) else 0L
  n_boot  <- if (!is.null(bootstrap) && is.matrix(bootstrap)) nrow(bootstrap) else 0L
  n_comm  <- if (!is.null(empirical)) {
    length(unique(stats::na.omit(as.integer(empirical))))
  } else {
    0L
  }

  stab_summary <- if (!is.null(stab) && length(stab)) {
    stats::quantile(
      stab,
      probs = c(0, 0.25, 0.5, 0.75, 1),
      na.rm = TRUE,
      names = FALSE
    )
  } else {
    rep(NA_real_, 5)
  }

  names(stab_summary) <- c("min", "q1", "median", "q3", "max")

  out <- list(
    n_nodes = n_nodes,
    n_bootstrap = n_boot,
    n_communities = n_comm,
    mean_stability = if (!is.null(stab) && length(stab)) mean(stab, na.rm = TRUE) else NA_real_,
    stability_summary = stab_summary,
    n_below_0_50 = if (!is.null(stab) && length(stab)) sum(stab < 0.50, na.rm = TRUE) else NA_integer_,
    n_below_0_70 = if (!is.null(stab) && length(stab)) sum(stab < 0.70, na.rm = TRUE) else NA_integer_
  )

  class(out) <- "summary.membershipStab"
  out
}

#' @export
print.summary.membershipStab <- function(x, digits = 3, ...) {
  if (!inherits(x, "summary.membershipStab")) {
    stop("`x` must be an object of class 'summary.membershipStab'.")
  }

  cat("Summary of MixMashNet membership stability\n")
  cat(strrep("=", 42), "\n", sep = "")
  cat("Nodes: ", x$n_nodes, "\n", sep = "")
  cat("Bootstrap replications: ", x$n_bootstrap, "\n", sep = "")
  cat("Communities: ", x$n_communities, "\n", sep = "")
  cat("Mean node stability: ", round(x$mean_stability, digits), "\n", sep = "")

  cat("\nNode stability distribution:\n")
  ss <- round(x$stability_summary, digits)
  cat("  Min:    ", ss["min"], "\n", sep = "")
  cat("  1st Qu.:", ss["q1"], "\n", sep = "")
  cat("  Median: ", ss["median"], "\n", sep = "")
  cat("  3rd Qu.:", ss["q3"], "\n", sep = "")
  cat("  Max:    ", ss["max"], "\n", sep = "")

  cat("\nNodes below stability thresholds:\n")
  cat("  < 0.50: ", x$n_below_0_50, "\n", sep = "")
  cat("  < 0.70: ", x$n_below_0_70, "\n", sep = "")

  invisible(x)
}

#' @export
plot.membershipStab <- function(
    x,
    title = "Node Stability by Community",
    cutoff = 0.7,
    ...
) {
  if (!inherits(x, "membershipStab")) {
    stop("`x` must be an object of class 'membershipStab'.")
  }

  membershipStab_plot(
    stab_obj = x,
    title = title,
    cutoff = cutoff
  )
}
