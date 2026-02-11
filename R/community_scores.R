#' Compute community scores from a fitted MixMashNet model
#'
#' @description
#' Computes subject-level community scores. Community scores are obtained as
#' weighted sums of the variables belonging to each detected community, where
#' weights correspond to the standardized community loadings estimated via
#' \code{EGAnet::net.loads} and stored in the fitted \code{mixMN_fit} object.
#' Scores are computed using the dataset provided via the \code{data} argument.
#' If \code{data = NULL}, the original dataset used to fit the model
#' (\code{fit$model$data}) is used by default.
#' Optionally, percentile bootstrap quantile regions for the community
#' scores can be computed if bootstrap community loadings are available in
#' \code{fit$community_loadings$boot}.
#'
#' @param fit A fitted object of class \code{c("mixmashnet","mixMN_fit", "multimixMN_fit")}
#'   returned by \code{mixMN()} or \code{multimixMN()}.
#' @param data Optional matrix/data.frame with variables in columns. If
#'   \code{NULL}, uses \code{fit$model$data}. Errors if both are \code{NULL}.
#' @param layer Optional. If fit is a multimixMN_fit, specify which layer to score (name or index).
#'   If NULL, scores are computed for all layers and returned as a named list.
#' @param scale Logical; if \code{TRUE} (default), z-standardize variables used
#'   for scoring, using the mean/SD computed from the dataset used for scoring.
#' @param quantile_level Optional numeric from 0 to 1, e.g. 0.95 or 0.99. If provided,
#'   percentile bootstrap quantile regions are computed for community scores
#'   (requires \code{fit$community_loadings$boot}).
#' @param return_quantile_region Logical; if \code{TRUE}, return quantile regions.
#' @param na_action Character. How to handle missing values in the scoring data:
#'   \code{"stop"} (default) stops if any missing value is present in the
#'   required variables; \code{"omit"} computes scores using row-wise omission
#'   within each community (i.e., uses available variables only, re-normalizing
#'   weights within community for that row).
#'
#' @return A list with class \code{c("mixmashnet","community_scores")} containing:
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{settings}}{List with \code{scale}, \code{quantile_level}, and \code{na_action}.}
#'   \item{\code{ids}}{Character vector of subject IDs (rownames of \code{data}).}
#'   \item{\code{communities}}{Character vector of community score names.}
#'   \item{\code{scores}}{Numeric matrix of scores (n × K).}
#'   \item{\code{quantile_region}}{If requested and available, a list with \code{lower} and \code{upper}
#'     matrices (n × K) for percentile bootstrap quantile regions; otherwise \code{NULL}.}
#'   \item{\code{details}}{List containing \code{nodes_used}, \code{loadings_true},
#'     \code{loadings_boot_available}, and scaling parameters (\code{center}, \code{scale}).}
#' }
#' If \code{fit} is a \code{mixMN_fit} (or a \code{multimixMN_fit} with \code{layer} specified),
#' returns a \code{c("mixmashnet","community_scores")} object.
#' If \code{fit} is a \code{multimixMN_fit} and \code{layer = NULL}, returns a named list
#' of \code{community_scores} objects (one per layer).
#'
#' @details
#' The function requires that \code{fit$community_loadings$true} exists and that
#' the input \code{data} contains all required variables in
#' \code{fit$community_loadings$nodes}. It errors otherwise.
#'
#' @references
#'
#' Christensen, A. P., Golino, H., Abad, F. J., & Garrido, L. E. (2025).
#' Revised network loadings. \emph{Behavior Research Methods}, 57(4), 114.
#' \doi{10.3758/s13428-025-02640-3}
#'
#' @export
community_scores <- function(
    fit,
    data = NULL,
    layer = NULL,
    scale = TRUE,
    quantile_level = 0.95,
    return_quantile_region = FALSE,
    na_action = c("stop", "omit")
) {
  na_action <- match.arg(na_action)

  if (!missing(quantile_level)) {
    return_quantile_region <- TRUE
  }

  # ---- MULTI: compute per layer ----
  if (is.null(fit) || !(inherits(fit, "mixMN_fit") || inherits(fit, "multimixMN_fit"))) {
    stop("`fit` must be a `mixMN_fit` or `multimixMN_fit` object.")
  }

  # ---- MULTI: compute per layer ----
  if (inherits(fit, "multimixMN_fit")) {

    # choose data: from multi by default
    if (is.null(data)) {
      data <- fit$model$data
      if (is.null(data)) {
        stop("`fit` is a multimixMN_fit but `fit$model$data` is NULL. Refit with save_data = TRUE or pass `data=`.")
      }
    }

    if (is.null(fit$layer_fits) || length(fit$layer_fits) == 0) {
      stop("`fit$layer_fits` is missing/empty: cannot compute per-layer scores.")
    }

    # if layer is NULL -> all layers
    if (is.null(layer)) {
      out_list <- lapply(names(fit$layer_fits), function(L) {
        community_scores(
          fit   = fit$layer_fits[[L]],
          data  = data,
          layer = NULL,
          scale = scale,
          quantile_level = quantile_level,
          return_quantile_region = return_quantile_region,
          na_action = na_action
        )
      })
      names(out_list) <- names(fit$layer_fits)
      return(out_list)
    }

    # otherwise: one layer
    layer_names <- names(fit$layer_fits)
    if (is.numeric(layer) && length(layer) == 1L) {
      if (layer < 1 || layer > length(layer_names)) stop("`layer` index is out of range.")
      layer <- layer_names[layer]
    } else {
      layer <- as.character(layer)[1]
      if (!layer %in% layer_names) {
        stop("`layer` not found in `fit$layer_fits`. Available layers: ", paste(layer_names, collapse = ", "))
      }
    }

    return(
      community_scores(
        fit   = fit$layer_fits[[layer]],
        data  = data,
        scale = scale,
        quantile_level = quantile_level,
        return_quantile_region = return_quantile_region,
        na_action = na_action
      )
    )
  }

  # ---- loadings ----
  L_true <- fit$community_loadings$true
  nodes  <- fit$community_loadings$nodes
  wc     <- fit$community_loadings$wc

  if (is.null(L_true) || is.null(nodes) || length(nodes) == 0) {
    stop("No community loadings found in `fit$community_loadings`. Did you run mixMN(compute_loadings = TRUE)?")
  }
  if (is.null(wc) || length(wc) != length(nodes)) {
    stop("`fit$community_loadings$wc` is missing or has wrong length: cannot zero cross-loadings.")
  }

  # Ensure dimnames
  if (is.null(rownames(L_true))) rownames(L_true) <- nodes
  if (is.null(colnames(L_true))) colnames(L_true) <- paste0("C", seq_len(ncol(L_true)))

  # helper: set cross-loadings to 0 using hard membership wc
  .zero_cross_loadings <- function(Lmat, nodes, wc) {
    Lmat <- Lmat[nodes, , drop = FALSE]
    K <- ncol(Lmat)

    wc_int <- as.integer(wc)
    if (any(is.na(wc_int))) stop("`wc` must be coercible to integers.")
    if (any(wc_int < 1 | wc_int > K)) {
      stop("`wc` contains community indices outside [1, K].")
    }

    Lhard <- matrix(0, nrow = nrow(Lmat), ncol = K,
                    dimnames = list(rownames(Lmat), colnames(Lmat)))

    for (j in seq_len(nrow(Lmat))) {
      k <- wc_int[j]
      Lhard[j, k] <- Lmat[j, k]  # keep only within-community loading
    }
    Lhard
  }

  # apply hardening to true loadings (cross-loadings -> 0)
  L_true <- .zero_cross_loadings(L_true, nodes = nodes, wc = wc)

  # ---- choose data ----
  used_fit_data <- FALSE
  if (is.null(data)) {
    data <- fit$model$data
    used_fit_data <- TRUE
    if (is.null(data)) {
      stop("No `data` provided and `fit$model$data` is NULL. Refit with save_data = TRUE or pass `data=`.")
    }
  }

  # ---- coerce and ids ----
  if (!is.data.frame(data) && !is.matrix(data)) data <- as.data.frame(data)
  if (is.null(colnames(data))) stop("`data` must have column names.")
  if (is.null(rownames(data))) rownames(data) <- sprintf("id_%d", seq_len(nrow(data)))
  ids <- rownames(data)

  # ---- variable checks ----
  missing_vars <- setdiff(nodes, colnames(data))
  if (length(missing_vars) > 0) {
    stop("`data` is missing required variables: ", paste(missing_vars, collapse = ", "))
  }

  X <- as.matrix(data[, nodes, drop = FALSE])

  # ---- NA handling ----
  if (na_action == "stop" && anyNA(X)) {
    stop("Missing values detected in required variables. Use `na_action = \"omit\"` if you want row-wise omission.")
  }

  # ---- scaling on the scoring data (default) ----
  center_vec <- rep(0, ncol(X)); names(center_vec) <- colnames(X)
  scale_vec  <- rep(1, ncol(X)); names(scale_vec)  <- colnames(X)

  if (isTRUE(scale)) {
    # compute mean/sd on the dataset used for scoring
    center_vec <- apply(X, 2, function(v) mean(v, na.rm = TRUE))
    scale_vec  <- apply(X, 2, function(v) stats::sd(v, na.rm = TRUE))
    scale_vec[is.na(scale_vec) | scale_vec == 0] <- 1
    X <- sweep(X, 2, center_vec, "-")
    X <- sweep(X, 2, scale_vec, "/")
  }

  # ---- score computation ----
  # Scores = X %*% L_true (n x p) %*% (p x K) -> (n x K)
  # We need L_true aligned to nodes
  L_true <- L_true[nodes, , drop = FALSE]

  compute_scores_matrix <- function(Xmat, Lmat) {
    # Xmat: n x p ; Lmat: p x K
    S <- Xmat %*% Lmat
    colnames(S) <- colnames(Lmat)
    S
  }

  # If NA and na_action = "omit": compute per community, per row with renormalized weights
  if (na_action == "omit" && anyNA(X)) {
    K <- ncol(L_true)
    S <- matrix(NA_real_, nrow = nrow(X), ncol = K,
                dimnames = list(ids, colnames(L_true)))

    for (k in seq_len(K)) {
      w <- L_true[, k]
      for (i in seq_len(nrow(X))) {
        xi <- X[i, ]
        ok <- !is.na(xi) & !is.na(w)
        if (!any(ok)) {
          S[i, k] <- NA_real_
        } else {
          # re-normalize weights to keep scale stable (optional; reasonable default)
          w_ok <- w[ok]
          # if all weights are 0, keep unnormalized
          denom <- sum(abs(w_ok))
          if (is.na(denom) || denom == 0) denom <- 1
          S[i, k] <- sum(xi[ok] * w_ok) / denom
        }
      }
    }
    scores <- S
    rownames(scores) <- ids
  } else {
    scores <- compute_scores_matrix(X, L_true)
    rownames(scores) <- ids
  }

  # ---- quantile regions (optional) ----
  quantile_region_out <- NULL
  if (isTRUE(return_quantile_region)) {
    if (!is.numeric(quantile_level) || length(quantile_level) != 1L ||
        is.na(quantile_level) || quantile_level <= 0 || quantile_level >= 1) {
      stop("`quantile_level` must be a single number strictly between 0 and 1 (e.g., 0.95).")
    }

    L_boot <- fit$community_loadings$boot
    if (is.null(L_boot) || length(L_boot) == 0) {
      stop("Quantile region requested but `fit$community_loadings$boot` is NULL/empty. Run mixMN with boot_what including \"loadings\".")
    }

    alpha <- 1 - quantile_level
    probs <- c(alpha/2, 1 - alpha/2)

    reps <- length(L_boot)
    n <- nrow(X)
    K <- ncol(L_true)

    # container to store boot scores as a list then stack

    S_boot <- array(NA_real_, dim = c(reps, n, K),
                    dimnames = list(NULL, ids, colnames(L_true)))
    for (r in seq_len(reps)) {
      Lr <- L_boot[[r]]
      if (is.null(Lr)) next

      if (is.null(colnames(Lr))) colnames(Lr) <- colnames(L_true)
      if (is.null(rownames(Lr))) rownames(Lr) <- nodes

      # harden bootstrap loadings: cross-loadings -> 0
      Lr <- .zero_cross_loadings(Lr, nodes = nodes, wc = wc)

      if (na_action == "omit" && anyNA(X)) {
        Sr <- matrix(NA_real_, nrow = n, ncol = K,
                     dimnames = list(ids, colnames(L_true)))

        for (k in seq_len(K)) {
          w <- Lr[, k]
          for (i in seq_len(n)) {
            xi <- X[i, ]
            ok <- !is.na(xi) & !is.na(w)
            if (!any(ok)) {
              Sr[i, k] <- NA_real_
            } else {
              w_ok <- w[ok]
              denom <- sum(abs(w_ok))
              if (is.na(denom) || denom == 0) denom <- 1
              Sr[i, k] <- sum(xi[ok] * w_ok) / denom
            }
          }
        }
      } else {
        Sr <- compute_scores_matrix(X, Lr)
        rownames(Sr) <- ids
      }

      S_boot[r, , ] <- Sr
    }

    # Quantiles over reps for each subject/community cell
    lower <- matrix(NA_real_, nrow = n, ncol = K, dimnames = list(ids, colnames(L_true)))
    upper <- matrix(NA_real_, nrow = n, ncol = K, dimnames = list(ids, colnames(L_true)))

    for (i in seq_len(n)) {
      for (k in seq_len(K)) {
        v <- S_boot[, i, k]
        if (all(is.na(v))) {
          lower[i, k] <- NA_real_
          upper[i, k] <- NA_real_
        } else {
          qs <- stats::quantile(v, probs = probs, na.rm = TRUE, names = FALSE)
          lower[i, k] <- qs[1]
          upper[i, k] <- qs[2]
        }
      }
    }

    quantile_region_out <- list(lower = lower, upper = upper, quantile_level = quantile_level)
  }

  out <- list(
    call = match.call(),
    settings = list(
      scale = isTRUE(scale),
      quantile_level = quantile_level,
      na_action = na_action
    ),
    ids = ids,
    communities = colnames(L_true),
    scores = scores,
    quantile_region = quantile_region_out,
    details = list(
      nodes_used = nodes,
      loadings_true = L_true,
      loadings_boot_available = !is.null(fit$community_loadings$boot) && length(fit$community_loadings$boot) > 0,
      center = center_vec,
      scale = scale_vec
    )
  )
  class(out) <- c("mixmashnet", "community_scores")
  return(out)
}
