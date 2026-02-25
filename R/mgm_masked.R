#' Masked MGM (modified from \pkg{mgm})
#'
#' @description
#' A modified version of the model-fitting routine from \pkg{mgm} that adds support
#' for a per-node predictor mask via \code{mask_list}. When provided, each node is
#' estimated using only the allowed predictors specified for that node. All other
#' functionality mirrors the original \pkg{mgm} nodewise estimation and graph
#' assembly.
#'
#' @param data Numeric matrix (n × p). No missing values allowed.
#' @param type Character vector of length \code{p} with variable types as in \pkg{mgm}
#'   (e.g., \code{"g"} for Gaussian, \code{"c"} for categorical, \code{"p"} for Poisson).
#' @param level Integer vector of length \code{p} with variable levels as in \pkg{mgm}.
#' @param regularize Logical; if \code{FALSE}, equivalent to no regularization
#'   (\code{lambdaSel = "EBIC"}, \code{lambdaSeq = 0}, \code{threshold = "none"}).
#' @param lambdaSeq,lambdaSel,lambdaFolds,lambdaGam Lambda grid/selection settings (see \pkg{mgm}).
#' @param alphaSeq,alphaSel,alphaFolds,alphaGam Mixing parameter grid/selection (see \pkg{mgm}).
#' @param k Interaction order (default \code{2}); same meaning as in \pkg{mgm}.
#' @param moderators Optional moderators specification (as in \pkg{mgm}).
#' @param ruleReg Regularization rule, default \code{"AND"} (as in \pkg{mgm}).
#' @param weights Observation weights (length \code{n}); internally normalized.
#' @param threshold Thresholding rule (default \code{"LW"}) as in \pkg{mgm}.
#' @param method Fitting backend (default \code{"glm"}) as in \pkg{mgm}.
#' @param binarySign Logical; if \code{TRUE}, store sign information.
#' @param scale Logical; if \code{TRUE}, standardize Gaussian variables.
#' @param verbatim,pbar,saveModels,saveData Overhead/UX flags as in \pkg{mgm}.
#' @param overparameterize Logical; use overparameterized design matrix (as in \pkg{mgm}).
#' @param thresholdCat Logical; categorical thresholding (as in \pkg{mgm}).
#' @param signInfo Logical; store sign information (as in \pkg{mgm}).
#' @param mask_list Optional list of length \code{p}; each element is an integer vector
#'   of column indices (1..p) indicating which predictors are allowed for the nodewise
#'   regression of the corresponding target node. If \code{NULL}, no masking is applied.
#' @param ... Further arguments forwarded consistently with \pkg{mgm}.
#'
#' @details
#' This function is adapted from the \pkg{mgm} workflow and preserves its nodewise
#' estimation, cross-validation / EBIC selection, thresholding, and graph assembly.
#' The only substantive extension is \code{mask_list}, which filters each node's design
#' matrix to the allowed predictors before fitting. This enables constrained structures
#' (e.g., multilayer masks) while keeping \pkg{mgm}’s estimation logic intact.
#'
#' @return An object of class \code{c("mgm","core")} with fields analogous to \pkg{mgm}:
#' \itemize{
#'   \item \code{$pairwise$}: \code{wadj}, \code{signs}, and nodewise variants.
#'   \item \code{$interactions$}, \code{$intercepts$}, \code{$nodemodels$}: as in \pkg{mgm}.
#'   \item \code{$call$}: input metadata (types, levels, thresholds, etc.).
#' }
#'
#' @seealso \code{\link[mgm]{mgm}} for the original implementation and argument semantics.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
#' @noRd
mgm_masked <- function(data, type, level,
                       regularize,
                       lambdaSeq, lambdaSel, lambdaFolds, lambdaGam,
                       alphaSeq, alphaSel, alphaFolds, alphaGam,
                       k,
                       moderators,
                       ruleReg,
                       weights,
                       threshold,
                       method,
                       binarySign,
                       scale,
                       verbatim,
                       pbar,
                       saveModels,
                       saveData,
                       overparameterize,
                       thresholdCat,
                       signInfo,
                       mask_list = NULL,
                       ...) {

  if (!("matrix" %in% class(data))) stop("Please provide data as matrix")
  if (any(is.na(data))) stop("No missing values permitted.")

  p <- ncol(data); n <- nrow(data)
  colnames(data)[1:p] <- paste("V", 1:p, '.', sep = "")

  args <- list(...)
  if (missing(regularize))        regularize <- TRUE
  if (missing(lambdaSeq))         lambdaSeq <- NULL
  if (missing(lambdaSel))         lambdaSel <- "CV"
  if (missing(lambdaFolds))       lambdaFolds <- 10
  if (missing(lambdaGam))         lambdaGam <- .25
  if (missing(alphaSeq))          alphaSeq <- 1
  if (missing(alphaSel))          alphaSel <- "CV"
  if (missing(alphaFolds))        alphaFolds <- 10
  if (missing(alphaGam))          alphaGam <- .25
  if (missing(k))                 k <- 2
  if (missing(moderators))        moderators <- NULL
  if (missing(ruleReg))           ruleReg <- "AND"
  if (missing(weights))           weights <- rep(1, n)
  if (missing(threshold))         threshold <- "LW"
  if (missing(method))            method <- "glm"
  if (missing(binarySign))        binarySign <- FALSE
  if (missing(scale))             scale <- TRUE
  if (missing(verbatim))          verbatim <- FALSE
  if (missing(pbar))              pbar <- TRUE
  if (missing(saveModels))        saveModels <- TRUE
  if (missing(saveData))          saveData <- FALSE
  if (missing(overparameterize))  overparameterize <- FALSE
  if (missing(signInfo))          signInfo <- TRUE
  if (missing(thresholdCat))      thresholdCat <- TRUE

  if (isTRUE(regularize == FALSE)) {
    lambdaSel <- "EBIC"
    lambdaSeq <- 0
    threshold <- "none"
  }

  if (!is.null(mask_list)) {
    if (!is.list(mask_list) || length(mask_list) != p)
      stop("mask_list must be a list of length p; each element is integer indices allowed for that target")
  }

  # REQUIREMENTS (no :::)
  .mgm_fun("glmnetRequirements")(data = data, type = type, weights = weights)

  d <- k - 1
  emp_lev <- rep(1, p)
  ind_cat_all <- which(type == "c")
  if (length(ind_cat_all) > 0) {
    for (i in seq_along(ind_cat_all)) {
      emp_lev[ind_cat_all][i] <- length(unique(data[, ind_cat_all[i]]))
    }
  }
  if (missing(level)) level <- emp_lev

  if (!missing(weights)) weights <- weights / max(weights)
  nadj <- sum(weights)

  ind_Gauss <- which(type == "g")
  if (scale && length(ind_Gauss) > 0) {
    for (i in ind_Gauss) data[, i] <- scale(data[, i])
  }

  unique_cats <- vector("list", length = p)
  for (i in seq_len(p)) if (type[i] == "c") unique_cats[[i]] <- unique(data[, i])

  data_df <- as.data.frame(data)
  for (i in which(type == "c")) data_df[, i] <- as.factor(data_df[, i])

  mgmobj <- list(
    call = NULL,
    pairwise = list(wadj = NULL, signs = NULL, edgecolor = NULL,
                    wadjNodewise = NULL, signsNodewise = NULL, edgecolorNodewise = NULL),
    interactions = list(indicator = NULL, weightsAgg = NULL, weights = NULL, signs = NULL),
    intercepts = NULL,
    nodemodels = list()
  )

  mgmobj$call <- list(
    data = if (isTRUE(saveData)) data_df else NULL,
    type = type,
    level = level,
    levelNames = NULL,
    regularize = regularize,
    lambdaSeq = lambdaSeq,
    lambdaSel = lambdaSel,
    lambdaFolds = lambdaFolds,
    lambdaGam = lambdaGam,
    alphaSeq = alphaSeq,
    alphaSel = alphaSel,
    alphaFolds = alphaFolds,
    alphaGam = alphaGam,
    k = k,
    moderators = moderators,
    ruleReg = ruleReg,
    weights = weights,
    threshold = threshold,
    method = method,
    binarySign = binarySign,
    scale = scale,
    verbatim = verbatim,
    pbar = pbar,
    saveModels = saveModels,
    saveData = saveData,
    overparameterize = overparameterize,
    thresholdCat = thresholdCat,
    signInfo = signInfo,
    npar = NULL,
    n = nrow(data_df),
    ind_cat = which(type == "c"),
    ind_binary = {
      indc <- which(type == "c")
      if (length(indc) == 0) logical(0) else {
        x <- logical(length(indc))
        for (i in seq_along(indc)) x[i] <- length(unique(data_df[, indc[i]])) == 2
        indc[x]
      }
    },
    unique_cats = {
      uc <- vector("list", length = ncol(data_df))
      for (i in seq_len(ncol(data_df))) if (type[i] == "c") uc[[i]] <- unique(data_df[, i])
      uc
    }
  )

  if (isTRUE(pbar)) pb <- txtProgressBar(min = 0, max = p, initial = 0, char = "-", style = 3)
  npar_standard <- rep(NA_integer_, p)

  .vars_in_term <- function(term) {
    parts <- unlist(strsplit(term, ":", fixed = TRUE))
    as.integer(sub("\\..*$", "", sub("^V", "", parts)))
  }
  .apply_mask <- function(X, v, allowed_preds_v) {
    allowed <- sort(unique(c(v, allowed_preds_v)))
    keep <- vapply(colnames(X), function(tt) {
      ids <- .vars_in_term(tt)
      all(ids %in% allowed)
    }, logical(1))
    X[, keep, drop = FALSE]
  }

  for (v in 1:p) {

    # Model matrix (no :::)
    X_standard <- X <- .mgm_fun("ModelMatrix_standard")(
      data = data_df, type = type, d = d, v = v, moderators = moderators
    )

    if (isTRUE(overparameterize)) {
      X <- .mgm_fun("ModelMatrix")(
        data = data_df, type = type, level = level, labels = colnames(data_df),
        d = d, moderators = moderators, v = v
      )
      X_standard <- X
    }

    if (!is.null(mask_list)) {
      X_standard <- .apply_mask(X_standard, v = v, allowed_preds_v = mask_list[[v]])
      X <- X_standard
    }
    npar_standard[v] <- ncol(X_standard)

    if (ncol(X) == 0) {
      mgmobj$nodemodels[[v]] <- list(
        fitobj = NULL, lambda = NA_real_, EBIC = NA_real_,
        parlist = NULL, modeltype = NULL
      )
      if (isTRUE(pbar)) setTxtProgressBar(pb, v)
      next
    }

    if (scale && any(type == "g")) {
      cn <- colnames(X)
      l_ints_split <- strsplit(cn, ":")
      ind_allc <- unlist(lapply(l_ints_split, function(x) {
        x2 <- sub("V", "", x)
        vars <- as.numeric(unlist(lapply(strsplit(x2, "[.]"), function(y) y[1])))
        all(type[vars] == "g")
      }))
      if (sum(ind_allc) > 0) {
        X[, ind_allc] <- apply(matrix(X[, ind_allc], ncol = sum(ind_allc)), 2, scale)
      }
    }

    y <- as.numeric(data_df[, v])
    n_alpha <- length(alphaSeq)

    if (alphaSel == "CV") {
      ind_cv <- sample(1:alphaFolds, size = n, replace = TRUE)
      v_mean_OOS_deviance <- rep(NA_real_, n_alpha)

      if (n_alpha > 1) {
        for (a in 1:n_alpha) {
          v_OOS_deviance <- rep(NA_real_, alphaFolds)
          for (fold in 1:alphaFolds) {
            train_X <- X[ind_cv != fold, , drop = FALSE]
            train_y <- y[ind_cv != fold]
            test_X  <- X[ind_cv == fold, , drop = FALSE]
            test_y  <- y[ind_cv == fold]
            n_train <- nrow(train_X); nadj_train <- sum(weights[ind_cv != fold])

            fit_fold <- .mgm_fun("nodeEst")(
              y = train_y, X = train_X,
              lambdaSeq = lambdaSeq, lambdaSel = lambdaSel,
              lambdaFolds = lambdaFolds, lambdaGam = lambdaGam,
              alpha = alphaSeq[a], weights = weights[ind_cv != fold],
              n = n_train, nadj = nadj_train, v = v, type = type,
              level = level, emp_lev = emp_lev,
              overparameterize = overparameterize, thresholdCat = thresholdCat
            )

            LL_model <- .mgm_fun("calcLL")(
              X = test_X, y = test_y, fit = fit_fold$fitobj, type = type, level = level,
              v = v, weights = weights[ind_cv == fold], lambda = fit_fold$lambda, LLtype = "model"
            )
            LL_sat <- .mgm_fun("calcLL")(
              X = test_X, y = test_y, fit = fit_fold$fitobj, type = type, level = level,
              v = v, weights = weights[ind_cv == fold], lambda = fit_fold$lambda, LLtype = "saturated"
            )
            v_OOS_deviance[fold] <- 2 * (LL_sat - LL_model)
          }
          v_mean_OOS_deviance[a] <- mean(v_OOS_deviance)
        }
        alpha_select <- alphaSeq[which.min(v_mean_OOS_deviance)]
      } else {
        alpha_select <- alphaSeq
      }

      model <- .mgm_fun("nodeEst")(
        y = y, X = X,
        lambdaSeq = lambdaSeq, lambdaSel = lambdaSel,
        lambdaFolds = lambdaFolds, lambdaGam = lambdaGam,
        alpha = alpha_select, weights = weights,
        n = n, nadj = nadj, v = v, type = type, level = level,
        emp_lev = emp_lev, overparameterize = overparameterize, thresholdCat = thresholdCat
      )
      mgmobj$nodemodels[[v]] <- model
    }

    if (alphaSel == "EBIC") {
      l_alphaModels <- vector("list", n_alpha)
      EBIC_Seq <- rep(NA_real_, n_alpha)
      for (a in 1:n_alpha) {
        fit_a <- .mgm_fun("nodeEst")(
          y = y, X = X,
          lambdaSeq = lambdaSeq, lambdaSel = lambdaSel,
          lambdaFolds = lambdaFolds, lambdaGam = lambdaGam,
          alpha = alphaSeq[a], weights = weights,
          n = n, nadj = nadj, v = v, type = type, level = level,
          emp_lev = emp_lev, overparameterize = overparameterize, thresholdCat = thresholdCat
        )
        l_alphaModels[[a]] <- fit_a
        EBIC_Seq[a] <- fit_a$EBIC
      }
      ind_min <- which.min(EBIC_Seq)
      mgmobj$nodemodels[[v]] <- l_alphaModels[[ind_min]]
    }

    if (isTRUE(pbar)) setTxtProgressBar(pb, v)
  }

  mgmobj$call$npar <- npar_standard

  levelNames <- vector("list", length = p)
  for (i in 1:p) {
    if (type[i] == "c") levelNames[[i]] <- sort(as.numeric(as.character(unique(data_df[, i])))) else levelNames[[i]] <- NA
  }
  mgmobj$call$levelNames <- levelNames

  mgmobj <- Reg2Graph_safe(mgmobj = mgmobj)
  if (!saveModels) mgmobj$nodemodels <- NULL

  class(mgmobj) <- c("mgm", "core")
  return(mgmobj)
}
