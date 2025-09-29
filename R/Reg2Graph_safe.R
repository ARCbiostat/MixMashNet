#' Safe graph assembly for mgm outputs (modified from \pkg{mgm})
#'
#' @description
#' Robust post-processing to assemble the pairwise (and nodewise) weighted adjacency
#' matrices from an \code{mgm}-style object. This variant tolerates masked designs
#' and empty parameter blocks by using conservative fallbacks and \code{na.rm = TRUE}
#' where appropriate. It preserves the structure of the original \pkg{mgm} output.
#'
#' @param mgmobj An object produced by a nodewise \pkg{mgm}-style routine.
#' @param thresholding Logical; if \code{TRUE}, applies the thresholding rule recorded in
#'   \code{mgmobj$call$threshold}.
#' @return The input object with populated \code{$pairwise} matrices and interaction lists.
#' @importFrom utils combn
#' @export
Reg2Graph_safe <- function(mgmobj, thresholding = TRUE) {

  # Prefer mgm's internal helpers if present; otherwise use safe fallbacks
  .FlagSymmetricFast <- function(x) {
    if (requireNamespace("mgm", quietly = TRUE) &&
        exists("FlagSymmetricFast", envir = asNamespace("mgm"), inherits = FALSE)) {
      return(.mgm_fun("FlagSymmetricFast")(x))
    } else {
      # Fallback: permutation-invariant row id (sort indices inside each row)
      keys <- apply(as.matrix(x), 1, function(r) paste(sort(as.numeric(r)), collapse = "_"))
      return(as.numeric(factor(keys)))
    }
  }

  .getSign_safe <- function(l_w_ind, l_w_par, type, set_signdefined, overparameterize, ord) {
    if (requireNamespace("mgm", quietly = TRUE) &&
        exists("getSign", envir = asNamespace("mgm"), inherits = FALSE)) {
      return(.mgm_fun("getSign")(l_w_ind, l_w_par, type, set_signdefined, overparameterize, ord))
    } else {
      # Conservative default when sign computation is unavailable
      return(list(voteSign = 0))
    }
  }

  # --- local helpers --------------------------------------------------------
  # stringr::str_count-compatible helper without adding a dependency
  if (!exists("str_count", mode = "function")) {
    str_count <- function(x, pattern) {
      sapply(regmatches(x, gregexpr(pattern, x, fixed = TRUE)), length)
    }
  }
  .safe_rbind <- function(lst) {
    lst <- Filter(function(el) !is.null(el) && length(el) > 0, lst)
    if (!length(lst)) return(matrix(numeric(0), nrow = 0))
    do.call(rbind, lst)
  }
  .safe_concat <- function(lst) {
    lst <- Filter(function(el) !is.null(el), lst)
    if (!length(lst)) return(list())
    do.call(c, lst)
  }
  .all_nonfinite <- function(x) {
    length(x) == 0L || all(!is.finite(x))
  }

  # ----- Basic info from model object ------
  p <- length(mgmobj$call$type)
  n <- mgmobj$call$n
  moderators <- mgmobj$call$moderators
  d <- mgmobj$call$k - 1
  if (!is.null(moderators)) d <- 2
  type <- mgmobj$call$type
  level <- mgmobj$call$level
  threshold <- mgmobj$call$threshold
  npar_standard <- mgmobj$call$npar
  binarySign <- mgmobj$call$binarySign
  ruleReg <- mgmobj$call$ruleReg
  overparameterize <- mgmobj$call$overparameterize
  ind_cat <- mgmobj$call$ind_cat
  ind_binary <- mgmobj$call$ind_binary

  # Determine moderation specification (if applicable)
  if (d < 2) {
    mSpec <- NULL
  } else {
    if (any(class(moderators) %in% c("integer","numeric"))) mSpec <- "vector"
    if (any(class(moderators) %in% "matrix")) mSpec <- "matrix"
    if (mgmobj$call$k > 2) mSpec <- "fullk"
  }

  # ----------------------------------------------------------------------------
  # -------------------- Fill glmnet output into list-structure ----------------
  # ----------------------------------------------------------------------------

  Pars_ind <- list()
  Pars_values <- list()
  mgmobj$intercepts <- vector('list', length = p)
  v_d_v <- rep(NA, p)

  for (v in 1:p) {

    model_obj <- mgmobj$nodemodels[[v]]$model
    predictor_set <- (1:p)[-v]

    if (!is.null(moderators)) {
      if (mSpec == "vector") ind_m <- TRUE
      if (mSpec == "matrix") {
        m_ind_m <- apply(moderators, 1, function(x) v %in% x)
        ind_m <- any(m_ind_m)
      }
      d_v <- if (ind_m) 2 else 1
    } else {
      ind_m <- FALSE
      d_v <- 1
    }
    if (mgmobj$call$k > 2) { d_v <- d; ind_m <- TRUE }

    v_Pars_ind <- vector('list', length = d_v)

    if (!ind_m) {
      v_Pars_ind[[1]] <- t(combn(predictor_set, 1))
    } else {
      if (mSpec == "vector") {
        if (v %in% moderators) {
          ind_mods <- t(combn((1:p)[-v], 2))
        } else {
          ind_mods <- expand.grid((1:p)[-v], moderators[moderators != v])
        }
        ind_mods <- ind_mods[!apply(ind_mods, 1, function(x) x[1] == x[2]), ]
        v_Pars_ind[[1]] <- matrix(predictor_set, ncol = 1)
        v_Pars_ind[[2]] <- ind_mods
      }
      if (mSpec == "matrix") {
        ind_v_inMod <- as.logical(apply(moderators, 1, function(x) v %in% x))
        ind_mods <- t(apply(matrix(moderators[ind_v_inMod], ncol = 3), 1, function(x) x[x != v]))
        v_Pars_ind[[1]] <- matrix(predictor_set, ncol = 1)
        v_Pars_ind[[2]] <- ind_mods
      }
      if (mSpec == "fullk") {
        for (ord in 1:d_v) v_Pars_ind[[ord]] <- t(combn(predictor_set, ord))
      }
    }
    for (ord in 1:d_v) v_Pars_ind[[ord]] <- matrix(as.matrix(v_Pars_ind[[ord]]), ncol = ord)
    no_interactions <- unlist(lapply(v_Pars_ind, nrow))

    v_Pars_values <- vector('list', length = d_v)
    for (ord in 1:d_v) v_Pars_values[[ord]] <- vector('list', length = no_interactions[ord])

    # ------ Fill (B) with parameter estimates ------
    if (type[v] == 'c') {
      n_cat <- level[v]
      for (cat in 1:n_cat) {
        model_obj_i <- as.numeric(model_obj[[cat]])
        interaction_names <- rownames(model_obj[[cat]])[-1]
        interaction_order <- str_count(interaction_names, ":") + 1
        mgmobj$intercepts[[v]][[cat]] <- model_obj_i[1]

        if (thresholding) {
          model_obj_i_ni <- model_obj_i[-1]
          if (threshold == 'LW') tau <- sqrt(d_v) * sqrt(sum(model_obj_i_ni^2)) * sqrt(log(npar_standard[v]) / n)
          if (threshold == 'HW') tau <- d_v * sqrt(log(npar_standard[v]) / n)
          if (threshold == 'none') tau <- 0
          model_obj_i[abs(model_obj_i) < tau] <- 0
          mgmobj$nodemodels[[v]]$tau <- tau
        }

        for (ord in 1:d_v) {
          part_ord <- model_obj_i[-1][ord == interaction_order]
          gotThemall <- rep(TRUE, length(part_ord))
          int_names_ord <- interaction_names[ord == interaction_order]
          if (no_interactions[ord] != nrow(v_Pars_ind[[ord]])) stop('Internal Error: parameter extraction type 1')

          find_int_dummy <- matrix(NA, nrow = length(int_names_ord), ncol = ord)
          for (t in 1:no_interactions[ord]) {
            for (cc in 1:ord) find_int_dummy[, cc] <- grepl(pattern = paste0('V', v_Pars_ind[[ord]][t, cc], '.'),
                                                            x = int_names_ord, fixed = TRUE)
            select_int <- rowSums(find_int_dummy) == ord
            parmat <- matrix(part_ord[select_int], ncol = 1)
            rownames(parmat) <- int_names_ord[select_int]
            v_Pars_values[[ord]][[t]][[cat]] <- parmat
            gotThemall[select_int == TRUE] <- FALSE
          }
          if (sum(gotThemall) > 0) stop('Internal Error: parameter extraction type 2')
        }
      }

    } else { # Continuous
      model_obj_i <- as.numeric(model_obj)
      interaction_names <- rownames(model_obj)[-1]
      interaction_order <- str_count(interaction_names, ":") + 1
      mgmobj$intercepts[[v]] <- model_obj_i[1]

      if (thresholding) {
        model_obj_i_ni <- model_obj_i[-1]
        if (threshold == 'LW') tau <- sqrt(d_v) * sqrt(sum(model_obj_i_ni^2)) * sqrt(log(npar_standard[v]) / n)
        if (threshold == 'HW') tau <- d_v * sqrt(log(npar_standard[v]) / n)
        if (threshold == 'none') tau <- 0
        model_obj_i[abs(model_obj_i) < tau] <- 0
        mgmobj$nodemodels[[v]]$tau <- tau
      }

      for (ord in 1:d_v) {
        if (no_interactions[ord] > 0) {
          part_ord <- model_obj_i[-1][ord == interaction_order]
          gotThemall <- rep(TRUE, length(part_ord))
          int_names_ord <- interaction_names[ord == interaction_order]
          if (no_interactions[ord] != nrow(v_Pars_ind[[ord]])) stop('Internal Error: parameter extraction 1')

          find_int_dummy <- matrix(NA, nrow = length(int_names_ord), ncol = ord)
          for (t in 1:no_interactions[ord]) {
            for (cc in 1:ord) find_int_dummy[, cc] <- grepl(pattern = paste0('V', v_Pars_ind[[ord]][t, cc], '.'),
                                                            x = int_names_ord, fixed = TRUE)
            select_int <- rowSums(find_int_dummy) == ord
            parmat <- matrix(part_ord[select_int], ncol = 1)
            rownames(parmat) <- int_names_ord[select_int]
            v_Pars_values[[ord]][[t]] <- parmat
            gotThemall[select_int == TRUE] <- FALSE
          }
          if (sum(gotThemall) > 0) stop('Internal Error: parameter extraction 2')
        }
      }
    }

    Pars_ind[[v]] <- v_Pars_ind
    Pars_values[[v]] <- v_Pars_values
    v_d_v[v] <- d_v
  }

  # ----------------------------------------------------------------------------
  # -------------------- Postprocess into (Factor) Graph -----------------------
  # ----------------------------------------------------------------------------

  Pars_ind_flip <- vector('list', length = d)
  dummy_list_flip <- vector('list', length = p)
  Pars_ind_flip <- lapply(Pars_ind_flip, function(x) dummy_list_flip)

  Pars_values_flip <- vector('list', length = d)
  Pars_values_flip <- lapply(Pars_values_flip, function(x) dummy_list_flip)

  for (v in 1:p) for (ord in 1:v_d_v[v]) {
    Pars_ind_part <- Pars_ind[[v]][[ord]]
    colnames(Pars_ind_part) <- NULL
    if (!is.null(Pars_ind_part)) {
      Pars_ind_part <- as.matrix(Pars_ind_part)
      Pars_ind[[v]][[ord]] <- cbind(rep(v, nrow(Pars_ind_part)), Pars_ind_part)
      Pars_ind_flip[[ord]][[v]] <- Pars_ind[[v]][[ord]]
      Pars_values_flip[[ord]][[v]] <- Pars_values[[v]][[ord]]
    } else {
      Pars_ind_flip[[ord]][[v]] <- NULL
      Pars_values_flip[[ord]][[v]] <- NULL
    }
  }

  # SAFE collapse
  Pars_ind_flip_red    <- lapply(Pars_ind_flip,    .safe_rbind)
  Pars_values_flip_red <- lapply(Pars_values_flip, .safe_concat)

  # Moderation?
  n_terms_d <- rep(NA, d)
  if (!is.null(moderators)) {
    n_terms_d[1] <- choose(p, 1+1)
    if (mSpec == "vector") {
      mod_terms <- expand.grid((1:p), (1:p), moderators)
      id_uni <- .FlagSymmetricFast(mod_terms)
      mod_terms2 <- mod_terms[!duplicated(id_uni), ]
      ind_diff <- as.numeric(apply(mod_terms2, 1, function(x) !any(duplicated(x))))
      mod_terms3 <- mod_terms2[ind_diff == 1, ]
    }
    if (mSpec == "matrix") {
      mod_terms3 <- mgmobj$call$moderators
    }
    n_terms_d[2] <- nrow(mod_terms3)
  } else {
    for (ord in 1:d) n_terms_d[ord] <- choose(p, ord+1)
  }

  l_factors <- lapply(1:d, function(ord) matrix(NA, nrow = n_terms_d[ord], ncol = ord+1))
  l_factor_par <- lapply(1:d, function(ord) vector('list', length = n_terms_d[ord]))
  l_sign_par <- lapply(1:d, function(ord) rep(NA, n_terms_d[ord]))
  l_factor_par_full <- l_factor_par_AggNodewise <- l_factor_par_SignNodewise <- l_factor_par

  if (binarySign) set_signdefined <- c(which(type == 'p'), which(type == 'g'), ind_cat[ind_binary])
  else           set_signdefined <- c(which(type == 'p'), which(type == 'g'))

  # Loop over order (ord = 1 for pairwise)
  for (ord in 1:d) {

    set_int_ord <- Pars_ind_flip_red[[ord]]
    if (is.null(set_int_ord) || length(set_int_ord) == 0) next
    row.names(set_int_ord) <- NULL

    ids <- .FlagSymmetricFast(x = set_int_ord)
    unique_set_int_ord <- cbind(set_int_ord, ids)[!duplicated(ids), ]
    unique_set_int_ord <- matrix(unique_set_int_ord, ncol = ord+1+1)
    n_unique <- if (length(unique_set_int_ord)) nrow(unique_set_int_ord) else 0
    if (n_unique == 0) next

    for (i in 1:n_unique) {

      l_w_ind <- lapply(1:(ord+1), function(j) set_int_ord[which(ids == i)[j], ])
      l_w_par <- lapply(1:(ord+1), function(j) Pars_values_flip_red[[ord]][[which(ids == i)[j]]])

      # Map regression -> edge parameter (mean of |beta|) with guards
      m_par_seq <- vapply(l_w_par, function(x) {
        x <- abs(unlist(x))
        if (.all_nonfinite(x)) return(NA_real_)
        mean(x, na.rm = TRUE)
      }, numeric(1))

      m_sign_seq <- vapply(l_w_par, function(x) {
        x <- unlist(x)
        if (length(x) > 1) 0 else if (length(x) == 1) sign(x) else NA_real_
      }, numeric(1))

      valid <- is.finite(m_par_seq)
      if (!any(valid)) {
        parcompr <- 0
        int_sign <- NA
      } else {
        m_par_mean <- mean(m_par_seq[valid], na.rm = TRUE)
        if (ruleReg == 'AND') parcompr <- m_par_mean * !(0 %in% m_par_seq[valid])
        if (ruleReg == 'OR')  parcompr <- m_par_mean

        if (!isTRUE(all.equal(m_par_mean, 0))) {
          pair <- l_w_ind[[1]]
          if (sum(!(pair %in% set_signdefined)) == 0) {
            sign_object <- .getSign_safe(
              l_w_ind, l_w_par, type, set_signdefined, overparameterize, ord
            )
            int_sign <- sign_object$voteSign
          } else {
            int_sign <- 0
          }
        } else {
          int_sign <- NA
        }
      }

      l_sign_par[[ord]][i] <- int_sign
      l_factors[[ord]][i, ] <- l_w_ind[[1]]

      for (i_ord in 1:(ord+1)) {
        l_factor_par_AggNodewise[[ord]][[i]][[i_ord]]  <- m_par_seq[i_ord]
        l_factor_par_SignNodewise[[ord]][[i]][[i_ord]] <- m_sign_seq[i_ord]
      }

      l_factor_par[[ord]][[i]]      <- parcompr
      l_factor_par_full[[ord]][[i]] <- l_w_par
    }
  }

  # -------------------- Weighted Adjacency (pairwise) -------------------

  l_factors_nz <- l_factors
  l_factor_par_nz <- l_factor_par
  l_factor_par_full_nz <- l_factor_par_full
  l_sign_par_nz <- l_sign_par

  for (ord in 1:d) {
    zero_indicator <- which(unlist(lapply(l_factor_par[[ord]], function(x) isTRUE(all.equal(x, 0)) )))
    if (length(zero_indicator) == 0) {
      l_factors_nz[[ord]] <- l_factors[[ord]]
      l_sign_par_nz[[ord]] <- l_sign_par[[ord]]
    } else {
      l_factors_nz[[ord]] <- l_factors[[ord]][-zero_indicator, , drop = FALSE]
      l_sign_par_nz[[ord]] <- l_sign_par[[ord]][-zero_indicator]
    }
    l_factor_par_nz[[ord]] <- l_factor_par_full_nz[[ord]] <- list()
    counter <- 1
    for (k in seq_along(l_factor_par[[ord]])) {
      if (!(k %in% zero_indicator)) {
        l_factor_par_nz[[ord]][[counter]]      <- l_factor_par[[ord]][[k]]
        l_factor_par_full_nz[[ord]][[counter]] <- l_factor_par_full[[ord]][[k]]
        counter <- counter + 1
      }
    }
  }

  mgmobj$interactions$indicator   <- l_factors_nz
  mgmobj$interactions$weightsAgg  <- l_factor_par_nz
  mgmobj$interactions$weights     <- l_factor_par_full_nz
  mgmobj$interactions$signs       <- l_sign_par_nz

  # ---------- p x p Matrix (pairwise) ----------
  m_signs <- matrix(NA, p, p)
  wadj <- matrix(0, p, p)

  edges_mat <- l_factors_nz[[1]]
  n_edges <- if (is.null(edges_mat) || length(edges_mat) == 0) 0 else nrow(as.matrix(edges_mat))

  if (n_edges > 0) {
    edges <- matrix(edges_mat, ncol = 2)
    for (i in 1:n_edges) {
      wadj[edges[i,1], edges[i,2]] <- wadj[edges[i,2], edges[i,1]] <- l_factor_par_nz[[1]][[i]]
      m_signs[edges[i,1], edges[i,2]] <- m_signs[edges[i,2], edges[i,1]] <- l_sign_par_nz[[1]][[i]]
    }
  }

  sign_colors <- matrix('darkgrey', p, p)
  sign_colors[m_signs == 1] <- 'darkgreen'
  sign_colors[m_signs == -1] <- 'red'

  sign_colors_cb <- matrix('darkgrey', p, p)
  sign_colors_cb[m_signs == 1] <- 'darkblue'
  sign_colors_cb[m_signs == -1] <- 'red'

  sign_ltys <- matrix(1, p, p)
  sign_ltys[m_signs == -1] <- 2

  mgmobj$pairwise$wadj        <- wadj
  mgmobj$pairwise$signs       <- m_signs
  mgmobj$pairwise$edgecolor   <- sign_colors
  mgmobj$pairwise$edgecolor_cb<- sign_colors_cb
  mgmobj$pairwise$edge_lty    <- sign_ltys

  # ---------- p x p Nodewise Matrix ----------
  m_wadj <- m_signs_dir <- matrix(0, p, p)

  ord <- 1
  set_int_ord <- Pars_ind_flip_red[[ord]]
  if (!is.null(set_int_ord) && length(set_int_ord) > 0) {
    row.names(set_int_ord) <- NULL
    ids <- .FlagSymmetricFast(x = set_int_ord)
    unique_set_int_ord <- cbind(set_int_ord, ids)[!duplicated(ids), ]
    unique_set_int_ord <- matrix(unique_set_int_ord, ncol = ord+1+1)
    n_edges <- if (length(unique_set_int_ord)) nrow(unique_set_int_ord) else 0

    if (n_edges > 0) {
      ED <- unique_set_int_ord
      for (i in 1:n_edges) {
        m_wadj[ED[i,1], ED[i,2]] <- l_factor_par_AggNodewise[[1]][[i]][[2]]
        m_wadj[ED[i,2], ED[i,1]] <- l_factor_par_AggNodewise[[1]][[i]][[1]]
        m_signs_dir[ED[i,1], ED[i,2]] <- l_factor_par_SignNodewise[[1]][[i]][[2]]
        m_signs_dir[ED[i,2], ED[i,1]] <- l_factor_par_SignNodewise[[1]][[i]][[1]]
      }
    }
  }

  sign_colors <- matrix('darkgrey', p, p)
  sign_colors[m_signs_dir == 1] <- 'darkgreen'
  sign_colors[m_signs_dir == -1] <- 'red'

  sign_colors_cb <- matrix('darkgrey', p, p)
  sign_colors_cb[m_signs_dir == 1] <- 'darkblue'
  sign_colors_cb[m_signs_dir == -1] <- 'red'

  sign_ltys <- matrix(1, p, p)
  sign_ltys[m_signs_dir == -1] <- 2

  mgmobj$pairwise$wadjNodewise         <- m_wadj
  mgmobj$pairwise$signsNodewise        <- m_signs_dir
  mgmobj$pairwise$edgecolorNodewise    <- sign_colors
  mgmobj$pairwise$edgecolor_cbNodewise <- sign_colors_cb
  mgmobj$pairwise$edge_ltyNodewise     <- sign_ltys

  return(mgmobj)
}
