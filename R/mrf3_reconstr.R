# -------------------------------------------------------------------------------------------------------------
# Reconstruction clustering method
# -------------------------------------------------------------------------------------------------------------

#' MRF unsupervised clustering -- Reconstruction method
#' @param rfit A model list of random forest models.
#' @param recon Optional reconstruction object returned by `get_reconstr_matrix()`.
#' If `NULL`, reconstruction is built from `rfit`.
#' @param k Pre-defined number of clusters. The default is selecting the optimal k by tuning method.
#' @param weights_cluster Logical; whether to cluster based on reconstructed weight affinity.
#' @param dat_use_to_cluster If clustering method selected as Reconstr, the user can choose the data to conduct clustering.
#' In weight-clustering mode, shared clustering uses `W_all`; this argument only
#' applies when `weights_cluster = FALSE`.
#' @param model_top_v Model-level top-v cutoff applied to each single-model
#' forest weight matrix before fusion.
#' @param recon_fusion Reconstruction fusion mode.
#' `"weighted"` (default) uses `connection_score`; `"uniform"` uses equal weights.
#' @param connection_score Optional directional score matrix from `find_connection(return_score = TRUE)`.
#' @param score_power Exponent applied to raw connection scores before normalization.
#' @param score_floor Non-negative floor applied to raw connection scores before normalization.
#' @param fallback_uniform Logical; whether to fallback to uniform averaging when weighted scores are unavailable.
#' @param fused_top_v Optional integer. If set, apply row-wise top-v truncation to fused weights
#' (`W_all`) after fusion. If `FALSE`, no fused top-v truncation is applied.
#' @param fused_row_normalize Logical; whether to row-normalize fused weights after optional truncation.
#' @param fused_keep_ties Logical; whether fused top-v truncation keeps ties at cutoff.
#' @param fusion_mode Clustering-affinity fusion mode.
#' `"average"` keeps current behavior; `"ao_reg"` uses alternating optimization with simplex-regularized weights.
#' @param gamma L2 regularization strength for the simplex weights in `"ao_reg"` mode.
#' @param alpha_init Optional initial fusion weights for `"ao_reg"` mode.
#' @param ao_max_iter Maximum alternating-optimization iterations.
#' @param ao_tol Convergence tolerance for `||alpha_t - alpha_{t-1}||_1`.
#' @param knn_q Optional kNN sparsification per base similarity (rows keep top-q).
#' @param hollow Logical; whether to zero out diagonal in base/fused similarities.
#' @param ao_symm Logical; whether to symmetrize kNN-sparsified similarities by `pmax(S, t(S))`.
#' @param ao_verbose Logical; whether to print AO iteration diagnostics.
#' @param ... Additional arguments passed to clustering backends.
#'
#' @return mrf3 clustering object
#' @export mrf3_reconstr
mrf3_reconstr <- function(recon = NULL,
                          rfit = NULL,
                          k = NULL,
                          weights_cluster = TRUE,
                          dat_use_to_cluster = "ALL",
                          model_top_v = 10,
                          recon_fusion = c("weighted", "uniform"),
                          global_fusion = c("average", "pmin"),
                          connection_score = NULL,
                          score_power = 1,
                          score_floor = 0,
                          fallback_uniform = TRUE,
                          fused_top_v = NULL,
                          fused_row_normalize = TRUE,
                          fused_keep_ties = TRUE,
                          fusion_mode = "average",
                          gamma = 0.1,
                          alpha_init = NULL,
                          ao_max_iter = 20,
                          ao_tol = 1e-4,
                          knn_q = NULL,
                          hollow = TRUE,
                          ao_symm = TRUE,
                          ao_verbose = FALSE,
                          ...){

  if (is.null(recon)) {
    if (is.null(rfit)) {
      stop("Either `recon` or `rfit` must be provided.")
    }
    recon <- get_reconstr_matrix(
      rfit = rfit,
      model_top_v = model_top_v,
      recon_fusion = recon_fusion,
      global_fusion = global_fusion,
      connection_score = connection_score,
      score_power = score_power,
      score_floor = score_floor,
      fallback_uniform = fallback_uniform,
      fused_top_v = fused_top_v,
      fused_row_normalize = fused_row_normalize,
      fused_keep_ties = fused_keep_ties
    )
  } else if (!is.list(recon) || is.null(recon$W) || is.null(recon$fused_mat)) {
    stop("`recon` must be a reconstruction object returned by `get_reconstr_matrix()`.")
  }

  embed <- NULL
  fusion_mode <- tolower(fusion_mode)
  ao_fit <- NULL
  k_tuned <- is.null(k)
  dat_used_actual <- dat_use_to_cluster

  if (weights_cluster) {
    if (fusion_mode == "ao_reg") {
      if (is.null(k)) {
        stop("For fusion_mode = 'ao_reg', please provide k.")
      }
      W_models <- recon$W$W_models
      ao_fit <- ao_fuse_similarity(
        W_models = W_models,
        k = k,
        gamma = gamma,
        alpha_init = alpha_init,
        max_iter = ao_max_iter,
        tol = ao_tol,
        knn_q = knn_q,
        hollow = hollow,
        symm = ao_symm,
        verbose = ao_verbose
      )
      dat <- ao_fit$S_fused
      embed <- ao_fit$U
      clm <- pam_cl(embed, k_tune = k, diss = FALSE)
      cl <- clm$cl
      method <- "AO_REG_PAM"
      dat_used_actual <- "W_models"
    } else {
      dat_ls <- recon$W
      dat <- dat_ls$W_all
      dat <- dat %*% t(dat)
      method <- "Spectral"
      dat_used_actual <- "ALL"
    }
  } else {
    dat_ls <- recon$fused_mat
    if (toupper(dat_use_to_cluster) == "ALL") {
      dat <- Reduce(cbind, dat_ls)
    } else {
      if (!dat_use_to_cluster %in% names(dat_ls)) {
        stop("`dat_use_to_cluster` is not in reconstructed matrices.")
      }
      dat <- dat_ls[[dat_use_to_cluster]]
    }
    method <- "PAM"
    dat_used_actual <- dat_use_to_cluster
  }

  sample_names <- recon$sample_names
  if (is.null(sample_names) && !is.null(rfit)) {
    sample_names <- rownames(rfit[[1]]$xvar)
  }
  if (!is.null(sample_names)) {
    rownames(dat) <- sample_names
    if (ncol(dat) == length(sample_names)) {
      colnames(dat) <- sample_names
    }
  }

  if (!(weights_cluster && fusion_mode == "ao_reg")) {
    if (is.null(k)) {
      message("Start tuning k step..")
      clm <- tune_k_clusters(dat, return_cluster = TRUE, method = method)
      cl <- clm$cl
      if (method == "Spectral") {
        embed <- clm$embed
      }
    } else {
      if (method == "PAM") {
        clm <- pam_cl(dat, k_tune = k, diss = FALSE)
      }
      if (method == "Spectral") {
        clm <- spectral_cl(dat, k_tune = k)
        embed <- clm$embed
      }
      cl <- clm$cl
    }
  }

  message("Done!")
  k_selected <- if (!is.null(clm$best_k) && length(clm$best_k) > 0L) {
    as.integer(clm$best_k)[1]
  } else {
    as.integer(length(unique(cl)))
  }
  out <- list(
    dat = dat,
    cl = cl,
    k = k_selected,
    k_tuned = k_tuned,
    cl_mod = clm,
    embed = embed,
    fused_mat = recon$fused_mat,
    W_mat = recon$W,
    ao_fit = ao_fit,
    fusion_mode = fusion_mode,
    recon_fusion_mode = recon$fusion_mode,
    model_score = recon$model_score,
    model_fusion_weights = recon$model_fusion_weights,
    dat_used = dat_used_actual,
    method = "Reconstr"
  )

  class(out) <- "reconstr"
  out
}


#' @param rfit A model list of random forest models.
#' @rdname mrf3_reconstr
get_reconstr_matrix <- function(rfit,
                                model_top_v = 10,
                                recon_fusion = c("weighted", "uniform"),
                                global_fusion = c("average", "pmin"),
                                connection_score = NULL,
                                score_power = 1,
                                score_floor = 0,
                                fallback_uniform = TRUE,
                                fused_top_v = NULL,
                                fused_row_normalize = TRUE,
                                fused_keep_ties = TRUE){

  if (!is.list(rfit) || length(rfit) == 0) {
    stop("`rfit` must be a non-empty list of RF models.")
  }
  if (is.null(names(rfit)) || any(names(rfit) == "")) {
    stop("`rfit` must be a named list. Model names should follow `response_predictor`.")
  }
  if (isTRUE(model_top_v) || isFALSE(model_top_v) ||
      !is.numeric(model_top_v) || length(model_top_v) != 1L ||
      !is.finite(model_top_v) || model_top_v <= 0) {
    stop("`model_top_v` must be a single positive integer.")
  }
  model_top_v <- as.integer(model_top_v)
  if (isTRUE(fused_top_v) || isFALSE(fused_top_v)) {
    fused_top_v <- NULL
  }
  if (!is.null(fused_top_v)) {
    if (!is.numeric(fused_top_v) || length(fused_top_v) != 1L ||
        !is.finite(fused_top_v) || fused_top_v <= 0) {
      stop("`fused_top_v` must be NULL or a single positive integer.")
    }
    fused_top_v <- as.integer(fused_top_v)
  }
  if (!is.numeric(score_power) || length(score_power) != 1L || !is.finite(score_power) || score_power <= 0) {
    stop("`score_power` must be a single positive numeric.")
  }
  if (!is.numeric(score_floor) || length(score_floor) != 1L || !is.finite(score_floor) || score_floor < 0) {
    stop("`score_floor` must be a single non-negative numeric.")
  }
  recon_fusion <- match.arg(recon_fusion)
  global_fusion <- match.arg(global_fusion)
  model_names <- names(rfit)
  n_models <- length(model_names)

  model_score <- rep(NA_real_, n_models)
  names(model_score) <- model_names
  weighted_sum <- NA_real_
  if (recon_fusion == "weighted") {
    model_score <- match_model_scores(model_names = model_names, connection_score = connection_score)
    model_score <- pmax(model_score, score_floor)^score_power
    weighted_sum <- sum(model_score, na.rm = TRUE)
    ## normalize_fusion_weights() handles fallback to uniform when scores are invalid
  } else {
    model_score[] <- 1
    weighted_sum <- sum(model_score)
  }

  global_alpha <- normalize_fusion_weights(
    model_score,
    fallback_uniform = fallback_uniform
  )
  names(global_alpha) <- model_names

  ## ══════════════════════════════════════════════════════════════
  ## Step 1: Extract and preprocess per-connection W_m matrices
  ## ══════════════════════════════════════════════════════════════
  fused_info <- vector("list", length(model_names))
  names(fused_info) <- model_names
  W_models <- vector("list", length(model_names))
  names(W_models) <- model_names

  for (mi in seq_along(model_names)) {
    m <- model_names[mi]
    mod <- rfit[[m]]
    fw <- mod$forest.wt
    if (is.null(fw)) {
      stop("Model `", m, "` does not contain `forest.wt`.")
    }

    W <- prepare_weight_matrix(
      W = fw,
      adjust = TRUE,
      top_v = model_top_v,
      row_normalize = TRUE,
      zero_diag = TRUE,
      keep_ties = TRUE
    )

    ## Compute reconstruction products (keep these, they are smaller)
    xvar <- mod$xvar
    yvar <- mod$yvar
    mod_names <- parse_model_pair(m)

    if (!is.null(yvar)) {
      out <- list(
        W %*% as.matrix(xvar),
        W %*% as.matrix(yvar)
      )
      if (length(mod_names) >= 2L) {
        names(out) <- rev(mod_names[1:2])
      } else {
        names(out) <- c("X", "Y")
      }
    } else {
      out <- list(W %*% as.matrix(xvar))
      names(out) <- mod_names[1]
    }

    W_models[[m]] <- W
    fused_info[[m]] <- list(model = m, mat = out, W = W)
  }

  ## Legacy reconstruction products (for backward compat / tSNE)
  mat_records <- list()
  rec_i <- 1L
  for (m in model_names) {
    mats <- fused_info[[m]]$mat
    for (dat_name in names(mats)) {
      mat_records[[rec_i]] <- list(
        dat_name = dat_name,
        model = m,
        mat = mats[[dat_name]]
      )
      rec_i <- rec_i + 1L
    }
  }

  dat_names <- unique(vapply(mat_records, `[[`, FUN.VALUE = character(1), "dat_name"))
  full_fused_mat <- list()
  block_model_weights <- list()

  for (d in dat_names) {
    idx <- which(vapply(mat_records, `[[`, FUN.VALUE = character(1), "dat_name") == d)
    d_models <- vapply(mat_records[idx], `[[`, FUN.VALUE = character(1), "model")
    d_mats <- lapply(mat_records[idx], `[[`, "mat")
    raw <- model_score[d_models]
    d_alpha <- normalize_fusion_weights(raw, fallback_uniform = fallback_uniform)
    names(d_alpha) <- d_models

    full_fused_mat[[d]] <- fuse_matrix_list(d_mats, d_alpha)
    block_model_weights[[d]] <- d_alpha
  }

  ## ══════════════════════════════════════════════════════════════
  ## Step 2: Per-response W^(k) — only connections where k is response
  ## W^(k) = sum_{i != k} alpha_{ki} * W_{ki},  normalised within
  ## the response-k subset.  Used for shared-specific decomposition.
  ## ══════════════════════════════════════════════════════════════
  W_per_response <- vector("list", length(dat_names))
  names(W_per_response) <- dat_names
  per_response_alpha <- vector("list", length(dat_names))
  names(per_response_alpha) <- dat_names

  for (d in dat_names) {
    resp_idx <- vapply(model_names, function(m) {
      pair <- parse_model_pair(m)
      identical(pair[1], d)
    }, logical(1))
    resp_models <- model_names[resp_idx]

    if (length(resp_models) == 0L) {
      warning(
        "No connection has block `", d, "` as response. ",
        "Falling back to uniform W for W_per_response.",
        call. = FALSE
      )
      ## Use uniform average of all W_m as fallback
      fallback_alpha <- rep(1 / length(model_names), length(model_names))
      names(fallback_alpha) <- model_names
      W_fb <- fuse_matrix_list(W_models, fallback_alpha)
      rs <- rowSums(W_fb); rs[rs <= 0 | !is.finite(rs)] <- 1
      W_per_response[[d]] <- W_fb / rs
      per_response_alpha[[d]] <- fallback_alpha
    } else {
      resp_scores <- model_score[resp_models]
      resp_scores[!is.finite(resp_scores)] <- 1.0
      resp_scores <- pmax(resp_scores, 0)
      resp_alpha <- normalize_fusion_weights(resp_scores, fallback_uniform = fallback_uniform)
      names(resp_alpha) <- resp_models
      per_response_alpha[[d]] <- resp_alpha

      W_k <- fuse_matrix_list(W_models[resp_models], resp_alpha)
      rs <- rowSums(W_k); rs[rs <= 0 | !is.finite(rs)] <- 1
      W_per_response[[d]] <- W_k / rs
    }
  }

  ## ══════════════════════════════════════════════════════════════
  ## Step 3: Global W_M* — built FROM per-response W^(k) matrices
  ## Two modes:
  ##   "average" : W_M* = (1/K) * sum_k W^(k)  + optional fused_top_v
  ##   "pmin"    : W_M*(i,j) = min_k W^(k)(i,j) + row-normalise
  ##               (cross-modal intersection; fused_top_v not needed)
  ## ══════════════════════════════════════════════════════════════
  K <- length(W_per_response)
  if (global_fusion == "pmin") {
    ## Element-wise minimum across all per-response matrices
    ## Samples are "shared-similar" only if ALL modalities agree
    W_all <- W_per_response[[1]]
    if (K > 1L) {
      for (ki in seq_len(K)[-1]) {
        W_all <- pmin(W_all, W_per_response[[ki]])
      }
    }
    ## Row-normalise to maintain linear smoother property
    rs <- rowSums(W_all); rs[rs <= 0 | !is.finite(rs)] <- 1
    W_all <- W_all / rs
  } else {
    ## "average": uniform average of per-response matrices
    W_all <- Reduce("+", W_per_response) / K
    ## Apply optional fused_top_v and row-normalisation
    W_all <- postprocess_fused_weight(
      W = W_all,
      top_v = fused_top_v,
      row_normalize = fused_row_normalize,
      keep_ties = fused_keep_ties
    )
  }

  list(
    fused_mat = full_fused_mat,
    W = list(
      W_all = W_all,
      W_models = W_models,
      W_per_response = W_per_response
    ),
    sample_names = rownames(rfit[[1]]$xvar),
    fusion_mode = if (recon_fusion == "weighted" && (!is.finite(weighted_sum) || weighted_sum <= 0)) {
      "uniform"
    } else {
      recon_fusion
    },
    global_fusion = global_fusion,
    model_score = model_score,
    model_fusion_weights = global_alpha,
    per_response_alpha = per_response_alpha,
    block_model_weights = block_model_weights
  )
}

postprocess_fused_weight <- function(W,
                                     top_v = NULL,
                                     row_normalize = TRUE,
                                     keep_ties = TRUE) {
  if (isFALSE(top_v)) {
    top_v <- NULL
  }
  out <- W
  if (!is.null(top_v)) {
    out <- truncate_top_v_rows(out, top_v = top_v, keep_ties = keep_ties)
  }
  if (isTRUE(row_normalize)) {
    out <- row_normalize_weights(out)
  }
  out
}

fuse_matrix_list <- function(mat_list, alpha) {
  if (length(mat_list) == 0L) {
    stop("Cannot fuse an empty matrix list.")
  }
  if (length(mat_list) == 1L) {
    return(mat_list[[1]])
  }
  if (length(alpha) != length(mat_list)) {
    stop("`alpha` length must match `mat_list` length.")
  }
  ## Sequential accumulation: only 2 matrices in memory at a time
  ## (the running sum + the current matrix being weighted).
  ## Previous implementation used Reduce("+", Map(...)) which

  ## materialised all weighted matrices simultaneously.
  out <- mat_list[[1]] * alpha[[1]]
  for (k in seq_along(mat_list)[-1]) {
    out <- out + mat_list[[k]] * alpha[[k]]
  }
  out
}

normalize_fusion_weights <- function(score, fallback_uniform = TRUE) {
  score <- as.numeric(score)
  score[!is.finite(score)] <- NA_real_
  score[score < 0] <- 0
  s <- sum(score, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    if (!fallback_uniform) {
      stop("No valid weighted-fusion score is available.")
    }
    return(rep(1 / length(score), length(score)))
  }
  score[is.na(score)] <- 0
  score / sum(score)
}

parse_model_pair <- function(model_name) {
  split <- as.character(stringr::str_split(model_name, "_", simplify = TRUE))
  split <- split[split != ""]
  if (length(split) == 0L) {
    return(model_name)
  }
  if (length(split) == 1L) {
    return(split[1])
  }
  split[1:2]
}

match_model_scores <- function(model_names, connection_score = NULL) {
  out <- rep(NA_real_, length(model_names))
  names(out) <- model_names
  if (is.null(connection_score)) {
    return(out)
  }
  if (!is.matrix(connection_score) || is.null(rownames(connection_score)) || is.null(colnames(connection_score))) {
    stop("`connection_score` must be a matrix with row/column names.")
  }

  for (i in seq_along(model_names)) {
    pair <- parse_model_pair(model_names[i])
    if (length(pair) < 2L) {
      out[i] <- NA_real_
      next
    }
    resp <- pair[1]
    pred <- pair[2]
    if (resp %in% rownames(connection_score) && pred %in% colnames(connection_score)) {
      out[i] <- connection_score[resp, pred]
    } else {
      out[i] <- NA_real_
    }
  }
  out
}

project_simplex <- function(v) {
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u) - 1
  rho <- max(which(u - cssv/seq_along(u) > 0))
  theta <- cssv[rho]/rho
  pmax(v - theta, 0)
}

knn_sparsify <- function(S, q, symm = TRUE) {
  n <- nrow(S)
  if (is.null(q) || q >= n) return(S)
  q <- max(1L, as.integer(q))
  out <- matrix(0, nrow = n, ncol = n, dimnames = dimnames(S))
  for (i in seq_len(n)) {
    idx <- order(S[i, ], decreasing = TRUE)[seq_len(q)]
    out[i, idx] <- S[i, idx]
  }
  if (symm) out <- pmax(out, t(out))
  out
}

build_similarity_laplacian <- function(W, knn_q = NULL, hollow = TRUE, symm = TRUE) {
  S <- W %*% t(W)
  if (hollow) diag(S) <- 0
  S <- knn_sparsify(S, q = knn_q, symm = symm)
  if (hollow) diag(S) <- 0
  D <- rowSums(S)
  L <- diag(D) - S
  list(S = S, L = L)
}

ao_fuse_similarity <- function(W_models, k, gamma = 0.1, alpha_init = NULL,
                               max_iter = 20, tol = 1e-4, knn_q = NULL,
                               hollow = TRUE, symm = TRUE, verbose = FALSE) {
  if (length(W_models) == 0L) stop("W_models is empty.")
  if (k < 1L) stop("k must be >= 1.")
  if (gamma <= 0) stop("gamma must be > 0.")

  SL <- lapply(W_models, build_similarity_laplacian, knn_q = knn_q, hollow = hollow, symm = symm)
  S_list <- lapply(SL, `[[`, "S")
  L_list <- lapply(SL, `[[`, "L")
  m <- length(L_list)

  if (is.null(alpha_init)) {
    alpha <- rep(1 / m, m)
  } else {
    if (length(alpha_init) != m) stop("alpha_init length must equal number of base similarities.")
    alpha <- project_simplex(as.numeric(alpha_init))
  }

  obj_path <- numeric(max_iter)
  U <- NULL
  for (iter in seq_len(max_iter)) {
    L_alpha <- Reduce("+", Map(`*`, alpha, L_list))
    eig <- eigen(L_alpha, symmetric = TRUE)
    ord <- order(eig$values, decreasing = FALSE)
    U <- eig$vectors[, ord[seq_len(k)], drop = FALSE]

    c_vec <- vapply(L_list, function(Lm) sum((Lm %*% U) * U), numeric(1))
    alpha_new <- project_simplex(-c_vec / (2 * gamma))
    obj <- sum(alpha_new * c_vec) + gamma * sum(alpha_new^2)
    obj_path[iter] <- obj

    if (verbose) {
      message(sprintf("AO iter %d: obj=%.6f, d_alpha=%.6e", iter, obj, sum(abs(alpha_new - alpha))))
    }
    if (sum(abs(alpha_new - alpha)) < tol) {
      alpha <- alpha_new
      obj_path <- obj_path[seq_len(iter)]
      break
    }
    alpha <- alpha_new
  }
  if (length(obj_path) == max_iter && obj_path[max_iter] == 0) {
    obj_path <- obj_path[obj_path != 0]
  }
  S_fused <- Reduce("+", Map(`*`, alpha, S_list))
  if (hollow) diag(S_fused) <- 0

  list(S_fused = S_fused, alpha = alpha, U = U, obj_path = obj_path)
}

resolve_top_v_values <- function(dat_input,
                                 mod_input,
                                 connection_input,
                                 connection_score,
                                 model_top_v_input,
                                 fused_top_v_input,
                                 disable_fused_top_v = FALSE,
                                 shared_k_for_tune = NULL,
                                 top_v_method = c("entropy_elbow", "neff"),
                                 neff_quantile = 0.5,
                                 model_top_v_tune_args = list(),
                                 fused_top_v_tune_args = list(),
                                 stage_prefix = "[Stage 3/5]",
                                 verbose = TRUE) {
  top_v_method <- match.arg(top_v_method)
  tuning <- list(
    model_top_v = NULL,
    fused_top_v = NULL
  )
  tune_mod <- list(
    mod = mod_input,
    connection = connection_input,
    connection_score = connection_score,
    recon = NULL
  )

  model_use <- model_top_v_input
  if (identical(model_use, Inf)) {
    # Inf = no truncation (use all neighbors)
    n_samples <- nrow(dat_input[[1]])
    model_use <- as.integer(n_samples)
    if (verbose) message("  model_top_v: all ", n_samples, " neighbors.")
  } else if (is.null(model_use)) {
    if (identical(top_v_method, "neff")) {
      # Fast path: select v from effective neighbourhood size
      if (verbose) message("Selecting model_top_v via effective neighbourhood size..")
      rfit <- if (!is.null(mod_input$mod)) mod_input$mod else mod_input
      # Compute adjusted weight matrices and pick v from their neff
      neff_vals <- vapply(rfit, function(m) {
        fw <- m$forest.wt
        if (is.null(fw)) return(NA_real_)
        W_adj <- prepare_weight_matrix(fw, adjust = TRUE, top_v = NULL,
                                       row_normalize = TRUE, zero_diag = TRUE)
        stats::median(effective_neighbourhood_size(W_adj), na.rm = TRUE)
      }, numeric(1))
      neff_vals <- neff_vals[is.finite(neff_vals)]
      if (length(neff_vals) == 0L) {
        model_use <- as.integer(nrow(dat_input[[1]]))
        if (verbose) message("  model_top_v: fallback to all ", model_use, " neighbors.")
      } else {
        model_use <- as.integer(ceiling(stats::quantile(neff_vals, probs = neff_quantile, na.rm = TRUE)))
        model_use <- max(model_use, 10L)
        if (verbose) message("  model_top_v = ", model_use, " (neff method, median neff across connections)")
      }
    } else {
      if (verbose) message("Tuning model_top_v..")
      tune_defaults <- list(
        dat.list = dat_input,
        mod = tune_mod,
        tmin = max(10L, as.integer(ceiling(0.08 * nrow(dat_input[[1]])))),
        object = "entropy_elbow"
      )
      if (is.null(model_top_v_tune_args$k) && !is.null(shared_k_for_tune)) {
        tune_defaults$k <- shared_k_for_tune
      }
      final_tune_args <- utils::modifyList(tune_defaults, model_top_v_tune_args)
      final_tune_args$object <- "entropy_elbow"
      tuning$model_top_v <- do.call(tune_model_top_v, final_tune_args)
      model_use <- pick_tuned_value(
        tb = tuning$model_top_v$tmax_tb,
        column = "model_top_v",
        allow_infinite = FALSE
      )
      if (verbose) message("  model_top_v = ", model_use)
    }
  } else {
    model_use <- as.integer(model_use)
  }

  fused_use <- if (disable_fused_top_v) NULL else fused_top_v_input
  if (identical(fused_use, Inf)) {
    # Inf = no truncation
    fused_use <- NULL
    if (verbose) message("  fused_top_v: no truncation.")
  } else if (is.null(fused_use) && !disable_fused_top_v) {
    if (identical(top_v_method, "neff")) {
      # Fast path: build fused weight matrix, then select v from neff
      if (verbose) message("Selecting fused_top_v via effective neighbourhood size..")
      rfit <- if (!is.null(mod_input$mod)) mod_input$mod else mod_input
      alpha <- compute_tune_model_alpha(
        model_names = names(rfit),
        connection_score = connection_score
      )
      cache_list <- build_model_weight_cache_list(rfit, vmax = model_use)
      W_fused <- build_fused_weight_from_cache(
        cache_list = cache_list,
        top_v = model_use,
        alpha = alpha,
        keep_ties = TRUE
      )
      W_fused <- postprocess_fused_weight(
        W = W_fused, top_v = NULL, row_normalize = TRUE, keep_ties = TRUE
      )
      fused_use <- select_top_v_neff(W_fused, quantile_prob = neff_quantile, min_v = 10L)
      # If neff is close to n, skip truncation
      if (fused_use >= as.integer(0.8 * ncol(W_fused))) {
        fused_use <- NULL
        if (verbose) message("  fused_top_v: no truncation (neff >= 80% of n).")
      } else {
        if (verbose) message("  fused_top_v = ", fused_use, " (neff method)")
      }
    } else {
      if (verbose) message("Tuning fused_top_v..")
      tune_defaults <- list(
        dat.list = dat_input,
        mod = tune_mod,
        model_top_v = model_use,
        vmin = min(10L, nrow(dat_input[[1]])),
        object = "entropy_elbow"
      )
      if (is.null(fused_top_v_tune_args$k) && !is.null(shared_k_for_tune)) {
        tune_defaults$k <- shared_k_for_tune
      }
      final_tune_args <- utils::modifyList(tune_defaults, fused_top_v_tune_args)
      final_tune_args$object <- "entropy_elbow"
      tuning$fused_top_v <- do.call(tune_fused_top_v, final_tune_args)
      fused_use <- pick_tuned_value(
        tb = tuning$fused_top_v$vtop_tb,
        column = "fused_top_v",
        allow_infinite = TRUE
      )
      if (is.null(fused_use)) {
        if (verbose) message("  fused_top_v: no truncation.")
      } else {
        if (verbose) message("  fused_top_v = ", fused_use)
      }
    }
  } else if (!is.null(fused_use)) {
    fused_use <- as.integer(fused_use)
  }

  list(
    model_top_v = model_use,
    fused_top_v = fused_use,
    tuning = tuning
  )
}
