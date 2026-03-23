#' Build Similarity Matrix From Forest Weights
#'
#' @param W A square forest-weight matrix.
#' @param similarity_type Similarity structure:
#' `"second"` uses `W %*% t(W)`;
#' `"first"` uses `(W + t(W))/2`.
#' @param zero_diag Logical; whether to zero out diagonal entries in similarity.
#' @param symm Logical; whether to force symmetry by averaging with transpose.
#'
#' @return A sample-by-sample similarity matrix.
build_similarity_from_weights <- function(W,
                                          similarity_type = c("second", "first"),
                                          zero_diag = TRUE,
                                          symm = TRUE) {
  W <- validate_weight_matrix(W)

  similarity_type <- match.arg(similarity_type)
  if (identical(similarity_type, "first")) {
    S <- (W + t(W)) / 2
  } else {
    S <- W %*% t(W)
    if (isTRUE(symm)) {
      S <- (S + t(S)) / 2
    }
  }
  if (isTRUE(zero_diag)) {
    diag(S) <- 0
  }
  S
}


#' Cluster A Similarity Matrix
#'
#' @param S A square similarity matrix.
#' @param k Optional integer. If `NULL`, the function tunes `k`.
#' @param method Clustering backend: `"Spectral"` or `"PAM"`.
#' @param tune_method Tuning criterion for PAM when `k` is `NULL`.
#' @param gap_w Weighting scheme for spectral eigengap when `k` is `NULL`.
#' @param ... Additional arguments passed to clustering backends.
#'
#' @return A list with `cl`, `cl_mod`, `k`, `method`, and `embed`.
cluster_similarity_matrix <- function(S,
                                      k = NULL,
                                      method = c("PAM", "Spectral"),
                                      tune_method = "silhouette",
                                      gap_w = "uniform",
                                      ...) {
  method <- match.arg(method)
  S <- as.matrix(S)
  if (!is.numeric(S) || nrow(S) != ncol(S)) {
    stop("`S` must be a square numeric matrix.")
  }

  if (is.null(k)) {
    clm <- tune_k_clusters(
      S,
      return_cluster = TRUE,
      method = method,
      tune_method = tune_method,
      gap_w = gap_w,
      prox = identical(method, "PAM"),
      ...
    )
  } else {
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2) {
      stop("`k` must be NULL or a single integer >= 2.")
    }
    k <- as.integer(k)
    if (identical(method, "Spectral")) {
      clm <- spectral_cl(S, k_tune = k, ...)
    } else {
      clm <- pam_cl(1 - S, k_tune = k, diss = TRUE, tune_method = tune_method, ...)
    }
  }

  list(
    cl = clm$cl,
    cl_mod = clm,
    k = clm$best_k,
    method = method,
    embed = if (!is.null(clm$embed)) clm$embed else NULL
  )
}


#' Cluster Shared Similarity From Reconstruction Weights
#'
#' @param recon Reconstruction output from `get_reconstr_matrix()`.
#' @param mode Shared clustering mode: `"average"` or `"ao_reg"`.
#' @param dat_use Shared weight source for `"average"` mode. Shared clustering
#' now always uses `W_all`; this argument is kept for compatibility.
#' @param k Optional number of clusters. Required for `"ao_reg"` mode.
#' @param method Clustering backend for `"average"` mode.
#' @param tune_method Tuning criterion for PAM in `"average"` mode when `k` is `NULL`.
#' @param gap_w Weighting scheme for spectral eigengap in `"average"` mode when `k` is `NULL`.
#' @param similarity_type Similarity structure used in `"average"` mode.
#' @param gamma L2 regularization strength for `"ao_reg"` simplex weights.
#' @param alpha_init Optional initial AO fusion weights.
#' @param ao_max_iter Maximum AO iterations.
#' @param ao_tol AO convergence tolerance.
#' @param knn_q Optional kNN sparsification per base similarity in AO mode.
#' @param hollow Logical; whether to zero diagonal in AO mode.
#' @param ao_symm Logical; whether to symmetrize kNN sparsification in AO mode.
#' @param ao_verbose Logical; whether to print AO diagnostics.
#' @param ... Additional arguments passed to clustering backends.
#'
#' @return A list with shared similarity, clustering result, and AO details (if used).
cluster_shared_similarity <- function(recon,
                                      mode = c("average", "ao_reg"),
                                      dat_use = "ALL",
                                      k = NULL,
                                      method = c("PAM", "Spectral"),
                                      tune_method = "silhouette",
                                      gap_w = "uniform",
                                      similarity_type = c("second", "first"),
                                      gamma = 0.1,
                                      alpha_init = NULL,
                                      ao_max_iter = 20,
                                      ao_tol = 1e-4,
                                      knn_q = NULL,
                                      hollow = TRUE,
                                      ao_symm = TRUE,
                                      ao_verbose = FALSE,
                                      ...) {
  if (!is.list(recon) || is.null(recon$W) || is.null(recon$W$W_all)) {
    stop("`recon` must be a reconstruction object returned by `get_reconstr_matrix()`.")
  }

  mode <- match.arg(mode)
  method <- match.arg(method)
  similarity_type <- match.arg(similarity_type)

  if (identical(mode, "ao_reg")) {
    if (is.null(k)) {
      stop("`k` is required when `mode = 'ao_reg'`.")
    }
    W_models <- recon$W$W_models
    if (is.null(W_models) || length(W_models) == 0L) {
      stop("`recon$W$W_models` is required for `mode = 'ao_reg'`.")
    }
    ao_fit <- ao_fuse_similarity(
      W_models = W_models,
      k = as.integer(k),
      gamma = gamma,
      alpha_init = alpha_init,
      max_iter = ao_max_iter,
      tol = ao_tol,
      knn_q = knn_q,
      hollow = hollow,
      symm = ao_symm,
      verbose = ao_verbose
    )
    clm <- pam_cl(ao_fit$U, k_tune = as.integer(k), diss = FALSE, ...)
    return(list(
      mode = mode,
      dat_used = "W_models",
      similarity = ao_fit$S_fused,
      similarity_type = "ao_reg_fused",
      cl = clm$cl,
      cl_mod = clm,
      k = clm$best_k,
      method = "AO_REG_PAM",
      embed = ao_fit$U,
      ao_fit = ao_fit
    ))
  }

  W <- recon$W$W_all

  S <- build_similarity_from_weights(
    W = W,
    similarity_type = similarity_type,
    zero_diag = TRUE,
    symm = TRUE
  )
  cl_out <- cluster_similarity_matrix(
    S = S,
    k = k,
    method = method,
    tune_method = tune_method,
    gap_w = gap_w,
    ...
  )

  c(
    list(
      mode = mode,
      dat_used = "ALL",
      similarity = S,
      similarity_type = similarity_type,
      ao_fit = NULL
    ),
    cl_out
  )
}


resolve_specific_k <- function(k, dat_names) {
  out <- vector("list", length(dat_names))
  names(out) <- dat_names

  if (is.null(k)) {
    return(out)
  }

  if (is.list(k)) {
    for (d in dat_names) {
      out[[d]] <- k[[d]]
    }
    return(out)
  }

  if (length(k) == 1L && is.numeric(k)) {
    for (d in dat_names) {
      out[[d]] <- as.integer(k)
    }
    return(out)
  }

  if (is.numeric(k) && !is.null(names(k))) {
    for (d in dat_names) {
      if (d %in% names(k)) {
        out[[d]] <- as.integer(k[[d]])
      }
    }
    return(out)
  }

  stop("`k` for specific clustering must be NULL, a single integer, a named numeric vector, or a named list.")
}


#' Cluster Specific Similarities Per Omics Block
#'
#' @param shared_specific Output from `get_shared_specific_weights()`.
#' @param k Optional specific-clustering `k` configuration.
#' `NULL` tunes each block; single integer uses same `k` for all blocks;
#' named numeric/list allows block-specific `k`.
#' @param method Clustering backend for each specific block.
#' Choose from `"Spectral"`, `"PAM"`, `"Proximity"`, and `"Enhanced_Proximity"`.
#' @param prox_method_cl Proximity clustering backend used when
#' `method` is `"Proximity"` or `"Enhanced_Proximity"`.
#' @param tune_method Tuning criterion for PAM when block-level `k` is `NULL`.
#' @param gap_w Weighting scheme for spectral eigengap when block-level `k` is `NULL`.
#' @param similarity_type Similarity structure for each specific block.
#' @param ... Additional arguments passed to clustering backends.
#'
#' @return A list with per-omics specific similarities and clustering outputs.
cluster_specific_similarity <- function(shared_specific,
                                        k = NULL,
                                        method = c("PAM", "Spectral", "Proximity", "Enhanced_Proximity"),
                                        prox_method_cl = c("PAM", "Spectral"),
                                        tune_method = "silhouette",
                                        gap_w = "uniform",
                                        similarity_type = c("second", "first"),
                                        ...) {
  method <- match.arg(method)
  prox_method_cl <- match.arg(prox_method_cl)
  similarity_type <- match.arg(similarity_type)
  if (!is.list(shared_specific) ||
      is.null(shared_specific$specific) ||
      is.null(shared_specific$specific$W)) {
    stop("`shared_specific` must be the output of `get_shared_specific_weights()`.")
  }

  W_spec <- shared_specific$specific$W
  if (!is.list(W_spec) || length(W_spec) == 0L) {
    stop("`shared_specific$specific$W` must be a non-empty list.")
  }
  dat_names <- names(W_spec)
  if (is.null(dat_names) || any(dat_names == "")) {
    stop("`shared_specific$specific$W` must be named.")
  }
  residual_mod <- shared_specific$specific$residual_mod

  k_map <- resolve_specific_k(k, dat_names = dat_names)
  out <- vector("list", length(dat_names))
  names(out) <- dat_names

  for (d in dat_names) {
    if (method %in% c("Proximity", "Enhanced_Proximity")) {
      if (is.null(residual_mod) || is.null(residual_mod[[d]])) {
        stop(
          "Specific proximity clustering requires `shared_specific$specific$residual_mod[[", d, "]]`."
        )
      }
      prox_fit <- mrf3_cl_prox(
        rfit = list(residual_mod[[d]]),
        k = k_map[[d]],
        enhanced = identical(method, "Enhanced_Proximity"),
        method_cl = prox_method_cl,
        ...
      )
      out[[d]] <- list(
        dat_used = d,
        similarity = prox_fit$dat,
        similarity_type = if (identical(method, "Enhanced_Proximity")) {
          "enhanced_proximity"
        } else {
          "proximity"
        },
        cl = prox_fit$cl,
        cl_mod = prox_fit$cl_mod,
        k = length(unique(prox_fit$cl)),
        method = method,
        embed = NULL
      )
    } else {
      S <- build_similarity_from_weights(
        W = W_spec[[d]],
        similarity_type = similarity_type,
        zero_diag = TRUE,
        symm = TRUE
      )
      cl_out <- cluster_similarity_matrix(
        S = S,
        k = k_map[[d]],
        method = method,
        tune_method = tune_method,
        gap_w = gap_w,
        ...
      )
      out[[d]] <- c(
        list(
          dat_used = d,
          similarity = S,
          similarity_type = similarity_type
        ),
        cl_out
      )
    }
  }

  list(
    method = method,
    similarity_type = if (method %in% c("Proximity", "Enhanced_Proximity")) {
      if (identical(method, "Enhanced_Proximity")) "enhanced_proximity" else "proximity"
    } else {
      similarity_type
    },
    by_omics = out
  )
}


#' Joint Shared/Specific Similarity Clustering Wrapper
#'
#' @param recon Reconstruction output from `get_reconstr_matrix()`.
#' @param shared_specific Output from `get_shared_specific_weights()`.
#' @param shared_mode Shared clustering mode passed to `cluster_shared_similarity()`.
#' @param shared_dat_use Shared weight source for `"average"` mode.
#' @param shared_k Optional `k` for shared clustering. Required for `"ao_reg"` mode.
#' @param shared_method Clustering backend for shared `"average"` mode.
#' @param shared_similarity_type Similarity structure for shared `"average"` mode.
#' @param specific_k Optional `k` configuration passed to `cluster_specific_similarity()`.
#' @param specific_method Clustering backend for specific clustering.
#' Choose from `"Spectral"`, `"PAM"`, `"Proximity"`, and `"Enhanced_Proximity"`.
#' @param specific_prox_method_cl Proximity clustering backend used when
#' `specific_method` is `"Proximity"` or `"Enhanced_Proximity"`.
#' @param specific_similarity_type Similarity structure for specific clustering.
#' @param tune_method Tuning criterion for PAM when tuning `k`.
#' @param gap_w Weighting scheme for spectral eigengap when tuning `k`.
#' @param gamma L2 regularization for shared AO fusion.
#' @param alpha_init Optional initial AO fusion weights.
#' @param ao_max_iter Maximum AO iterations.
#' @param ao_tol AO convergence tolerance.
#' @param knn_q Optional kNN sparsification for shared AO fusion.
#' @param hollow Logical; whether to zero diagonal in shared AO fusion.
#' @param ao_symm Logical; whether to symmetrize kNN sparsification in AO fusion.
#' @param ao_verbose Logical; whether to print AO diagnostics.
#' @param ... Additional arguments passed to clustering backends.
#'
#' @return A list with `shared` and `specific` clustering outputs.
mrf3_specific_clustering <- function(recon,
                                   shared_specific,
                                   shared_mode = c("average", "ao_reg"),
                                   shared_dat_use = "ALL",
                                   shared_k = NULL,
                                   shared_method = c("PAM", "Spectral"),
                                   shared_similarity_type = c("second", "first"),
                                   specific_k = NULL,
                                   specific_method = c("PAM", "Spectral", "Proximity", "Enhanced_Proximity"),
                                   specific_prox_method_cl = c("PAM", "Spectral"),
                                   specific_similarity_type = c("second", "first"),
                                   tune_method = "silhouette",
                                   gap_w = "uniform",
                                   gamma = 0.1,
                                   alpha_init = NULL,
                                   ao_max_iter = 20,
                                   ao_tol = 1e-4,
                                   knn_q = NULL,
                                   hollow = TRUE,
                                   ao_symm = TRUE,
                                   ao_verbose = FALSE,
                                   ...) {
  shared_mode <- match.arg(shared_mode)
  shared_method <- match.arg(shared_method)
  shared_similarity_type <- match.arg(shared_similarity_type)
  specific_method <- match.arg(specific_method)
  specific_prox_method_cl <- match.arg(specific_prox_method_cl)
  specific_similarity_type <- match.arg(specific_similarity_type)

  shared_out <- cluster_shared_similarity(
    recon = recon,
    mode = shared_mode,
    dat_use = shared_dat_use,
    k = shared_k,
    method = shared_method,
    similarity_type = shared_similarity_type,
    tune_method = tune_method,
    gap_w = gap_w,
    gamma = gamma,
    alpha_init = alpha_init,
    ao_max_iter = ao_max_iter,
    ao_tol = ao_tol,
    knn_q = knn_q,
    hollow = hollow,
    ao_symm = ao_symm,
    ao_verbose = ao_verbose,
    ...
  )

  specific_out <- cluster_specific_similarity(
    shared_specific = shared_specific,
    k = specific_k,
    method = specific_method,
    prox_method_cl = specific_prox_method_cl,
    similarity_type = specific_similarity_type,
    tune_method = tune_method,
    gap_w = gap_w,
    ...
  )

  list(
    shared = shared_out,
    specific = specific_out
  )
}
