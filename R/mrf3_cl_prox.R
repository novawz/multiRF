# -------------------------------------------------------------------------------------------------------------
# Proximity clustering method
# -------------------------------------------------------------------------------------------------------------

reconstruct_plain_proximity <- function(mod,
                                        parallel = TRUE,
                                        cores = NULL,
                                        model_index = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package `Matrix` is required to reconstruct proximity.", call. = FALSE)
  }

  model_label <- if (is.null(model_index)) "" else paste0(" for model ", model_index)
  membership <- mod$membership
  if (!is.matrix(membership) || !is.numeric(membership) ||
      nrow(membership) < 1L || ncol(membership) < 1L) {
    stop(
      "Cannot reconstruct proximity", model_label,
      ": `membership` must be a non-empty numeric matrix.",
      call. = FALSE
    )
  }
  if (any(!is.finite(membership))) {
    stop(
      "Cannot reconstruct proximity", model_label,
      ": `membership` contains missing or non-finite values.",
      call. = FALSE
    )
  }
  if (any(membership < 0) || any(membership != round(membership))) {
    stop(
      "Cannot reconstruct proximity", model_label,
      ": `membership` must contain non-negative integer terminal-node IDs.",
      call. = FALSE
    )
  }
  if (!is.null(mod$ntree)) {
    if (length(mod$ntree) != 1L || !is.finite(mod$ntree) ||
        mod$ntree < 1 || mod$ntree != as.integer(mod$ntree)) {
      stop(
        "Cannot reconstruct proximity", model_label,
        ": `ntree` must be a single positive integer.",
        call. = FALSE
      )
    }
    if (ncol(membership) != as.integer(mod$ntree)) {
      stop(
        "Cannot reconstruct proximity", model_label,
        ": `membership` column count does not match `ntree`.",
        call. = FALSE
      )
    }
  }
  if (!is.null(mod$xvar) && nrow(mod$xvar) != nrow(membership)) {
    stop(
      "Cannot reconstruct proximity", model_label,
      ": `membership` and `xvar` have different sample counts.",
      call. = FALSE
    )
  }

  sample_names <- rownames(membership)
  if (is.null(sample_names) && !is.null(mod$xvar)) {
    sample_names <- rownames(mod$xvar)
  }
  if (is.null(sample_names)) {
    sample_names <- as.character(seq_len(nrow(membership)))
  }
  if (length(sample_names) != nrow(membership) || anyNA(sample_names) ||
      any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
    stop(
      "Cannot reconstruct proximity", model_label,
      ": sample names must be complete and unique.",
      call. = FALSE
    )
  }
  xvar_names <- if (is.null(mod$xvar)) NULL else rownames(mod$xvar)
  if (!is.null(xvar_names) && !identical(sample_names, xvar_names)) {
    stop(
      "Cannot reconstruct proximity", model_label,
      ": `membership` and `xvar` sample names are not aligned.",
      call. = FALSE
    )
  }

  force(membership)
  force(sample_names)
  n <- nrow(membership)
  tree_ids <- seq_len(ncol(membership))

  sum_tree_chunk <- function(ids) {
    mem <- membership[, ids, drop = FALSE]
    n_chunk_trees <- ncol(mem)
    keys <- paste(
      rep(seq_len(n_chunk_trees), each = n),
      as.vector(mem),
      sep = "\r"
    )
    leaf_group <- match(keys, unique(keys))
    indicator <- Matrix::sparseMatrix(
      i = rep.int(seq_len(n), times = n_chunk_trees),
      j = leaf_group,
      x = 1,
      dims = c(n, max(leaf_group))
    )
    Matrix::tcrossprod(indicator)
  }

  if (isTRUE(parallel) && length(tree_ids) > 1L && !is.null(cores)) {
    worker_count <- cores
    worker_count <- min(
      sanitize_mc_cores(worker_count, fallback = 1L),
      length(tree_ids)
    )
    chunk_id <- ceiling(seq_along(tree_ids) * worker_count / length(tree_ids))
    tree_chunks <- split(tree_ids, chunk_id)
    chunk_sums <- parallel_lapply_psock(
      tree_chunks,
      sum_tree_chunk,
      cores = worker_count
    )
    prox <- Reduce("+", chunk_sums)
  } else {
    prox <- sum_tree_chunk(tree_ids)
  }

  prox <- prox / length(tree_ids)
  dimnames(prox) <- list(sample_names, sample_names)
  prox
}

plain_proximity_for_model <- function(mod,
                                      parallel = TRUE,
                                      cores = NULL,
                                      model_index = NULL) {
  proximity_mode <- if (is.null(mod$proximity.mode)) {
    NA_character_
  } else {
    as.character(mod$proximity.mode)[1L]
  }
  prox <- mod$proximity

  if (identical(proximity_mode, "none") || is.null(prox) || identical(prox, FALSE)) {
    return(reconstruct_plain_proximity(
      mod,
      parallel = parallel,
      cores = cores,
      model_index = model_index
    ))
  }

  model_label <- if (is.null(model_index)) "" else paste0(" for model ", model_index)
  if (is.null(prox) || length(dim(prox)) != 2L || nrow(prox) != ncol(prox) ||
      nrow(prox) < 1L) {
    stop(
      "Invalid proximity", model_label, ": expected a non-empty square matrix.",
      call. = FALSE
    )
  }
  if (any(!is.finite(prox))) {
    stop(
      "Invalid proximity", model_label, ": matrix contains non-finite values.",
      call. = FALSE
    )
  }
  if (any(prox < 0)) {
    stop(
      "Invalid proximity", model_label, ": matrix contains negative values.",
      call. = FALSE
    )
  }
  if (!isTRUE(isSymmetric(prox, tol = 1e-12))) {
    stop(
      "Invalid proximity", model_label, ": matrix must be symmetric.",
      call. = FALSE
    )
  }
  d_prox <- if (inherits(prox, "Matrix")) Matrix::diag(prox) else diag(prox)
  if (any(!is.finite(d_prox)) || any(d_prox <= 0)) {
    stop(
      "Invalid proximity", model_label,
      ": diagonal must be positive and finite.",
      call. = FALSE
    )
  }

  expected_names <- if (!is.null(mod$xvar)) rownames(mod$xvar) else NULL
  if (is.null(expected_names) && !is.null(mod$membership)) {
    expected_names <- rownames(mod$membership)
  }
  if (!is.null(expected_names)) {
    if (!identical(rownames(prox), expected_names) ||
        !identical(colnames(prox), expected_names)) {
      stop(
        "Invalid proximity", model_label,
        ": matrix sample names are not aligned with the fitted model.",
        call. = FALSE
      )
    }
  }
  prox
}

combine_proximity_matrices <- function(prox) {
  if (!is.list(prox) || length(prox) < 1L) {
    stop("No proximity matrices were supplied.", call. = FALSE)
  }

  first_dim <- dim(prox[[1L]])
  if (length(first_dim) != 2L || first_dim[1L] < 1L ||
      first_dim[1L] != first_dim[2L]) {
    stop("Proximity matrix for model 1 must be a non-empty square matrix.", call. = FALSE)
  }
  ref_n <- nrow(prox[[1L]])
  ref_names <- rownames(prox[[1L]])
  for (i in seq_along(prox)) {
    current <- prox[[i]]
    if (length(dim(current)) != 2L || nrow(current) != ref_n ||
        ncol(current) != ref_n) {
      stop(
        "Proximity matrix for model ", i,
        " does not match the reference dimensions.",
        call. = FALSE
      )
    }

    current_rows <- rownames(current)
    current_cols <- colnames(current)
    if (is.null(ref_names)) {
      if (!is.null(current_rows) || !is.null(current_cols)) {
        stop(
          "Proximity matrix for model ", i,
          " has sample names but the reference matrix does not.",
          call. = FALSE
        )
      }
    } else {
      if (is.null(current_rows) || is.null(current_cols) ||
          anyDuplicated(current_rows) || anyDuplicated(current_cols) ||
          !setequal(current_rows, ref_names) || !setequal(current_cols, ref_names)) {
        stop(
          "Proximity matrix for model ", i,
          " has sample names that do not match the reference matrix.",
          call. = FALSE
        )
      }
      if (!identical(current_rows, ref_names) ||
          !identical(current_cols, ref_names)) {
        current <- current[ref_names, ref_names, drop = FALSE]
        prox[[i]] <- current
      }
    }
  }

  Reduce("+", prox)
}

#' MRF unsupervised clustering -- Proximity method
#'
#' @param rfit A model list of random forest models.
#' @param k Pre-defined number of clusters. By default, the optimal value is
#'   selected by the requested tuning method.
#' @param enhanced Logical; whether to calculate enhanced proximity.
#' @param size_min Minimum terminal-node size when computing enhanced proximity.
#' @param use Which importance side to use in tree traversal, `"X"` or `"Y"`.
#' @param symm Logical; whether to symmetrize directional proximities.
#' @param leaf_embed_dim Dimension used for low-dimensional leaf embedding.
#' @param merge_quantile Quantile threshold for sibling-pair merges in enhanced
#'   mode. `0.9` keeps the top 10\% highest sibling correlations at each
#'   iteration. Used only when `merge_mode = "hard"`.
#' @param merge_mode Enhanced-proximity merge mode: `"soft"` (default) or
#'   `"hard"`. `"hard"` performs iterative sibling merges; `"soft"` directly
#'   builds per-tree similarity with same-leaf = 1 and sibling-leaf =
#'   `sibling_gamma * f(corr)`. `"soft"` is generally more stable.
#' @param sibling_gamma Multiplicative weight `gamma` used in soft merge mode.
#' @param sibling_fun Correlation transform `f(corr)` used in soft merge mode.
#'   One of `"constant"`, `"positive"`, `"shift01"`, `"abs"`, or `"identity"`.
#' @param sibling_cap Logical; whether to cap soft sibling similarities at 1.
#' @param hard_prox_mode Proximity construction used after hard merge:
#'   `"soft_enhanced"` uses same-enhanced-leaf = 1 and sibling-enhanced-leaf =
#'   `sibling_gamma * f(corr)`, while `"binary"` uses same-enhanced-leaf only.
#' @param parallel Logical; whether to allow forest-level computation in fresh
#'   PSOCK worker processes. Parallel execution also requires an explicit
#'   `cores` value.
#' @param sparse Logical; whether to sparsify the enhanced proximity matrix.
#' @param method_cl Clustering backend (`"PAM"` or `"Spectral"`).
#' @param tune_method Tuning criterion for PAM when `k` is `NULL`
#'   (`"silhouette"` or `"ratio"`).
#' @param gap_w Weighting scheme for spectral eigengap when `k` is `NULL`
#'   (`"uniform"` or `"log"`).
#' @param cores Number of CPU cores used by parallel steps. Default `NULL`
#'   keeps matrix-heavy forest reductions serial to limit peak memory.
#' @param ... Additional arguments passed to downstream clustering helpers.
#'
#' @return An `mrf3` clustering object.
#' @export
mrf3_cl_prox <- function(rfit, k = NULL,
                         enhanced = TRUE,
                         size_min = 5, use = "X", symm = TRUE,
                         leaf_embed_dim = 10,
                         merge_quantile = 0.9,
                         merge_mode = c("soft", "hard"),
                         sibling_gamma = 0.5,
                         sibling_fun = c("constant", "positive", "shift01", "abs", "identity"),
                         sibling_cap = TRUE,
                         hard_prox_mode = c("soft_enhanced", "binary"),
                         parallel = TRUE,
                         sparse = FALSE,
                         method_cl = "PAM",
                         tune_method = "silhouette",
                         gap_w = "uniform",
                         cores = NULL,
                         ...){
  dot_args <- list(...)
  if ("calc_imp" %in% names(dot_args)) {
    warning("`calc_imp` has been removed from proximity clustering and is ignored.", call. = FALSE)
  }
  if ("thres" %in% names(dot_args)) {
    warning("`thres` is no longer used in proximity clustering and is ignored.", call. = FALSE)
  }
  merge_mode <- match.arg(merge_mode)
  sibling_fun <- match.arg(sibling_fun)
  hard_prox_mode <- match.arg(hard_prox_mode)
  method_cl <- match.arg(method_cl, c("PAM", "Spectral"))
  if (isTRUE(enhanced) && identical(merge_mode, "hard")) {
    warning(
      "`merge_mode = 'hard'` is experimental and can over-merge leaves. ",
      "Use `merge_mode = 'soft'` for more stable enhanced proximity.",
      call. = FALSE
    )
  }
  if (!is.numeric(leaf_embed_dim) || length(leaf_embed_dim) != 1L ||
      !is.finite(leaf_embed_dim) || leaf_embed_dim < 1) {
    stop("`leaf_embed_dim` must be a single positive numeric value.")
  }
  if (!is.numeric(merge_quantile) || length(merge_quantile) != 1L ||
      !is.finite(merge_quantile) || merge_quantile < 0 || merge_quantile > 1) {
    stop("`merge_quantile` must be a single numeric value in [0, 1].")
  }
  if (!is.numeric(sibling_gamma) || length(sibling_gamma) != 1L ||
      !is.finite(sibling_gamma) || sibling_gamma < 0) {
    stop("`sibling_gamma` must be a single non-negative numeric value.")
  }

  if(length(rfit) == 1){
    if(is.null(rfit[[1]]$yvar)) symm <- FALSE
  }
  if(enhanced){
    ## For sub-MRF models that already carry a pre-computed enhanced_prox
    ## matrix (computed inside fit_sub_mrf with full-data embeddings),
    ## use it directly instead of calling cl_forest() which requires
    ## rfsrc internals ($forest$nativeArray) that sub-models don't retain.
    has_precomputed <- vapply(rfit, function(r) !is.null(r$enhanced_prox), logical(1))
    is_sub_mrf <- vapply(rfit, function(r) !is.null(r$sub_mrf_info), logical(1))

    if (all(has_precomputed)) {
      message("[mrf3_cl_prox] Using pre-computed C++ enhanced proximity (fast path)")
      cl_mod <- NULL
      prox <- purrr::map(rfit, "enhanced_prox")
      prox <- combine_proximity_matrices(prox)
    } else if (any(is_sub_mrf) && !all(has_precomputed)) {
      ## Sub-MRF models without pre-computed enhanced_prox: fall back to
      ## plain proximity (enhanced prox was not requested at fit time).
      warning(
        "Sub-MRF models lack pre-computed enhanced_prox. ",
        "Set enhanced = TRUE in sub_mrf_args to enable. Falling back to plain proximity.",
        call. = FALSE
      )
      cl_mod <- NULL
      prox <- lapply(seq_along(rfit), function(i) {
        plain_proximity_for_model(
          rfit[[i]],
          parallel = parallel,
          cores = cores,
          model_index = i
        )
      })
      prox <- combine_proximity_matrices(prox)
    } else {
      message("[mrf3_cl_prox] enhanced_prox not pre-computed (",
              sum(has_precomputed), "/", length(has_precomputed),
              " models have it). Falling back to R-level cl_forest().")
      cl_mod <- plyr::llply(
        rfit,
        .fun = function(r){
          cl_forest(r,
                    size_min = size_min,
                    parallel = parallel,
                    use = use,
                    symm = symm,
                    leaf_embed_dim = leaf_embed_dim,
                    merge_quantile = merge_quantile,
                    merge_mode = merge_mode,
                    sibling_gamma = sibling_gamma,
                    sibling_fun = sibling_fun,
                    sibling_cap = sibling_cap,
                    hard_prox_mode = hard_prox_mode,
                    cores = cores,
                    ...)
        }
      )

      prox <- purrr::map(cl_mod, "prox")
      prox <- combine_proximity_matrices(prox)
    }
  } else {

    cl_mod <- NULL
    proximity_mode <- vapply(
      rfit,
      function(r) {
        if (is.null(r$proximity.mode)) NA_character_ else as.character(r$proximity.mode)[1L]
      },
      character(1)
    )
    proximity_missing <- vapply(
      rfit,
      function(r) is.null(r$proximity) || identical(r$proximity, FALSE),
      logical(1)
    )
    needs_reconstruction <- !isTRUE(enhanced) & (
      (!is.na(proximity_mode) & proximity_mode == "none") | proximity_missing
    )
    if (any(needs_reconstruction)) {
      message(
        "[mrf3_cl_prox] Reconstructing ordinary proximity from tree membership for ",
        sum(needs_reconstruction), "/", length(rfit),
        " models without stored proximity."
      )
    }
    prox <- lapply(seq_along(rfit), function(i) {
      plain_proximity_for_model(
        rfit[[i]],
        parallel = parallel,
        cores = cores,
        model_index = i
      )
    })

    prox <- combine_proximity_matrices(prox)

  }

  if (inherits(prox, "Matrix")) {
    prox <- as.matrix(prox)
  }

  d_prox <- diag(prox)
  if (any(!is.finite(d_prox)) || any(d_prox <= 0)) {
    stop(
      "Proximity matrix diagonal must be positive and finite for normalization.",
      call. = FALSE
    )
  }
  prox <- prox / sqrt(d_prox %o% d_prox)
  
  # Make sparse proximity
  if(enhanced & sparse){
    diag(prox) <- 0
    prox[prox < estimate_density_mode(prox)] <- 0
    diag(prox) <- 1
  }
  
  
  nm <- rownames(rfit[[1]]$xvar)
  if (!is.null(nm)) rownames(prox) <- colnames(prox) <- nm

  if(method_cl == "PAM") {
    # Keep similarity here and let PAM consume dissimilarity (1 - prox) explicitly.
    p <- prox
  }
  if(method_cl == "Spectral"){
    p <- prox
    diag(p) <- 0
  }
  k_tuned <- is.null(k)
  if(is.null(k)){
    message("Start tuning k step..")

    if(method_cl == "PAM"){
      k_fit <- tune_k_clusters(p, return_cluster = TRUE, method = method_cl,
                               tune_method = tune_method, gap_w = gap_w, prox = TRUE)
    } else {
      k_fit <- tune_k_clusters(p, return_cluster = TRUE, method = method_cl,
                               tune_method = tune_method, gap_w = gap_w)
    }
    cl <- k_fit$cl
    k_selected <- as.integer(k_fit$best_k)[1]

  } else {
    if(method_cl == "PAM"){
      cl_fit <- pam_cl(1 - p, k_tune = k, diss = TRUE)
    } 
    if(method_cl == "Spectral"){
      cl_fit <- spectral_cl(p, k_tune = k)
    }
    cl <- cl_fit$cl
    k_selected <- as.integer(k)[1]
  }

  message("Done!")
  out <- list(dat = prox,
              cl = cl,
              k = k_selected,
              k_tuned = k_tuned,
              cl_mod = cl_mod,
              enhanced = enhanced,
              method = "Proximity")
  class(out) <- "prox"

  return(out)
}

estimate_density_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
