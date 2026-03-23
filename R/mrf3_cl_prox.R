#' MRF unsupervised clustering -- Proximity method
#' @param rfit A model list of random forest models
#' @param k Pre-defined number of clusters. The default is selecting the optimal k by tuning method
#' @param enhanced A logical parameter that determines whether enhanced proximity is calculated when selecting Proximity as clustering method.
#' The default is True.
#' @param size_min Minimum terminal-node size when computing enhanced proximity.
#' @param use Which importance side to use in tree traversal, `"X"` or `"Y"`.
#' @param symm Logical; whether to symmetrize directional proximities.
#' @param leaf_embed_dim Dimension used for low-dimensional leaf embedding.
#' @param merge_quantile Quantile threshold for sibling-pair merges in enhanced mode.
#' `0.9` means keep top 10\% highest sibling correlations each iteration.
#' Used only when `merge_mode = "hard"`.
#' @param merge_mode Enhanced-proximity merge mode: `"soft"` (default) or `"hard"`.
#' `"hard"` performs iterative sibling merges; `"soft"` directly builds per-tree
#' similarity with same-leaf = 1 and sibling-leaf = `sibling_gamma * f(corr)`.
#' In practice, `"soft"` is generally more stable.
#' @param sibling_gamma Multiplicative weight `gamma` used in soft merge mode.
#' @param sibling_fun Correlation transform `f(corr)` used in soft merge mode.
#' One of `"constant"`, `"positive"`, `"shift01"`, `"abs"`, `"identity"`.
#' @param sibling_cap Logical; whether to cap soft sibling similarities at 1.
#' @param hard_prox_mode Proximity construction used after hard merge:
#' `"soft_enhanced"` uses same-enhanced-leaf = 1 and sibling-enhanced-leaf
#' = `sibling_gamma * f(corr)`, while `"binary"` uses same-enhanced-leaf only.
#' @param parallel Logical; whether to parallelize forest-level computation.
#' @param sparse Logical; whether to sparsify the enhanced proximity matrix.
#' @param method_cl Clustering backend (`"PAM"` or `"Spectral"`).
#' @param cores Number of CPU cores used by parallel steps.
#' @param ... Additional arguments passed to downstream clustering helpers.
#'
#' @return mrf3 clustering object
#' @export
#'
# -------------------------------------------------------------------------------------------------------------
# Proximity clustering method
# -------------------------------------------------------------------------------------------------------------

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
      prox <- Reduce("+", prox)
    } else if (any(is_sub_mrf) && !all(has_precomputed)) {
      ## Sub-MRF models without pre-computed enhanced_prox: fall back to
      ## plain proximity (enhanced prox was not requested at fit time).
      warning(
        "Sub-MRF models lack pre-computed enhanced_prox. ",
        "Set enhanced = TRUE in sub_mrf_args to enable. Falling back to plain proximity.",
        call. = FALSE
      )
      cl_mod <- NULL
      prox <- purrr::map(rfit, "proximity")
      prox <- Reduce("+", prox)
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
      prox <- Reduce("+", prox)
    }
  } else {

    cl_mod <- NULL
    prox <- purrr::map(rfit, "proximity")

    prox <- Reduce("+", prox)

  }

  if (inherits(prox, "sparseMatrix")) {
    prox <- as.matrix(prox)
  }

  num_dim <- diag(prox)[1]
  prox <- prox/num_dim
  
  # Make sparse proximity
  if(enhanced & sparse){
    diag(prox) <- 0
    prox[prox < estimate_density_mode(prox)] <- 0
    diag(prox) <- 1
  }
  
  
  rownames(prox) <- colnames(prox) <- rownames(rfit[[1]]$xvar)

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
      k_fit <- tune_k_clusters(p, return_cluster = TRUE, method = method_cl, prox = TRUE)
    } else {
      k_fit <- tune_k_clusters(p, return_cluster = TRUE, method = method_cl)
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

get_prox <-  function(class_mem){

  if(length(unique(class_mem)) == 1){

    class_new <- matrix(0, nrow = length(class_mem), ncol = length(class_mem))

  } else {

    class_new <- data.frame(class = as.factor(class_mem))
    one_hot <-  model.matrix(~class + 0, class_new)
    class_new <- one_hot %*% t(one_hot)

  }

  return(class_new)
}

estimate_density_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
