#' Compute OOB forest weights from membership and inbag info
#'
#' Given an rfsrc model fitted with `membership = TRUE`, reconstructs
#' the out-of-bag forest weight matrix.  For each sample i, only trees
#' where i is OOB contribute to its weight row.
#'
#' @param mod  A fitted rfsrc model with `membership` and `inbag` slots.
#' @return An n x n OOB forest-weight matrix.
#' @keywords internal
compute_oob_fw <- function(mod) {
  membership <- mod$membership   # n x ntree
  inbag      <- mod$inbag        # n x ntree
  n     <- nrow(membership)
  ntree <- ncol(membership)

  fw_oob    <- matrix(0, n, n)
  oob_count <- integer(n)

  for (t in seq_len(ntree)) {
    nodes    <- membership[, t]
    oob_mask <- inbag[, t] == 0

    ## Split samples by terminal node
    node_ids <- unique(nodes)
    for (nd in node_ids) {
      members <- which(nodes == nd)
      ns <- length(members)
      oob_in_node <- members[oob_mask[members]]
      if (length(oob_in_node) > 0L) {
        fw_oob[oob_in_node, members] <- fw_oob[oob_in_node, members] + 1 / ns
      }
    }
    oob_count <- oob_count + as.integer(oob_mask)
  }

  ## Normalize: each row by its number of OOB trees
  nonzero <- oob_count > 0L
  if (any(nonzero)) {
    fw_oob[nonzero, ] <- fw_oob[nonzero, , drop = FALSE] / oob_count[nonzero]
  }

  fw_oob
}


#' Fit an ensemble of sub-sampled MRF models for a single connection
#'
#' When the number of response and/or predictor features is large, fitting a
#' single MRF with all features is computationally expensive.
#' `fit_sub_mrf()` addresses this by repeatedly sub-sampling features from both
#' the response and predictor blocks of a **single connection**, fitting a
#' smaller MRF each time, and averaging the resulting forest-weight and
#' proximity matrices.
#'
#' @param X  Predictor data frame (samples x features) — a single omics block.
#' @param Y  Response data frame (samples x features) — a single omics block.
#' @param n_sub  Integer; number of sub-MRF replicates to fit.
#' @param frac_response  Fraction of response columns to sample per replicate
#'   (default 0.1).  Convergence analysis shows symmetric low fractions
#'   (0.1--0.3 for both response and predictor) yield the best
#'   approximation of the full forest-weight matrix.
#' @param frac_predictor  Fraction of predictor columns to sample per replicate
#'   (default 0.1).  Avoid setting this to 1.0 when `frac_response` is
#'   small, as the asymmetric case leads to poor convergence (each
#'   sub-model overfits to its response subset).
#' @param ntree_per_sub  Number of trees per sub-MRF (default 50).
#' @param mtry Number of candidate X variables per split. Passed through to
#'   `fit_forest()`.
#' @param ytry  Number of candidate Y variables per split. Default `0L` means
#'   the C++ engine uses `min(qy, ceiling(p/3))`. A positive integer overrides this.
#' @param min_response  Minimum number of response columns per sub-MRF.
#' @param min_predictor  Minimum number of predictor columns per sub-MRF.
#' @param enhanced  Logical; if `TRUE`, compute soft enhanced proximity
#'   using full-data sample embeddings within each sub-model.  The result
#'   is stored as `$enhanced_prox` in the output so that
#'   `mrf3_cl_prox()` can use it directly without re-traversing trees.
#'   Default `FALSE` because the extra computation is non-trivial.
#' @param compute_imd Logical; whether to pre-compute IMD weights across sub-models.
#' @param imd_args Named list of arguments passed to `get_imp_forest()` when
#'   `compute_imd = TRUE`.
#' @param seed  Base random seed.
#' @param parallel  Logical; if `TRUE`, use [parallel::mclapply()].
#' @param cores  Number of cores when `parallel = TRUE`.
#' @param verbose  Logical; print progress messages.
#' @param ...  Additional arguments forwarded to `fit_forest()`.
#'
#' @return An rfsrc-like list with:
#' \describe{
#'   \item{forest.wt}{Averaged n x n forest-weight matrix (all samples).}
#'   \item{forest.wt.oob}{Averaged n x n OOB forest-weight matrix.}
#'   \item{proximity}{Averaged n x n proximity matrix.}
#'   \item{xvar}{Full predictor data (all columns of `X`).}
#'   \item{yvar}{Full response data (all columns of `Y`).}
#'   \item{ntree}{Total tree count across all sub-MRFs.}
#'   \item{sub_mrf_info}{List with coverage and timing metadata.}
#' }
#'
#' @details
#' Both `forest.wt` and `proximity` are n x n matrices that do NOT depend
#' on feature dimension.  Averaging across sub-models with different feature
#' subsets produces a smoothed estimate.
#'
#' `forest.wt.oob` is computed from each sub-model's `membership` and
#' `inbag` matrices: for sample i, only trees where i is out-of-bag
#' contribute to its weight row.  This provides a regularized version
#' suitable for connection selection via `find_connection()`. Because this
#' requires `inbag`, sub-MRF currently uses the `randomForestSRC` fallback
#' path even when the package default engine is native.
#'
fit_sub_mrf <- function(X, Y,
                        n_sub = 15L,
                        frac_response = 0.2,
                        frac_predictor = 0.2,
                        ntree_per_sub = 25L,
                        mtry = NULL,
                        ytry = NULL,
                        min_response = 20L,
                        min_predictor = 30L,
                        enhanced = FALSE,
                        compute_imd = FALSE,
                        imd_args = list(),
                        seed = 529L,
                        parallel = FALSE,
                        cores = NULL,
                        verbose = TRUE,
                        ...) {

  ## ---- input checks --------------------------------------------------------
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  stopifnot(nrow(X) == nrow(Y))

  n  <- nrow(X)
  pX <- ncol(X)
  pY <- ncol(Y)

  n_resp <- max(min_response, ceiling(pY * frac_response))
  n_resp <- min(n_resp, pY)

  n_pred <- max(min_predictor, ceiling(pX * frac_predictor))
  n_pred <- min(n_pred, pX)

  ## ---- coverage tracker ----------------------------------------------------
  cov_resp <- setNames(integer(pY), colnames(Y))
  cov_pred <- setNames(integer(pX), colnames(X))

  ## ---- pre-compute leaf embedding from full data (once) --------------------
  ## Only needed when enhanced proximity is requested.  Embeddings are
  ## computed from the original full X/Y so that sibling-leaf correlations
  ## are evaluated in the full-feature space, not the sub-sampled one.
  full_embed <- NULL
  if (isTRUE(enhanced)) {
    if (verbose) message("  Pre-computing full-data sample embeddings for enhanced proximity...")
    full_embed <- list(
      X = build_sample_embedding(X, leaf_embed_dim = 10L)
    )
    if (!is.null(Y) && ncol(Y) > 0L) {
      full_embed$Y <- build_sample_embedding(Y, leaf_embed_dim = 10L)
    }
  }

  ## ---- worker function -----------------------------------------------------
  fit_one <- function(b) {
    set.seed(seed + b)

    y_idx <- sort(sample.int(pY, size = n_resp))
    Y_sub <- Y[, y_idx, drop = FALSE]

    x_idx <- sort(sample.int(pX, size = n_pred))
    X_sub <- X[, x_idx, drop = FALSE]

    mod <- fit_forest(
      X = X_sub,
      Y = Y_sub,
      ntree = ntree_per_sub,
      forest.wt = "all",
      proximity = "all",
      mtry = mtry,
      ytry = ytry,
      seed = seed + b,
      ...
    )

    ## Compute OOB forest weights from membership
    fw_oob <- compute_oob_fw(mod)

    ## Compute soft enhanced proximity using full-data embeddings
    ## (only when enhanced = TRUE; otherwise skip to save time)
    enh_prox <- NULL
    if (isTRUE(enhanced) && !is.null(full_embed)) {
      enh <- tryCatch(
        cl_forest(
          mod,
          merge_mode     = "soft",
          symm           = TRUE,
          parallel       = FALSE,
          sample_embed_list = full_embed
        ),
        error = function(e) NULL
      )
      enh_prox <- if (!is.null(enh)) enh$prox else mod$proximity
    }

    res <- list(
      forest.wt     = mod$forest.wt,
      forest.wt.oob = fw_oob,
      proximity     = mod$proximity,
      enhanced_prox = enh_prox,
      y_cols = colnames(Y_sub),
      x_cols = colnames(X_sub)
    )

    ## Compute IMD on this sub-model before it is discarded
    if (isTRUE(compute_imd)) {
      imd_defaults <- list(
        mod = mod, parallel = FALSE, robust = FALSE,
        calc = "Both", normalized = FALSE, use_depth = FALSE,
        weighted = FALSE, ytry = NULL, seed = seed + b
      )
      imd_call <- utils::modifyList(imd_defaults, imd_args)
      imd_out <- tryCatch(
        do.call(get_imp_forest, imd_call),
        error = function(e) NULL
      )
      if (!is.null(imd_out)) {
        res$imd_imp <- imd_out$imp_ls
      }
    }

    res
  }

  ## ---- fit all sub-MRFs ---------------------------------------------------
  if (verbose) message(sprintf(
    "Fitting %d sub-MRFs (resp: %d/%d, pred: %d/%d, trees: %d each)...",
    n_sub, n_resp, pY, n_pred, pX, ntree_per_sub
  ))

  t0 <- proc.time()[3]

  cores <- sanitize_mc_cores(cores = cores, fallback = 1L)

  if (parallel && cores > 1L) {
    results <- parallel::mclapply(
      seq_len(n_sub),
      fit_one,
      mc.cores = cores,
      mc.set.seed = FALSE
    )
  } else {
    results <- lapply(seq_len(n_sub), function(b) {
      if (verbose && b %% 5 == 0) message(sprintf("  sub-MRF %d / %d", b, n_sub))
      fit_one(b)
    })
  }

  elapsed <- proc.time()[3] - t0
  if (verbose) message(sprintf("  Done in %.1f sec", elapsed))

  ## ---- aggregate n x n matrices -------------------------------------------
  fw_sum     <- matrix(0, n, n)
  fw_oob_sum <- matrix(0, n, n)
  prox_sum   <- matrix(0, n, n)
  enh_prox_sum <- if (isTRUE(enhanced)) matrix(0, n, n) else NULL

  ## IMD accumulators: named numeric vectors (full feature space)
  imd_X_sum   <- NULL
  imd_Y_sum   <- NULL
  imd_X_count <- NULL
  imd_Y_count <- NULL
  if (isTRUE(compute_imd)) {
    imd_X_sum   <- setNames(numeric(pX), colnames(X))
    imd_Y_sum   <- setNames(numeric(pY), colnames(Y))
    imd_X_count <- setNames(integer(pX), colnames(X))
    imd_Y_count <- setNames(integer(pY), colnames(Y))
  }

  for (b in seq_len(n_sub)) {
    r <- results[[b]]
    fw_sum     <- fw_sum     + r$forest.wt
    fw_oob_sum <- fw_oob_sum + r$forest.wt.oob
    prox_sum   <- prox_sum   + r$proximity
    if (!is.null(enh_prox_sum) && !is.null(r$enhanced_prox)) {
      enh_prox_sum <- enh_prox_sum + r$enhanced_prox
    }

    for (col in r$y_cols) cov_resp[col] <- cov_resp[col] + 1L
    for (col in r$x_cols) cov_pred[col] <- cov_pred[col] + 1L

    ## Accumulate IMD weights (each sub-model covers a feature subset)
    if (isTRUE(compute_imd) && !is.null(r$imd_imp)) {
      if (!is.null(r$imd_imp$X)) {
        x_names <- names(r$imd_imp$X)
        common_x <- intersect(x_names, names(imd_X_sum))
        imd_X_sum[common_x]   <- imd_X_sum[common_x]   + r$imd_imp$X[common_x]
        imd_X_count[common_x] <- imd_X_count[common_x] + 1L
      }
      if (!is.null(r$imd_imp$Y)) {
        y_names <- names(r$imd_imp$Y)
        common_y <- intersect(y_names, names(imd_Y_sum))
        imd_Y_sum[common_y]   <- imd_Y_sum[common_y]   + r$imd_imp$Y[common_y]
        imd_Y_count[common_y] <- imd_Y_count[common_y] + 1L
      }
    }

    results[b] <- list(NULL)
  }

  fw_avg       <- fw_sum     / n_sub
  fw_oob_avg   <- fw_oob_sum / n_sub
  prox_avg     <- prox_sum   / n_sub
  enh_prox_avg <- if (!is.null(enh_prox_sum)) enh_prox_sum / n_sub else NULL

  ## ---- build rfsrc-compatible output ---------------------------------------
  rn <- rownames(X)
  if (!is.null(rn)) {
    rownames(fw_avg) <- colnames(fw_avg) <- rn
    rownames(fw_oob_avg) <- colnames(fw_oob_avg) <- rn
    rownames(prox_avg) <- colnames(prox_avg) <- rn
    if (!is.null(enh_prox_avg)) {
      rownames(enh_prox_avg) <- colnames(enh_prox_avg) <- rn
    }
  }

  ## ---- average IMD weights -------------------------------------------------
  imd_weights <- NULL
  if (isTRUE(compute_imd) && !is.null(imd_X_sum)) {
    ## Average only over sub-models that included each feature
    safe_div <- function(s, cnt) {
      out <- ifelse(cnt > 0L, s / cnt, 0)
      out
    }
    imd_weights <- list(
      X = safe_div(imd_X_sum, imd_X_count),
      Y = safe_div(imd_Y_sum, imd_Y_count)
    )
    if (verbose) message("  Pre-computed IMD weights from sub-models.")
  }

  out <- list(
    forest.wt     = fw_avg,
    forest.wt.oob = fw_oob_avg,
    proximity     = prox_avg,
    enhanced_prox = enh_prox_avg,
    xvar          = X,
    yvar          = Y,
    ntree         = as.integer(n_sub * ntree_per_sub),
    imd_weights   = imd_weights,
    sub_mrf_info  = list(
      n_sub      = n_sub,
      ntree_per_sub = ntree_per_sub,
      frac_response  = frac_response,
      frac_predictor = frac_predictor,
      ytry = ytry,
      col_coverage_response  = cov_resp,
      col_coverage_predictor = cov_pred,
      elapsed = elapsed
    )
  )

  out
}


#' Fit sub-MRF ensemble for a multi-omics connection list
#'
#' Higher-level wrapper that matches the interface of `fit_multi_forest()`.
#' For each directed connection (response_block -> predictor_block),
#' features are sub-sampled from both blocks independently, and the
#' resulting n x n matrices are averaged across replicates.
#'
#' @param dat.list  Named list of data frames (samples x features).
#' @param connect_list  A list of connections, each a character vector of
#'   length 2: `c(response_name, predictor_name)`.  Same format as used
#'   by `fit_multi_forest()`.
#' @param ytry  Number of candidate Y variables per split
#'   (0 = `min(qy, ceiling(p/3))` default).
#' @param min_response_for_sub Minimum response-block size below which all
#'   response features are used instead of sub-sampling.
#' @param min_predictor_for_sub Minimum predictor-block size below which all
#'   predictor features are used instead of sub-sampling.
#' @param ntree_full Deprecated compatibility argument. Ignored.
#' @inheritParams fit_sub_mrf
#' @param ...  Passed to `fit_sub_mrf()` and then `fit_forest()`.
#'
#' @return A **named list** of rfsrc-compatible objects, one per connection.
#'   Names follow the `"response_predictor"` convention used by
#'   `fit_multi_forest()`.  Each object contains `forest.wt`, `forest.wt.oob`,
#'   `proximity`, `xvar`, `yvar`.
#' @details
#' Sub-MRF uses the default engine (native C++) for its internal fits.
#' OOB forest-weight reconstruction uses the `inbag` matrix returned by
#' the native engine.
fit_sub_multi_rfsrc <- function(dat.list,
                                connect_list,
                                n_sub = 15L,
                                frac_response = 0.2,
                                frac_predictor = 0.2,
                                ntree_per_sub = 25L,
                                mtry = NULL,
                                ytry = NULL,
                                min_response = 20L,
                                min_predictor = 30L,
                                min_response_for_sub = 500L,
                                min_predictor_for_sub = 500L,
                                ntree_full = 300L,
                                enhanced = FALSE,
                                compute_imd = FALSE,
                                imd_args = list(),
                                seed = 529L,
                                parallel = FALSE,
                                cores = 2L,
                                verbose = TRUE,
                                ...) {

  mod_l <- plyr::llply(
    connect_list,
    .fun = function(d) {
      resp_name <- d[1]
      pred_name <- d[2]

      Y <- dat.list[[resp_name]]
      X <- dat.list[[pred_name]]

      if (verbose) message(sprintf(
        "\n--- Connection: %s -> %s (resp: %d, pred: %d) ---",
        resp_name, pred_name, ncol(Y), ncol(X)
      ))

      ## Per-dimension decision: if a block is too small for sub-sampling,
      ## use all its features (frac = 1) instead of falling back entirely.
      ## This way the other (large) block still benefits from sub-MRF.
      frac_resp_use <- frac_response
      frac_pred_use <- frac_predictor
      resp_small <- min_response_for_sub > 0L && ncol(Y) < min_response_for_sub
      pred_small <- min_predictor_for_sub > 0L && ncol(X) < min_predictor_for_sub

      if (resp_small) frac_resp_use <- 1.0
      if (pred_small) frac_pred_use <- 1.0

      ## Guard against pathological asymmetry: frac_pred ≈ 1 with low
      ## frac_resp on a LARGE predictor block causes poor convergence.
      ## This only triggers if the user manually sets min_predictor_for_sub
      ## very high or forces frac_predictor = 1.0.
      if (frac_pred_use >= 0.9 && frac_resp_use <= 0.3 && ncol(X) >= 500L) {
        if (verbose) message(
          "  WARNING: frac_predictor ~1.0 with low frac_response on a large ",
          "predictor block (p=", ncol(X), ") leads to poor convergence. ",
          "Consider lowering frac_predictor or raising frac_response."
        )
      }

      ## If BOTH are small, full forest is simpler and equivalent
      if (resp_small && pred_small) {
        if (verbose) message(sprintf(
          "  resp (%d) and pred (%d) both small -> using full forest", ncol(Y), ncol(X)
        ))
        mod <- fit_forest(X, Y, mtry = mtry, ytry = ytry, ntree = ntree_full,
                         seed = seed, forest.wt = "all", ...)
        return(mod)
      }

      if (resp_small || pred_small) {
        reasons <- character(0)
        if (resp_small) reasons <- c(reasons, sprintf("resp %d < %d -> frac=1.0", ncol(Y), min_response_for_sub))
        if (pred_small) reasons <- c(reasons, sprintf("pred %d < %d -> frac=1.0", ncol(X), min_predictor_for_sub))
        if (verbose) message("  ", paste(reasons, collapse = "; "))
      }

      fit_sub_mrf(
        X = X, Y = Y,
        n_sub = n_sub,
        frac_response = frac_resp_use,
        frac_predictor = frac_pred_use,
        ntree_per_sub = ntree_per_sub,
        mtry = mtry,
        ytry = ytry,
        min_response = min_response,
        min_predictor = min_predictor,
        enhanced = enhanced,
        compute_imd = compute_imd,
        imd_args = imd_args,
        seed = seed,
        parallel = parallel,
        cores = cores,
        verbose = verbose,
        ...
      )
    }
  )

  names(mod_l) <- purrr::map_chr(connect_list, ~paste0(., collapse = "_"))
  mod_l
}
