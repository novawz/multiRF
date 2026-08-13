.combine_specific_imd_per_tree <- function(matrices, feature_names, block) {
  matrices <- Filter(Negate(is.null), matrices)
  if (!length(matrices)) return(NULL)
  aligned <- lapply(seq_along(matrices), function(i) {
    mat <- as.matrix(matrices[[i]])
    rn <- rownames(mat)
    if (is.null(rn) || anyNA(rn) || any(!nzchar(rn)) || anyDuplicated(rn)) {
      stop("Per-tree specific IMD for block `", block, "` in consensus run ",
           i, " must have unique, non-empty feature row names.")
    }
    missing_features <- setdiff(feature_names, rn)
    extra_features <- setdiff(rn, feature_names)
    if (length(missing_features) || length(extra_features)) {
      stop("Per-tree specific IMD feature mismatch for block `", block,
           "` in consensus run ", i, ".")
    }
    mat <- mat[feature_names, , drop = FALSE]
    if (any(!is.finite(mat))) {
      stop("Per-tree specific IMD for block `", block,
           "` contains non-finite values in consensus run ", i, ".")
    }
    mat
  })
  do.call(cbind, aligned)
}

#' Build Shared/Specific Weights From Reconstruction
#'
#' @param dat.list A named list of omics matrices (samples in rows, features in columns).
#' @param recon Reconstruction output from `get_reconstr_matrix()`.
#' @param per_response_recon Logical; when \code{TRUE} (default), the shared
#' reconstruction for each block \code{k} uses only the forest weight matrices
#' from connections where block \code{k} is the \emph{response}.
#' Specifically, `W^(k) = sum_m alpha_m W_m` (modularity-weighted) over
#' models \code{m} whose response is block \code{k}, and then
#' `X_hat^(k) = W^(k) X^(k)`.
#' This ensures that each block's residual captures genuinely block-specific
#' variation rather than reconstruction artefacts from connections where
#' the block served as predictor.
#' When \code{FALSE}, reverts to the legacy behaviour that uses the global
#' fused reconstruction (\code{recon$fused_mat}) or \code{W_{all} \%*\% X}.
#' @param specific_top_v Optional integer. If set, each row of fused specific
#' residual-forest weights keeps top-v entries.
#' @param specific_keep_ties Logical; whether specific top-v truncation keeps ties at cutoff.
#' @param specific_row_normalize Logical; whether to row-normalize fused specific
#' weights after optional truncation.
#' @param specific_seed Integer seed for the residual unsupervised forests.
#' When \code{NULL} (default), falls back to 529 with a warning.
#' \code{mrf3_fit()} automatically injects the pipeline seed here, so
#' users calling \code{mrf3()} or \code{mrf3_fit()} need not set this.
#' Direct callers of \code{get_shared_specific_weights()} should pass
#' an explicit positive seed for reproducibility.
#' @param specific_n_consensus Integer; number of consensus runs for the
#' residual unsupervised forests. When \code{> 1}, fits the residual RF
#' \code{specific_n_consensus} times with seeds
#' \code{specific_seed, specific_seed + 1, ...} and averages the forest
#' weights and IMD across runs. This reduces variance in the specific
#' signal at the cost of proportionally more computation. Default \code{1L}
#' (single run, no consensus averaging).
#' @param specific_ntree,specific_samptype,specific_ytry,specific_proximity RF
#' structure inherited from the main forest workflow for residual
#' unsupervised forests.
#' @param specific_nsplit,specific_nthread Candidate cutpoint count and native
#' thread count inherited from the main workflow.
#' @param specific_nodesize,specific_max_depth Residual-tree stopping settings.
#' @param specific_forest_wt Forest-weight mode for residual forests. The
#' bioRxiv clustering definition uses `"all"`.
#' @param ... Deprecated/ignored arguments kept for compatibility.
#'
#' @return A list with `shared` and `specific` components:
#' - `shared$W_all`: shared fused weights (global, for shared clustering).
#' - `shared$W_per_response`: named list of per-response fused weight matrices
#'   (only present when `per_response_recon = TRUE`).
#' - `specific$residual`: residual omics matrices
#'   `R = X - X_pred`.
#' - `specific$predicted`: predicted omics matrices
#'   `X_pred` from per-response (or global) reconstruction.
#' - `specific$residual_mod`: unsupervised RF models fitted on residual matrices.
#' - `specific$W`: specific residual weights from residual RF models
#'   after adjusted/truncate/row-sum-normalize.
#' - `specific$imd`: named list of per-block variable-level IMD weights
#'   (named numeric vectors) from unsupervised RF on residuals.
get_shared_specific_weights <- function(dat.list,
                                        recon,
                                        per_response_recon = TRUE,
                                        specific_top_v = NULL,
                                        specific_keep_ties = TRUE,
                                        specific_row_normalize = TRUE,
                                        specific_seed = NULL,
                                        specific_n_consensus = 1L,
                                        specific_ntree = 500L,
                                        specific_samptype = c("swor", "swr"),
                                        specific_ytry = NULL,
                                        specific_proximity = "all",
                                        specific_nsplit = 10L,
                                        specific_nthread = getOption("multiRF.nthread", 0L),
                                        specific_nodesize = 3L,
                                        specific_max_depth = 0L,
                                        specific_forest_wt = "all",
                                        ...) {

  dot_args <- list(...)
  deprecated_args <- intersect(
    names(dot_args),
    c("mod_list", "connection_score", "specific_model_top_v", "score_power",
      "score_floor", "fallback_uniform", "specific_consensus_seed")
  )
  if (length(deprecated_args) > 0L) {
    warning(
      "Ignoring deprecated arguments in `get_shared_specific_weights()`: ",
      paste(deprecated_args, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.list(dat.list) || length(dat.list) == 0) {
    stop("`dat.list` must be a non-empty named list of omics matrices.")
  }
  if (is.null(names(dat.list)) || any(names(dat.list) == "")) {
    stop("`dat.list` must be named.")
  }
  if (!is.list(recon) || is.null(recon$W) || is.null(recon$W$W_all)) {
    stop("`recon` must be a reconstruction object returned by `get_reconstr_matrix()`.")
  }

  dat_names <- names(dat.list)
  specific_samptype <- match.arg(specific_samptype)
  specific_proximity <- match.arg(
    specific_proximity, c("all", "inbag", "oob", "none")
  )
  specific_forest_wt <- match.arg(
    specific_forest_wt, c("all", "inbag", "oob")
  )
  n_cons <- max(as.integer(specific_n_consensus), 1L)
  if (is.null(specific_seed)) {
    specific_seed <- 529L
    warning(
      "`specific_seed` not set; defaulting to 529. ",
      "Pass `specific_seed` explicitly for reproducible specific weights.",
      call. = FALSE
    )
  }
  spec_seed <- as.integer(specific_seed)

  predicted <- vector("list", length(dat_names))
  residual <- vector("list", length(dat_names))
  names(predicted) <- names(residual) <- dat_names

  specific_W <- vector("list", length(dat_names))
  specific_imd <- vector("list", length(dat_names))
  specific_imd_per_tree <- vector("list", length(dat_names))
  residual_mod <- vector("list", length(dat_names))
  names(specific_W) <- names(specific_imd) <- names(specific_imd_per_tree) <- names(residual_mod) <- dat_names

  # --- Per-response weight matrices W^(k) -----------------------------------
  # Pre-computed in get_reconstr_matrix() and stored in recon$W$W_per_response.
  # Each W^(k) = sum_{i != k} alpha_{ki} W_{ki}, normalised within the
  # response-k subset and row-normalised, so X_hat = W^(k) X^(k) is a proper
  # linear smoother using only forests trained to predict block k.
  W_per_response <- if (per_response_recon) recon$W$W_per_response else NULL

  for (d in dat_names) {
    X <- as.matrix(dat.list[[d]])
    X_hat <- NULL

    if (per_response_recon && !is.null(W_per_response[[d]])) {
      # Per-response reconstruction: X_hat = W^(d) %*% X
      W_d <- as.matrix(W_per_response[[d]])
      if (nrow(W_d) != nrow(X) || ncol(W_d) != nrow(X)) {
        stop("Dimension mismatch for per-response W of block `", d,
             "`: expected n x n with n = ", nrow(X), ".")
      }
      X_hat <- W_d %*% X
    } else if (!is.null(recon$fused_mat) && !is.null(recon$fused_mat[[d]])) {
      # Legacy: use global fused reconstruction
      X_hat <- as.matrix(recon$fused_mat[[d]])

      if (!is.null(rownames(X)) && !is.null(rownames(X_hat)) &&
          setequal(rownames(X), rownames(X_hat))) {
        X_hat <- X_hat[rownames(X), , drop = FALSE]
      }
      if (!is.null(colnames(X)) && !is.null(colnames(X_hat)) &&
          setequal(colnames(X), colnames(X_hat))) {
        X_hat <- X_hat[, colnames(X), drop = FALSE]
      }

      if (!all(dim(X_hat) == dim(X))) {
        warning(
          "Dimension mismatch for fused reconstruction in block `", d,
          "`. Falling back to weight-based prediction `W_shared %*% X`.",
          call. = FALSE
        )
        X_hat <- NULL
      }
    }

    if (is.null(X_hat)) {
      # Final fallback: use shared fused weight matrix to predict X
      W_shared <- as.matrix(recon$W$W_all)
      if (nrow(W_shared) != nrow(X) || ncol(W_shared) != nrow(X)) {
        stop("Dimension mismatch for `", d, "`: shared weight must be n x n with n = nrow(dat.list[[d]]).")
      }
      X_hat <- W_shared %*% X
    }

    R <- X - X_hat
    rownames(X_hat) <- rownames(R) <- rownames(X)
    colnames(X_hat) <- colnames(R) <- colnames(X)
    predicted[[d]] <- X_hat
    residual[[d]] <- R

    # Inherit the main forest's structural settings; enhanced proximity
    # remains an optional extension forwarded through `...`.
    fit_extra <- c(
      list(
        ntree = as.integer(specific_ntree),
        samptype = specific_samptype,
        ytry = specific_ytry,
        proximity = specific_proximity,
        nsplit = as.integer(specific_nsplit),
        nthread = as.integer(specific_nthread),
        nodesize = as.integer(specific_nodesize),
        max_depth = as.integer(specific_max_depth),
        forest.wt = specific_forest_wt
      ),
      dot_args[intersect(
        names(dot_args), c("enhanced_prox", "sibling_gamma", "leaf_embed_dim")
      )]
    )

    # ---- Fit residual unsupervised RF (with seed and optional consensus) ----
    if (n_cons == 1L) {
      # Single run: pass seed directly to fit_forest for determinism
      r_mod <- do.call(fit_forest, c(
        list(
          X = as.data.frame(R, check.names = FALSE),
          Y = NULL,
          type = "unsupervised",
          seed = spec_seed
        ),
        fit_extra
      ))
      if (is.null(r_mod$forest.wt)) {
        stop("Residual unsupervised model for `", d, "` does not contain `forest.wt`.")
      }
      residual_mod[[d]] <- r_mod
      specific_W[[d]] <- prepare_weight_matrix(
        W = r_mod$forest.wt,
        adjust = TRUE,
        top_v = specific_top_v,
        row_normalize = specific_row_normalize,
        zero_diag = TRUE,
        keep_ties = specific_keep_ties
      )
      if (!is.null(r_mod$imd_weights) && !is.null(r_mod$imd_weights$X)) {
        specific_imd[[d]] <- r_mod$imd_weights$X
      } else {
        specific_imd[[d]] <- setNames(rep(1.0 / ncol(X), ncol(X)), colnames(X))
      }
      if (!is.null(r_mod$imd_weights_per_tree) && !is.null(r_mod$imd_weights_per_tree$X)) {
        specific_imd_per_tree[[d]] <- r_mod$imd_weights_per_tree$X
      }

    } else {
      # Consensus path: average forest weights and IMD across n_cons runs
      wt_accum <- NULL
      imd_accum <- NULL
      imd_per_tree_runs <- vector("list", n_cons)
      last_mod <- NULL

      for (ci in seq_len(n_cons)) {
        # Use sequential seeds: spec_seed, spec_seed+1, ...
        ci_seed <- spec_seed + ci - 1L
        r_mod_ci <- do.call(fit_forest, c(
          list(
            X = as.data.frame(R, check.names = FALSE),
            Y = NULL,
            type = "unsupervised",
            seed = ci_seed
          ),
          fit_extra
        ))
        if (is.null(r_mod_ci$forest.wt)) {
          stop("Residual unsupervised model (consensus run ", ci,
               ") for `", d, "` does not contain `forest.wt`.")
        }

        wt_i <- as.matrix(r_mod_ci$forest.wt)
        if (is.null(wt_accum)) {
          wt_accum <- wt_i
        } else {
          wt_accum <- wt_accum + wt_i
        }

        if (!is.null(r_mod_ci$imd_weights) && !is.null(r_mod_ci$imd_weights$X)) {
          imd_i <- as.numeric(r_mod_ci$imd_weights$X)
          names(imd_i) <- names(r_mod_ci$imd_weights$X)
          imd_names <- names(imd_i)
          if (is.null(imd_names) || anyNA(imd_names) || any(!nzchar(imd_names)) ||
              anyDuplicated(imd_names) ||
              !setequal(imd_names, colnames(X))) {
            stop("Specific IMD feature mismatch for block `", d,
                 "` in consensus run ", ci, ".")
          }
          imd_i <- imd_i[colnames(X)]
          if (any(!is.finite(imd_i))) {
            stop("Specific IMD for block `", d,
                 "` contains non-finite values in consensus run ", ci, ".")
          }
          if (is.null(imd_accum)) {
            imd_accum <- imd_i
          } else {
            imd_accum <- imd_accum + imd_i
          }
        }
        if (!is.null(r_mod_ci$imd_weights_per_tree) &&
            !is.null(r_mod_ci$imd_weights_per_tree$X)) {
          imd_per_tree_runs[[ci]] <- r_mod_ci$imd_weights_per_tree$X
        }
        last_mod <- r_mod_ci
      }

      # Average
      wt_avg <- wt_accum / n_cons
      residual_mod[[d]] <- last_mod
      specific_W[[d]] <- prepare_weight_matrix(
        W = wt_avg,
        adjust = TRUE,
        top_v = specific_top_v,
        row_normalize = specific_row_normalize,
        zero_diag = TRUE,
        keep_ties = specific_keep_ties
      )
      if (!is.null(imd_accum)) {
        imd_avg <- imd_accum / n_cons
        specific_imd[[d]] <- imd_avg
      } else {
        specific_imd[[d]] <- setNames(rep(1.0 / ncol(X), ncol(X)), colnames(X))
      }
      # Preserve all consensus runs. With the same tree count per run, the row
      # means of this combined matrix equal the consensus mean IMD rather than
      # silently reverting to the final run.
      specific_imd_per_tree[[d]] <- .combine_specific_imd_per_tree(
        imd_per_tree_runs,
        feature_names = colnames(X),
        block = d
      )
    }
  }

  shared_out <- list(
    source = if (per_response_recon) "per_response" else "W_all",
    W_all = recon$W$W_all
  )
  if (!is.null(W_per_response)) {
    shared_out$W_per_response <- W_per_response
  }

  list(
    shared = shared_out,
    specific = list(
      residual = residual,
      predicted = predicted,
      residual_mod = residual_mod,
      W = specific_W,
      imd = specific_imd,
      imd_per_tree = specific_imd_per_tree
    )
  )
}
