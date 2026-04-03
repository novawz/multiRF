#' Build Shared/Specific Weights From Reconstruction
#'
#' @param dat.list A named list of omics matrices (samples in rows, features in columns).
#' @param recon Reconstruction output from `get_reconstr_matrix()`.
#' @param per_response_recon Logical; when \code{TRUE} (default), the shared
#' reconstruction for each block \code{k} uses only the forest weight matrices
#' from connections where block \code{k} is the \emph{response}.
#' Specifically, \code{W^{(k)} = \sum \alpha_m W_m} (modularity-weighted) over
#' models \code{m} whose response is block \code{k}, and then
#' \code{\hat{X}^{(k)} = W^{(k)} X^{(k)}}.
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
#' When positive, guarantees deterministic specific signal across runs.
#' When negative (default \code{-10}), each call uses a random device seed,
#' which causes non-reproducible specific weights.
#' The recommended practice is to pass the same seed used for the main
#' pipeline (e.g., \code{specific_seed = seed}).
#' @param specific_n_consensus Integer; number of consensus runs for the
#' residual unsupervised forests. When \code{> 1}, fits the residual RF
#' \code{specific_n_consensus} times with seeds
#' \code{specific_seed, specific_seed + 1, ...} and averages the forest
#' weights and IMD across runs. This reduces variance in the specific
#' signal at the cost of proportionally more computation. Default \code{1L}
#' (single run, no consensus averaging).
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
                                        specific_seed = -10L,
                                        specific_n_consensus = 1L,
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
  n_cons <- max(as.integer(specific_n_consensus), 1L)
  spec_seed <- as.integer(specific_seed)

  predicted <- vector("list", length(dat_names))
  residual <- vector("list", length(dat_names))
  names(predicted) <- names(residual) <- dat_names

  specific_W <- vector("list", length(dat_names))
  specific_imd <- vector("list", length(dat_names))
  residual_mod <- vector("list", length(dat_names))
  names(specific_W) <- names(specific_imd) <- names(residual_mod) <- dat_names

  # --- Build per-response weight matrices W^(k) ----------------------------
  # For each block k, W^(k) is the modularity-weighted average of W_m's
  # from connections where k is the response.  This ensures that the shared
  # reconstruction X_hat^(k) = W^(k) X^(k) only uses weights optimised for
  # predicting block k, so the residual is a clean specific signal.
  W_per_response <- NULL
  if (per_response_recon) {
    W_models <- recon$W$W_models
    model_names <- names(W_models)
    model_score <- recon$model_score
    if (is.null(model_score)) {
      model_score <- setNames(rep(1.0, length(model_names)), model_names)
    }

    W_per_response <- vector("list", length(dat_names))
    names(W_per_response) <- dat_names

    for (d in dat_names) {
      # Find models where d is the response (first element of model name)
      resp_idx <- vapply(model_names, function(m) {
        pair <- parse_model_pair(m)
        identical(pair[1], d)
      }, logical(1))

      resp_models <- model_names[resp_idx]

      if (length(resp_models) == 0L) {
        # Fallback: no connection has d as response — use global W_all
        warning(
          "No connection has block `", d, "` as response. ",
          "Falling back to global W_all for reconstruction.",
          call. = FALSE
        )
        W_per_response[[d]] <- recon$W$W_all
      } else {
        # Modularity-weighted average of response-side W_m's
        resp_scores <- model_score[resp_models]
        resp_scores[!is.finite(resp_scores)] <- 1.0
        resp_scores <- pmax(resp_scores, 0)
        alpha_sum <- sum(resp_scores)
        if (alpha_sum <= 0) {
          resp_alpha <- rep(1.0 / length(resp_models), length(resp_models))
        } else {
          resp_alpha <- resp_scores / alpha_sum
        }
        names(resp_alpha) <- resp_models

        W_per_response[[d]] <- fuse_matrix_list(
          W_models[resp_models], resp_alpha
        )
        # Row-normalize the per-response fused weight
        rs <- rowSums(W_per_response[[d]])
        rs[rs <= 0 | !is.finite(rs)] <- 1
        W_per_response[[d]] <- W_per_response[[d]] / rs
      }
    }
  }

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

    # Forward enhanced_prox / sibling_gamma / leaf_embed_dim from ... if present
    fit_extra <- dot_args[intersect(names(dot_args), c("enhanced_prox", "sibling_gamma", "leaf_embed_dim"))]

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

    } else {
      # Consensus path: average forest weights and IMD across n_cons runs
      wt_accum <- NULL
      imd_accum <- NULL
      last_mod <- NULL

      for (ci in seq_len(n_cons)) {
        # Use sequential seeds: spec_seed, spec_seed+1, ...
        # If spec_seed < 0, use ci as positive seed to ensure reproducibility
        ci_seed <- if (spec_seed > 0L) (spec_seed + ci - 1L) else ci
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
          imd_i <- r_mod_ci$imd_weights$X
          if (is.null(imd_accum)) {
            imd_accum <- imd_i
          } else {
            imd_accum <- imd_accum + imd_i
          }
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
        # Re-normalize to sum to 1
        imd_avg <- imd_avg / sum(imd_avg)
        specific_imd[[d]] <- imd_avg
      } else {
        specific_imd[[d]] <- setNames(rep(1.0 / ncol(X), ncol(X)), colnames(X))
      }
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
      imd = specific_imd
    )
  )
}
