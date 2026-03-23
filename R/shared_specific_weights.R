#' Build Shared/Specific Weights From Reconstruction
#'
#' @param dat.list A named list of omics matrices (samples in rows, features in columns).
#' @param recon Reconstruction output from `get_reconstr_matrix()`.
#' @param specific_top_v Optional integer. If set, each row of fused specific
#' residual-forest weights keeps top-v entries.
#' @param specific_keep_ties Logical; whether specific top-v truncation keeps ties at cutoff.
#' @param specific_row_normalize Logical; whether to row-normalize fused specific
#' weights after optional truncation.
#' @param ... Deprecated/ignored arguments kept for compatibility.
#'
#' @return A list with `shared` and `specific` components:
#' - `shared$W_all`: shared fused weights.
#' - `specific$residual`: residual omics matrices
#'   `R = X - X_pred`.
#' - `specific$predicted`: predicted omics matrices
#'   `X_pred` from shared reconstruction.
#' - `specific$residual_mod`: unsupervised RF models fitted on residual matrices.
#' - `specific$W`: specific residual weights from residual RF models
#'   after adjusted/truncate/row-sum-normalize.
get_shared_specific_weights <- function(dat.list,
                                        recon,
                                        specific_top_v = NULL,
                                        specific_keep_ties = TRUE,
                                        specific_row_normalize = TRUE,
                                        ...) {

  dot_args <- list(...)
  deprecated_args <- intersect(
    names(dot_args),
    c("mod_list", "connection_score", "specific_model_top_v", "score_power", "score_floor", "fallback_uniform")
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

  predicted <- vector("list", length(dat_names))
  residual <- vector("list", length(dat_names))
  names(predicted) <- names(residual) <- dat_names

  specific_W <- vector("list", length(dat_names))
  residual_mod <- vector("list", length(dat_names))
  names(specific_W) <- names(residual_mod) <- dat_names

  for (d in dat_names) {
    X <- as.matrix(dat.list[[d]])
    X_hat <- NULL

    if (!is.null(recon$fused_mat) && !is.null(recon$fused_mat[[d]])) {
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
    r_mod <- do.call(fit_forest, c(
      list(
        X = as.data.frame(R, check.names = FALSE),
        Y = NULL,
        type = "unsupervised"
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
  }

  list(
    shared = list(
      source = "W_all",
      W_all = recon$W$W_all
    ),
    specific = list(
      residual = residual,
      predicted = predicted,
      residual_mod = residual_mod,
      W = specific_W
    )
  )
}
