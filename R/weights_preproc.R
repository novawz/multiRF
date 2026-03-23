#' Preprocess a weight matrix for downstream steps
#'
#' @param W A square numeric weight matrix.
#' @param adjust Logical; whether to apply adjusted-weight scaling.
#' @param top_v Optional integer. If set, each row keeps only the top-v entries.
#' @param row_normalize Logical; whether to row-normalize after adjustment/truncation.
#' @param zero_diag Logical; whether to set diagonal to zero.
#' @param eps Small positive value used to avoid division by zero.
#' @param keep_ties Logical; whether top-v truncation keeps ties at the cutoff.
#'
#' @return A processed square numeric matrix.
prepare_weight_matrix <- function(W,
                                  adjust = TRUE,
                                  top_v = NULL,
                                  row_normalize = TRUE,
                                  zero_diag = TRUE,
                                  eps = 1e-8,
                                  keep_ties = TRUE) {
  W <- validate_weight_matrix(W)

  if (adjust) {
    W <- adjust_weight_matrix(W, zero_diag = zero_diag, eps = eps)
  } else if (zero_diag) {
    diag(W) <- 0
  }

  W <- truncate_top_v_rows(W, top_v = top_v, keep_ties = keep_ties)

  if (row_normalize) {
    W <- row_normalize_weights(W, eps = eps)
  }

  W
}


#' @rdname prepare_weight_matrix
#' @param W_list A named list of square numeric weight matrices.
#' @param ... Additional arguments passed to `prepare_weight_matrix()`.
#' @return A list of processed weight matrices.
prepare_weight_list <- function(W_list, ...) {
  if (!is.list(W_list) || length(W_list) == 0) {
    stop("`W_list` must be a non-empty list of matrices.")
  }
  lapply(W_list, function(W) prepare_weight_matrix(W, ...))
}


#' @rdname prepare_weight_matrix
#' @param W A square numeric weight matrix.
#' @return Adjusted weight matrix using row-wise scaling by `1 - diag(W)`.
adjust_weight_matrix <- function(W, zero_diag = TRUE, eps = 1e-8) {
  W <- validate_weight_matrix(W)
  d <- 1 - diag(W)
  d[!is.finite(d)] <- eps
  d[d < eps] <- eps
  W_adj <- W / d
  if (zero_diag) {
    diag(W_adj) <- 0
  }
  W_adj
}


#' @rdname prepare_weight_matrix
#' @param top_v Optional integer. If set, each row keeps only the top-v entries.
#' @param keep_ties Logical; whether top-v truncation keeps ties at the cutoff.
#' @return Top-v truncated matrix (row-wise).
truncate_top_v_rows <- function(W, top_v = NULL, keep_ties = TRUE) {
  W <- validate_weight_matrix(W)
  n <- nrow(W)

  if (is.null(top_v)) {
    return(W)
  }
  top_v <- as.integer(top_v)
  if (!is.finite(top_v) || length(top_v) != 1L) {
    stop("`top_v` must be a single finite integer or NULL.")
  }
  if (top_v <= 0L) {
    out <- matrix(0, nrow = n, ncol = ncol(W), dimnames = dimnames(W))
    return(out)
  }
  if (top_v >= ncol(W)) {
    return(W)
  }

  out <- W
  for (i in seq_len(n)) {
    row <- W[i, ]
    row[!is.finite(row)] <- -Inf
    if (keep_ties) {
      thr <- sort(row, decreasing = TRUE)[top_v]
      keep <- which(row >= thr)
    } else {
      keep <- order(row, decreasing = TRUE)[seq_len(top_v)]
    }
    drop <- setdiff(seq_len(ncol(W)), keep)
    if (length(drop) > 0) {
      out[i, drop] <- 0
    }
  }
  out
}


#' @rdname prepare_weight_matrix
#' @return Row-normalized weight matrix.
row_normalize_weights <- function(W, eps = 1e-12) {
  W <- validate_weight_matrix(W)

  d <- rowSums(W)

  d[!is.finite(d)] <- 0
  out <- W
  ok <- abs(d) > eps
  if (any(ok)) {
    out[ok, ] <- W[ok, , drop = FALSE] / d[ok]
  }
  if (any(!ok)) {
    out[!ok, ] <- 0
  }
  out
}


validate_weight_matrix <- function(W) {
  W <- as.matrix(W)
  if (!is.numeric(W)) {
    stop("`W` must be numeric.")
  }
  if (length(dim(W)) != 2L) {
    stop("`W` must be a 2D matrix.")
  }
  if (nrow(W) == 0 || ncol(W) == 0) {
    stop("`W` cannot be empty.")
  }
  if (nrow(W) != ncol(W)) {
    stop("`W` must be square.")
  }
  W
}
