#' Compute Shared Fraction Per Omics Block
#'
#' Shared fraction is computed by the signal-based definition:
#' `1 - (||R||_F^2 / ||X||_F^2)`,
#' where `R = X - X_pred`.
#'
#' @param dat.list Optional named list of original omics matrices (`X`).
#' @param residual Optional named list of residual matrices (`R`).
#' @param shared_specific Optional output from `get_shared_specific_weights()`.
#' If provided, `residual` and (when needed) reconstructed `dat.list` are
#' extracted from it.
#' @param eps Small positive threshold to avoid division by zero.
#'
#' @return A data.frame with one row per omics block:
#' `data`, `frob_shared`, `frob_specific`, `specific_ratio`, `shared_frac`,
#' and `method`. For backward compatibility, it also includes
#' `frob_data`, `frob_residual`, and `residual_ratio`.
#' @export
get_shared_frac <- function(dat.list = NULL,
                            residual = NULL,
                            shared_specific = NULL,
                            eps = 1e-12) {
  if (is.null(residual) && !is.null(shared_specific)) {
    if (!is.list(shared_specific) ||
        is.null(shared_specific$specific) ||
        is.null(shared_specific$specific$residual)) {
      stop("`shared_specific` must be output from `get_shared_specific_weights()`.")
    }
    residual <- shared_specific$specific$residual
  }

  if (is.null(residual)) {
    stop("Please provide `residual`, or provide `shared_specific` containing residuals.")
  }
  if (!is.list(residual) || length(residual) == 0L) {
    stop("`residual` must be a non-empty named list of matrices.")
  }
  if (is.null(names(residual)) || any(names(residual) == "")) {
    stop("`residual` must be named.")
  }

  if (is.null(dat.list) && !is.null(shared_specific)) {
    pred <- shared_specific$specific$predicted
    if (is.null(pred)) {
      stop("`shared_specific$specific$predicted` is missing. Please provide `dat.list`.")
    }
    dat.list <- lapply(names(residual), function(d) {
      if (is.null(pred[[d]])) {
        stop("Missing predicted matrix for block `", d, "` in `shared_specific`.")
      }
      as.matrix(residual[[d]]) + as.matrix(pred[[d]])
    })
    names(dat.list) <- names(residual)
  }

  if (is.null(dat.list)) {
    stop("Please provide `dat.list`, or provide `shared_specific` with predicted matrices.")
  }
  if (!is.list(dat.list) || length(dat.list) == 0L) {
    stop("`dat.list` must be a non-empty named list of matrices.")
  }
  if (is.null(names(dat.list)) || any(names(dat.list) == "")) {
    stop("`dat.list` must be named.")
  }

  dat_names <- intersect(names(dat.list), names(residual))
  if (length(dat_names) == 0L) {
    stop("No overlapping block names between `dat.list` and `residual`.")
  }

  out <- lapply(dat_names, function(d) {
    X <- as.matrix(dat.list[[d]])
    R <- as.matrix(residual[[d]])

    if (!all(dim(X) == dim(R))) {
      stop("Dimension mismatch in block `", d, "`: `dat.list[[d]]` and `residual[[d]]` must have identical dimensions.")
    }

    ss_data <- sum(X^2)
    ss_residual <- sum(R^2)
    frob_data <- sqrt(ss_data)
    frob_residual <- sqrt(ss_residual)

    if (!is.finite(frob_data) || !is.finite(frob_residual)) {
      stop("Non-finite values detected in block `", d, "` while computing Frobenius norm.")
    }

    if (frob_data <= eps) {
      residual_ratio <- NA_real_
      shared_frac <- NA_real_
      specific_ratio <- NA_real_
    } else {
      residual_ratio <- ss_residual / ss_data
      shared_frac <- 1 - residual_ratio
      specific_ratio <- residual_ratio
    }

    data.frame(
      data = d,
      frob_shared = NA_real_,
      frob_specific = frob_residual,
      specific_ratio = specific_ratio,
      frob_data = frob_data,
      frob_residual = frob_residual,
      residual_ratio = residual_ratio,
      shared_frac = shared_frac,
      method = "signal_frobenius",
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}
