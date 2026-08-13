# -------------------------------------------------------------------------------------------------------------
#' Tune the number of clusters
#'
#' @param x A proximity/similarity matrix (or embedding input) for clustering.
#' @param return_cluster Logical; return the complete tuning result instead of
#'   only the selected number of clusters.
#' @param plot_k Logical; draw the tuning criterion over the evaluated values
#'   of `k`.
#' @param method Clustering backend: `"PAM"` or `"Spectral"`.
#' @param tune_method Tuning criterion for PAM (`"silhouette"` or `"ratio"`).
#' @param gap_w Weighting scheme for spectral eigengap (`"uniform"` or `"log"`).
#' @param prox Logical; whether `x` is a proximity matrix (converted to a
#'   dissimilarity matrix for PAM).
#' @param k_tune Candidate cluster counts to evaluate.
#' @param d Optional degree vector used by spectral clustering.
#' @param ... Additional arguments passed to lower-level clustering functions.
#'
#' @return The selected number of clusters, or the complete clustering result
#'   when `return_cluster = TRUE`.
#' @inheritParams cluster::pam
#' @export
tune_k_clusters <- function(x, ...) {
  UseMethod("tune_k_clusters")
}

.tune_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.tune_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    stop(
      "`", name, "` must be one of: ",
      paste(sprintf('"%s"', choices), collapse = ", "), ".",
      call. = FALSE
    )
  }
  x
}

.tune_n_obs <- function(x) {
  n <- if (inherits(x, "dist")) {
    attr(x, "Size")
  } else if (length(dim(x)) == 2L) {
    nrow(x)
  } else {
    NULL
  }

  if (is.null(n) || length(n) != 1L || !is.finite(n)) {
    stop("`x` must be a matrix, data frame, or `dist` object.", call. = FALSE)
  }
  n <- as.integer(n)
  if (n < 3L) {
    stop("At least three samples are required to tune the number of clusters.", call. = FALSE)
  }
  n
}

.tune_k_grid <- function(k_tune, n) {
  if (!is.numeric(k_tune) || length(k_tune) == 0L || any(!is.finite(k_tune))) {
    stop("`k_tune` must be a non-empty finite numeric vector.", call. = FALSE)
  }
  if (any(abs(k_tune - round(k_tune)) > sqrt(.Machine$double.eps))) {
    stop("Every value in `k_tune` must be an integer.", call. = FALSE)
  }

  k_grid <- sort(unique(as.integer(round(k_tune))))
  k_grid <- k_grid[k_grid >= 2L & k_grid < n]
  if (length(k_grid) == 0L) {
    stop("`k_tune` must contain at least one integer in [2, n - 1].", call. = FALSE)
  }
  k_grid
}

.tune_plot_series <- function(x, y, best_k, ylab, xlab = "k",
                              mark_selected = TRUE,
                              selected_label = NULL) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) {
    stop("Internal tuning plot data have inconsistent lengths.", call. = FALSE)
  }

  if (length(x) == 0L || !any(is.finite(x) & is.finite(y))) {
    graphics::plot.new()
    graphics::title(xlab = xlab, ylab = ylab)
    graphics::text(0.5, 0.5, "Not enough candidate values")
    return(invisible(NULL))
  }

  selected <- isTRUE(mark_selected) & x == best_k
  point_col <- ifelse(selected, "#E28E2C", "#3182BD")
  graphics::plot(
    x, y,
    type = "b", pch = 19, lwd = 1, cex = 0.85,
    ylab = ylab, xlab = xlab,
    col = "#3182BD", bg = point_col,
    bty = "l", las = 1,
    panel.first = graphics::abline(h = pretty(y), col = "#D8D8D8", lwd = 0.5)
  )
  graphics::points(x, y, pch = 21, bg = point_col, col = "white", cex = 1.05)
  if (isTRUE(mark_selected) && best_k %in% x) {
    if (is.null(selected_label)) selected_label <- paste0("Selected k = ", best_k)
    graphics::abline(v = best_k, col = "#767676", lty = 3, lwd = 0.8)
    graphics::legend(
      "topright", legend = selected_label,
      pch = 21, pt.bg = "#E28E2C", col = "white",
      bty = "n", cex = 0.75
    )
  }
  invisible(NULL)
}

.tune_plot_result <- function(cl, method, tune_method) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (identical(method, "PAM")) {
    if (identical(tune_method, "ratio")) {
      graphics::par(mfrow = c(1L, 2L))
      .tune_plot_series(
        cl$diff_k, cl$diff, cl$diff_k[which.max(cl$diff)],
        ylab = "Ratio difference",
        selected_label = paste0("Selected k = ", cl$best_k)
      )
      .tune_plot_series(
        cl$k_eval, cl$sil, cl$best_k,
        ylab = "Ratio"
      )
    } else {
      .tune_plot_series(
        cl$k_eval, cl$sil, cl$best_k,
        ylab = "Silhouette"
      )
    }
  } else {
    graphics::par(mfrow = c(1L, 2L))
    .tune_plot_series(
      cl$k_tune, cl$diff_e, cl$best_k,
      ylab = "Eigenvalue difference"
    )
    .tune_plot_series(
      seq_along(cl$ev), cl$ev, cl$best_k,
      ylab = "Eigenvalue", xlab = "Eigenvalue index",
      mark_selected = FALSE
    )
  }

  invisible(NULL)
}

#' @rdname tune_k_clusters
#' @method tune_k_clusters default
#' @export
tune_k_clusters.default <- function(x, return_cluster = FALSE, plot_k = FALSE,
                                    method = "Spectral", tune_method = "silhouette",
                                    gap_w = "uniform", prox = FALSE, ...) {
  return_cluster <- .tune_scalar_logical(return_cluster, "return_cluster")
  plot_k <- .tune_scalar_logical(plot_k, "plot_k")
  prox <- .tune_scalar_logical(prox, "prox")
  method <- .tune_choice(method, c("Spectral", "PAM"), "method")
  tune_method <- .tune_choice(tune_method, c("silhouette", "ratio"), "tune_method")
  gap_w <- .tune_choice(gap_w, c("uniform", "log"), "gap_w")

  n_obs <- .tune_n_obs(x)
  dots <- list(...)
  supplied_k <- dots$k_tune
  dots$k_tune <- NULL
  k_tune <- if (is.null(supplied_k)) {
    seq.int(2L, min(12L, n_obs - 1L))
  } else {
    .tune_k_grid(supplied_k, n_obs)
  }

  if (identical(method, "Spectral")) {
    cl <- do.call(
      spectral_cl,
      c(list(x = x, k_tune = k_tune, gap_w = gap_w), dots)
    )
  } else {
    x_use <- x
    diss <- FALSE
    if (prox) {
      x_use <- as.matrix(x_use)
      if (!is.numeric(x_use) || nrow(x_use) != ncol(x_use)) {
        stop("When `prox = TRUE`, `x` must be a square numeric matrix.", call. = FALSE)
      }
      if (any(!is.finite(x_use))) {
        stop("When `prox = TRUE`, `x` must contain only finite values.", call. = FALSE)
      }
      x_use <- (x_use + t(x_use)) / 2
      x_use <- 1 - x_use
      diag(x_use) <- 0
      diss <- TRUE
    }
    cl <- do.call(
      pam_cl,
      c(
        list(
          x = x_use, k_tune = k_tune, diss = diss,
          tune_method = tune_method
        ),
        dots
      )
    )
  }

  if (plot_k) {
    .tune_plot_result(cl, method = method, tune_method = tune_method)
  }

  if (return_cluster) cl else cl$best_k
}

#' @rdname tune_k_clusters
#' @method tune_k_clusters mrf3
#' @export
tune_k_clusters.mrf3 <- function(x, return_cluster = FALSE, plot_k = FALSE,
                                 method = "Spectral", tune_method = "silhouette",
                                 gap_w = "uniform", prox = FALSE, ...) {
  if (!inherits(x, "mrf3") || !is.list(x) || is.null(x[["dat"]])) {
    stop("`x` must be an `mrf3` object containing a `dat` matrix.", call. = FALSE)
  }

  tune_k_clusters.default(
    x = x[["dat"]],
    return_cluster = return_cluster,
    plot_k = plot_k,
    method = method,
    tune_method = tune_method,
    gap_w = gap_w,
    prox = prox,
    ...
  )
}

#' @rdname tune_k_clusters
#' @export
spectral_cl <- function(x, k_tune = seq(2, 12, by = 1), gap_w = "uniform", d = NULL, ...) {
  gap_w <- .tune_choice(gap_w, c("uniform", "log"), "gap_w")
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) != ncol(x)) {
    stop("`x` must be a square numeric matrix.", call. = FALSE)
  }
  n <- nrow(x)
  if (n < 3L) {
    stop("At least three samples are required for spectral clustering.", call. = FALSE)
  }
  k_grid <- .tune_k_grid(k_tune, n)

  x[!is.finite(x)] <- 0
  x <- (x + t(x)) / 2
  x[x < 0] <- 0

  if (is.null(d)) {
    d <- rowSums(x)
  } else if (!is.numeric(d) || length(d) != n || any(!is.finite(d))) {
    stop("`d` must be NULL or a finite numeric vector with one value per sample.", call. = FALSE)
  }
  d <- as.numeric(d)
  eps <- 1e-8
  d[d <= eps] <- eps
  inv_sqrt_d <- 1 / sqrt(d)
  l <- diag(1, n) - diag(inv_sqrt_d) %*% x %*% diag(inv_sqrt_d)
  l[!is.finite(l)] <- 0
  l <- (l + t(l)) / 2

  e <- tryCatch(eigen(l, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e)) {
    pam_x <- 1 - x
    pam_x[pam_x < 0] <- 0
    diag(pam_x) <- 0
    pam_fit <- pam_cl(
      pam_x,
      k_tune = k_grid,
      diss = TRUE,
      tune_method = "silhouette",
      ...
    )
    return(list(
      best_k = pam_fit$best_k,
      cl = pam_fit$cl,
      diff_e = rep(NA_real_, length(k_grid)),
      ev = NA_real_,
      embed = NULL,
      obj = pam_fit$obj,
      sil = NA_real_,
      k_tune = k_grid
    ))
  }

  eigenvectors <- e$vectors
  lambda <- rev(e$values)
  w <- if (identical(gap_w, "log")) log(k_grid) else rep(1, length(k_grid))
  raw_gap <- lambda[k_grid + 1L] - lambda[k_grid]
  diff_e <- raw_gap * w
  k <- k_grid[which.max(diff_e)]
  eigenvalues <- lambda[seq_len(min(n, max(k_grid) + 1L))]

  mat <- as.matrix(eigenvectors[, (n - k + 1L):n, drop = FALSE])
  mat_norm <- mat / pmax(sqrt(rowSums(mat^2)), .Machine$double.eps)
  cl <- cluster::pam(mat, k = k, ...)

  list(
    best_k = k,
    cl = cl$cluster,
    diff_e = diff_e,
    ev = eigenvalues,
    embed = mat_norm,
    obj = cl$objective[2],
    sil = cl$silinfo$avg.width,
    k_tune = k_grid
  )
}

#' @rdname tune_k_clusters
#' @export
pam_cl <- function(x, k_tune = seq(2, 9, by = 1), diss = TRUE,
                   tune_method = "silhouette", ...) {
  diss <- .tune_scalar_logical(diss, "diss")
  tune_method <- .tune_choice(tune_method, c("silhouette", "ratio"), "tune_method")

  n <- if (inherits(x, "dist")) {
    attr(x, "Size")
  } else if (length(dim(x)) == 2L) {
    nrow(x)
  } else {
    NULL
  }
  if (is.null(n) || !is.finite(n)) {
    stop("`x` must be a matrix, data frame, or `dist` object.", call. = FALSE)
  }
  n <- as.integer(n)
  if (n < 3L) {
    stop("At least three samples are required for PAM clustering.", call. = FALSE)
  }
  k_grid <- .tune_k_grid(k_tune, n)

  diss_use <- diss
  if (diss_use) {
    if (!inherits(x, "dist")) {
      x <- as.matrix(x)
      if (!is.numeric(x) || nrow(x) != ncol(x)) {
        stop("When `diss = TRUE`, `x` must be a square numeric matrix or a `dist` object.", call. = FALSE)
      }
      if (any(!is.finite(x))) {
        stop("When `diss = TRUE`, `x` must contain only finite values.", call. = FALSE)
      }
      x <- (x + t(x)) / 2
      diag(x) <- 0
    }
    x_pam <- x
  } else {
    x <- as.matrix(x)
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop("When `diss = FALSE`, `x` must be a finite numeric data matrix.", call. = FALSE)
    }
    x_pam <- x
    # Preserve the historical tuning path for multiple candidates while
    # leaving a fixed-k PAM call on the original feature matrix.
    if (length(k_grid) > 1L) {
      x_pam <- as.matrix(stats::dist(x_pam))
      diss_use <- TRUE
    }
  }

  k_eval <- k_grid
  if (identical(tune_method, "ratio")) {
    extra_k <- max(k_grid) + 1L
    if (extra_k < n) {
      k_eval <- c(k_grid, extra_k)
    }
    ratio_dissimilarity <- if (diss_use) {
      as.matrix(x_pam)
    } else {
      as.matrix(stats::dist(x_pam))
    }
    similarity <- 1 - ratio_dissimilarity
    dist_all <- mean(similarity)
  }

  scores <- numeric(length(k_eval))
  for (ii in seq_along(k_eval)) {
    k <- k_eval[ii]
    cl_fit <- cluster::pam(x_pam, k = k, diss = diss_use, ...)
    if (identical(tune_method, "silhouette")) {
      scores[ii] <- cl_fit$silinfo$avg.width
    } else {
      cl_unique <- unique(cl_fit$cluster)
      class_dist <- vapply(
        cl_unique,
        function(i) mean(similarity[cl_fit$cluster == i, cl_fit$cluster != i, drop = FALSE]),
        numeric(1)
      )
      scores[ii] <- mean(class_dist) / dist_all
    }
  }

  if (identical(tune_method, "silhouette")) {
    k <- k_eval[which.max(scores)]
    diff_s <- numeric(0)
    diff_k <- integer(0)
  } else {
    diff_s <- abs(diff(scores)) * log(k_eval[-length(k_eval)])
    diff_k <- k_eval[-length(k_eval)]
    if (length(diff_s) == 0L || !any(is.finite(diff_s))) {
      k <- k_grid[1L]
    } else {
      k <- diff_k[which.max(diff_s)]
    }
  }

  cl <- cluster::pam(x_pam, k = k, diss = diss_use, ...)

  list(
    best_k = k,
    cl = cl$cluster,
    sil = scores,
    diff = diff_s,
    obj = cl$objective[2],
    k_tune = k_grid,
    k_eval = k_eval,
    diff_k = diff_k
  )
}
