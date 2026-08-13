#' compute tSNE embedding for MRF similarity matrix
#' @param mod mrf3 model
#' @param dims Number of dims
#' @param max_iter integer; Number of iterations (default: 1000)
#' @param learning_rate numeric; Learning rate (default: 200.0)
#' @param verbose Logical; whether to print optimization progress.
#' @param tol Convergence tolerance.
#' @param patience Early-stopping patience measured in iterations.
#' @param seed Optional integer seed used only to initialize the embedding.
#' @return Two dimensional tSNE embedding
#'
#' @export mrf3_tsne
#'

mrf3_tsne <- function(mod, dims = 2, max_iter = 1000,
                      learning_rate = 200, verbose = FALSE,
                      tol = 1e-05, patience = 10, seed = NULL) {
  if (!is.list(mod) || is.null(mod$dat)) {
    stop("`mod` must be a list containing a similarity matrix in `mod$dat`.",
         call. = FALSE)
  }

  tsne(
    X = mod$dat,
    dims = dims,
    max_iter = max_iter,
    learning_rate = learning_rate,
    verbose = verbose,
    patience = patience,
    tol = tol,
    seed = seed
  )
}

.tsne_integer_scalar <- function(x, name, minimum = 1) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop("`", name, "` must be a finite integer >= ", minimum, ".",
         call. = FALSE)
  }
  as.integer(x)
}

.prepare_tsne_probabilities <- function(X) {
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`mod$dat` must be a numeric matrix.", call. = FALSE)
  }
  if (length(dim(X)) != 2L || nrow(X) != ncol(X)) {
    stop("`mod$dat` must be a square matrix.", call. = FALSE)
  }
  if (nrow(X) < 2L) {
    stop("`mod$dat` must contain at least two observations.", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("`mod$dat` must contain only finite values.", call. = FALSE)
  }
  if (any(X < 0)) {
    stop("`mod$dat` must contain only non-negative similarities.",
         call. = FALSE)
  }

  storage.mode(X) <- "double"
  X <- (X + t(X)) / 2
  diag(X) <- 0
  total_mass <- sum(X)
  if (!is.finite(total_mass) || total_mass <= 0) {
    stop("`mod$dat` must have positive off-diagonal similarity mass.",
         call. = FALSE)
  }
  X / total_mass
}

.initialize_tsne_embedding <- function(n, dims, seed) {
  if (is.null(seed)) {
    return(matrix(stats::rnorm(n * dims), nrow = n, ncol = dims))
  }

  seed <- .tsne_integer_scalar(seed, "seed", minimum = 0)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  matrix(stats::rnorm(n * dims), nrow = n, ncol = dims)
}

tsne <- function(X, dims = 2, max_iter = 1000, learning_rate = 200,
                 verbose = FALSE, patience = 10, tol = 1e-05,
                 seed = NULL) {
  X <- .prepare_tsne_probabilities(X)
  dims <- .tsne_integer_scalar(dims, "dims")
  max_iter <- .tsne_integer_scalar(max_iter, "max_iter")
  patience <- .tsne_integer_scalar(patience, "patience")
  if (length(learning_rate) != 1L || !is.numeric(learning_rate) ||
      !is.finite(learning_rate) || learning_rate <= 0) {
    stop("`learning_rate` must be a finite positive number.", call. = FALSE)
  }
  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a finite non-negative number.", call. = FALSE)
  }
  if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(seed)) {
    .tsne_integer_scalar(seed, "seed", minimum = 0)
  }

  n <- nrow(X)
  Y <- .initialize_tsne_embedding(n, dims, seed)
  Y <- sweep(Y, 2L, colMeans(Y), FUN = "-")
  cost_history <- rep(NA_real_, max_iter)
  patience_count <- 0L
  converged <- FALSE

  # Faster squared Euclidean distance matrix via matrix algebra.
  sq_dist_mat <- function(A) {
    rs <- rowSums(A * A)
    D <- outer(rs, rs, "+") - 2 * tcrossprod(A)
    pmax(D, 0)
  }

  iterations <- max_iter
  for (iter in seq_len(max_iter)) {
    grad_Y <- tsne_cost_gradient_cpp(Y, X)
    if (any(!is.finite(grad_Y))) {
      stop("The tSNE gradient became non-finite at iteration ", iter, ".",
           call. = FALSE)
    }

    Y <- Y - learning_rate * grad_Y
    if (any(!is.finite(Y))) {
      stop("The tSNE embedding became non-finite at iteration ", iter, ".",
           call. = FALSE)
    }
    Y <- sweep(Y, 2L, colMeans(Y), FUN = "-")

    Q <- 1 / (1 + sq_dist_mat(Y))
    diag(Q) <- 0
    sum_Q <- sum(Q)
    if (!is.finite(sum_Q) || sum_Q <= 0) {
      stop("The low-dimensional tSNE probabilities became invalid at iteration ",
           iter, ".", call. = FALSE)
    }
    Q <- Q / sum_Q

    cost <- sum(X * log((X + .Machine$double.eps) /
                        (Q + .Machine$double.eps)))
    if (!is.finite(cost)) {
      stop("The tSNE cost became non-finite at iteration ", iter, ".",
           call. = FALSE)
    }
    cost_history[iter] <- cost

    if (iter > 1L) {
      cost_change <- abs(cost_history[iter] - cost_history[iter - 1L])
      if (cost_change <= tol) {
        patience_count <- patience_count + 1L
      } else {
        patience_count <- 0L
      }

      if (patience_count >= patience) {
        converged <- TRUE
        iterations <- iter
        if (verbose) {
          cat("Converged at iteration", iter, "with cost", cost, "\n")
        }
        break
      }
    }

    if (verbose && iter %% 10L == 0L) {
      cat("Iteration", iter, "Cost", cost, "\n")
    }
  }

  colnames(Y) <- paste0("tSNE", seq_len(ncol(Y)))
  attr(Y, "converged") <- converged
  attr(Y, "iterations") <- iterations
  attr(Y, "cost_history") <- cost_history[seq_len(iterations)]
  Y
}
