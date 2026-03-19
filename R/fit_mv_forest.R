#' Fit a multivariate regression forest (native C++ engine)
#'
#' Drop-in replacement for `fit_rfsrc()` when `type = "regression"`.
#' Returns an object with the same interface (`$forest.wt`, `$proximity`,
#' `$membership`, `$xvar`, `$yvar`, `$ntree`, `$xvar.names`) so downstream
#' code in multiRF works without changes.
#'
#' @param X Data frame or matrix of predictor features (n x px).
#' @param Y Data frame or matrix of response features (n x qy).
#' @param ntree Number of trees.
#' @param mtry Number of candidate X variables per split. Default: sqrt(px).
#' @param ytry Number of candidate Y variables per split. Default: qy/3.
#' @param nodesize Minimum terminal node size. Default: 5.
#' @param max_depth Maximum tree depth (0 = unlimited). Default: 0.
#' @param seed Random seed.
#' @return A list compatible with rfsrc output, containing:
#'   \describe{
#'     \item{forest.wt}{n x n forest weight matrix}
#'     \item{proximity}{n x n proximity matrix}
#'     \item{membership}{n x ntree terminal node membership}
#'     \item{xvar}{predictor data frame}
#'     \item{yvar}{response data frame}
#'     \item{xvar.names}{character vector of predictor names}
#'     \item{ntree}{number of trees}
#'     \item{tree_info}{per-tree structure for IMD}
#'   }
#' Fit an unsupervised random forest (native C++ engine)
#'
#' Emulates rfsrc unsupervised mode: at each tree, randomly partition the
#' columns of `X` into pseudo-predictor and pseudo-response halves, then fit
#' a multivariate regression tree.  Only `forest.wt` and `proximity` are
#' meaningful; there is no real Y.
#'
#' @param X Data frame or matrix of features (n x p).
#' @param ntree Number of trees.
#' @param ytry Number of candidate pseudo-Y columns per split.
#' @param nodesize Minimum terminal node size.
#' @param max_depth Maximum tree depth (0 = unlimited).
#' @param seed Random seed.
#' @return A list with `forest.wt`, `proximity`, `membership`, `xvar`,
#'   `xvar.names`, `ntree`, and `engine = "multiRF"`.
#' @keywords internal
fit_mv_forest_unsup <- function(X, ntree = 500L, ytry = 10L,
                                 nodesize = 5L, max_depth = 0L, seed = -1L,
                                 nthread = getOption("multiRF.nthread", 0L)) {

  X <- as.data.frame(X, check.names = FALSE)
  n <- nrow(X)
  all_names <- colnames(X)
  X_mat <- as.matrix(X)

  # All-C++ unsupervised forest (random column splits done in C++)
  res <- fit_mv_forest_unsup_cpp(
    data = X_mat,
    ntree = as.integer(ntree),
    ytry = as.integer(ytry),
    nodesize_min = as.integer(nodesize),
    max_depth = as.integer(max_depth),
    seed = as.integer(seed),
    nthread = as.integer(nthread)
  )

  sample_names <- rownames(X)
  if (is.null(sample_names)) sample_names <- paste0("S", seq_len(n))
  rownames(res$forest.wt) <- colnames(res$forest.wt) <- sample_names
  rownames(res$proximity) <- colnames(res$proximity) <- sample_names
  rownames(res$membership) <- sample_names

  # Remap membership: 0-indexed node index -> sequential leaf ID (DFS order)
  mem <- res$membership
  tree_info <- res$tree_info
  for (t in seq_len(ntree)) {
    ti <- tree_info[[t]]
    n_nodes <- length(ti$split_var)
    is_leaf <- ti$is_leaf
    dfs_order <- integer(n_nodes)
    stack <- 1L; idx <- 0L
    while (length(stack) > 0L) {
      node <- stack[length(stack)]; stack <- stack[-length(stack)]
      idx <- idx + 1L; dfs_order[idx] <- node
      r <- ti$right[node]; l <- ti$left[node]
      if (r >= 0) stack <- c(stack, r + 1L)
      if (l >= 0) stack <- c(stack, l + 1L)
    }
    leaf_map <- integer(n_nodes); leaf_id <- 0L
    for (ni in dfs_order) {
      if (is_leaf[ni]) { leaf_id <- leaf_id + 1L; leaf_map[ni] <- leaf_id }
    }
    mem[, t] <- leaf_map[mem[, t] + 1L]
  }

  out <- list(
    forest.wt  = res$forest.wt,
    proximity  = res$proximity,
    membership = mem,
    xvar       = X,
    yvar       = X,
    xvar.names = all_names,
    ntree      = as.integer(ntree),
    tree_info  = tree_info,
    n          = n,
    engine     = "multiRF"
  )
  out$forest <- list(nativeArray = NULL)
  out$node.stats <- NULL
  out
}

fit_mv_forest <- function(X, Y, ntree = 500L, mtry = 0L, ytry = 0L,
                           nodesize = 5L, max_depth = 0L, seed = -1L,
                           nthread = getOption("multiRF.nthread", 0L)) {

  X <- as.data.frame(X, check.names = FALSE)
  Y <- as.data.frame(Y, check.names = FALSE)

  stopifnot(nrow(X) == nrow(Y))

  X_mat <- as.matrix(X)
  Y_mat <- as.matrix(Y)

  # C++ engine
  res <- fit_mv_forest_cpp(
    X = X_mat,
    Y = Y_mat,
    ntree = as.integer(ntree),
    mtry = as.integer(mtry),
    ytry = as.integer(ytry),
    nodesize_min = as.integer(nodesize),
    max_depth = as.integer(max_depth),
    nthread = as.integer(nthread),
    seed = as.integer(seed)
  )

  # Set row/col names
  sample_names <- rownames(X)
  if (is.null(sample_names)) sample_names <- paste0("S", seq_len(nrow(X)))
  rownames(res$forest.wt) <- colnames(res$forest.wt) <- sample_names
  rownames(res$proximity)  <- colnames(res$proximity)  <- sample_names
  rownames(res$membership) <- sample_names

  # Remap membership from 0-indexed C++ node indices to sequential
  # 1-indexed leaf IDs so that IMD code's mem_id matching works correctly.
  # In get_tree_net, leaves are listed in DFS (pre-order) traversal and
  # get var.tip.id = var_count = 1, 2, 3, ... in that order.
  mem <- res$membership
  tree_info <- res$tree_info
  for (t in seq_len(ntree)) {
    ti <- tree_info[[t]]
    n_nodes <- length(ti$split_var)
    is_leaf <- ti$split_var < 0
    # Get DFS order (same as .bfs_to_preorder in mrf3-util.R)
    dfs_order <- integer(n_nodes)
    stack <- 1L
    idx <- 0L
    while (length(stack) > 0L) {
      node <- stack[length(stack)]
      stack <- stack[-length(stack)]
      idx <- idx + 1L
      dfs_order[idx] <- node
      r <- ti$right[node]
      l <- ti$left[node]
      if (r >= 0) stack <- c(stack, r + 1L)
      if (l >= 0) stack <- c(stack, l + 1L)
    }
    # Assign sequential leaf IDs in DFS order
    leaf_map <- integer(n_nodes)
    leaf_id <- 0L
    for (ni in dfs_order) {
      if (is_leaf[ni]) {
        leaf_id <- leaf_id + 1L
        leaf_map[ni] <- leaf_id
      }
    }
    # Remap membership: mem values are 0-indexed node indices
    mem[, t] <- leaf_map[mem[, t] + 1L]  # +1 for 0-index to 1-index
  }

  # Build rfsrc-compatible output
  inbag_out <- res$inbag
  rownames(inbag_out) <- sample_names

  # IMD weights computed during tree building (zero-cost)
  imd_x <- as.numeric(res$imd_x)
  imd_y <- as.numeric(res$imd_y)
  names(imd_x) <- colnames(X)
  names(imd_y) <- colnames(Y)
  imd_weights <- list(X = imd_x, Y = imd_y)

  # Per-tree IMD distributions (for method="test" in mrf3_vs)
  imd_x_pt <- res$imd_x_per_tree  # px x ntree matrix
  imd_y_pt <- res$imd_y_per_tree  # qy x ntree matrix
  rownames(imd_x_pt) <- colnames(X)
  rownames(imd_y_pt) <- colnames(Y)
  imd_weights_per_tree <- list(X = imd_x_pt, Y = imd_y_pt)

  # Pairwise X-Y co-occurrence matrix (px x qy) for pairwise_imd
  pairwise_xy <- res$pairwise_xy
  rownames(pairwise_xy) <- colnames(X)
  colnames(pairwise_xy) <- colnames(Y)

  out <- list(
    forest.wt  = res$forest.wt,
    proximity  = res$proximity,
    membership = mem,
    inbag      = inbag_out,
    imd_weights = imd_weights,
    imd_weights_per_tree = imd_weights_per_tree,
    pairwise_xy = pairwise_xy,
    xvar       = X,
    yvar       = Y,
    xvar.names = colnames(X),
    ntree      = as.integer(ntree),
    tree_info  = tree_info,
    n          = nrow(X),
    engine     = "multiRF"
  )

  # Compatibility: create minimal forest$nativeArray stub
  # so existing get_tree_net fallback path doesn't crash
  out$forest <- list(nativeArray = NULL)
  out$node.stats <- NULL

  out
}

#' Fit a classification forest (native engine)
#'
#' Uses the native multivariate regression forest on a one-hot encoded response,
#' then reconstructs training-set class probabilities and labels.
#'
#' @keywords internal
fit_class_forest <- function(X, Y, ntree = 500L, mtry = 0L,
                             nodesize = 5L, max_depth = 0L, seed = -1L,
                             nthread = getOption("multiRF.nthread", 0L)) {

  X <- as.data.frame(X, check.names = FALSE)

  if (is.data.frame(Y)) {
    if (ncol(Y) != 1L) {
      stop("Native classification currently expects a single response column.")
    }
    Y <- Y[[1]]
  }

  y_fac <- as.factor(Y)
  if (nlevels(y_fac) < 2L) {
    stop("Classification requires at least two response classes.")
  }

  y_mat <- stats::model.matrix(~ y_fac - 1L)
  colnames(y_mat) <- levels(y_fac)
  fit <- fit_mv_forest(
    X = X,
    Y = y_mat,
    ntree = ntree,
    mtry = mtry,
    ytry = ncol(y_mat),
    nodesize = nodesize,
    max_depth = max_depth,
    seed = seed,
    nthread = nthread
  )

  fw <- fit$forest.wt
  rs <- rowSums(fw)
  prob <- matrix(0, nrow = nrow(fw), ncol = ncol(y_mat))
  colnames(prob) <- colnames(y_mat)
  rownames(prob) <- rownames(X)
  ok <- rs > 0
  if (any(ok)) {
    prob[ok, ] <- (fw[ok, , drop = FALSE] / rs[ok]) %*% y_mat
  }
  if (any(!ok)) {
    prob[!ok, ] <- y_mat[!ok, , drop = FALSE]
  }

  pred_idx <- max.col(prob, ties.method = "first")
  pred_class <- factor(colnames(prob)[pred_idx], levels = levels(y_fac))
  names(pred_class) <- rownames(X)

  fit$yvar <- data.frame(Y = y_fac, stringsAsFactors = FALSE)
  rownames(fit$yvar) <- rownames(X)
  fit$predicted <- pred_class
  fit$class.prob <- prob
  fit$err.rate <- mean(pred_class != y_fac)
  class(fit) <- c("multiRF_native", "grow", "class+")

  fit
}
