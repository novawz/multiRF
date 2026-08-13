#' Resolve a tuning parameter (mtry or ytry) from integer, formula string, or NULL
#'
#' Accepts an integer (used as-is), a formula string like `"sqrt(p)"`,
#' `"p/3"`, `"p/2"`, or `NULL` (returns `default`).
#' In the formula, `p` refers to the number of columns (px for mtry, qy for ytry).
#'
#' @param value Integer, character formula, or `NULL`.
#' @param p Number of columns to substitute into the formula.
#' @param default Value to return when `value` is `NULL`.
#' @param name Parameter name for error messages (e.g. "mtry", "ytry").
#' @return An integer value (>= 1).
#' @keywords internal
resolve_param <- function(value, p, default, name = "param") {
  if (is.null(value)) return(max(1L, as.integer(default)))
  if (is.numeric(value)) {
    if (length(value) != 1L || !is.finite(value) || value != floor(value)) {
      stop("`", name, "` must be one finite integer.")
    }
    return(max(1L, as.integer(value)))
  }
  if (is.character(value)) {
    if (length(value) != 1L || is.na(value)) {
      stop("`", name, "` formula must be one non-missing string.")
    }
    expr <- tryCatch(parse(text = value)[[1]], error = function(e) NULL)
    if (is.null(expr)) stop("Cannot parse ", name, " formula: ", value)
    val <- eval(expr, envir = list(p = p))
    if (length(val) != 1L || !is.numeric(val) || !is.finite(val)) {
      stop("`", name, "` formula must evaluate to one finite number.")
    }
    return(max(1L, as.integer(ceiling(val))))
  }
  stop("`", name, "` must be NULL, an integer, or a formula string like \"sqrt(p)\" or \"p/3\".")
}

.native_integer_scalar <- function(value, label, lower = 0L) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value != floor(value) || value < lower || value > .Machine$integer.max) {
    qualifier <- if (lower == 1L) "positive" else "non-negative"
    stop("`", label, "` must be one ", qualifier, " integer.")
  }
  as.integer(value)
}

.native_nonnegative_scalar <- function(value, label) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value < 0) {
    stop("`", label, "` must be one finite non-negative number.")
  }
  as.double(value)
}

.native_seed_scalar <- function(value) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value != floor(value) || abs(value) > .Machine$integer.max) {
    stop("`seed` must be one finite integer.")
  }
  as.integer(value)
}

.native_logical_scalar <- function(value, label) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop("`", label, "` must be TRUE or FALSE.")
  }
  value
}

.native_forest_wt_mode <- function(forest.wt) {
  forest.wt <- match.arg(forest.wt, c("all", "inbag", "oob"))
  list(name = forest.wt, code = match(forest.wt, c("all", "inbag", "oob")) - 1L)
}

.as_native_numeric_frame <- function(x, label) {
  x <- as.data.frame(x, check.names = FALSE)
  if (nrow(x) < 2L || ncol(x) < 1L) {
    stop("`", label, "` must contain at least two samples and one feature.")
  }
  if (!all(vapply(x, is.numeric, logical(1)))) {
    stop("All columns in `", label, "` must be numeric for the native engine.")
  }
  if (any(!is.finite(as.matrix(x)))) {
    stop("`", label, "` contains NA, NaN, or infinite values.")
  }
  if (anyNA(rownames(x)) || any(rownames(x) == "") || anyDuplicated(rownames(x))) {
    stop("`", label, "` must have unique, non-missing sample names.")
  }
  x
}

.resolve_variable_weights <- function(wt, feature_names, label) {
  if (is.null(wt)) return(NULL)
  wt_names <- names(wt)
  if (!is.null(wt_names)) {
    if (anyDuplicated(wt_names) || !setequal(wt_names, feature_names)) {
      stop("Named `", label, "` must contain each corresponding feature exactly once.")
    }
    wt <- wt[match(feature_names, wt_names)]
  }
  wt <- as.numeric(wt)
  if (length(wt) != length(feature_names)) {
    stop("`", label, "` must have length ", length(feature_names), ".")
  }
  if (any(!is.finite(wt)) || any(wt < 0) || sum(wt) <= 0) {
    stop("`", label, "` must be finite, non-negative, and contain a positive value.")
  }
  wt / sum(wt)
}

#' Fit an unsupervised random forest (native C++ engine)
#'
#' Emulates rfsrc unsupervised mode: at each tree, randomly partition the
#' columns of `X` into pseudo-predictor and pseudo-response halves, then fit
#' a multivariate regression tree.  Only `forest.wt` and `proximity` are
#' meaningful; there is no real Y.
#'
#' @param X Data frame or matrix of features (n x p).
#' @param ntree Number of trees.
#' @param ytry Number of candidate pseudo-Y columns per split. Default `NULL` = 15.
#' @param forest.wt Forest-weight mode: `"all"` uses every sample in a
#' terminal node, `"inbag"` uses bootstrap donors, and `"oob"` uses only
#' trees in which the target sample is out of bag.
#' @param proximity Proximity output mode.
#' @param nsplit Number of candidate cutpoints per pseudo-predictor; `0` scans
#' all distinct cutpoints.
#' @param nodesize Minimum terminal node size.
#' @param max_depth Maximum tree depth (0 = unlimited).
#' @param samptype Sampling scheme: `"swor"` or `"swr"`.
#' @param nthread Number of OpenMP threads used by the native engine.
#' @param enhanced_prox Logical; whether to compute enhanced proximity.
#' @param sibling_gamma Strength of the sibling-leaf correction used by
#'   enhanced proximity.
#' @param leaf_embed_dim Embedding dimension used by enhanced proximity.
#' @param seed Random seed.
#' @return A list with `forest.wt`, `proximity`, `membership`, `xvar`,
#'   `xvar.names`, `ntree`, and `engine = "multiRF"`.
#' @keywords internal
fit_mv_forest_unsup <- function(X, ntree = 500L, ytry = NULL, nsplit = 10L,
                                 forest.wt = "all",
                                 proximity = c("all", "inbag", "oob", "none"),
                                 nodesize = 3L, max_depth = 0L, seed = 529L,
                                 samptype = c("swor", "swr"),
                                 nthread = getOption("multiRF.nthread", 0L),
                                 enhanced_prox = FALSE,
                                 sibling_gamma = 0.5,
                                 leaf_embed_dim = 10L) {

  X <- .as_native_numeric_frame(X, "X")
  ntree <- .native_integer_scalar(ntree, "ntree", 1L)
  nodesize <- .native_integer_scalar(nodesize, "nodesize", 1L)
  max_depth <- .native_integer_scalar(max_depth, "max_depth", 0L)
  nthread <- .native_integer_scalar(nthread, "nthread", 0L)
  seed <- .native_seed_scalar(seed)
  enhanced_prox <- .native_logical_scalar(enhanced_prox, "enhanced_prox")
  leaf_embed_dim <- .native_integer_scalar(leaf_embed_dim, "leaf_embed_dim", 1L)
  sibling_gamma <- .native_nonnegative_scalar(sibling_gamma, "sibling_gamma")
  n <- nrow(X)
  all_names <- colnames(X)
  X_mat <- as.matrix(X)
  proximity <- match.arg(proximity)
  prox_mode <- if (proximity == "none") -1L else match(proximity, c("all", "inbag", "oob")) - 1L
  fw_mode <- .native_forest_wt_mode(forest.wt)
  nsplit <- .native_integer_scalar(nsplit, "nsplit", 0L)

  # Default ytry = 15 for unsupervised; accept string formulas like "sqrt(p)"
  p_unsup <- ncol(X_mat)
  ytry_int <- min(
    p_unsup,
    resolve_param(ytry, p = p_unsup, default = 15L, name = "ytry")
  )

  # Map samptype string to integer: 0 = swor, 1 = swr
  samptype <- match.arg(samptype)
  samptype_int <- if (samptype == "swr") 1L else 0L

  # Build embedding for enhanced proximity (PCA on X, computed once)
  embed_mat <- NULL
  if (isTRUE(enhanced_prox)) {
    embed_k <- max(1L, min(as.integer(leaf_embed_dim),
                            ncol(X_mat) - 1L,
                            nrow(X_mat) - 1L))
    embed_mat <- tryCatch({
      pc <- stats::prcomp(X_mat, center = TRUE, scale. = TRUE, rank. = embed_k)
      pc$x[, seq_len(min(embed_k, ncol(pc$x))), drop = FALSE]
    }, error = function(e) {
      scale(X_mat)[, seq_len(embed_k), drop = FALSE]
    })
  }

  # All-C++ unsupervised forest (random column splits done in C++)
  res <- fit_mv_forest_unsup_cpp(
    data = X_mat,
    ntree = as.integer(ntree),
    ytry = ytry_int,
    nodesize_min = as.integer(nodesize),
    max_depth = as.integer(max_depth),
    seed = as.integer(seed),
    nthread = as.integer(nthread),
    samptype = samptype_int,
    prox_mode = as.integer(prox_mode),
    embed = embed_mat,
    sibling_gamma = as.double(sibling_gamma),
    enhanced_prox_mode = if (isTRUE(enhanced_prox)) 1L else 0L,
    forest_wt_mode = as.integer(fw_mode$code),
    nsplit = nsplit
  )

  sample_names <- rownames(X)
  if (is.null(sample_names)) sample_names <- paste0("S", seq_len(n))
  rownames(res$forest.wt) <- colnames(res$forest.wt) <- sample_names
  rownames(res$proximity) <- colnames(res$proximity) <- sample_names
  if (isTRUE(enhanced_prox) && !is.null(res$enhanced_prox) && nrow(res$enhanced_prox) == n) {
    rownames(res$enhanced_prox) <- colnames(res$enhanced_prox) <- sample_names
  }
  rownames(res$membership) <- sample_names
  inbag_out <- res$inbag
  rownames(inbag_out) <- sample_names

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

  eprox_out <- if (isTRUE(enhanced_prox) && !is.null(res$enhanced_prox) &&
                    nrow(res$enhanced_prox) == n) res$enhanced_prox else NULL

  # Raw IMD-X weights (per-tree inverse minimal depth, averaged over trees)
  imd_x <- as.numeric(res$imd_x)
  names(imd_x) <- all_names
  imd_weights <- list(X = imd_x)

  imd_x_pt <- res$imd_x_per_tree  # p x ntree matrix
  rownames(imd_x_pt) <- all_names
  imd_weights_per_tree <- list(X = imd_x_pt)

  out <- list(
    forest.wt  = res$forest.wt,
    proximity  = res$proximity,
    enhanced_prox = eprox_out,
    membership = mem,
    inbag      = inbag_out,
    xvar       = X,
    yvar       = X,
    xvar.names = all_names,
    ntree      = as.integer(ntree),
    ytry       = as.integer(ytry_int),
    nsplit     = nsplit,
    nodesize   = as.integer(nodesize),
    max_depth  = as.integer(max_depth),
    nthread    = as.integer(nthread),
    samptype   = samptype,
    seed       = as.integer(seed),
    forest.wt.mode = fw_mode$name,
    proximity.mode = proximity,
    tree_info  = tree_info,
    n          = n,
    engine     = "multiRF",
    imd_weights = imd_weights,
    imd_weights_per_tree = imd_weights_per_tree,
    structural_params = list(
      ntree = as.integer(ntree), ytry = as.integer(ytry_int), nsplit = nsplit,
      nodesize = as.integer(nodesize), max_depth = as.integer(max_depth),
      nthread = as.integer(nthread), samptype = samptype,
      seed = as.integer(seed), forest.wt = fw_mode$name, proximity = proximity
    )
  )
  out$forest <- list(nativeArray = NULL)
  out$node.stats <- NULL
  out
}

#' Fit a multivariate regression forest (native C++ engine)
#'
#' Drop-in replacement for `fit_forest()` when `type = "regression"`.
#' Returns an object with the same interface (`$forest.wt`, `$proximity`,
#' `$membership`, `$xvar`, `$yvar`, `$ntree`, `$xvar.names`) so downstream
#' code in multiRF works without changes.
#'
#' @param X Data frame or matrix of predictor features (n x px).
#' @param Y Data frame or matrix of response features (n x qy).
#' @param ntree Number of trees.
#' @param mtry Number of candidate X variables per split. Accepts an integer,
#'   a formula string (`"sqrt(p)"`, `"p/3"`, `"p/2"`), or `NULL` for the
#'   default (`ceiling(px/3)` for regression, `ceiling(sqrt(px))` for
#'   classification). In formulas, `p` is the number of predictor columns.
#' @param ytry Number of candidate Y variables per split. Accepts an integer,
#'   a formula string (`"sqrt(p)"`, `"p/3"`), or `NULL` for the default
#'   (`ceiling(qy/3)`). In formulas, `p` is the number of Y columns.
#' @param nsplit Number of candidate numeric cutpoints per split variable.
#' @param forest.wt Forest-weight mode: `"all"`, `"inbag"`, or `"oob"`.
#' @param proximity Proximity output mode.
#' @param nodesize Minimum terminal node size. Default: 5.
#' @param max_depth Maximum tree depth (0 = unlimited). Default: 0.
#' @param samptype Sampling scheme: `"swor"` or `"swr"`.
#' @param nthread Number of OpenMP threads used by the native engine.
#' @param enhanced_prox Logical; whether to compute enhanced proximity.
#' @param sibling_gamma Strength of the sibling-leaf correction used by
#'   enhanced proximity.
#' @param leaf_embed_dim Embedding dimension used by enhanced proximity.
#' @param seed Random seed.
#' @param xvar.wt Optional non-negative predictor sampling weights.
#' @param yvar.wt Optional non-negative response sampling weights.
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
#' @keywords internal
fit_mv_forest <- function(X, Y, ntree = 500L,
                           mtry = NULL, ytry = NULL, nsplit = 10L,
                           forest.wt = "all",
                           proximity = c("all", "inbag", "oob", "none"),
                           nodesize = 5L, max_depth = 0L, seed = 529L,
                           samptype = c("swor", "swr"),
                           nthread = getOption("multiRF.nthread", 0L),
                           xvar.wt = NULL, yvar.wt = NULL,
                           enhanced_prox = FALSE,
                           sibling_gamma = 0.5,
                           leaf_embed_dim = 10L) {

  X <- .as_native_numeric_frame(X, "X")
  Y <- .as_native_numeric_frame(Y, "Y")
  ntree <- .native_integer_scalar(ntree, "ntree", 1L)
  nodesize <- .native_integer_scalar(nodesize, "nodesize", 1L)
  max_depth <- .native_integer_scalar(max_depth, "max_depth", 0L)
  nthread <- .native_integer_scalar(nthread, "nthread", 0L)
  seed <- .native_seed_scalar(seed)
  enhanced_prox <- .native_logical_scalar(enhanced_prox, "enhanced_prox")
  leaf_embed_dim <- .native_integer_scalar(leaf_embed_dim, "leaf_embed_dim", 1L)
  sibling_gamma <- .native_nonnegative_scalar(sibling_gamma, "sibling_gamma")
  proximity <- match.arg(proximity)
  # -1 = skip proximity, 0 = all, 1 = inbag, 2 = oob
  prox_mode <- if (proximity == "none") -1L else match(proximity, c("all", "inbag", "oob")) - 1L
  fw_mode <- .native_forest_wt_mode(forest.wt)

  if (nrow(X) != nrow(Y) || !identical(rownames(X), rownames(Y))) {
    stop("`X` and `Y` must contain identical samples in identical row order.")
  }

  X_mat <- as.matrix(X)
  Y_mat <- as.matrix(Y)

  # Resolve mtry and ytry: accept integer, formula string ("sqrt(p)", "p/3"), or NULL
  px <- ncol(X_mat)
  qy <- ncol(Y_mat)
  mtry <- min(px, resolve_param(mtry, p = px, default = ceiling(px / 3), name = "mtry"))
  # Default ytry = ceiling(qy/3), analogous to mtry's ceiling(px/3) heuristic.
  ytry <- min(qy, resolve_param(ytry, p = qy, default = ceiling(qy / 3), name = "ytry"))
  nsplit <- .native_integer_scalar(nsplit, "nsplit", 0L)
  xvar.wt <- .resolve_variable_weights(xvar.wt, colnames(X), "xvar.wt")
  yvar.wt <- .resolve_variable_weights(yvar.wt, colnames(Y), "yvar.wt")

  # Map samptype string to integer: 0 = swor, 1 = swr
  samptype <- match.arg(samptype)
  samptype_int <- if (samptype == "swr") 1L else 0L

  # Build predictor and response PCA embeddings separately. Enhanced proximity
  # averages their sibling-leaf Spearman correlations, matching the legacy R
  # `build_embedding_list(symm = TRUE)` implementation.
  embed_mat <- NULL
  embed_split <- 0L
  if (isTRUE(enhanced_prox)) {
    build_embed <- function(dat) {
      embed_k <- max(1L, min(
        as.integer(leaf_embed_dim), ncol(dat), nrow(dat) - 1L
      ))
      tryCatch({
        pc <- stats::prcomp(
          dat, center = TRUE, scale. = TRUE, rank. = embed_k
        )
        pc$x[, seq_len(min(embed_k, ncol(pc$x))), drop = FALSE]
      }, error = function(e) {
        scale(dat)[, seq_len(embed_k), drop = FALSE]
      })
    }
    embed_x <- build_embed(X_mat)
    embed_y <- build_embed(Y_mat)
    embed_split <- ncol(embed_x)
    embed_mat <- cbind(embed_x, embed_y)
  }

  # C++ engine
  res <- fit_mv_forest_cpp(
    X = X_mat,
    Y = Y_mat,
    ntree = as.integer(ntree),
    mtry = as.integer(mtry),
    ytry = as.integer(ytry),
    nsplit = as.integer(nsplit),
    nodesize_min = as.integer(nodesize),
    max_depth = as.integer(max_depth),
    nthread = as.integer(nthread),
    seed = as.integer(seed),
    samptype = samptype_int,
    prox_mode = as.integer(prox_mode),
    embed = embed_mat,
    embed_split = as.integer(embed_split),
    sibling_gamma = as.double(sibling_gamma),
    enhanced_prox_mode = if (isTRUE(enhanced_prox)) 1L else 0L,
    forest_wt_mode = as.integer(fw_mode$code),
    xvar_wt = xvar.wt,
    yvar_wt = yvar.wt
  )

  # Set row/col names
  sample_names <- rownames(X)
  if (is.null(sample_names)) sample_names <- paste0("S", seq_len(nrow(X)))
  rownames(res$forest.wt) <- colnames(res$forest.wt) <- sample_names
  rownames(res$proximity)  <- colnames(res$proximity)  <- sample_names
  if (isTRUE(enhanced_prox) && !is.null(res$enhanced_prox) && nrow(res$enhanced_prox) == nrow(X)) {
    rownames(res$enhanced_prox) <- colnames(res$enhanced_prox) <- sample_names
  }
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

  # Enhanced proximity (NULL when not requested — zero overhead)
  eprox_out <- if (isTRUE(enhanced_prox) && !is.null(res$enhanced_prox) &&
                    nrow(res$enhanced_prox) == nrow(X)) res$enhanced_prox else NULL

  out <- list(
    forest.wt  = res$forest.wt,
    proximity  = res$proximity,
    enhanced_prox = eprox_out,
    membership = mem,
    inbag      = inbag_out,
    imd_weights = imd_weights,
    imd_weights_per_tree = imd_weights_per_tree,
    pairwise_xy = pairwise_xy,
    xvar       = X,
    yvar       = Y,
    xvar.names = colnames(X),
    ntree      = as.integer(ntree),
    mtry       = as.integer(mtry),
    ytry       = as.integer(ytry),
    nsplit     = as.integer(nsplit),
    nodesize   = as.integer(nodesize),
    max_depth  = as.integer(max_depth),
    nthread    = as.integer(nthread),
    samptype   = samptype,
    seed       = as.integer(seed),
    forest.wt.mode = fw_mode$name,
    proximity.mode = proximity,
    xvar.wt    = xvar.wt,
    yvar.wt    = yvar.wt,
    tree_info  = tree_info,
    n          = nrow(X),
    engine     = "multiRF",
    structural_params = list(
      ntree = as.integer(ntree), mtry = as.integer(mtry),
      ytry = as.integer(ytry), nsplit = as.integer(nsplit),
      nodesize = as.integer(nodesize), max_depth = as.integer(max_depth),
      nthread = as.integer(nthread), samptype = samptype,
      seed = as.integer(seed), forest.wt = fw_mode$name, proximity = proximity,
      xvar.wt = xvar.wt, yvar.wt = yvar.wt
    )
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
fit_class_forest <- function(X, Y, ntree = 500L, mtry = NULL, nsplit = 10L,
                             forest.wt = "all",
                             proximity = c("all", "inbag", "oob", "none"),
                             nodesize = 1L, max_depth = 0L, seed = 529L,
                             samptype = c("swor", "swr"),
                             nthread = getOption("multiRF.nthread", 0L),
                             xvar.wt = NULL) {

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
  proximity <- match.arg(proximity)

  # Default mtry for classification: ceiling(sqrt(p)), matching rfsrc class/class+
  mtry <- resolve_param(mtry, p = ncol(X), default = ceiling(sqrt(ncol(X))), name = "mtry")

  y_mat <- stats::model.matrix(~ y_fac - 1L)
  colnames(y_mat) <- levels(y_fac)
  rownames(y_mat) <- rownames(X)
  fit <- fit_mv_forest(
    X = X,
    Y = y_mat,
    ntree = ntree,
    mtry = mtry,
    ytry = ncol(y_mat),
    nsplit = nsplit,
    forest.wt = forest.wt,
    proximity = proximity,
    nodesize = nodesize,
    max_depth = max_depth,
    seed = seed,
    samptype = match.arg(samptype),
    nthread = nthread,
    xvar.wt = xvar.wt
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
