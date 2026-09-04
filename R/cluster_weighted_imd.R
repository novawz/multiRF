#' Compute cluster-weighted IMD from pre-computed per-node split scores
#'
#' Aggregates the per-node split statistics stored on every fitted forest
#' (`tree_info[[t]]$imd_x_score` and
#' `tree_info[[t]]$imd_y_stats`) into cluster-specific variable weights,
#' without refitting the forest and without re-traversing trees per sample in
#' R.
#'
#' @section Estimator:
#' For cluster \eqn{c} and predictor \eqn{j}, the weight is the
#' cluster-occupancy-weighted mean of node split scores:
#' \deqn{w_j(c) = \frac{\sum_{v : \mathrm{split}(v) = j} s_v \, n_c(v)}
#'                     {\sum_{v : \mathrm{split}(v) = j} n_c(v)}}
#' where \eqn{s_v} is the (standardized multivariate) split score recorded at
#' internal node \eqn{v} during fitting and \eqn{n_c(v)} is the number of
#' cluster-\eqn{c} samples whose root-to-leaf path passes through \eqn{v}.
#' The response-side weight replaces \eqn{s_v} with the per-response component
#' `imd_y_stats[j, v]` (nodes where that component is zero are excluded from
#' both numerator and denominator).
#'
#' This is a *split-score-weighted* descriptive estimator. It is **not** the
#' Eq. 6-8 inverse-minimal-depth (IMD) statistic computed by
#' `get_multi_weights()`: depth-based IMD scores a variable by how close to
#' the root it splits, while this estimator scores it by the size of the
#' variance reduction at its split nodes, localized to the cluster through
#' path occupancy. Because the node scores \eqn{s_v} are global statistics
#' (computed from all inbag samples during fitting) and occupancy of the
#' high-score shallow nodes is similar for every cluster, the weights stay
#' close to the global importance ranking: in benchmarks (synthetic with
#' planted cluster-specific drivers, and TCGA-BRCA subtypes) per-cluster
#' weights correlate at Spearman 0.9-0.99 *across* clusters, versus 0.6-0.8
#' for the depth-based path, and agree with the depth-based per-cluster
#' weights at Spearman ~0.4-0.6 (top-20 overlap ~50-75%). Treat the output
#' as a fast, deterministic descriptive screen with a mild cluster tilt --
#' not as a substitute for the depth-based cluster IMD.
#'
#' @param mod A model fitted with `engine = "native"` (with
#'   `$tree_info`, `$membership`, `$xvar`, `$yvar`). `$membership` is expected
#'   in the remapped form stored by `fit_mv_forest()`: sequential 1-based DFS
#'   leaf IDs per tree (this function inverts that remap internally to recover
#'   raw node indices).
#' @param cluster Named character/factor vector of cluster labels for ALL
#'   samples. When named, it is aligned to `rownames(mod$xvar)`.
#' @param normalized Logical; L2-normalize the output weights per side.
#' @return A named list: one element per cluster label, each containing
#'   `list(X = named_vector, Y = named_vector)` of cluster-weighted scores.
#'   `Y` has length zero for unsupervised (single-block) models.
#' @keywords internal
cluster_weighted_imd <- function(mod, cluster, normalized = TRUE) {

  tree_info <- mod$tree_info
  membership <- mod$membership  # n x ntree, 1-based DFS leaf IDs (remapped in fit_mv_forest)
  if (!is.list(tree_info) || length(tree_info) == 0L) {
    stop("`mod` has no `tree_info`; cluster_weighted_imd() requires `engine = 'native'`.")
  }
  if (is.null(tree_info[[1L]]$imd_x_score)) {
    stop("`mod$tree_info` does not carry `imd_x_score`; refit with `engine = 'native'`.")
  }
  if (is.null(membership)) {
    stop("`mod` has no `membership` matrix.")
  }

  ntree <- mod$ntree
  n <- nrow(mod$xvar)
  px <- ncol(mod$xvar)
  x_names <- colnames(mod$xvar)
  qy <- if (!is.null(mod$yvar) && !identical(mod$xvar, mod$yvar)) ncol(mod$yvar) else 0L
  y_names <- if (qy > 0) colnames(mod$yvar) else character(0)

  # Align cluster labels to sample order
  sample_ids <- rownames(mod$xvar)
  if (!is.null(names(cluster)) && !is.null(sample_ids)) {
    if (!all(sample_ids %in% names(cluster))) {
      stop("`cluster` names do not cover the model's sample names.")
    }
    cluster <- cluster[sample_ids]
  }
  stopifnot(length(cluster) == n)
  cluster <- as.character(cluster)
  labs <- unique(cluster[!is.na(cluster)])
  if (length(labs) == 0L) {
    stop("`cluster` has no non-NA labels.")
  }
  nlab <- length(labs)
  cl_idx <- match(cluster, labs)  # NA for NA labels

  # Accumulators: rows = variables, cols = clusters
  wx <- matrix(0, px, nlab, dimnames = list(x_names, labs))
  cx <- matrix(0, px, nlab, dimnames = list(x_names, labs))
  wy <- matrix(0, qy, nlab, dimnames = list(y_names, labs))
  cy <- matrix(0, qy, nlab, dimnames = list(y_names, labs))

  for (t in seq_len(ntree)) {
    ti <- tree_info[[t]]
    n_nodes <- length(ti$split_var)
    split_var <- ti$split_var     # -1 for leaf, 0-indexed X column for internal
    left_child <- ti$left         # -1 if no child (0-indexed otherwise)
    right_child <- ti$right

    # Internal nodes in node-array order: this is exactly the column order of
    # `imd_y_stats` (the C++ export walks nodes in array order).
    internal <- which(split_var >= 0L)
    if (length(internal) == 0L) next  # single-leaf tree: nothing to attribute

    # `membership` stores sequential 1-based DFS leaf IDs (remapped in
    # fit_mv_forest), so rebuild the per-tree DFS leaf order and invert it
    # back to 1-indexed raw node indices.
    dfs_order <- integer(n_nodes)
    stack <- 1L
    d_idx <- 0L
    while (length(stack) > 0L) {
      node <- stack[length(stack)]
      stack <- stack[-length(stack)]
      d_idx <- d_idx + 1L
      dfs_order[d_idx] <- node
      r <- right_child[node]
      l <- left_child[node]
      if (r >= 0) stack <- c(stack, r + 1L)
      if (l >= 0) stack <- c(stack, l + 1L)
    }
    leaf_nodes <- dfs_order[split_var[dfs_order] < 0L]  # leaf_nodes[k] = node index of leaf ID k

    mem_t <- membership[, t]
    ok <- !is.na(mem_t) & mem_t >= 1L & mem_t <= length(leaf_nodes) & !is.na(cl_idx)
    if (!any(ok)) next

    # Per-node cluster occupancy: scatter leaf counts, then propagate
    # bottom-up (children before parents, guaranteed by decreasing depth).
    cnt <- matrix(0, n_nodes, nlab)
    leaf_tab <- table(
      factor(leaf_nodes[mem_t[ok]], levels = seq_len(n_nodes)),
      factor(cl_idx[ok], levels = seq_len(nlab))
    )
    cnt[] <- as.numeric(leaf_tab)
    for (v in order(ti$depth, decreasing = TRUE)) {
      if (split_var[v] >= 0L) {
        cnt[v, ] <- cnt[left_child[v] + 1L, ] + cnt[right_child[v] + 1L, ]
      }
    }
    cnt_int <- cnt[internal, , drop = FALSE]  # n_internal x nlab

    # X side: group node contributions by split variable
    score_int <- ti$imd_x_score[internal]
    grp <- split_var[internal] + 1L  # 1-indexed X column
    wx_add <- rowsum(cnt_int * score_int, group = grp)
    cx_add <- rowsum(cnt_int, group = grp)
    rows <- as.integer(rownames(wx_add))
    wx[rows, ] <- wx[rows, ] + wx_add
    cx[rows, ] <- cx[rows, ] + cx_add

    # Y side: per-response components at each internal node (qy x n_internal)
    imd_y_mat <- ti$imd_y_stats
    if (qy > 0 && !is.null(imd_y_mat) && nrow(imd_y_mat) == qy &&
        ncol(imd_y_mat) == length(internal)) {
      pos <- imd_y_mat > 0
      wy <- wy + (imd_y_mat * pos) %*% cnt_int
      cy <- cy + pos %*% cnt_int
    }
  }

  # Normalize: divide by count, then optionally L2-normalize
  out <- lapply(seq_len(nlab), function(li) {
    wxl <- ifelse(cx[, li] > 0, wx[, li] / cx[, li], 0)
    names(wxl) <- x_names
    wyl <- if (qy > 0) {
      v <- ifelse(cy[, li] > 0, wy[, li] / cy[, li], 0)
      names(v) <- y_names
      v
    } else {
      numeric(0)
    }

    if (normalized) {
      norm_x <- sqrt(sum(wxl^2))
      if (is.finite(norm_x) && norm_x > 0) wxl <- wxl / norm_x
      if (qy > 0) {
        norm_y <- sqrt(sum(wyl^2))
        if (is.finite(norm_y) && norm_y > 0) wyl <- wyl / norm_y
      }
    }

    list(X = wxl, Y = wyl)
  })
  names(out) <- labs
  out
}
