#' Compute cluster-weighted IMD from pre-computed per-node split scores
#'
#' For each tree, traces which internal nodes each sample passes through
#' (root to leaf path). The IMD contribution of variable \eqn{j} at node
#' \eqn{v} is weighted by the fraction of the cluster's samples that
#' pass through \eqn{v}. This gives cluster-specific importance without
#' refitting the forest.
#'
#' @param mod A single fitted forest model (with `$tree_info`, `$membership`,
#'   `$xvar`, `$yvar`).
#' @param cluster Named character/factor vector of cluster labels for ALL samples.
#' @param normalized Logical; L2-normalize the output weights per block.
#' @return A named list: one element per cluster label, each containing
#'   `list(X = named_vector, Y = named_vector)` of cluster-weighted IMD.
#' @keywords internal
cluster_weighted_imd <- function(mod, cluster, normalized = TRUE) {

  tree_info <- mod$tree_info
  membership <- mod$membership  # n x ntree, 0-indexed node IDs (before remap)
  ntree <- mod$ntree
  n <- nrow(mod$xvar)
  px <- ncol(mod$xvar)
  x_names <- colnames(mod$xvar)
  qy <- if (!is.null(mod$yvar) && !identical(mod$xvar, mod$yvar)) ncol(mod$yvar) else 0L
  y_names <- if (qy > 0) colnames(mod$yvar) else character(0)

  # Align cluster labels to sample order
  sample_ids <- rownames(mod$xvar)
  if (!is.null(names(cluster)) && !is.null(sample_ids)) {
    cluster <- cluster[sample_ids]
  }
  stopifnot(length(cluster) == n)
  cluster <- as.character(cluster)
  labs <- unique(cluster[!is.na(cluster)])

  # For each tree, trace each sample's root-to-leaf path to find
  # which internal nodes it visits.
  # node_visits[[t]] is an n-length list of integer vectors (visited node indices)
  # But building full path lists for all samples is expensive.
  # Instead: for each tree + each internal node, count how many samples
  # from each cluster pass through it. We can do this by noting that
  # nodesize = sum of all leaf sizes in the subtree rooted at that node.
  # Equivalently: a sample visits node v iff its leaf is a descendant of v.

  # Precompute: for each tree, build a mapping from leaf_id -> set of
  # ancestor internal nodes. Then for each sample, look up its leaf
  # and add its cluster to all ancestors.

  # Initialize accumulators per cluster
  imd_x <- lapply(labs, function(lb) setNames(numeric(px), x_names))
  names(imd_x) <- labs
  imd_x_cnt <- lapply(labs, function(lb) setNames(integer(px), x_names))
  names(imd_x_cnt) <- labs

  imd_y <- if (qy > 0) {
    lapply(labs, function(lb) setNames(numeric(qy), y_names))
  } else {
    lapply(labs, function(lb) numeric(0))
  }
  names(imd_y) <- labs
  imd_y_cnt <- if (qy > 0) {
    lapply(labs, function(lb) setNames(integer(qy), y_names))
  } else {
    lapply(labs, function(lb) integer(0))
  }
  names(imd_y_cnt) <- labs

  # Total samples per cluster (for normalization)
  n_per_cluster <- table(cluster)

  for (t in seq_len(ntree)) {
    ti <- tree_info[[t]]
    n_nodes <- length(ti$split_var)
    split_var <- ti$split_var     # -1 for leaf, 0-indexed X column for internal
    left_child <- ti$left         # -1 if no child
    right_child <- ti$right
    imd_x_score <- ti$imd_x_score  # per-node X split score
    imd_y_mat <- ti$imd_y_stats    # qy x n_internal matrix (or NULL)

    # Build leaf -> ancestor map via DFS
    # For each node, record the path of internal ancestors
    # ancestors[[node_1indexed]] = vector of 1-indexed internal node indices
    ancestors <- vector("list", n_nodes)
    ancestors[[1]] <- integer(0)  # root has no ancestors above it

    # BFS to propagate ancestor lists
    queue <- 1L
    while (length(queue) > 0) {
      node <- queue[1]
      queue <- queue[-1]

      # Current node's ancestor list: all ancestors of this node + this node if internal
      my_ancestors <- ancestors[[node]]
      if (split_var[node] >= 0) {
        # This is an internal node: children inherit my ancestors + me
        child_ancestors <- c(my_ancestors, node)
      } else {
        child_ancestors <- my_ancestors
      }

      l <- left_child[node]
      r <- right_child[node]
      if (l >= 0) {
        ancestors[[l + 1L]] <- child_ancestors  # +1 for 0-index to 1-index
        queue <- c(queue, l + 1L)
      }
      if (r >= 0) {
        ancestors[[r + 1L]] <- child_ancestors
        queue <- c(queue, r + 1L)
      }
    }

    # Get leaf assignments for all samples in this tree
    mem_t <- membership[, t]  # 0-indexed node IDs

    # Build internal node index mapping for imd_y_mat columns
    internal_col <- integer(n_nodes)
    col_idx <- 0L
    for (ni in seq_len(n_nodes)) {
      if (split_var[ni] >= 0) {
        col_idx <- col_idx + 1L
        internal_col[ni] <- col_idx
      }
    }

    # Count samples per cluster at each internal node
    # node_cluster_count[[node_1indexed]][cluster_label] = count
    # For efficiency, accumulate directly into IMD
    for (i in seq_len(n)) {
      cl <- cluster[i]
      if (is.na(cl)) next
      leaf_node <- mem_t[i] + 1L  # 0-indexed to 1-indexed
      if (leaf_node < 1 || leaf_node > n_nodes) next
      anc <- ancestors[[leaf_node]]
      if (length(anc) == 0) next

      for (v in anc) {
        xv <- split_var[v]  # 0-indexed X column
        if (xv < 0) next
        xv1 <- xv + 1L  # 1-indexed

        # Weight = 1 / total nodesize at this node (so cluster fraction emerges)
        # Actually: accumulate raw, normalize later
        score <- imd_x_score[v]
        imd_x[[cl]][xv1] <- imd_x[[cl]][xv1] + score
        imd_x_cnt[[cl]][xv1] <- imd_x_cnt[[cl]][xv1] + 1L

        # Y side
        if (qy > 0 && !is.null(imd_y_mat) && nrow(imd_y_mat) == qy) {
          col <- internal_col[v]
          if (col > 0 && col <= ncol(imd_y_mat)) {
            for (j in seq_len(qy)) {
              ys <- imd_y_mat[j, col]
              if (ys > 0) {
                imd_y[[cl]][j] <- imd_y[[cl]][j] + ys
                imd_y_cnt[[cl]][j] <- imd_y_cnt[[cl]][j] + 1L
              }
            }
          }
        }
      }
    }
  }

  # Normalize: divide by count, then optionally L2-normalize
  out <- lapply(labs, function(lb) {
    wx <- imd_x[[lb]]
    cx <- imd_x_cnt[[lb]]
    wx <- ifelse(cx > 0, wx / cx, 0)

    wy <- if (qy > 0) {
      cy <- imd_y_cnt[[lb]]
      ifelse(cy > 0, imd_y[[lb]] / cy, 0)
    } else {
      numeric(0)
    }

    if (normalized) {
      norm_x <- sqrt(sum(wx^2))
      if (is.finite(norm_x) && norm_x > 0) wx <- wx / norm_x
      if (qy > 0) {
        norm_y <- sqrt(sum(wy^2))
        if (is.finite(norm_y) && norm_y > 0) wy <- wy / norm_y
      }
    }

    list(X = wx, Y = wy)
  })
  names(out) <- labs
  out
}
