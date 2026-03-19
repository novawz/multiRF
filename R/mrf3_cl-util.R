## Utility functions for mrf3_cl.R

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) {
    return(0)
  }
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(0)
  }
  out <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  if (!is.finite(out)) 0 else out
}

build_sample_embedding <- function(dat, leaf_embed_dim = 10L) {
  dat <- as.matrix(dat)
  if (is.null(rownames(dat))) {
    rownames(dat) <- as.character(seq_len(nrow(dat)))
  }

  k <- max(
    1L,
    min(
      as.integer(leaf_embed_dim),
      as.integer(ncol(dat)),
      max(1L, as.integer(nrow(dat) - 1L))
    )
  )

  emb <- tryCatch(
    {
      pc <- stats::prcomp(dat, center = TRUE, scale. = TRUE, rank. = k)
      pc$x[, seq_len(min(k, ncol(pc$x))), drop = FALSE]
    },
    error = function(e) NULL
  )

  if (is.null(emb) || !is.matrix(emb) || ncol(emb) == 0L) {
    emb <- scale(dat)
    if (is.null(dim(emb))) {
      emb <- matrix(emb, ncol = 1)
    }
    emb <- emb[, seq_len(min(ncol(emb), k)), drop = FALSE]
  }

  rownames(emb) <- rownames(dat)
  emb
}

build_embedding_list <- function(mod, symm = TRUE, use = "X", leaf_embed_dim = 10L) {
  if (symm) {
    sample_embed_list <- list(
      X = build_sample_embedding(mod$xvar, leaf_embed_dim = leaf_embed_dim)
    )
    if (!is.null(mod$yvar)) {
      sample_embed_list$Y <- build_sample_embedding(mod$yvar, leaf_embed_dim = leaf_embed_dim)
    }
    return(sample_embed_list)
  }

  dat_use <- if (use == "Y" && !is.null(mod$yvar)) mod$yvar else mod$xvar
  nm <- if (use == "Y" && !is.null(mod$yvar)) "Y" else "X"
  out <- list()
  out[[nm]] <- build_sample_embedding(dat_use, leaf_embed_dim = leaf_embed_dim)
  out
}

get_leaf_centroid_embedding <- function(sample_embed, tree.membership) {
  ids <- sort(unique(tree.membership))
  out <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = ncol(sample_embed),
    dimnames = list(as.character(ids), colnames(sample_embed))
  )

  for (i in seq_along(ids)) {
    idx <- which(tree.membership == ids[i])
    out[i, ] <- colMeans(sample_embed[idx, , drop = FALSE], na.rm = TRUE)
  }

  out
}

# Get Spearman correlation of one sibling-leaf pair based on low-dimensional embeddings.
get_corr <- function(mem_pair, centroid_list) {
  corr_vals <- vapply(
    centroid_list,
    FUN.VALUE = numeric(1),
    FUN = function(centroid) {
      id1 <- as.character(mem_pair[1])
      id2 <- as.character(mem_pair[2])
      if (!(id1 %in% rownames(centroid)) || !(id2 %in% rownames(centroid))) {
        return(NA_real_)
      }
      safe_spearman(
        as.numeric(centroid[id1, , drop = TRUE]),
        as.numeric(centroid[id2, , drop = TRUE])
      )
    }
  )

  corr_vals <- corr_vals[is.finite(corr_vals)]
  if (length(corr_vals) == 0L) {
    return(NA_real_)
  }
  mean(corr_vals)
}

# Compute correlations for sibling leaf candidates only.
get_leaf_corr <- function(tree.membership, net, sample_embed_list) {
  if (is.null(net$corr)) {
    net$corr <- 0
  }
  if (is.null(net$is_leaf)) {
    net$is_leaf <- ifelse(grepl("leaf", net$to), 1, 0)
  }
  if (is.null(net$leaf_used)) {
    net$leaf_used <- 0
  }
  if (is.null(net$mem_id)) {
    net$mem_id <- ifelse(net$is_leaf == 1, net$to_id, 0)
  }

  sibling_parent <- dplyr::filter(.data = net, is_leaf == 1)
  sibling_parent <- dplyr::summarise(.data = sibling_parent, n = n(), .by = from)
  sibling_parent <- dplyr::filter(.data = sibling_parent, n == 2)

  if (nrow(sibling_parent) == 0L) {
    return(list(
      net = net,
      leaf_select = net[0, , drop = FALSE],
      parent_tbl = data.frame(
        parent = character(0),
        corr = numeric(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  leaf_select <- dplyr::filter(net, from %in% sibling_parent$from)
  centroid_list <- purrr::map(
    sample_embed_list,
    ~get_leaf_centroid_embedding(.x, tree.membership = tree.membership)
  )

  parent_ids <- unique(leaf_select$from)
  corr_leaf <- vapply(
    parent_ids,
    FUN.VALUE = numeric(1),
    FUN = function(pid) {
      mem_pair <- unique(leaf_select$mem_id[leaf_select$from == pid])
      if (length(mem_pair) != 2L) {
        return(NA_real_)
      }
      get_corr(mem_pair = mem_pair, centroid_list = centroid_list)
    }
  )

  net$leaf_used[match(leaf_select$to, net$to)] <- 1
  parent_idx <- match(parent_ids, net$to)
  keep <- which(!is.na(parent_idx))
  net$corr[parent_idx[keep]] <- corr_leaf[keep]

  parent_tbl <- data.frame(
    parent = as.character(parent_ids),
    corr = as.numeric(corr_leaf),
    stringsAsFactors = FALSE
  )

  list(
    net = net,
    leaf_select = leaf_select,
    parent_tbl = parent_tbl
  )
}

transform_sibling_corr <- function(corr, sibling_fun = "positive") {
  if (!is.finite(corr)) {
    return(0)
  }
  switch(
    sibling_fun,
    constant = 1,
    positive = max(corr, 0),
    shift01 = (corr + 1) / 2,
    abs = abs(corr),
    identity = corr,
    max(corr, 0)
  )
}

build_sibling_edges <- function(leaf_select,
                                parent_tbl,
                                sibling_gamma = 0.5,
                                sibling_fun = "positive",
                                sibling_cap = TRUE) {
  if (nrow(parent_tbl) == 0L || nrow(leaf_select) == 0L) {
    return(data.frame(
      leaf1 = character(0),
      leaf2 = character(0),
      weight = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  edge_ls <- lapply(seq_len(nrow(parent_tbl)), function(i) {
    pid <- parent_tbl$parent[i]
    corr <- parent_tbl$corr[i]
    pair_rows <- leaf_select[leaf_select$from == pid, , drop = FALSE]
    pair_mem <- unique(pair_rows$mem_id)
    if (length(pair_mem) != 2L) {
      return(NULL)
    }

    w <- sibling_gamma * transform_sibling_corr(corr, sibling_fun = sibling_fun)
    if (!is.finite(w)) {
      return(NULL)
    }
    w <- max(0, w)
    if (isTRUE(sibling_cap)) {
      w <- min(1, w)
    }
    if (w <= 0) {
      return(NULL)
    }

    data.frame(
      leaf1 = as.character(pair_mem[1]),
      leaf2 = as.character(pair_mem[2]),
      weight = w,
      stringsAsFactors = FALSE
    )
  })

  edge_ls <- Filter(Negate(is.null), edge_ls)
  if (length(edge_ls) == 0L) {
    return(data.frame(
      leaf1 = character(0),
      leaf2 = character(0),
      weight = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, edge_ls)
}

build_sparse_leaf_prox <- function(tree.membership, sibling_edges = NULL) {
  if (is.null(names(tree.membership))) {
    names(tree.membership) <- as.character(seq_along(tree.membership))
  }

  n <- length(tree.membership)
  idx_by_leaf <- split(seq_len(n), as.character(tree.membership))

  i_all <- integer(0)
  j_all <- integer(0)
  x_all <- numeric(0)

  # Same enhanced leaf edges (weight = 1).
  for (idx in idx_by_leaf) {
    m <- length(idx)
    if (m == 0L) {
      next
    }
    i_all <- c(i_all, rep(idx, each = m))
    j_all <- c(j_all, rep(idx, times = m))
    x_all <- c(x_all, rep.int(1, m * m))
  }

  # Neighbor enhanced leaf edges (weight from sibling edges).
  if (!is.null(sibling_edges) && nrow(sibling_edges) > 0L) {
    for (k in seq_len(nrow(sibling_edges))) {
      leaf1 <- sibling_edges$leaf1[k]
      leaf2 <- sibling_edges$leaf2[k]
      w <- sibling_edges$weight[k]
      idx1 <- idx_by_leaf[[as.character(leaf1)]]
      idx2 <- idx_by_leaf[[as.character(leaf2)]]
      if (is.null(idx1) || is.null(idx2) || length(idx1) == 0L || length(idx2) == 0L) {
        next
      }
      i_all <- c(i_all, rep(idx1, each = length(idx2)), rep(idx2, each = length(idx1)))
      j_all <- c(j_all, rep(idx2, times = length(idx1)), rep(idx1, times = length(idx2)))
      x_all <- c(x_all, rep(w, 2L * length(idx1) * length(idx2)))
    }
  }

  prox <- Matrix::sparseMatrix(
    i = i_all,
    j = j_all,
    x = x_all,
    dims = c(n, n),
    dimnames = list(names(tree.membership), names(tree.membership)),
    giveCsparse = TRUE
  )
  if (length(prox@x) > 0L) {
    prox@x <- pmin(prox@x, 1)
  }
  prox
}

build_soft_tree_prox <- function(tree.membership,
                                 leaf_select,
                                 parent_tbl,
                                 sibling_gamma = 0.5,
                                 sibling_fun = "positive",
                                 sibling_cap = TRUE) {
  sibling_edges <- build_sibling_edges(
    leaf_select = leaf_select,
    parent_tbl = parent_tbl,
    sibling_gamma = sibling_gamma,
    sibling_fun = sibling_fun,
    sibling_cap = sibling_cap
  )
  build_sparse_leaf_prox(tree.membership = tree.membership, sibling_edges = sibling_edges)
}

# Merge sibling pairs conservatively: keep only candidates that satisfy
# the old "mem_drop" style condition before applying optional quantile pruning.
update_tree_leaf <- function(tree.membership,
                             net,
                             sample_embed_list,
                             merge_quantile = 0.9,
                             size_min = 10) {
  update <- get_leaf_corr(
    tree.membership = tree.membership,
    net = net,
    sample_embed_list = sample_embed_list
  )

  updated_net <- update$net
  parent_tbl <- update$parent_tbl

  if (nrow(parent_tbl) == 0L) {
    return(list(net = updated_net, mem = tree.membership, merged = FALSE))
  }

  parent_tbl <- parent_tbl[is.finite(parent_tbl$corr), , drop = FALSE]
  if (nrow(parent_tbl) == 0L) {
    return(list(net = updated_net, mem = tree.membership, merged = FALSE))
  }

  leaf_select <- update$leaf_select
  merge_tbl <- lapply(seq_len(nrow(parent_tbl)), function(i) {
    pid <- as.character(parent_tbl$parent[i])
    new_corr <- as.numeric(parent_tbl$corr[i])
    pair_rows <- leaf_select[leaf_select$from == pid, , drop = FALSE]
    pair_mem <- unique(pair_rows$mem_id)
    if (length(pair_mem) != 2L) {
      return(NULL)
    }

    child_corr <- pair_rows$corr
    child_corr <- child_corr[is.finite(child_corr)]
    corr_max <- if (length(child_corr) > 0L) max(child_corr) else 0

    node_sz <- pair_rows$nodesize
    node_sz <- node_sz[is.finite(node_sz)]
    nodesize_min <- if (length(node_sz) > 0L) min(node_sz) else 0

    mem_drop <- ((new_corr < corr_max) || (corr_max < 0)) && (nodesize_min > size_min)
    if (!isTRUE(mem_drop)) {
      return(NULL)
    }

    data.frame(
      parent = pid,
      corr = new_corr,
      corr_max = corr_max,
      nodesize_min = nodesize_min,
      stringsAsFactors = FALSE
    )
  })
  merge_tbl <- Filter(Negate(is.null), merge_tbl)
  if (length(merge_tbl) == 0L) {
    return(list(net = updated_net, mem = tree.membership, merged = FALSE))
  }
  merge_tbl <- do.call(rbind, merge_tbl)

  # Optional pruning among valid mem_drop candidates.
  merge_quantile <- min(max(as.numeric(merge_quantile), 0), 1)
  if (nrow(merge_tbl) > 1L) {
    corr_cut <- stats::quantile(
      merge_tbl$corr,
      probs = merge_quantile,
      na.rm = TRUE,
      names = FALSE,
      type = 8
    )
    merge_tbl <- merge_tbl[merge_tbl$corr >= corr_cut, , drop = FALSE]
  }
  if (nrow(merge_tbl) == 0L) {
    return(list(net = updated_net, mem = tree.membership, merged = FALSE))
  }

  mem_new <- tree.membership

  for (pid in merge_tbl$parent) {
    pair_rows <- leaf_select[leaf_select$from == pid, , drop = FALSE]
    pair_mem <- unique(pair_rows$mem_id)
    if (length(pair_mem) != 2L) {
      next
    }

    keep_id <- max(pair_mem)
    drop_id <- min(pair_mem)

    if (!any(mem_new == drop_id)) {
      next
    }

    mem_new[mem_new == drop_id] <- keep_id
    updated_net$mem_id[updated_net$mem_id == drop_id] <- keep_id

    p_idx <- which(updated_net$to == pid)
    if (length(p_idx) > 0L) {
      updated_net$is_leaf[p_idx] <- 1
      updated_net$mem_id[p_idx] <- keep_id
    }

    child_idx <- match(pair_rows$to, updated_net$to)
    child_idx <- child_idx[!is.na(child_idx)]
    if (length(child_idx) > 0L) {
      updated_net$is_leaf[child_idx] <- 0
    }
  }

  merged <- !identical(unname(mem_new), unname(tree.membership))
  list(net = updated_net, mem = mem_new, merged = merged)
}

# Update terminal nodes from bottom to top
update_iter_cl <- function(mod,
                           tree.id,
                           size_min = 10,
                           symm = TRUE,
                           use = "X",
                           leaf_embed_dim = 10,
                           sample_embed_list = NULL,
                           merge_quantile = 0.9,
                           merge_mode = c("soft", "hard"),
                           sibling_gamma = 0.5,
                           sibling_fun = c("constant", "positive", "shift01", "abs", "identity"),
                           sibling_cap = TRUE,
                           hard_prox_mode = c("soft_enhanced", "binary")) {
  merge_mode <- match.arg(merge_mode)
  sibling_fun <- match.arg(sibling_fun)
  hard_prox_mode <- match.arg(hard_prox_mode)

  net <- get_tree_net(mod = mod, tree.id = tree.id)
  mem <- mod$membership[, tree.id]
  names(mem) <- rownames(mod$xvar)
  if (is.null(sample_embed_list)) {
    sample_embed_list <- build_embedding_list(
      mod = mod,
      symm = symm,
      use = use,
      leaf_embed_dim = leaf_embed_dim
    )
  }

  if (merge_mode == "soft") {
    if (is.null(net$is_leaf)) {
      net$is_leaf <- ifelse(grepl("leaf", net$to), 1, 0)
    }
    if (is.null(net$mem_id)) {
      net$mem_id <- ifelse(net$is_leaf == 1, net$to_id, 0)
    }

    leaf_info <- get_leaf_corr(
      tree.membership = mem,
      net = net,
      sample_embed_list = sample_embed_list
    )
    prox <- build_soft_tree_prox(
      tree.membership = mem,
      leaf_select = leaf_info$leaf_select,
      parent_tbl = leaf_info$parent_tbl,
      sibling_gamma = sibling_gamma,
      sibling_fun = sibling_fun,
      sibling_cap = sibling_cap
    )

    class_d <- mem[match(rownames(mod$xvar), names(mem))]
    if (anyNA(class_d)) {
      class_d <- as.numeric(mem)
      names(class_d) <- rownames(mod$xvar)
    }
    return(list(class_mem = class_d, prox = prox))
  }

  max_iter <- nrow(net) + 5L
  for (iter in seq_len(max_iter)) {
    update_ls <- update_tree_leaf(
      tree.membership = mem,
      net = net,
      sample_embed_list = sample_embed_list,
      merge_quantile = merge_quantile,
      size_min = size_min
    )

    net <- update_ls$net
    if (!isTRUE(update_ls$merged)) {
      break
    }
    mem <- update_ls$mem
  }

  class_d <- mem[match(rownames(mod$xvar), names(mem))]
  if (anyNA(class_d)) {
    class_d <- as.numeric(mem)
    names(class_d) <- rownames(mod$xvar)
  }

  if (hard_prox_mode == "binary") {
    prox <- build_sparse_leaf_prox(tree.membership = class_d, sibling_edges = NULL)
  } else {
    leaf_info <- get_leaf_corr(
      tree.membership = mem,
      net = net,
      sample_embed_list = sample_embed_list
    )
    prox <- build_soft_tree_prox(
      tree.membership = mem,
      leaf_select = leaf_info$leaf_select,
      parent_tbl = leaf_info$parent_tbl,
      sibling_gamma = sibling_gamma,
      sibling_fun = sibling_fun,
      sibling_cap = sibling_cap
    )
  }
  list(class_mem = class_d, prox = prox)
}

# Create new membership and new proximity matrix of forest
cl_forest <- function(mod,
                      size_min = 5,
                      use = "X",
                      symm = FALSE,
                      leaf_embed_dim = 10,
                      merge_quantile = 0.9,
                      merge_mode = c("soft", "hard"),
                      sibling_gamma = 0.5,
                      sibling_fun = c("constant", "positive", "shift01", "abs", "identity"),
                      sibling_cap = TRUE,
                      hard_prox_mode = c("soft_enhanced", "binary"),
                      parallel = TRUE,
                      cores = NULL,
                      sample_embed_list = NULL,
                      ...) {
  merge_mode <- match.arg(merge_mode)
  sibling_fun <- match.arg(sibling_fun)
  hard_prox_mode <- match.arg(hard_prox_mode)
  nt <- mod$ntree

  if (parallel) {
    cores <- sanitize_mc_cores(cores = cores, fallback = 1L)
    if (Sys.info()["sysname"] == "Windows") {
      cluster <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cluster)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
    } else {
      doParallel::registerDoParallel(cores)
    }
  }
  `%myinfix%` <- ifelse(parallel, `%dopar%`, `%do%`)
  if (is.null(sample_embed_list)) {
    sample_embed_list <- build_embedding_list(
      mod = mod,
      symm = symm,
      use = use,
      leaf_embed_dim = leaf_embed_dim
    )
  }

  combine_tree_stats <- function(lhs, rhs) {
    lhs$prox <- lhs$prox + rhs$prox
    lhs$n_trees <- lhs$n_trees + rhs$n_trees
    lhs
  }

  forest_stat <- foreach(
    t = seq_len(nt),
    .errorhandling = "remove",
    .combine = combine_tree_stats
  ) %myinfix% {
    one_tree <- update_iter_cl(
      mod,
      tree.id = t,
      size_min = size_min,
      use = use,
      symm = symm,
      leaf_embed_dim = leaf_embed_dim,
      sample_embed_list = sample_embed_list,
      merge_quantile = merge_quantile,
      merge_mode = merge_mode,
      sibling_gamma = sibling_gamma,
      sibling_fun = sibling_fun,
      sibling_cap = sibling_cap,
      hard_prox_mode = hard_prox_mode
    )
    list(
      prox = one_tree$prox,
      n_trees = 1L
    )
  }

  if (is.null(forest_stat) || is.null(forest_stat$n_trees) || forest_stat$n_trees < 1L) {
    stop("No valid trees were processed in `cl_forest()`.")
  }

  nt_eff <- as.integer(forest_stat$n_trees)
  prox <- forest_stat$prox / nt_eff

  list(
    prox = prox
  )
}
