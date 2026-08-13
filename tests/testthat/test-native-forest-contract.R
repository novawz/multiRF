native_supervised <- function(...) {
  env <- if (exists("fit_mv_forest_cpp", envir = .GlobalEnv, inherits = FALSE)) {
    .GlobalEnv
  } else {
    asNamespace("multiRF")
  }
  get("fit_mv_forest_cpp", envir = env)(...)
}

native_unsupervised <- function(...) {
  env <- if (exists("fit_mv_forest_unsup_cpp", envir = .GlobalEnv, inherits = FALSE)) {
    .GlobalEnv
  } else {
    asNamespace("multiRF")
  }
  get("fit_mv_forest_unsup_cpp", envir = env)(...)
}

reconstruct_forest_weights <- function(membership, inbag, mode) {
  n <- nrow(membership)
  ntree <- ncol(membership)
  total <- matrix(0, n, n)
  denom <- numeric(n)
  for (tt in seq_len(ntree)) {
    leaf <- membership[, tt]
    freq <- inbag[, tt]
    for (ii in seq_len(n)) {
      if (identical(mode, "oob") && freq[ii] > 0L) next
      donors <- which(leaf == leaf[ii])
      donor_weight <- if (identical(mode, "all")) {
        rep(1, length(donors))
      } else {
        freq[donors]
      }
      keep <- donor_weight > 0
      donors <- donors[keep]
      donor_weight <- donor_weight[keep]
      if (!length(donors)) next
      total[ii, donors] <- total[ii, donors] + donor_weight / sum(donor_weight)
      denom[ii] <- denom[ii] + 1
    }
  }
  out <- matrix(NA_real_, n, n)
  keep <- denom > 0
  out[keep, ] <- total[keep, , drop = FALSE] / denom[keep]
  out
}

expected_imd_x <- function(tree_info, p) {
  out <- numeric(p)
  internal <- which(!tree_info$is_leaf)
  if (!length(internal)) return(out)
  for (jj in seq_len(p)) {
    at <- internal[tree_info$split_var[internal] == jj - 1L]
    if (length(at)) out[jj] <- 1 / (min(tree_info$depth[at]) + 1)
  }
  out
}

reconstruct_grouped_enhanced_proximity <- function(
    membership, tree_info, embed, embed_split, sibling_gamma) {
  n <- nrow(membership)
  ntree <- ncol(membership)
  out <- matrix(0, n, n)

  safe_spearman <- function(a, b) {
    if (length(a) < 2L || stats::sd(a) == 0 || stats::sd(b) == 0) return(0)
    value <- suppressWarnings(stats::cor(a, b, method = "spearman"))
    if (is.finite(value)) value else 0
  }

  for (tt in seq_len(ntree)) {
    leaf <- membership[, tt]
    for (node in unique(leaf)) {
      idx <- which(leaf == node)
      out[idx, idx] <- out[idx, idx] + 1
    }

    ti <- tree_info[[tt]]
    internal <- which(!ti$is_leaf)
    for (node in internal) {
      left <- ti$left[node]
      right <- ti$right[node]
      if (left < 0L || right < 0L ||
          !ti$is_leaf[left + 1L] || !ti$is_leaf[right + 1L]) next
      idx_left <- which(leaf == left)
      idx_right <- which(leaf == right)
      if (!length(idx_left) || !length(idx_right)) next

      cent_left <- colMeans(embed[idx_left, , drop = FALSE])
      cent_right <- colMeans(embed[idx_right, , drop = FALSE])
      corr_x <- safe_spearman(
        cent_left[seq_len(embed_split)], cent_right[seq_len(embed_split)]
      )
      y_idx <- seq.int(embed_split + 1L, ncol(embed))
      corr_y <- safe_spearman(cent_left[y_idx], cent_right[y_idx])
      weight <- min(1, sibling_gamma * max((corr_x + corr_y) / 2, 0))
      if (weight > 0) {
        out[idx_left, idx_right] <- out[idx_left, idx_right] + weight
        out[idx_right, idx_left] <- out[idx_right, idx_left] + weight
      }
    }
  }
  out / ntree
}

test_that("SWR split statistics retain bootstrap multiplicity", {
  set.seed(2801)
  n <- 48L
  X <- matrix(rnorm(n * 4L), n, 4L)
  Y <- cbind(
    1.7 * X[, 1L] + rnorm(n, sd = 0.2),
    -1.2 * X[, 2L] + rnorm(n, sd = 0.3),
    0.8 * X[, 3L] + rnorm(n, sd = 0.4)
  )

  fit <- native_supervised(
    X, Y, ntree = 1L, mtry = ncol(X), ytry = ncol(Y), nsplit = 0L,
    nodesize_min = 2L, max_depth = 1L, seed = 77L, nthread = 1L,
    samptype = 1L, prox_mode = -1L
  )
  freq <- fit$inbag[, 1L]
  expect_equal(sum(freq), n)
  expect_equal(fit$tree_info[[1L]]$nodesize[1L], n)
  expect_true(any(freq > 1L))

  total <- sum(freq)
  y_mean <- colSums(Y * freq) / total
  y_centered <- sweep(Y, 2L, y_mean, "-")
  y_sd <- sqrt(colSums(y_centered^2 * freq) / total)
  y_std <- sweep(y_centered, 2L, y_sd, "/")

  candidates <- list()
  for (xj in seq_len(ncol(X))) {
    observed <- sort(unique(X[freq > 0L, xj]))
    for (cut in head(observed, -1L)) {
      left_freq <- freq * (X[, xj] <= cut)
      n_left <- sum(left_freq)
      n_right <- total - n_left
      s_left <- colSums(y_std * left_freq)
      component <- s_left^2 / n_left + s_left^2 / n_right
      candidates[[length(candidates) + 1L]] <- list(
        x = xj - 1L, cut = cut, score = mean(component), component = component
      )
    }
  }
  score <- vapply(candidates, `[[`, numeric(1L), "score")
  best <- candidates[[which.max(score)]]
  root <- fit$tree_info[[1L]]

  expect_equal(root$imd_x_score[1L], best$score, tolerance = 1e-10)
  expect_equal(root$split_var[1L], best$x)
  expect_equal(root$split_val[1L], best$cut, tolerance = 1e-12)
  expect_equal(root$imd_y_stats[, 1L], best$component, tolerance = 1e-10)
  expect_equal(root$msrv[1L], which.max(best$component) - 1L)
})

test_that("all, inbag, and oob forest-weight modes implement their definitions", {
  set.seed(2802)
  X <- matrix(rnorm(36L * 4L), 36L, 4L)
  Y <- matrix(rnorm(36L * 3L), 36L, 3L)

  fits <- lapply(0:2, function(mode) {
    native_supervised(
      X, Y, ntree = 12L, mtry = 4L, ytry = 3L, nsplit = 4L,
      nodesize_min = 3L, max_depth = 3L, seed = 529L, nthread = 1L,
      samptype = 1L, prox_mode = -1L, forest_wt_mode = mode
    )
  })
  expect_identical(fits[[1L]]$membership, fits[[2L]]$membership)
  expect_identical(fits[[1L]]$inbag, fits[[3L]]$inbag)

  for (kk in seq_along(fits)) {
    mode <- c("all", "inbag", "oob")[[kk]]
    expected <- reconstruct_forest_weights(
      fits[[kk]]$membership, fits[[kk]]$inbag, mode
    )
    expect_equal(fits[[kk]]$forest.wt, expected, tolerance = 1e-12)
  }
  expect_gt(max(abs(fits[[1L]]$forest.wt - fits[[2L]]$forest.wt)), 0)

  unsup <- lapply(0:2, function(mode) {
    native_unsupervised(
      X, ntree = 12L, ytry = 3L, nodesize_min = 3L, max_depth = 3L,
      seed = 529L, nthread = 1L, samptype = 1L, prox_mode = -1L,
      forest_wt_mode = mode, nsplit = 4L
    )
  })
  expect_identical(unsup[[1L]]$membership, unsup[[2L]]$membership)
  expect_identical(unsup[[1L]]$inbag, unsup[[3L]]$inbag)
  for (kk in seq_along(unsup)) {
    expected <- reconstruct_forest_weights(
      unsup[[kk]]$membership, unsup[[kk]]$inbag,
      c("all", "inbag", "oob")[[kk]]
    )
    expect_equal(unsup[[kk]]$forest.wt, expected, tolerance = 1e-12)
  }
})

test_that("native IMD is inverse minimal depth with one MSRV per split", {
  set.seed(2803)
  X <- matrix(rnorm(45L * 5L), 45L, 5L)
  Y <- cbind(X[, 1L] + rnorm(45L, sd = 0.1),
             X[, 2L] + rnorm(45L, sd = 0.2),
             matrix(rnorm(45L * 2L), 45L, 2L))
  fit <- native_supervised(
    X, Y, ntree = 9L, mtry = 5L, ytry = 4L, nsplit = 0L,
    nodesize_min = 3L, max_depth = 3L, seed = 529L, nthread = 1L,
    samptype = 0L, prox_mode = -1L
  )

  for (tt in seq_len(fit$ntree)) {
    ti <- fit$tree_info[[tt]]
    expect_equal(fit$imd_x_per_tree[, tt], expected_imd_x(ti, ncol(X)))
    internal <- which(!ti$is_leaf)
    expect_true(all(!is.na(ti$msrv[internal])))
    expect_true(all(ti$msrv[internal] >= 0L & ti$msrv[internal] < ncol(Y)))
    expected_y <- numeric(ncol(Y))
    for (jj in seq_len(ncol(Y))) {
      at <- internal[ti$msrv[internal] == jj - 1L]
      if (length(at)) expected_y[jj] <- 1 / (min(ti$depth[at]) + 1)
    }
    expect_equal(fit$imd_y_per_tree[, tt], expected_y)
  }
  expect_equal(fit$imd_x, rowMeans(fit$imd_x_per_tree))
  expect_equal(fit$imd_y, rowMeans(fit$imd_y_per_tree))
  expect_true(all(fit$imd_x >= 0 & fit$imd_x <= 1))
  expect_true(all(fit$imd_y >= 0 & fit$imd_y <= 1))

  unsup <- native_unsupervised(
    X, ntree = 9L, ytry = 3L, nodesize_min = 3L, max_depth = 3L,
    seed = 529L, nthread = 1L, samptype = 1L, prox_mode = -1L,
    nsplit = 0L
  )
  for (tt in seq_len(unsup$ntree)) {
    ti <- unsup$tree_info[[tt]]
    expect_equal(unsup$imd_x_per_tree[, tt], expected_imd_x(ti, ncol(X)))
    internal <- which(!ti$is_leaf)
    if (length(internal)) {
      expect_equal(
        ti$nodesize[ti$left[internal] + 1L] + ti$nodesize[ti$right[internal] + 1L],
        ti$nodesize[internal]
      )
    }
  }
  expect_equal(unsup$imd_x, rowMeans(unsup$imd_x_per_tree))
  expect_true(all(unsup$imd_x >= 0 & unsup$imd_x <= 1))
})

test_that("weighted candidate sampling excludes zero-weight variables", {
  set.seed(2804)
  X <- matrix(rnorm(40L * 4L), 40L, 4L)
  Y <- matrix(rnorm(40L * 3L), 40L, 3L)
  fit <- native_supervised(
    X, Y, ntree = 8L, mtry = 4L, ytry = 3L, nsplit = 0L,
    nodesize_min = 2L, max_depth = 1L, seed = 529L, nthread = 1L,
    samptype = 0L, prox_mode = -1L,
    xvar_wt = c(0, 1, 0, 0), yvar_wt = c(0, 0, 1)
  )
  expect_true(all(vapply(fit$tree_info, function(ti) ti$split_var[1L] == 1L,
                         logical(1L))))
  expect_true(all(vapply(fit$tree_info, function(ti) ti$msrv[1L] == 2L,
                         logical(1L))))
})

test_that("enhanced proximity averages separate X and Y embedding correlations", {
  set.seed(2806)
  n <- 42L
  X <- matrix(rnorm(n * 5L), n, 5L)
  Y <- cbind(
    X[, 1L] - 0.6 * X[, 2L] + rnorm(n, sd = 0.3),
    X[, 3L] + rnorm(n, sd = 0.4),
    matrix(rnorm(n * 2L), n, 2L)
  )
  embed_x <- cbind(X[, 1L], X[, 2L]^2, X[, 4L] - X[, 5L])
  embed_y <- cbind(Y[, 1L], -Y[, 2L], Y[, 3L] + Y[, 4L])
  embed <- cbind(embed_x, embed_y)
  gamma <- 0.65

  fit <- native_supervised(
    X, Y, ntree = 11L, mtry = 5L, ytry = 4L, nsplit = 4L,
    nodesize_min = 2L, max_depth = 4L, seed = 529L, nthread = 1L,
    samptype = 0L, prox_mode = -1L, embed = embed,
    embed_split = ncol(embed_x), sibling_gamma = gamma,
    enhanced_prox_mode = 1L
  )
  expected <- reconstruct_grouped_enhanced_proximity(
    fit$membership, fit$tree_info, embed, ncol(embed_x), gamma
  )
  expect_equal(fit$enhanced_prox, expected, tolerance = 1e-12)
})

test_that("node standardization is invariant to nonzero response rescaling", {
  set.seed(2805)
  n <- 50L
  X <- matrix(rnorm(n * 5L), n, 5L)
  Y <- cbind(
    1.8 * X[, 1L] + rnorm(n, sd = 0.25),
    -1.1 * X[, 2L] + rnorm(n, sd = 0.35),
    rnorm(n)
  )
  args <- list(
    X = X, ntree = 31L, mtry = 5L, ytry = 3L, nsplit = 0L,
    nodesize_min = 3L, max_depth = 3L, seed = 529L, nthread = 1L,
    samptype = 0L, prox_mode = -1L
  )
  ordinary <- do.call(native_supervised, c(args, list(Y = Y)))
  small <- do.call(native_supervised, c(args, list(Y = Y * 1e-5)))

  expect_identical(ordinary$membership, small$membership)
  expect_identical(ordinary$inbag, small$inbag)
  expect_identical(ordinary$forest.wt, small$forest.wt)
  expect_identical(
    lapply(ordinary$tree_info, function(x) x[c("split_var", "split_val", "left", "right")]),
    lapply(small$tree_info, function(x) x[c("split_var", "split_val", "left", "right")])
  )
})

test_that("native entry points reject invalid dimensions, values, and controls", {
  X <- matrix(seq_len(24L), nrow = 8L, ncol = 3L)
  Y <- matrix(seq_len(16L), nrow = 8L, ncol = 2L)

  expect_error(native_supervised(X, Y, ntree = 0L), "ntree must be a positive")
  expect_error(native_unsupervised(X, ntree = 0L), "ntree must be a positive")
  expect_error(native_supervised(X, Y[-1L, , drop = FALSE]), "nrow\\(X\\)")

  X_bad <- X
  X_bad[1L, 1L] <- Inf
  expect_error(native_supervised(X_bad, Y), "X must contain only finite")
  expect_error(native_unsupervised(X_bad), "data must contain only finite")
  Y_bad <- Y
  Y_bad[1L, 1L] <- NA_real_
  expect_error(native_supervised(X, Y_bad), "Y must contain only finite")

  expect_error(
    native_supervised(
      X, Y, embed = matrix(0, nrow = 7L, ncol = 2L),
      enhanced_prox_mode = 1L
    ),
    "embed must have nrow\\(X\\) rows"
  )
  expect_error(
    native_supervised(
      X, Y, embed = matrix(0, nrow = 8L, ncol = 2L),
      embed_split = -1L, enhanced_prox_mode = 1L
    ),
    "embed_split must be a non-negative"
  )
  expect_error(
    native_supervised(
      X, Y, embed = matrix(0, nrow = 8L, ncol = 2L),
      embed_split = 2L, enhanced_prox_mode = 1L
    ),
    "embed_split must be smaller"
  )
  expect_error(
    native_unsupervised(
      X, embed = matrix(0, nrow = 7L, ncol = 2L),
      enhanced_prox_mode = 1L
    ),
    "embed must have nrow\\(data\\) rows"
  )
  expect_error(native_supervised(X, Y, xvar_wt = c(0, 0, 0)), "positive finite mass")
  expect_error(native_supervised(X, Y, nsplit = -1L), "nsplit")
  expect_error(native_supervised(X, Y, seed = NA_integer_), "seed")
  expect_error(native_supervised(X, Y, prox_mode = 3L), "prox_mode")
  expect_error(native_unsupervised(X, samptype = 2L), "samptype")
})
