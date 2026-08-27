# Tests for cluster_weighted_imd() and the cluster_imd() node-score fast path.

# Shared small supervised fixture -----------------------------------------

.cwi_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(4201)
    n <- 90L
    px <- 8L
    qy <- 5L
    X <- matrix(rnorm(n * px), n, px,
                dimnames = list(sprintf("s%03d", seq_len(n)), paste0("x", seq_len(px))))
    Y <- cbind(
      1.6 * X[, 1L] + rnorm(n, sd = 0.3),
      -1.2 * X[, 2L] + rnorm(n, sd = 0.3),
      X[, 3L] + rnorm(n, sd = 0.5),
      rnorm(n),
      rnorm(n)
    )
    dimnames(Y) <- list(rownames(X), paste0("y", seq_len(qy)))
    mod <- fit_forest(X, Y, ntree = 25L, nodesize = 5L, seed = 7L)
    cl <- stats::setNames(rep(c("A", "B", "C"), length.out = n), rownames(X))
    cache <<- list(X = X, Y = Y, mod = mod, cl = cl)
    cache
  }
})

# Brute-force reference: walk every sample down every tree by the split
# rules alone (never touching `membership`), accumulating node statistics.
# This independently validates the DFS leaf-ID inversion of `membership`.
.cwi_reference <- function(mod, cluster, normalized = TRUE) {
  n <- nrow(mod$xvar)
  px <- ncol(mod$xvar)
  x_names <- colnames(mod$xvar)
  qy <- if (!is.null(mod$yvar) && !identical(mod$xvar, mod$yvar)) ncol(mod$yvar) else 0L
  y_names <- if (qy > 0) colnames(mod$yvar) else character(0)
  cluster <- as.character(cluster[rownames(mod$xvar)])
  labs <- unique(cluster[!is.na(cluster)])
  acc <- function(p, nm) stats::setNames(numeric(p), nm)
  wx <- lapply(labs, function(l) acc(px, x_names)); names(wx) <- labs
  cx <- lapply(labs, function(l) acc(px, x_names)); names(cx) <- labs
  wy <- lapply(labs, function(l) acc(qy, y_names)); names(wy) <- labs
  cy <- lapply(labs, function(l) acc(qy, y_names)); names(cy) <- labs
  Xm <- as.matrix(mod$xvar)
  for (t in seq_len(mod$ntree)) {
    ti <- mod$tree_info[[t]]
    n_nodes <- length(ti$split_var)
    internal_col <- integer(n_nodes)
    cidx <- 0L
    for (ni in seq_len(n_nodes)) {
      if (ti$split_var[ni] >= 0) { cidx <- cidx + 1L; internal_col[ni] <- cidx }
    }
    for (i in seq_len(n)) {
      cl_i <- cluster[i]
      if (is.na(cl_i)) next
      node <- 1L
      repeat {
        sv <- ti$split_var[node]
        if (sv < 0) break
        wx[[cl_i]][sv + 1L] <- wx[[cl_i]][sv + 1L] + ti$imd_x_score[node]
        cx[[cl_i]][sv + 1L] <- cx[[cl_i]][sv + 1L] + 1
        if (qy > 0 && !is.null(ti$imd_y_stats) && nrow(ti$imd_y_stats) == qy) {
          ys <- ti$imd_y_stats[, internal_col[node]]
          pos <- ys > 0
          wy[[cl_i]][pos] <- wy[[cl_i]][pos] + ys[pos]
          cy[[cl_i]][pos] <- cy[[cl_i]][pos] + 1
        }
        node <- if (Xm[i, sv + 1L] <= ti$split_val[node]) {
          ti$left[node] + 1L
        } else {
          ti$right[node] + 1L
        }
      }
    }
  }
  out <- lapply(labs, function(l) {
    ox <- ifelse(cx[[l]] > 0, wx[[l]] / cx[[l]], 0)
    oy <- if (qy > 0) ifelse(cy[[l]] > 0, wy[[l]] / cy[[l]], 0) else numeric(0)
    if (normalized) {
      nx <- sqrt(sum(ox^2)); if (is.finite(nx) && nx > 0) ox <- ox / nx
      if (qy > 0) { ny <- sqrt(sum(oy^2)); if (is.finite(ny) && ny > 0) oy <- oy / ny }
    }
    list(X = ox, Y = oy)
  })
  names(out) <- labs
  out
}

test_that("cluster_weighted_imd matches hand-computed weights on a depth-1 tree", {
  fx <- .cwi_fixture()
  mod1 <- fit_forest(fx$X, fx$Y, ntree = 1L, max_depth = 1L,
                     nodesize = 2L, seed = 31L)
  ti <- mod1$tree_info[[1L]]
  expect_length(ti$split_var, 3L)          # root + 2 leaves
  expect_gte(ti$split_var[1L], 0L)

  out <- cluster_weighted_imd(mod1, fx$cl, normalized = FALSE)
  root_var <- colnames(fx$X)[ti$split_var[1L] + 1L]
  ys <- ti$imd_y_stats[, 1L]
  expected_y <- ifelse(ys > 0, ys, 0)
  names(expected_y) <- colnames(fx$Y)
  for (lb in c("A", "B", "C")) {
    # every sample passes the root and nothing else: the weight for the
    # root's split variable is exactly the root's split score
    expect_equal(unname(out[[lb]]$X[root_var]), ti$imd_x_score[1L])
    expect_true(all(out[[lb]]$X[setdiff(names(out[[lb]]$X), root_var)] == 0))
    expect_equal(out[[lb]]$Y, expected_y, tolerance = 1e-12)
  }
})

test_that("cluster_weighted_imd equals a rule-traversal reference (remap correctness)", {
  fx <- .cwi_fixture()
  ref <- .cwi_reference(fx$mod, fx$cl, normalized = FALSE)
  got <- cluster_weighted_imd(fx$mod, fx$cl, normalized = FALSE)
  expect_identical(names(got), names(ref))
  for (lb in names(ref)) {
    expect_equal(got[[lb]]$X, ref[[lb]]$X, tolerance = 1e-10)
    expect_equal(got[[lb]]$Y, ref[[lb]]$Y, tolerance = 1e-10)
  }

  # normalized variant and NA labels are handled identically
  cl_na <- fx$cl
  cl_na[1:4] <- NA
  refn <- .cwi_reference(fx$mod, cl_na, normalized = TRUE)
  gotn <- cluster_weighted_imd(fx$mod, cl_na, normalized = TRUE)
  for (lb in names(refn)) {
    expect_equal(gotn[[lb]]$X, refn[[lb]]$X, tolerance = 1e-10)
    expect_equal(gotn[[lb]]$Y, refn[[lb]]$Y, tolerance = 1e-10)
  }
})

test_that("cluster_weighted_imd is deterministic with sane, aligned output", {
  fx <- .cwi_fixture()
  a <- cluster_weighted_imd(fx$mod, fx$cl, normalized = TRUE)
  b <- cluster_weighted_imd(fx$mod, fx$cl, normalized = TRUE)
  expect_identical(a, b)

  # name-based alignment: a shuffled cluster vector gives the same answer
  shuffled <- fx$cl[sample(names(fx$cl))]
  expect_identical(cluster_weighted_imd(fx$mod, shuffled, normalized = TRUE), a)

  for (lb in names(a)) {
    expect_identical(names(a[[lb]]$X), colnames(fx$X))
    expect_identical(names(a[[lb]]$Y), colnames(fx$Y))
    expect_true(all(is.finite(a[[lb]]$X)))
    expect_true(all(is.finite(a[[lb]]$Y)))
    expect_gt(stats::sd(a[[lb]]$X), 0)     # non-degenerate
    expect_equal(sqrt(sum(a[[lb]]$X^2)), 1, tolerance = 1e-8)
    expect_equal(sqrt(sum(a[[lb]]$Y^2)), 1, tolerance = 1e-8)
  }
})

test_that("cluster_weighted_imd handles unsupervised models and rejects unsupported ones", {
  fx <- .cwi_fixture()
  umod <- fit_forest(fx$X, NULL, type = "unsupervised", ntree = 10L,
                     nodesize = 5L, seed = 13L)
  out <- cluster_weighted_imd(umod, fx$cl, normalized = FALSE)
  expect_length(out$A$Y, 0L)
  expect_identical(names(out$A$X), colnames(fx$X))
  ref <- .cwi_reference(umod, fx$cl, normalized = FALSE)
  expect_equal(out$A$X, ref$A$X, tolerance = 1e-10)

  bad <- fx$mod
  bad$tree_info <- NULL
  expect_error(cluster_weighted_imd(bad, fx$cl), "tree_info")
  bad2 <- fx$mod
  bad2$tree_info <- lapply(bad2$tree_info, function(ti) { ti$imd_x_score <- NULL; ti })
  expect_error(cluster_weighted_imd(bad2, fx$cl), "imd_x_score")
})

# cluster_imd() wiring -----------------------------------------------------

.cwi_mrf3_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(90210)
    n <- 90L
    cl <- rep(c("c1", "c2", "c3"), each = n / 3L)
    A <- matrix(rnorm(n * 30L), n, 30L,
                dimnames = list(sprintf("t%03d", seq_len(n)), paste0("a", 1:30)))
    B <- matrix(rnorm(n * 20L), n, 20L,
                dimnames = list(rownames(A), paste0("b", 1:20)))
    B[, 1:5] <- B[, 1:5] + 2.0 * A[, 1:5]   # strong shared signal
    dat <- list(A = A, B = B)
    mods <- fit_multi_forest(dat, ntree = 60L, seed = 17L)
    mrf <- structure(
      list(mod = mods, connection = enumerate_connections(names(dat)),
           dat.list = dat, cl = stats::setNames(cl, rownames(A)),
           ntree = 60L, ytry = NULL),
      class = "mrf3"
    )
    cache <<- mrf
    cache
  }
})

test_that("cluster_imd fast path runs, is deterministic, and reports its path", {
  mrf <- .cwi_mrf3_fixture()
  res <- cluster_imd(mrf, use_node_imd = TRUE, parallel = FALSE)
  expect_identical(res$params$imd_path, "node_score")
  expect_true(all(res$summary$status == "ok"))
  res2 <- cluster_imd(mrf, use_node_imd = TRUE, parallel = FALSE)
  expect_identical(res$by_cluster, res2$by_cluster)

  for (lb in names(res$by_cluster)) {
    wl <- res$by_cluster[[lb]]$imd$weight_list
    expect_named(wl, c("A", "B"))
    for (bn in names(wl)) {
      v <- wl[[bn]]
      expect_identical(names(v), colnames(mrf$dat.list[[bn]]))
      expect_true(all(is.finite(v)))
      expect_gt(stats::sd(v), 0)
    }
  }
})

test_that("cluster_imd fast and slow paths agree on strong-signal rankings", {
  mrf <- .cwi_mrf3_fixture()
  fast <- cluster_imd(mrf, use_node_imd = TRUE, parallel = FALSE)
  slow <- suppressMessages(cluster_imd(mrf, use_node_imd = FALSE, parallel = FALSE))
  expect_identical(slow$params$imd_path, "depth_traversal")

  for (lb in names(fast$by_cluster)) {
    for (bn in c("A", "B")) {
      fv <- fast$by_cluster[[lb]]$imd$weight_list[[bn]]
      sv <- slow$by_cluster[[lb]]$imd$weight_list[[bn]][names(fv)]
      # deliberately different estimators: expect positive rank agreement,
      # not equality (measured Spearman ~0.5-0.8 on this fixture)
      expect_gt(suppressWarnings(cor(fv, sv, method = "spearman")), 0.25)
      # the planted strong signal variables must surface in both top-10s
      top_f <- names(sort(fv, decreasing = TRUE))[1:10]
      top_s <- names(sort(sv, decreasing = TRUE))[1:10]
      expect_gte(length(intersect(top_f, top_s)), 4L)
    }
  }
})

test_that("use_node_imd = FALSE is honored even when node stats are available", {
  mrf <- .cwi_mrf3_fixture()
  res <- suppressMessages(cluster_imd(mrf, use_node_imd = FALSE, parallel = FALSE))
  expect_identical(res$params$imd_path, "depth_traversal")
  # default is the depth-based estimator
  res_def <- suppressMessages(cluster_imd(mrf, parallel = FALSE))
  expect_identical(res_def$params$imd_path, "depth_traversal")
})

test_that("fast path falls back or errors cleanly when unavailable or blocked", {
  mrf <- .cwi_mrf3_fixture()

  # models without node stats: auto falls back, forced errors
  mrf_bare <- mrf
  mrf_bare$mod <- lapply(mrf$mod, function(m) { m$tree_info <- NULL; m })
  res_auto <- suppressMessages(cluster_imd(mrf_bare, use_node_imd = NULL, parallel = FALSE))
  expect_identical(res_auto$params$imd_path, "depth_traversal")
  expect_error(
    cluster_imd(mrf_bare, use_node_imd = TRUE, parallel = FALSE),
    "use_node_imd"
  )

  # auto mode with an option the fast path cannot honor: message + slow path
  expect_message(
    res_blocked <- cluster_imd(mrf, use_node_imd = NULL,
                               imd_args = list(use_depth = TRUE),
                               parallel = FALSE),
    "cannot honor"
  )
  expect_identical(res_blocked$params$imd_path, "depth_traversal")

  # forcing the fast path with a blocked option is an error
  expect_error(
    cluster_imd(mrf, use_node_imd = TRUE, imd_args = list(use_depth = TRUE),
                parallel = FALSE),
    "cannot honor"
  )

  # auto mode with nothing blocking selects the fast path
  res_ok <- cluster_imd(mrf, use_node_imd = NULL, parallel = FALSE)
  expect_identical(res_ok$params$imd_path, "node_score")
})
