# Regression tests locking in the verified review fixes.
# Each block adapts one of the behavioral validation scripts used to confirm
# the fixes; data are tiny and seeds fixed so the whole file runs in seconds.

# ---- shared fixtures --------------------------------------------------------

# Two separable groups across two named blocks (A: 6 features, B: 5 features).
review_dat <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(42)
      n_per <- 20L
      n <- 2L * n_per
      grp <- rep(c(1, 2), each = n_per)
      mkblock <- function(p, shift, prefix) {
        m <- matrix(rnorm(n * p), nrow = n)
        m[grp == 2, seq_len(ceiling(p / 2))] <-
          m[grp == 2, seq_len(ceiling(p / 2))] + shift
        rownames(m) <- paste0("S", seq_len(n))
        colnames(m) <- paste0(prefix, "_f", seq_len(p))
        m
      }
      cache <<- list(
        dat = list(A = mkblock(6L, 3, "A"), B = mkblock(5L, 3, "B")),
        truth = setNames(grp, paste0("S", seq_len(n)))
      )
    }
    cache
  }
})

# Evaluate `expr` silently (progress cat()s and messages), keeping the value.
review_quiet <- function(expr) {
  res <- NULL
  utils::capture.output(res <- suppressMessages(expr), type = "output")
  res
}

review_fit_basic <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- review_quiet(mrf3_fit(
        review_dat()$dat, ntree = 30, top_v = 10, filter_mode = "none",
        clustering_args = list(shared_k = 2, specific_k = 2),
        run_imd = TRUE, verbose = FALSE, seed = 529
      ))
    }
    cache
  }
})

review_fit_robust <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- review_quiet(mrf3_fit(
        review_dat()$dat, ntree = 30, top_v = 10, filter_mode = "none",
        clustering_args = list(shared_k = 2, specific_k = 2),
        run_imd = TRUE, run_variable_selection = TRUE,
        run_robust_clustering = TRUE, verbose = FALSE, seed = 529
      ))
    }
    cache
  }
})

test_that("a manual one-way connection is preserved through the full workflow", {
  manual <- list(c("A", "B"))
  fit <- review_quiet(suppressWarnings(mrf3_fit(
    review_dat()$dat,
    ntree = 3,
    connect_list = manual,
    top_v = 10,
    filter_mode = "none",
    clustering_args = list(shared_k = 2, specific_k = 2),
    shared_specific_args = list(
      specific_ntree = 3,
      specific_nthread = 1,
      specific_proximity = "none"
    ),
    nthread = 1,
    filter_verbose = FALSE,
    verbose = FALSE,
    seed = 529
  )))

  expect_identical(fit$connection, manual)
  expect_identical(
    fit$reconstruction$block_weight_source,
    c(A = "response", B = "predictor")
  )
  expect_identical(names(fit$reconstruction$W$W_by_block), c("A", "B"))
  expect_identical(
    fit$shared$weights$block_weight_source,
    fit$reconstruction$block_weight_source
  )
  expect_identical(summary(fit)$metadata$block_names, c("A", "B"))
})

# ---- edge-ratio invariant in build_tree_network_cpp -------------------------

test_that("build_tree_network_cpp edge equals parent/child nodesize ratio", {
  # Hand-built preorder tree:
  # R(100) -> A(60) -> {L1(30), L2(30)};  R -> B(40) -> {L3(20), L4(20)}
  net <- multiRF:::build_tree_network_cpp(
    c("R_1", "A_1", "<leaf>_1", "<leaf>_2", "B_1", "<leaf>_3", "<leaf>_4"),
    c(5L, 6L, 1L, 2L, 7L, 3L, 4L),
    c(100L, 60L, 30L, 30L, 40L, 20L, 20L),
    c(0, 1, 2, 2, 1, 2, 2),
    c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE)
  )
  net <- data.frame(net, stringsAsFactors = FALSE)

  expect_identical(net$from, c("R_1", "A_1", "A_1", "R_1", "B_1", "B_1"))
  expect_identical(net$to,
                   c("A_1", "<leaf>_1", "<leaf>_2", "B_1", "<leaf>_3", "<leaf>_4"))
  expect_equal(net$edge, c(100 / 60, 2, 2, 2.5, 2, 2), tolerance = 1e-12)

  # The previously buggy case: an edge emitted after walking back up to an
  # ancestor (R -> B) must use R's size (100/40 = 2.5), not the sibling
  # subtree's (pre-fix value was 0.75).
  expect_equal(net$edge[net$from == "R_1" & net$to == "B_1"], 2.5)
})

test_that("edge-ratio invariant holds across a fitted forest", {
  set.seed(42)
  n <- 60; p <- 5
  X <- as.data.frame(matrix(rnorm(n * p), n, p))
  colnames(X) <- paste0("x", seq_len(p))
  Y <- data.frame(y1 = X$x1 + rnorm(n, sd = 0.3),
                  y2 = X$x2 + rnorm(n, sd = 0.3))
  mod <- fit_forest(X, Y, ntree = 10, nodesize = 5, seed = 11)
  tree_dfs <- multiRF:::.prep_tree_dfs(mod)

  n_edges <- 0L
  n_second_child <- 0L
  for (t in seq_len(mod$ntree)) {
    nt <- multiRF:::get_tree_net(mod, tree.id = t, tree_dfs = tree_dfs)
    if (nrow(nt) == 0L) next

    # Each node's own size is the `nodesize` of the edge where it is the
    # child; the root's size comes from preorder row 1.
    size_of <- setNames(as.numeric(nt$nodesize), nt$to)
    root_nm <- setdiff(unique(nt$from), nt$to)
    expect_length(root_nm, 1L)
    size_of[root_nm] <- as.numeric(tree_dfs[[as.character(t)]]$nodeSZ[1])

    expected <- unname(size_of[nt$from]) / as.numeric(nt$nodesize)
    expect_equal(nt$edge, expected, tolerance = 1e-9)

    # Binary-tree sanity: parent size == sum of children's sizes.
    child_sums <- tapply(as.numeric(nt$nodesize), nt$from, sum)
    expect_equal(as.numeric(child_sums), unname(size_of[names(child_sums)]),
                 tolerance = 1e-9)

    n_edges <- n_edges + nrow(nt)
    n_second_child <- n_second_child + sum(duplicated(nt$from))
  }

  expect_gt(n_edges, 50)          # invariant checked on a nontrivial forest
  expect_gt(n_second_child, 10)   # the buggy "second child" path exercised
})

# ---- spectral_cl clusters the row-normalized embedding ----------------------

test_that("spectral_cl clusters the row-normalized embedding and recovers groups", {
  set.seed(1)
  n_per <- 20
  X <- rbind(matrix(rnorm(n_per * 2, mean = 0), ncol = 2),
             matrix(rnorm(n_per * 2, mean = 8), ncol = 2))
  truth <- rep(1:2, each = n_per)

  D <- as.matrix(dist(X))
  S <- exp(-D^2 / (2 * 2^2))
  res <- spectral_cl(S, k_tune = 2:5)

  expect_identical(as.integer(res$best_k), 2L)
  expect_equal(cluster_ari(res$cl, truth), 1)

  # Returned embedding rows are unit L2 norm (row-normalized).
  expect_identical(nrow(res$embed), nrow(S))
  expect_equal(unname(sqrt(rowSums(res$embed^2))), rep(1, nrow(S)),
               tolerance = 1e-8)

  # Labels must come from clustering the normalized matrix: rerunning pam on
  # res$embed (deterministic) reproduces res$cl exactly.
  cl2 <- cluster::pam(res$embed, k = res$best_k)$cluster
  expect_identical(as.integer(unname(res$cl)), as.integer(unname(cl2)))
})

# ---- mrf3_stability with all-default arguments ------------------------------

test_that("mrf3_stability defaults warn about missing robust branch and still run", {
  fit <- review_fit_basic()
  expect_null(fit$robust_clusters)

  warns <- character(0)
  res <- review_quiet(withCallingHandlers(
    mrf3_stability(fit),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ))

  # Pre-fix behavior was stop(); fixed behavior warns and skips the branch.
  expect_true(any(grepl("robust_clustering.*unavailable.*skipping", warns)))
  expect_identical(names(res$by_branch), "specific_shared")

  mt <- res$by_branch$specific_shared$metrics
  expect_s3_class(mt, "data.frame")
  expect_identical(nrow(mt), 50L)  # default n_rep

  ok <- mt$status == "ok"
  expect_gt(sum(ok), 0L)
  expect_false(all(is.na(mt$ari_vs_base[ok])))
  expect_false(all(is.na(mt$jaccard_vs_base[ok])))
  expect_false(all(is.na(mt$nmi_vs_base[ok])))
  expect_true(all(mt$ari_vs_base[ok] >= -1 - 1e-8 &
                    mt$ari_vs_base[ok] <= 1 + 1e-8, na.rm = TRUE))

  s <- res$summary
  expect_true(is.finite(s$ari_vs_base_mean[s$branch == "specific_shared"]))
})

# ---- mrf3_stability robust branch base labels -------------------------------

test_that("mrf3_stability robust branch uses shared labels of length n as base", {
  fit <- review_fit_robust()
  expect_true(is.list(fit$robust_clusters))
  expect_false(is.null(fit$robust_clusters$shared))

  n <- length(multiRF:::mrf3_sample_names(fit))
  expect_identical(n, 40L)
  k_expected <- length(unique(as.character(fit$robust_clusters$shared)))

  res <- review_quiet(suppressWarnings(
    mrf3_stability(fit, branches = "robust_clustering",
                   n_rep = 8L, seed = 1L, verbose = FALSE)
  ))
  br <- res$by_branch$robust_clustering
  expect_false(is.null(br))

  # Pre-fix code used the whole robust_clusters LIST as base labels, forcing
  # k = length(list) and garbage metrics; fixed code extracts $shared.
  expect_identical(br$settings$k, k_expected)

  cl_reps <- Filter(Negate(is.null), br$cluster_by_rep)
  expect_gt(length(cl_reps), 0L)
  expect_true(all(vapply(cl_reps, length, integer(1)) == n))
  expect_true(all(vapply(
    cl_reps, function(z) identical(names(z), paste0("S", 1:40)), logical(1)
  )))

  mt <- br$metrics
  ok <- mt$status == "ok"
  expect_gt(sum(ok), 0L)
  expect_true(all(is.finite(mt$ari_vs_base[ok])))
  expect_true(all(mt$ari_vs_base[ok] >= -1 - 1e-8 &
                    mt$ari_vs_base[ok] <= 1 + 1e-8))
  # Separable data: agreement with a sane base labeling must be well above 0.
  expect_gt(mean(mt$ari_vs_base[ok]), 0.5)

  expect_true(is.matrix(br$coassign))
  expect_identical(dim(br$coassign), c(n, n))
})

# ---- get_multi_weights(weighted = TRUE) -------------------------------------

test_that("get_multi_weights(weighted = TRUE) runs on the slow traversal path", {
  d <- review_dat()
  mods <- review_quiet(fit_multi_forest(d$dat, ntree = 15, seed = 11))

  # Force the post-hoc tree-traversal path (models without precomputed IMD),
  # which previously crashed multiplying by a NULL edge column.
  mods_slow <- lapply(mods, function(m) {
    m$imd_weights <- NULL
    m$imd_weights_per_tree <- NULL
    m
  })
  names(mods_slow) <- names(mods)

  res <- multiRF:::get_multi_weights(mods_slow, d$dat, weighted = TRUE)
  wl <- res$weight_list
  expect_setequal(names(wl), c("A", "B"))
  for (b in c("A", "B")) {
    expect_length(wl[[b]], ncol(d$dat[[b]]))
    expect_setequal(names(wl[[b]]), colnames(d$dat[[b]]))
    expect_true(all(is.finite(wl[[b]])))
    expect_true(any(wl[[b]] != 0))
  }

  # weighted = TRUE must actually apply the edge factor.
  res0 <- multiRF:::get_multi_weights(mods_slow, d$dat, weighted = FALSE)
  expect_false(isTRUE(all.equal(res$weight_list, res0$weight_list,
                                tolerance = 1e-12)))
})

# ---- pairwise_imd guards ----------------------------------------------------

test_that("pairwise_imd rejects duplicated cross-block feature names", {
  d <- review_dat()
  dat_dup <- d$dat
  colnames(dat_dup$B)[1] <- colnames(dat_dup$A)[1]  # "A_f1" now in both blocks

  fit_dup <- review_quiet(mrf3_fit(
    dat_dup, ntree = 25, top_v = 10, filter_mode = "none",
    clustering_args = list(shared_k = 2, specific_k = 2),
    run_imd = TRUE, verbose = FALSE, seed = 529
  ))

  err <- tryCatch(review_quiet(suppressWarnings(pairwise_imd(fit_dup))),
                  error = function(e) e)
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "unique across blocks", fixed = TRUE)
  expect_match(conditionMessage(err), "A_f1", fixed = TRUE)
})

test_that("pairwise_imd rejects a single-block connection informatively", {
  fit_self <- review_fit_basic()
  fit_self$connection <- c(fit_self$connection, list("A"))

  err <- tryCatch(review_quiet(suppressWarnings(pairwise_imd(fit_self))),
                  error = function(e) e)
  expect_s3_class(err, "error")
  expect_false(grepl("subscript out of bounds", conditionMessage(err),
                     fixed = TRUE))
  expect_match(conditionMessage(err), "two-block connections", fixed = TRUE)

  # A clean two-block fit still works.
  res_ok <- review_quiet(suppressWarnings(pairwise_imd(review_fit_basic())))
  expect_true(is.matrix(res_ok$adj_var_mat))
})

test_that("pairwise_imd accepts a bare mrf3 object (slow path)", {
  # Regression: pairwise_imd_extract_object() stored the IMD weights under
  # `$weights` while the slow path reads `mod$imd`, so every bare `mrf3`
  # object errored with "requires `net`, `connection`, and `imd`".
  bare <- structure(
    list(
      # Connection A-B: one tree whose traversal links A-side responses
      # (Y_id) to B-side split variables (from, suffixed with the node id).
      net = list(
        list(
          list(
            Y_id = c("A_f1", "A_f2"),
            from = c("B_f1_1", "B_f1_2"),
            inv_d = c(0.5, 0.25)
          )
        )
      ),
      connection = list(c("A", "B")),
      imd = list(
        A = c(A_f1 = 0.5, A_f2 = 0.2),
        B = c(B_f1 = 0.7, B_f2 = 0)
      ),
      ntree = 1L
    ),
    class = "mrf3"
  )

  # Default feature source falls back from "selected" with a warning; the
  # non-zero-weight filter then drops B_f2.
  expect_warning(res <- review_quiet(pairwise_imd(bare)),
                 "falling back to non-zero IMD weights")
  expect_s3_class(res, "pairwise_imd_analysis")
  expect_identical(sort(rownames(res$adj_var_mat)), c("A_f1", "A_f2", "B_f1"))
  expect_equal(res$adj_var_mat["A_f1", "B_f1"], 0.5)
  expect_equal(res$adj_var_mat["A_f2", "B_f1"], 0.25)
  expect_equal(res$adj_var_mat, t(res$adj_var_mat))
  expect_equal(res$adj_dat_mat,
               matrix(1, 2, 2, dimnames = list(c("A", "B"), c("A", "B"))))
  expect_identical(res$var_use, list(A = c("A_f1", "A_f2"), B = "B_f1"))

  # feature_source = "all" also reads mod$imd and keeps zero-weight features.
  res_all <- review_quiet(pairwise_imd(bare, feature_source = "all"))
  expect_identical(sort(rownames(res_all$adj_var_mat)),
                   c("A_f1", "A_f2", "B_f1", "B_f2"))
})

# ---- ARI degenerate partitions ----------------------------------------------

test_that("cluster_ari scores identical degenerate partitions as 1", {
  # Zero-denominator cases: both single-cluster, both all-singletons.
  expect_equal(cluster_ari(rep(1L, 5), rep(1L, 5)), 1)
  expect_equal(cluster_ari(1:5, 1:5), 1)

  # Label values must not matter, only partition structure.
  expect_equal(cluster_ari(rep("x", 4), rep(7L, 4)), 1)
  expect_equal(cluster_ari(c(10, 20, 30), c("a", "b", "c")), 1)

  # Non-degenerate cases unaffected.
  expect_equal(cluster_ari(c(1, 1, 2, 2), c(1, 1, 2, 2)), 1)
  expect_equal(cluster_ari(rep(1L, 4), c(1L, 1L, 2L, 2L)), 0)
  # Mixing the two degenerate shapes is non-degenerate and disagrees.
  expect_equal(cluster_ari(rep(1L, 4), 1:4), 0)
})

# ---- compute_oob_fw self-leakage --------------------------------------------

test_that("compute_oob_fw has zero diagonal and renormalized rows", {
  # n = 6, ntree = 3; sample 5 is inbag in every tree (never OOB).
  membership <- cbind(c(1, 1, 2, 2, 3, 3),
                      c(1, 2, 1, 2, 1, 2),
                      c(1, 1, 1, 2, 2, 2))
  inbag <- cbind(c(0, 1, 1, 0, 1, 1),
                 c(1, 0, 1, 1, 1, 1),
                 c(1, 1, 0, 0, 1, 1))

  fw <- multiRF:::compute_oob_fw(list(membership = membership, inbag = inbag))

  expect_identical(max(abs(diag(fw))), 0)                 # no self-leakage
  rs <- rowSums(fw)
  expect_true(all(abs(rs - 1) < 1e-10 | abs(rs) < 1e-10)) # rows sum to 1 or 0
  expect_gte(min(fw), 0)
  expect_lt(abs(rs[5]), 1e-12)                            # never-OOB row is 0
  expect_gt(sum(abs(rs - 1) < 1e-10), 0)

  # Self-weight is redistributed (renormalized), not just zeroed: row 1's
  # whole unit mass sits on non-self donors.
  expect_equal(sum(fw[1, -1]), 1, tolerance = 1e-12)
})

# ---- RNG-state preservation -------------------------------------------------

test_that("prepare_tune_inputs preserves the caller's RNG state", {
  d <- review_dat()
  set.seed(123)
  expected_draws <- runif(5)

  set.seed(123)
  seed_before <- get(".Random.seed", envir = globalenv())
  out <- review_quiet(multiRF:::prepare_tune_inputs(
    dat.list = list(a = as.data.frame(d$dat$A), b = as.data.frame(d$dat$B)),
    mod = NULL, sample_n = 10, seed = 99
  ))
  expect_identical(nrow(out$dat.list[[1]]), 10L)

  # The internal seed = 99 must not disturb the caller's RNG stream.
  expect_identical(get(".Random.seed", envir = globalenv()), seed_before)
  expect_identical(runif(5), expected_draws)
})

# ---- unsupervised forest reproducibility ------------------------------------

test_that("unsupervised C++ forest is seed-reproducible", {
  set.seed(77)
  n <- 30; p <- 6
  X <- rbind(matrix(rnorm(15 * p), ncol = p),
             matrix(rnorm(15 * p, mean = 3), ncol = p))
  rownames(X) <- paste0("s", seq_len(n))
  colnames(X) <- paste0("v", seq_len(p))
  X <- as.data.frame(X)

  run_forest <- function(seed) {
    fit_forest(X = X, type = "unsupervised", ntree = 50, seed = seed,
               forest.wt = "all", proximity = "all")
  }
  f1 <- run_forest(1L)
  f2 <- run_forest(1L)
  f3 <- run_forest(2L)

  expect_identical(f1$engine, "multiRF")  # native C++ path exercised

  expect_identical(f1$membership, f2$membership)
  expect_identical(f1$proximity, f2$proximity)
  expect_identical(f1$forest.wt, f2$forest.wt)

  expect_false(identical(f1$membership, f3$membership))
  expect_false(identical(f1$proximity, f3$proximity))
})
