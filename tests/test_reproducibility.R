#!/usr/bin/env Rscript
# Minimal reproducibility test for multiRF
# Run: Rscript tests/test_reproducibility.R

library(multiRF)

# Internal C++ wrappers are not exported; bind them explicitly for this
# diagnostics script.
fit_mv_forest_cpp <- multiRF:::fit_mv_forest_cpp
fit_mv_forest_unsup_cpp <- multiRF:::fit_mv_forest_unsup_cpp
compute_oob_forest_wt_cpp <- multiRF:::compute_oob_forest_wt_cpp

cat("=== multiRF Reproducibility Test ===\n\n")

# Small synthetic data
set.seed(42)
n <- 50; p1 <- 10; p2 <- 8
X1 <- matrix(rnorm(n * p1), n, p1, dimnames = list(paste0("s", 1:n), paste0("x", 1:p1)))
X2 <- matrix(rnorm(n * p2), n, p2, dimnames = list(paste0("s", 1:n), paste0("y", 1:p2)))

# ── Test 1: Raw C++ supervised forest ──
cat("--- Test 1: C++ supervised forest (fit_mv_forest_cpp) ---\n")
X_mat <- as.matrix(data.frame(X1))
Y_mat <- as.matrix(data.frame(X2))
res1a <- fit_mv_forest_cpp(X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 1L)
res1b <- fit_mv_forest_cpp(X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 1L)
fw_diff_1 <- max(abs(res1a$forest.wt - res1b$forest.wt), na.rm = TRUE)
cat("  forest.wt max diff (1 thread):", fw_diff_1, "\n")
cat("  PASS:", fw_diff_1 == 0, "\n\n")

# ── Test 1b: Multi-thread ──
cat("--- Test 1b: C++ supervised forest (multi-thread) ---\n")
res1c <- fit_mv_forest_cpp(X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 4L)
res1d <- fit_mv_forest_cpp(X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 4L)
fw_diff_1b <- max(abs(res1c$forest.wt - res1d$forest.wt), na.rm = TRUE)
cat("  forest.wt max diff (4 threads):", fw_diff_1b, "\n")
cat("  PASS:", fw_diff_1b == 0, "\n\n")

# ── Test 1c: 1 thread vs 4 threads ──
cat("--- Test 1c: 1 thread vs 4 threads ---\n")
fw_diff_1c <- max(abs(res1a$forest.wt - res1c$forest.wt), na.rm = TRUE)
cat("  forest.wt max diff (1 vs 4 threads):", fw_diff_1c, "\n")
cat("  PASS:", fw_diff_1c == 0, "\n\n")

# ── Test 2: Raw C++ unsupervised forest ──
cat("--- Test 2: C++ unsupervised forest (fit_mv_forest_unsup_cpp) ---\n")
data_mat <- as.matrix(data.frame(X1))
res2a <- fit_mv_forest_unsup_cpp(data_mat, ntree = 50L, seed = 529L, nthread = 1L)
res2b <- fit_mv_forest_unsup_cpp(data_mat, ntree = 50L, seed = 529L, nthread = 1L)
fw_diff_2 <- max(abs(res2a$forest.wt - res2b$forest.wt), na.rm = TRUE)
cat("  forest.wt max diff (1 thread):", fw_diff_2, "\n")
cat("  PASS:", fw_diff_2 == 0, "\n\n")

# ── Test 2b: Unsupervised multi-thread ──
cat("--- Test 2b: C++ unsupervised forest (multi-thread) ---\n")
res2c <- fit_mv_forest_unsup_cpp(data_mat, ntree = 50L, seed = 529L, nthread = 4L)
res2d <- fit_mv_forest_unsup_cpp(data_mat, ntree = 50L, seed = 529L, nthread = 4L)
fw_diff_2b <- max(abs(res2c$forest.wt - res2d$forest.wt), na.rm = TRUE)
cat("  forest.wt max diff (4 threads):", fw_diff_2b, "\n")
cat("  PASS:", fw_diff_2b == 0, "\n\n")

# ── Test 3: R-level fit_forest wrapper ──
cat("--- Test 3: R fit_forest (supervised) ---\n")
X1_df <- as.data.frame(X1)
X2_df <- as.data.frame(X2)
r3a <- fit_forest(X = X1_df, Y = X2_df, ntree = 50, seed = 529L)
r3b <- fit_forest(X = X1_df, Y = X2_df, ntree = 50, seed = 529L)
fw_diff_3 <- max(abs(r3a$forest.wt - r3b$forest.wt))
cat("  forest.wt max diff:", fw_diff_3, "\n")
cat("  PASS:", fw_diff_3 == 0, "\n\n")

# ── Test 4: R-level fit_forest (unsupervised) ──
cat("--- Test 4: R fit_forest (unsupervised) ---\n")
r4a <- fit_forest(X = X1_df, Y = NULL, type = "unsupervised", ntree = 50, seed = 529L)
r4b <- fit_forest(X = X1_df, Y = NULL, type = "unsupervised", ntree = 50, seed = 529L)
fw_diff_4 <- max(abs(r4a$forest.wt - r4b$forest.wt))
cat("  forest.wt max diff:", fw_diff_4, "\n")
cat("  PASS:", fw_diff_4 == 0, "\n\n")

# ── Test 5: OOB forest weight ──
cat("--- Test 5: OOB forest weight (compute_oob_forest_wt_cpp) ---\n")
mem1 <- res1a$membership
inb1 <- res1a$inbag
oob_a <- compute_oob_forest_wt_cpp(mem1, inb1)
oob_b <- compute_oob_forest_wt_cpp(mem1, inb1)
oob_diff <- max(abs(oob_a - oob_b))
cat("  OOB wt max diff:", oob_diff, "\n")
cat("  PASS:", oob_diff == 0, "\n\n")

# ── Test 6: Full mrf3_fit pipeline ──
cat("--- Test 6: Full mrf3_fit pipeline ---\n")
dat_list <- list(omics1 = X1_df, omics2 = X2_df)
m6a <- mrf3_fit(dat_list, ntree = 50, seed = 529L, verbose = FALSE)
m6b <- mrf3_fit(dat_list, ntree = 50, seed = 529L, verbose = FALSE)

# Compare shared similarity matrix
if (!is.null(m6a$shared_clustering$S) && !is.null(m6b$shared_clustering$S)) {
  S_diff <- max(abs(m6a$shared_clustering$S - m6b$shared_clustering$S), na.rm = TRUE)
  cat("  Shared S max diff:", S_diff, "\n")
} else {
  S_diff <- NA_real_
  cat("  Shared S max diff: NA (shared_clustering$S unavailable)\n")
}

# Compare shared clusters
cl_a <- m6a$clusters
cl_b <- m6b$clusters
ari <- if (requireNamespace("mclust", quietly = TRUE)) {
  mclust::adjustedRandIndex(cl_a, cl_b)
} else {
  NA
}
cat("  Shared cluster ARI:", ari, "\n")
cat("  Clusters identical:", identical(as.integer(cl_a), as.integer(cl_b)), "\n")

# Compare model_top_v
cat("  model_top_v a:", m6a$top_v$model_top_v, "\n")
cat("  model_top_v b:", m6b$top_v$model_top_v, "\n")
cat("  model_top_v identical:", identical(m6a$top_v$model_top_v, m6b$top_v$model_top_v), "\n\n")

# Compare specific
if (!is.null(m6a$specific_clustering$S) && !is.null(m6b$specific_clustering$S)) {
  for (nm in names(m6a$specific_clustering$S)) {
    sp_diff <- max(abs(m6a$specific_clustering$S[[nm]] - m6b$specific_clustering$S[[nm]]), na.rm = TRUE)
    cat("  Specific S[", nm, "] max diff:", sp_diff, "\n")
  }
}

cat("\n=== DONE ===\n")
