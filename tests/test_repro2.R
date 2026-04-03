#!/usr/bin/env Rscript
# Pinpoint where mrf3_fit diverges — layer by layer
library(multiRF)

fit_mv_forest_cpp       <- multiRF:::fit_mv_forest_cpp
fit_mv_forest_unsup_cpp <- multiRF:::fit_mv_forest_unsup_cpp

set.seed(42)
n <- 50; p1 <- 10; p2 <- 8
X1 <- data.frame(matrix(rnorm(n*p1), n, p1)); names(X1) <- paste0("x",1:p1)
X2 <- data.frame(matrix(rnorm(n*p2), n, p2)); names(X2) <- paste0("y",1:p2)
rownames(X1) <- rownames(X2) <- paste0("s",1:n)
dat_list <- list(omics1 = X1, omics2 = X2)

cat("=== Layer-by-layer mrf3_fit reproducibility ===\n\n")

# Run mrf3_fit twice
m1 <- mrf3_fit(dat_list, ntree = 50, seed = 529L, verbose = FALSE)
m2 <- mrf3_fit(dat_list, ntree = 50, seed = 529L, verbose = FALSE)

# ── Layer 1: model_top_v ──
cat("1. model_top_v:  a =", m1$model_top_v, " b =", m2$model_top_v,
    " identical:", identical(m1$model_top_v, m2$model_top_v), "\n")

# ── Layer 2: fused_top_v ──
cat("2. fused_top_v:  a =", m1$fused_top_v, " b =", m2$fused_top_v,
    " identical:", identical(m1$fused_top_v, m2$fused_top_v), "\n")

# ── Layer 3: Reconstruction W_all ──
W1 <- m1$reconstruction$W$W_all
W2 <- m2$reconstruction$W$W_all
if (!is.null(W1) && !is.null(W2)) {
  cat("3. W_all max diff:", max(abs(W1 - W2)), "\n")
} else {
  cat("3. W_all: NULL\n")
}

# ── Layer 4: Reconstruction W_models (per connection) ──
wm1 <- m1$reconstruction$W$W_models
wm2 <- m2$reconstruction$W$W_models
if (!is.null(wm1) && !is.null(wm2)) {
  for (nm in names(wm1)) {
    d <- max(abs(wm1[[nm]] - wm2[[nm]]))
    cat("4. W_models[", nm, "] max diff:", d, "\n")
  }
} else {
  cat("4. W_models: NULL\n")
}

# ── Layer 5: Shared similarity ──
sc1 <- m1$shared$clustering
sc2 <- m2$shared$clustering
if (!is.null(sc1$S) && !is.null(sc2$S)) {
  cat("5. Shared S max diff:", max(abs(sc1$S - sc2$S)), "\n")
} else {
  cat("5. Shared S: not in clustering; checking shared$weights$W_all...\n")
  sw1 <- m1$shared$weights$W_all
  sw2 <- m2$shared$weights$W_all
  if (!is.null(sw1) && !is.null(sw2)) {
    cat("   shared$weights$W_all max diff:", max(abs(sw1 - sw2)), "\n")
  }
}

# ── Layer 6: Shared clusters ──
cl1 <- m1$clusters; cl2 <- m2$clusters
if (is.list(cl1) && !is.null(names(cl1))) {
  # New list structure: shared + per-block specific
  cat("6a. Shared clusters identical:",
      identical(as.integer(cl1$shared), as.integer(cl2$shared)), "\n")
  if (!identical(as.integer(cl1$shared), as.integer(cl2$shared))) {
    cat("   cl1$shared:", head(cl1$shared, 20), "\n")
    cat("   cl2$shared:", head(cl2$shared, 20), "\n")
  }
  for (nm in setdiff(names(cl1), "shared")) {
    cat("6b. Specific cluster[", nm, "] identical:",
        identical(as.integer(cl1[[nm]]), as.integer(cl2[[nm]])), "\n")
  }
} else {
  # Legacy vector structure
  cat("6. Shared clusters identical:", identical(as.integer(cl1), as.integer(cl2)), "\n")
  if (!identical(as.integer(cl1), as.integer(cl2))) {
    cat("   cl1:", head(cl1, 20), "\n")
    cat("   cl2:", head(cl2, 20), "\n")
  }
}

# ── Layer 7: Specific weights ──
sp1 <- m1$specific$weights
sp2 <- m2$specific$weights
if (!is.null(sp1$W) && !is.null(sp2$W)) {
  for (nm in names(sp1$W)) {
    d <- max(abs(sp1$W[[nm]] - sp2$W[[nm]]))
    cat("7. Specific W[", nm, "] max diff:", d, "\n")
  }
}
if (!is.null(sp1$residual) && !is.null(sp2$residual)) {
  for (nm in names(sp1$residual)) {
    d <- max(abs(sp1$residual[[nm]] - sp2$residual[[nm]]))
    cat("7b. Specific residual[", nm, "] max diff:", d, "\n")
  }
}

# ── Layer 8: Specific clusters ──
spc1 <- m1$specific$clustering
spc2 <- m2$specific$clustering
if (!is.null(spc1) && !is.null(spc2)) {
  for (nm in names(spc1)) {
    c1 <- spc1[[nm]]; c2 <- spc2[[nm]]
    # Extract $cl if present
    if (is.list(c1) && !is.null(c1$cl)) c1 <- c1$cl
    if (is.list(c2) && !is.null(c2$cl)) c2 <- c2$cl
    # Only compare numeric/integer vectors (skip metadata like method, similarity_type)
    if (is.null(c1) || is.null(c2) || is.list(c1) || is.list(c2)) next
    if (!is.atomic(c1) || !is.atomic(c2)) next
    if (is.character(c1) || is.character(c2)) next
    cat("8. Specific cluster[", nm, "] identical:",
        identical(as.integer(c1), as.integer(c2)), "\n")
  }
}

# ── Layer 9: Print all top-level names for reference ──
cat("\n--- m1 top-level names ---\n")
cat(paste(names(m1), collapse=", "), "\n")

cat("\n=== DONE ===\n")
