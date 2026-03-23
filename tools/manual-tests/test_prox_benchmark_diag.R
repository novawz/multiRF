# test_prox_benchmark_diag.R
# ─────────────────────────────────────────────────────────────────
# Diagnose why native vs rfsrc proximity clustering differs.
# Uses the SAME data setup as the multiRF_prox benchmark.
# ─────────────────────────────────────────────────────────────────

library(multiRF)

# ── Same data as benchmark ──
sim <- simulate_intersim_shared_specific(
  n = 500,
  Z_prop = c(0.1, 0.2, 0.3, 0.4),
  delta_methyl = 0.5,
  delta_expr = 0.5,
  delta_protein = 0.5,
  p_DMP = 0.2,
  p_spec_feat = 0.2,
  U_nclass = 2,
  U_within_Z = TRUE,
  return_workflow_data = TRUE,
  keep_sim_raw = FALSE,
  seed = 20260320L
)
dat_list <- sim$dat.list
truth_z <- as.character(sim$truth)

cat("═══════════════════════════════════════════════════════\n")
cat("  Data: ", length(dat_list), " omics blocks, n =", nrow(dat_list[[1]]), "\n")
cat("  True Z clusters: ", paste(table(truth_z), collapse=", "), "\n")
cat("═══════════════════════════════════════════════════════\n\n")

# ── Fit forests with BOTH engines ──
ntree <- 300L
seed_fit <- 529L

cat("Fitting native forests...\n")
options(multiRF.engine = "native")
t0 <- proc.time()
native_mods <- fit_multi_forest(dat_list, ntree = ntree, seed = seed_fit)
t_native <- (proc.time() - t0)[3]

cat("Fitting rfsrc forests...\n")
options(multiRF.engine = "randomForestSRC")
t0 <- proc.time()
rfsrc_mods <- fit_multi_forest(dat_list, ntree = ntree, seed = seed_fit)
t_rfsrc <- (proc.time() - t0)[3]

cat(sprintf("  Native: %.1fs, rfsrc: %.1fs\n\n", t_native, t_rfsrc))

# ── Compare proximity per forest ──
cat("── Per-forest proximity comparison ──\n")
for (nm in names(native_mods)) {
  pn <- native_mods[[nm]]$proximity
  pr <- rfsrc_mods[[nm]]$proximity

  # Basic checks
  cat(sprintf("  %s:\n", nm))
  cat(sprintf("    dims: native %dx%d, rfsrc %dx%d\n",
      nrow(pn), ncol(pn), nrow(pr), ncol(pr)))
  cat(sprintf("    diag: native [%.4f, %.4f], rfsrc [%.4f, %.4f]\n",
      min(diag(pn)), max(diag(pn)), min(diag(pr)), max(diag(pr))))
  cat(sprintf("    range: native [%.4f, %.4f], rfsrc [%.4f, %.4f]\n",
      min(pn), max(pn), min(pr), max(pr)))
  cat(sprintf("    any NA: native %s, rfsrc %s\n",
      any(is.na(pn)), any(is.na(pr))))

  # Correlation
  if (nrow(pn) == nrow(pr) && ncol(pn) == ncol(pr)) {
    cc <- cor(as.numeric(pn), as.numeric(pr))
    cat(sprintf("    cor: %.4f\n", cc))

    # Within vs between cluster separation
    same_z <- outer(truth_z, truth_z, "==")
    within_n <- mean(pn[same_z])
    between_n <- mean(pn[!same_z])
    within_r <- mean(pr[same_z])
    between_r <- mean(pr[!same_z])
    cat(sprintf("    native  within=%.4f between=%.4f ratio=%.2f\n",
        within_n, between_n, within_n / between_n))
    cat(sprintf("    rfsrc   within=%.4f between=%.4f ratio=%.2f\n",
        within_r, between_r, within_r / between_r))
  }
  cat("\n")
}

# ── Combined proximity ──
cat("── Combined proximity ──\n")
combined_n <- Reduce("+", lapply(native_mods, "[[", "proximity"))
combined_r <- Reduce("+", lapply(rfsrc_mods, "[[", "proximity"))

# Normalize same way as mrf3_cl_prox
num_dim_n <- diag(combined_n)[1]
num_dim_r <- diag(combined_r)[1]
cat(sprintf("  num_dim: native=%.4f, rfsrc=%.4f\n", num_dim_n, num_dim_r))
combined_n <- combined_n / num_dim_n
combined_r <- combined_r / num_dim_r

cc <- cor(as.numeric(combined_n), as.numeric(combined_r))
cat(sprintf("  Combined cor: %.4f\n", cc))

same_z <- outer(truth_z, truth_z, "==")
within_n <- mean(combined_n[same_z])
between_n <- mean(combined_n[!same_z])
within_r <- mean(combined_r[same_z])
between_r <- mean(combined_r[!same_z])
cat(sprintf("  native  within=%.4f between=%.4f ratio=%.2f\n",
    within_n, between_n, within_n / between_n))
cat(sprintf("  rfsrc   within=%.4f between=%.4f ratio=%.2f\n",
    within_r, between_r, within_r / between_r))

# ── Cluster both and compare ──
cat("\n── Clustering comparison ──\n")
do_pam <- function(prox_mat, k) {
  d <- 1 - prox_mat
  d <- (d + t(d)) / 2
  diag(d) <- 0
  cl <- cluster::pam(as.dist(d), k = k)$clustering
  cl
}
k <- length(unique(truth_z))
cl_n <- do_pam(combined_n, k)
cl_r <- do_pam(combined_r, k)

ari_n <- multiRF:::cluster_metrics(as.character(cl_n), truth_z)$ari[1]
ari_r <- multiRF:::cluster_metrics(as.character(cl_r), truth_z)$ari[1]
cat(sprintf("  Native ARI_Z: %.4f\n", ari_n))
cat(sprintf("  rfsrc  ARI_Z: %.4f\n", ari_r))

# ── Check leaf statistics ──
cat("\n── Leaf statistics (first forest) ──\n")
nm1 <- names(native_mods)[1]
mem_n <- native_mods[[nm1]]$membership
mem_r <- rfsrc_mods[[nm1]]$membership
cat(sprintf("  Native: mean leaves=%.1f\n",
    mean(apply(mem_n, 2, function(x) length(unique(x))))))
cat(sprintf("  rfsrc:  mean leaves=%.1f\n",
    mean(apply(mem_r, 2, function(x) length(unique(x))))))
