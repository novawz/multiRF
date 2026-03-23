## ============================================================================
## Convergence analysis: how many sub-MRFs are needed?
##
## Sweeps n_sub from 3 to 50, tracking:
##   - Correlation of proximity / forest.wt with full model
##   - Clustering quality (NMI / ARI)
##   - Cumulative runtime
##
## Also computes a "total effective trees" metric for comparison with full MRF.
## ============================================================================

pkg_dir <- if (file.exists("DESCRIPTION")) {
  "."
} else if (file.exists("multiRF/DESCRIPTION")) {
  "multiRF"
} else if (file.exists("../DESCRIPTION")) {
  ".."
} else {
  stop("Cannot find multiRF package root.")
}
devtools::load_all(pkg_dir, quiet = TRUE)
library(randomForestSRC)
library(cluster)

cat("=== Sub-MRF Convergence Analysis ===\n\n")

## ---------- helpers ----------------------------------------------------------
fw_similarity_cluster <- function(fw, k, top_v = 10) {
  W <- prepare_weight_matrix(
    W = fw, adjust = TRUE, top_v = top_v,
    row_normalize = TRUE, zero_diag = TRUE, keep_ties = TRUE
  )
  S <- W %*% t(W)
  diag(S) <- 0
  cl_fit <- spectral_cl(S, k_tune = k)
  cl_fit$cl
}

## ---- simulate data ---------------------------------------------------------
set.seed(42)
n <- 300
k <- 4
cluster_true <- rep(1:k, each = n / k)

pA <- 1500
sig_A <- 80
mu_A <- matrix(rnorm(k * pA, sd = 0.1), nrow = k)
for (j in seq_len(k)) {
  idx <- ((j - 1) * sig_A + 1):(j * sig_A)
  mu_A[j, idx] <- mu_A[j, idx] + 2
}
A <- mu_A[cluster_true, ] + matrix(rnorm(n * pA), nrow = n)
colnames(A) <- paste0("mRNA_", seq_len(pA))

pB <- 1500
sig_B <- 60
mu_B <- matrix(rnorm(k * pB, sd = 0.1), nrow = k)
for (j in seq_len(k)) {
  idx <- ((j - 1) * sig_B + 1):(j * sig_B)
  mu_B[j, idx] <- mu_B[j, idx] + 1.5
}
B <- mu_B[cluster_true, ] + matrix(rnorm(n * pB), nrow = n)
colnames(B) <- paste0("CpG_", seq_len(pB))

dat_list <- list(mRNA = as.data.frame(A), CpG = as.data.frame(B))
connect_list <- list(c("mRNA", "CpG"))

cat(sprintf("Data: n=%d, mRNA=%d, CpG=%d, k=%d\n", n, pA, pB, k))
cat(sprintf("Connection: mRNA -> CpG\n\n"))

## ---- full MRF reference ----------------------------------------------------
cat("Fitting full MRF (ntree=500)...\n")
t_full <- system.time({
  full_mods <- fit_multi_forest(
    dat.list = dat_list,
    connect_list = connect_list,
    ntree = 500,
    forest.wt = "all",
    proximity = "all",
    yprob = 0.5,
    seed = 42
  )
})
cat(sprintf("  Full MRF runtime: %.1f sec\n\n", t_full[3]))

full_mod <- full_mods[[1]]
full_prox <- as.vector(full_mod$proximity)
full_fw   <- as.vector(full_mod$forest.wt)

## ---- sweep n_sub -----------------------------------------------------------
ntree_per_sub <- 50
frac_r <- 0.15
frac_p <- 0.3

n_sub_values <- c(3, 5, 8, 10, 15, 20, 30, 50)

results <- data.frame(
  n_sub = integer(),
  total_trees = integer(),
  time_sec = numeric(),
  prox_cor = numeric(),
  fw_cor = numeric(),
  prox_nmi = numeric(),
  fw_nmi = numeric(),
  resp_cov_min = integer(),
  pred_cov_min = integer(),
  stringsAsFactors = FALSE
)

cat(sprintf("Sweep: ntree_per_sub=%d, frac_r=%.2f, frac_p=%.1f\n", ntree_per_sub, frac_r, frac_p))
cat(sprintf("%-8s %-12s %-10s %-12s %-12s %-10s %-10s\n",
            "n_sub", "total_trees", "time(s)", "prox_cor", "fw_cor", "prox_NMI", "fw_NMI"))
cat(strrep("-", 82), "\n")

for (ns in n_sub_values) {
  t_sub <- system.time({
    sub_mods <- fit_sub_multi_rfsrc(
      dat.list = dat_list,
      connect_list = connect_list,
      n_sub = ns,
      frac_response = frac_r,
      frac_predictor = frac_p,
      ntree_per_sub = ntree_per_sub,
      yprob = 0.5,
      seed = 42,
      verbose = FALSE
    )
  })

  sub_mod <- sub_mods[[1]]

  prox_cor <- cor(full_prox, as.vector(sub_mod$proximity))
  fw_cor   <- cor(full_fw, as.vector(sub_mod$forest.wt))

  d_sub <- as.dist(1 - sub_mod$proximity)
  pam_sub <- pam(d_sub, k = k)
  prox_nmi <- cluster_nmi(pam_sub$clustering, cluster_true)

  cl_fw <- fw_similarity_cluster(sub_mod$forest.wt, k = k, top_v = 10)
  fw_nmi <- cluster_nmi(cl_fw, cluster_true)

  info <- sub_mod$sub_mrf_info
  total_trees <- ns * ntree_per_sub

  row <- data.frame(
    n_sub = ns,
    total_trees = total_trees,
    time_sec = round(t_sub[3], 1),
    prox_cor = round(prox_cor, 4),
    fw_cor = round(fw_cor, 4),
    prox_nmi = round(prox_nmi, 3),
    fw_nmi = round(fw_nmi, 3),
    resp_cov_min = min(info$col_coverage_response),
    pred_cov_min = min(info$col_coverage_predictor),
    stringsAsFactors = FALSE
  )
  results <- rbind(results, row)

  cat(sprintf("%-8d %-12d %-10.1f %-12.4f %-12.4f %-10.3f %-10.3f\n",
              ns, total_trees, t_sub[3], prox_cor, fw_cor, prox_nmi, fw_nmi))
}

cat("\n")
cat(sprintf("Full MRF reference: 500 trees, %.1f sec\n", t_full[3]))

## ---- marginal gain analysis ------------------------------------------------
cat("\n=== Marginal Gain Analysis ===\n")
cat("(Incremental improvement in correlation per additional sub-MRF)\n\n")

for (i in 2:nrow(results)) {
  delta_n <- results$n_sub[i] - results$n_sub[i - 1]
  delta_prox <- results$prox_cor[i] - results$prox_cor[i - 1]
  delta_fw   <- results$fw_cor[i] - results$fw_cor[i - 1]
  cat(sprintf("  n_sub %2d -> %2d: delta_prox_cor = %+.4f, delta_fw_cor = %+.4f  (per sub: %+.5f, %+.5f)\n",
              results$n_sub[i - 1], results$n_sub[i],
              delta_prox, delta_fw,
              delta_prox / delta_n, delta_fw / delta_n))
}

## ---- comparison: total trees needed ----------------------------------------
cat("\n=== Effective Tree Comparison ===\n")
cat("How many sub-MRF total trees match full MRF correlation?\n\n")

target_cors <- c(0.95, 0.96, 0.97, 0.98)
for (tc in target_cors) {
  idx_prox <- which(results$prox_cor >= tc)
  idx_fw   <- which(results$fw_cor >= tc)
  if (length(idx_prox) > 0) {
    cat(sprintf("  prox_cor >= %.2f: n_sub >= %d (total_trees = %d, %.1f sec)\n",
                tc, results$n_sub[idx_prox[1]], results$total_trees[idx_prox[1]],
                results$time_sec[idx_prox[1]]))
  } else {
    cat(sprintf("  prox_cor >= %.2f: not reached (max = %.4f)\n", tc, max(results$prox_cor)))
  }
}

cat("\n=== Convergence Test Complete ===\n")
