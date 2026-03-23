## ============================================================================
## Test: sub-MRF ensemble vs full MRF (per-connection architecture)
##
## Each connection is a directed pair: response_block -> predictor_block
## e.g., mRNA -> CpG, mRNA -> miRNA (matching fit_multi_forest design)
##
## Usage:  Rscript tests/test_sub_mrf.R   (from multiRF/ or paper2_vs_cov/)
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

cat("=== Sub-MRF Ensemble Test (Per-Connection Architecture) ===\n\n")

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

eval_cl <- function(pred, ref) {
  c(NMI = cluster_nmi(pred, ref), ARI = cluster_ari(pred, ref))
}

## ---- simulate data (HNSC-like scale) ---------------------------------------
set.seed(42)
n <- 500
k <- 4
cluster_true <- rep(1:k, each = n / k)
n <- length(cluster_true)

## mRNA-like (response): p = 2000
pA <- 2000
sig_A <- 100
mu_A <- matrix(rnorm(k * pA, sd = 0.1), nrow = k)
for (j in seq_len(k)) {
  idx <- ((j - 1) * sig_A + 1):(j * sig_A)
  mu_A[j, idx] <- mu_A[j, idx] + 2
}
A <- mu_A[cluster_true, ] + matrix(rnorm(n * pA), nrow = n)
colnames(A) <- paste0("mRNA_", seq_len(pA))

## methylation-like (predictor): p = 2000
pB <- 2000
sig_B <- 80
mu_B <- matrix(rnorm(k * pB, sd = 0.1), nrow = k)
for (j in seq_len(k)) {
  idx <- ((j - 1) * sig_B + 1):(j * sig_B)
  mu_B[j, idx] <- mu_B[j, idx] + 1.5
}
B <- mu_B[cluster_true, ] + matrix(rnorm(n * pB), nrow = n)
colnames(B) <- paste0("CpG_", seq_len(pB))

## miRNA-like (predictor): p = 300
pC <- 300
sig_C <- 30
mu_C <- matrix(rnorm(k * pC, sd = 0.1), nrow = k)
for (j in seq_len(k)) {
  idx <- ((j - 1) * sig_C + 1):(j * sig_C)
  mu_C[j, idx] <- mu_C[j, idx] + 1.8
}
C <- mu_C[cluster_true, ] + matrix(rnorm(n * pC), nrow = n)
colnames(C) <- paste0("miRNA_", seq_len(pC))

dat_list <- list(mRNA = as.data.frame(A), CpG = as.data.frame(B), miRNA = as.data.frame(C))

## Define connections: response -> predictor (same as mrf3_init would produce)
connect_list <- list(
  c("mRNA", "CpG"),
  c("mRNA", "miRNA")
)

cat(sprintf("Data: n=%d, mRNA=%d (resp), CpG=%d, miRNA=%d (pred), k=%d\n",
            n, pA, pB, pC, k))
cat(sprintf("Connections: %s\n\n",
            paste(sapply(connect_list, paste, collapse="->"), collapse=", ")))

## ---- 1) Full MRF per connection --------------------------------------------
cat("=== Full MRF (per-connection, ntree=500) ===\n")
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
cat(sprintf("  Total runtime: %.1f sec\n", t_full[3]))

## Evaluate each connection
for (conn_name in names(full_mods)) {
  mod <- full_mods[[conn_name]]
  cat(sprintf("\n  Connection: %s\n", conn_name))

  d_prox <- as.dist(1 - mod$proximity)
  pam_fit <- pam(d_prox, k = k)
  m_prox <- eval_cl(pam_fit$clustering, cluster_true)
  cat(sprintf("    [Proximity PAM]   NMI=%.3f  ARI=%.3f\n", m_prox["NMI"], m_prox["ARI"]))

  cl_fw <- fw_similarity_cluster(mod$forest.wt, k = k, top_v = 10)
  m_fw <- eval_cl(cl_fw, cluster_true)
  cat(sprintf("    [Forest.wt Sim]   NMI=%.3f  ARI=%.3f\n", m_fw["NMI"], m_fw["ARI"]))
}

## ---- 2) Sub-MRF settings ---------------------------------------------------
settings <- list(
  list(n_sub = 10, frac_r = 0.15, frac_p = 0.3, ntree = 50, yprob = 0.5, label = "conservative"),
  list(n_sub = 10, frac_r = 0.15, frac_p = 0.3, ntree = 50, yprob = 0.1, label = "low ytry"),
  list(n_sub = 20, frac_r = 0.1,  frac_p = 0.3, ntree = 25, yprob = 0.5, label = "aggressive subset")
)

for (s in settings) {
  total_trees <- s$n_sub * s$ntree * length(connect_list)
  cat(sprintf("\n=== Sub-MRF \"%s\" (n_sub=%d, frac_r=%.2f, frac_p=%.1f, ntree=%d, yprob=%.1f) ===\n",
              s$label, s$n_sub, s$frac_r, s$frac_p, s$ntree, s$yprob))
  cat(sprintf("  Total trees across all connections: %d\n", total_trees))

  t_sub <- system.time({
    sub_mods <- fit_sub_multi_rfsrc(
      dat.list = dat_list,
      connect_list = connect_list,
      n_sub = s$n_sub,
      frac_response = s$frac_r,
      frac_predictor = s$frac_p,
      ntree_per_sub = s$ntree,
      yprob = s$yprob,
      seed = 42,
      verbose = TRUE
    )
  })
  cat(sprintf("  Total runtime: %.1f sec  (Speedup: %.1fx)\n", t_sub[3], t_full[3] / t_sub[3]))

  ## Evaluate each connection
  for (conn_name in names(sub_mods)) {
    sub_mod <- sub_mods[[conn_name]]
    full_mod <- full_mods[[conn_name]]
    cat(sprintf("\n  Connection: %s\n", conn_name))

    ## proximity PAM
    d_sub <- as.dist(1 - sub_mod$proximity)
    pam_sub <- pam(d_sub, k = k)
    m_prox_sub <- eval_cl(pam_sub$clustering, cluster_true)
    cat(sprintf("    [Proximity PAM]   NMI=%.3f  ARI=%.3f\n",
                m_prox_sub["NMI"], m_prox_sub["ARI"]))

    ## forest.wt similarity
    cl_fw_sub <- fw_similarity_cluster(sub_mod$forest.wt, k = k, top_v = 10)
    m_fw_sub <- eval_cl(cl_fw_sub, cluster_true)
    cat(sprintf("    [Forest.wt Sim]   NMI=%.3f  ARI=%.3f\n",
                m_fw_sub["NMI"], m_fw_sub["ARI"]))

    ## correlations with full model
    prox_cor <- cor(as.vector(full_mod$proximity), as.vector(sub_mod$proximity))
    fw_cor   <- cor(as.vector(full_mod$forest.wt), as.vector(sub_mod$forest.wt))
    cat(sprintf("    Proximity cor: %.4f  |  Forest.wt cor: %.4f\n", prox_cor, fw_cor))

    ## coverage
    info <- sub_mod$sub_mrf_info
    cat(sprintf("    Resp coverage: min=%d mean=%.1f max=%d\n",
                min(info$col_coverage_response), mean(info$col_coverage_response),
                max(info$col_coverage_response)))
    cat(sprintf("    Pred coverage: min=%d mean=%.1f max=%d\n",
                min(info$col_coverage_predictor), mean(info$col_coverage_predictor),
                max(info$col_coverage_predictor)))
  }
}

## ---- 3) Sub-MRF with enhanced proximity ------------------------------------
cat("\n=== Sub-MRF with Enhanced Proximity ===\n")

t_enh <- system.time({
  sub_enh_mods <- fit_sub_multi_rfsrc(
    dat.list = dat_list,
    connect_list = connect_list,
    n_sub = 10,
    frac_response = 0.15,
    frac_predictor = 0.3,
    ntree_per_sub = 50,
    enhanced = TRUE,
    yprob = 0.5,
    seed = 42,
    verbose = TRUE
  )
})
cat(sprintf("  Total runtime: %.1f sec\n", t_enh[3]))

## Check enhanced_prox exists and is valid
for (conn_name in names(sub_enh_mods)) {
  mod_enh <- sub_enh_mods[[conn_name]]
  cat(sprintf("\n  Connection: %s\n", conn_name))

  has_enh <- !is.null(mod_enh$enhanced_prox)
  cat(sprintf("    enhanced_prox present: %s\n", has_enh))
  if (has_enh) {
    cat(sprintf("    enhanced_prox dim: %d x %d\n",
                nrow(mod_enh$enhanced_prox), ncol(mod_enh$enhanced_prox)))
    cat(sprintf("    enhanced_prox range: [%.4f, %.4f]\n",
                min(mod_enh$enhanced_prox), max(mod_enh$enhanced_prox)))

    ## enhanced prox PAM
    d_enh <- as.dist(1 - mod_enh$enhanced_prox)
    pam_enh <- pam(d_enh, k = k)
    m_enh <- eval_cl(pam_enh$clustering, cluster_true)
    cat(sprintf("    [Enhanced PAM]    NMI=%.3f  ARI=%.3f\n",
                m_enh["NMI"], m_enh["ARI"]))

    ## plain prox PAM for comparison
    d_plain <- as.dist(1 - mod_enh$proximity)
    pam_plain <- pam(d_plain, k = k)
    m_plain <- eval_cl(pam_plain$clustering, cluster_true)
    cat(sprintf("    [Plain Prox PAM]  NMI=%.3f  ARI=%.3f\n",
                m_plain["NMI"], m_plain["ARI"]))

    ## correlation between enhanced and plain
    enh_cor <- cor(as.vector(mod_enh$enhanced_prox), as.vector(mod_enh$proximity))
    cat(sprintf("    Cor(enhanced, plain): %.4f\n", enh_cor))
  }
}

## ---- 4) Test mrf3_cl_prox dispatch with pre-computed enhanced_prox ---------
cat("\n=== mrf3_cl_prox dispatch test ===\n")

## Should use pre-computed enhanced_prox (no cl_forest call)
tryCatch({
  cl_res <- mrf3_cl_prox(
    rfit = sub_enh_mods,
    k = k,
    enhanced = TRUE,
    method_cl = "PAM"
  )
  m_cl <- eval_cl(cl_res$cl, cluster_true)
  cat(sprintf("  mrf3_cl_prox(enhanced=TRUE, sub_mrf w/ precomputed):  NMI=%.3f  ARI=%.3f\n",
              m_cl["NMI"], m_cl["ARI"]))
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
})

## Should fallback with warning when enhanced_prox not available
cat("\n  Testing fallback (sub-MRF without enhanced_prox)...\n")
sub_no_enh <- sub_mods  # from settings[[1]] above, fitted without enhanced=TRUE
tryCatch({
  withCallingHandlers(
    {
      cl_res2 <- mrf3_cl_prox(
        rfit = sub_no_enh,
        k = k,
        enhanced = TRUE,
        method_cl = "PAM"
      )
      m_cl2 <- eval_cl(cl_res2$cl, cluster_true)
      cat(sprintf("  Fallback result:  NMI=%.3f  ARI=%.3f\n",
                  m_cl2["NMI"], m_cl2["ARI"]))
    },
    warning = function(w) {
      cat(sprintf("  Expected warning: %s\n", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
}, error = function(e) {
  cat(sprintf("  ERROR (unexpected): %s\n", conditionMessage(e)))
})

cat("\n=== Test Complete ===\n")
