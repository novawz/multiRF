## ============================================================================
## Debug: why does cluster_imd fail on full-MRF models?
##
## Usage: Run in R after preparing the HNSC pipeline state, OR standalone
##        with simulated data.
##
##   Rscript multiRF/tests/test_cluster_imd_debug.R
## ============================================================================

pkg_dir <- {
  if (file.exists("DESCRIPTION")) "."
  else if (file.exists("multiRF/DESCRIPTION")) "multiRF"
  else if (file.exists("../DESCRIPTION")) ".."
  else stop("Cannot find multiRF package root.")
}
devtools::load_all(pkg_dir, quiet = TRUE)
library(randomForestSRC)

cat("=== cluster_imd Debug Test ===\n\n")

## ---- Simulate small data ------------------------------------------------
set.seed(42)
n <- 200
k <- 3
cluster_true <- rep(1:k, each = n / k)
n <- length(cluster_true)  # ensure exact

pA <- 100; pB <- 80
A <- matrix(rnorm(n * pA), n, pA)
B <- matrix(rnorm(n * pB), n, pB)
colnames(A) <- paste0("gene_", seq_len(pA))
colnames(B) <- paste0("cpg_", seq_len(pB))
rownames(A) <- rownames(B) <- paste0("S", seq_len(n))

## Add some signal
for (j in 1:k) {
  idx <- which(cluster_true == j)
  sig_cols <- ((j-1)*10 + 1):(j*10)
  A[idx, sig_cols] <- A[idx, sig_cols] + 2
  B[idx, sig_cols[1:8]] <- B[idx, sig_cols[1:8]] + 1.5
}

dat_list <- list(gene = as.data.frame(A), cpg = as.data.frame(B))

## ---- Fit full MRF -------------------------------------------------------
cat("Fitting full MRF (ntree=100)...\n")
init_mod <- mrf3_init(
  dat.list = dat_list,
  ntree = 100,
  scale = TRUE,
  return_data = TRUE,
  seed = 42
)

cat("  class(init_mod):", class(init_mod), "\n")
cat("  connections:", length(init_mod$connection), "\n")
for (cn in names(init_mod$mod)) {
  m <- init_mod$mod[[cn]]
  cat(sprintf("  Model '%s': membership=%s, forest=%s, node.stats=%s, xvar=%dx%d, yvar=%s\n",
              cn,
              if (!is.null(m$membership)) paste(dim(m$membership), collapse="x") else "NULL",
              if (!is.null(m$forest$nativeArray)) paste(dim(m$forest$nativeArray), collapse="x") else "NULL",
              if (!is.null(m$node.stats)) paste(dim(m$node.stats), collapse="x") else "NULL",
              nrow(m$xvar), ncol(m$xvar),
              if (!is.null(m$yvar)) paste(nrow(m$yvar), "x", ncol(m$yvar)) else "NULL"))
}

## ---- Test cluster_imd ---------------------------------------------------
cat("\n--- Testing cluster_imd ---\n")
clusters <- as.character(cluster_true)
names(clusters) <- rownames(A)

imd_out <- tryCatch({
  cluster_imd(
    x = init_mod,
    cluster = clusters,
    dat.list = dat_list,
    min_cluster_size = 10L,
    keep_model = TRUE,
    keep_data = TRUE,
    seed = 42
  )
}, error = function(e) {
  cat("  TOP-LEVEL ERROR:", conditionMessage(e), "\n")
  cat("  Traceback:\n")
  traceback()
  NULL
})

if (!is.null(imd_out)) {
  cat("\n  Summary:\n")
  print(imd_out$summary)

  ## Check individual cluster failures
  for (i in seq_len(nrow(imd_out$summary))) {
    s <- imd_out$summary[i,]
    if (s$status != "ok") {
      cat(sprintf("\n  Cluster %s FAILED: %s\n", s$cluster, s$message))
    } else {
      cat(sprintf("\n  Cluster %s OK (n=%d)\n", s$cluster, s$n))
    }
  }
}

## ---- Manual test: subset + get_multi_weights on cluster 1 ---------------
cat("\n--- Manual single-cluster test ---\n")
idx <- which(cluster_true == 1)
dat_sub <- lapply(dat_list, function(d) d[idx, , drop = FALSE])
names(dat_sub) <- names(dat_list)

conn <- init_mod$connection[[1]]
conn_name <- paste(conn, collapse = "_")
mod_orig <- init_mod$mod[[conn_name]]

cat(sprintf("  Testing connection '%s' on cluster 1 (n=%d)\n", conn_name, length(idx)))

## Step 1: subset model
cat("  Step 1: Subsetting model...\n")
mod_sub <- tryCatch({
  ms <- mod_orig
  ms$membership <- mod_orig$membership[idx, , drop = FALSE]

  x_name <- conn[2]
  x_block <- dat_sub[[x_name]]
  x_cols <- colnames(mod_orig$xvar)
  x_block <- x_block[, x_cols, drop = FALSE]
  ms$xvar <- as.data.frame(x_block, check.names = FALSE)

  y_name <- conn[1]
  y_block <- dat_sub[[y_name]]
  y_cols <- colnames(mod_orig$yvar)
  y_block <- y_block[, y_cols, drop = FALSE]
  ms$yvar <- as.data.frame(y_block, check.names = FALSE)

  cat(sprintf("    xvar: %dx%d, yvar: %dx%d, membership: %s\n",
              nrow(ms$xvar), ncol(ms$xvar),
              nrow(ms$yvar), ncol(ms$yvar),
              paste(dim(ms$membership), collapse="x")))
  ms
}, error = function(e) {
  cat("    SUBSET ERROR:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(mod_sub)) {
  ## Step 2: get_tree_net on first tree
  cat("  Step 2: get_tree_net(tree.id=1)...\n")
  net <- tryCatch({
    get_tree_net(mod_sub, tree.id = 1)
  }, error = function(e) {
    cat("    get_tree_net ERROR:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(net)) {
    cat(sprintf("    net has %d edges\n", nrow(net)))
  }

  ## Step 3: get_multi_weights
  cat("  Step 3: get_multi_weights...\n")
  mod_list_sub <- list(mod_sub)
  names(mod_list_sub) <- conn_name

  imd_single <- tryCatch({
    get_multi_weights(
      mod_list = mod_list_sub,
      dat.list = dat_sub,
      yprob = 0.5,
      parallel = FALSE,
      normalized = TRUE,
      seed = 42
    )
  }, error = function(e) {
    cat("    get_multi_weights ERROR:", conditionMessage(e), "\n")
    cat("    Error class:", class(e), "\n")
    NULL
  })

  if (!is.null(imd_single)) {
    cat("    SUCCESS! Weight names:", paste(names(imd_single$weight_list), collapse=", "), "\n")
  }
}

cat("\n=== Debug Test Complete ===\n")
