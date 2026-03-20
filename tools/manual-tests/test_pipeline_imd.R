## ============================================================================
## Test: Full pipeline under different sub-MRF / IMD scenarios
##
## Covers:
##   1. Full forest + IMD (baseline)
##   2. Sub-MRF + compute_imd = TRUE (pre-computed IMD)
##   3. Sub-MRF + compute_imd = FALSE + run_imd = TRUE (should warn & skip)
##   4. Sub-MRF with per-dimension threshold (mixed full/sub per connection)
##   5. mrf3_fit end-to-end: full forest
##   6. mrf3_fit end-to-end: sub-MRF with auto compute_imd injection
##   7. Sub-MRF IMD weights sanity check (coverage, non-zero, correlation)
##   8. Single connection fit_sub_mrf with compute_imd
##   9. Full forest + variable selection
##  10. Sub-MRF + variable selection (auto fallback to mixture)
##  11. Full forest + robust clustering
##  12. Sub-MRF + robust clustering
##  13. Variable selection agreement (sub-MRF vs full forest)
##
## Usage:  Rscript tests/test_pipeline_imd.R   (from multiRF/ or paper2_vs_cov/)
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

## ---- global helpers ---------------------------------------------------------
pass <- function(msg) cat(sprintf("  PASS: %s\n", msg))
fail <- function(msg) { cat(sprintf("  FAIL: %s\n", msg)); FAILURES <<- FAILURES + 1L }
check <- function(cond, msg) if (isTRUE(cond)) pass(msg) else fail(msg)
FAILURES <- 0L
TOTAL_START <- proc.time()

## ---- report helpers ---------------------------------------------------------
## Print weight distribution summary for a named numeric vector
report_weights <- function(w, label = "") {
  if (is.null(w) || length(w) == 0) {
    cat(sprintf("    [%s] (empty or NULL)\n", label))
    return(invisible(NULL))
  }
  n_zero <- sum(w == 0)
  n_pos  <- sum(w > 0)
  cat(sprintf("    [%s] len=%d | non-zero=%d (%.1f%%) | range=[%.4f, %.4f] | mean=%.4f | sd=%.4f\n",
              label, length(w), n_pos, 100 * n_pos / length(w),
              min(w), max(w), mean(w), sd(w)))
}

## Print top-N features by weight
report_top <- function(w, label = "", top_n = 5) {
  if (is.null(w) || length(w) == 0) return(invisible(NULL))
  w_sorted <- sort(w, decreasing = TRUE)
  top <- head(w_sorted, top_n)
  feat_str <- paste(sprintf("%s=%.4f", names(top), top), collapse = ", ")
  cat(sprintf("    [%s] top-%d: %s\n", label, top_n, feat_str))
}

## Print model structure summary
report_model_structure <- function(mod, conn_name = "") {
  cat(sprintf("    [%s] class=%s\n", conn_name, paste(class(mod), collapse = ",")))
  fields <- names(mod)
  cat(sprintf("    [%s] fields: %s\n", conn_name, paste(fields, collapse = ", ")))
  is_sub <- !is.null(mod$sub_mrf_info)
  is_full <- !is.null(mod$membership)
  cat(sprintf("    [%s] sub_mrf=%s, full_forest=%s\n", conn_name, is_sub, is_full))
  if (is_sub) {
    info <- mod$sub_mrf_info
    cat(sprintf("    [%s] n_sub=%s, frac_response=%s, frac_predictor=%s\n",
                conn_name,
                if (!is.null(info$n_sub)) info$n_sub else "?",
                if (!is.null(info$frac_response)) info$frac_response else "?",
                if (!is.null(info$frac_predictor)) info$frac_predictor else "?"))
    if (!is.null(info$col_coverage_response)) {
      cov_r <- info$col_coverage_response
      cat(sprintf("    [%s] response coverage: %d/%d features covered (%.1f%%)\n",
                  conn_name, sum(cov_r > 0), length(cov_r), 100 * mean(cov_r > 0)))
    }
    if (!is.null(info$col_coverage_predictor)) {
      cov_p <- info$col_coverage_predictor
      cat(sprintf("    [%s] predictor coverage: %d/%d features covered (%.1f%%)\n",
                  conn_name, sum(cov_p > 0), length(cov_p), 100 * mean(cov_p > 0)))
    }
  }
  if (is_full) {
    cat(sprintf("    [%s] ntree=%s, n=%s, xvar.ncol=%s, yvar.ncol=%s\n",
                conn_name,
                if (!is.null(mod$ntree)) mod$ntree else "?",
                if (!is.null(mod$n)) mod$n else "?",
                if (!is.null(mod$xvar)) ncol(mod$xvar) else "?",
                if (!is.null(mod$yvar)) ncol(mod$yvar) else "?"))
  }
  if (!is.null(mod$forest.wt)) {
    cat(sprintf("    [%s] forest.wt: %s\n", conn_name, paste(dim(mod$forest.wt), collapse = "x")))
  }
  if (!is.null(mod$proximity)) {
    cat(sprintf("    [%s] proximity: %s\n", conn_name, paste(dim(mod$proximity), collapse = "x")))
  }
}

## Print mrf3_fit structure
report_fit_structure <- function(fit, label = "") {
  cat(sprintf("  [%s] mrf3_fit structure:\n", label))
  cat(sprintf("    fields: %s\n", paste(names(fit), collapse = ", ")))
  if (!is.null(fit$models)) {
    cat(sprintf("    n_models=%d, connections: %s\n",
                length(fit$models),
                paste(names(fit$models), collapse = ", ")))
    for (cn in names(fit$models)) {
      m <- fit$models[[cn]]
      is_sub <- !is.null(m$sub_mrf_info)
      cat(sprintf("    model[%s]: sub_mrf=%s\n", cn, is_sub))
    }
  }
  if (!is.null(fit$weights)) {
    cat(sprintf("    weights blocks: %s\n", paste(names(fit$weights), collapse = ", ")))
    for (bn in names(fit$weights)) {
      report_weights(fit$weights[[bn]], label = bn)
    }
  }
  if (!is.null(fit$weights_init)) {
    cat(sprintf("    weights_init: present (per-tree distributions)\n"))
  } else {
    cat(sprintf("    weights_init: NULL (no per-tree distributions)\n"))
  }
  if (!is.null(fit$selected_vars)) {
    cat(sprintf("    selected_vars blocks: %s\n", paste(names(fit$selected_vars), collapse = ", ")))
    for (bn in names(fit$selected_vars)) {
      sv <- fit$selected_vars[[bn]]
      cat(sprintf("      %s: %d vars selected\n", bn, length(sv)))
    }
  }
  if (!is.null(fit$robust_clusters)) {
    tbl <- table(fit$robust_clusters)
    cat(sprintf("    robust_clusters: %d labels, sizes: %s\n",
                length(tbl), paste(sprintf("%s=%d", names(tbl), tbl), collapse = ", ")))
  }
  if (!is.null(fit$config)) {
    cat(sprintf("    config: ntree=%s, yprob=%s\n",
                if (!is.null(fit$config$ntree)) fit$config$ntree else "?",
                if (!is.null(fit$config$yprob)) fit$config$yprob else "?"))
  }
}

## Elapsed time helper
elapsed <- function(t0) {
  dt <- proc.time() - t0
  sprintf("%.1fs", dt["elapsed"])
}

## ---- simulate data ----------------------------------------------------------
cat("=== Pipeline IMD Test Suite ===\n\n")
cat("Setting up simulated data...\n")
set.seed(42)
n <- 150
k <- 3
cluster_true <- rep(1:k, each = n / k)
n <- length(cluster_true)

## Block A: large response (p=500) -- will use sub-MRF
pA <- 500
A <- matrix(rnorm(n * pA), n, pA)
for (j in 1:k) {
  idx <- which(cluster_true == j)
  sig <- ((j - 1) * 30 + 1):(j * 30)
  A[idx, sig] <- A[idx, sig] + 2.5
}
colnames(A) <- paste0("geneA_", seq_len(pA))
rownames(A) <- paste0("S", seq_len(n))

## Block B: medium predictor (p=200)
pB <- 200
B <- matrix(rnorm(n * pB), n, pB)
for (j in 1:k) {
  idx <- which(cluster_true == j)
  sig <- ((j - 1) * 15 + 1):(j * 15)
  B[idx, sig] <- B[idx, sig] + 2.0
}
colnames(B) <- paste0("cpg_", seq_len(pB))
rownames(B) <- paste0("S", seq_len(n))

## Block C: small predictor (p=50) -- below threshold for sub
pC <- 50
C <- matrix(rnorm(n * pC), n, pC)
for (j in 1:k) {
  idx <- which(cluster_true == j)
  sig <- ((j - 1) * 5 + 1):(j * 5)
  C[idx, sig] <- C[idx, sig] + 2.2
}
colnames(C) <- paste0("miR_", seq_len(pC))
rownames(C) <- paste0("S", seq_len(n))

dat_list <- list(
  geneA = as.data.frame(A),
  cpg   = as.data.frame(B),
  miR   = as.data.frame(C)
)

connect_list <- list(
  c("geneA", "cpg"),
  c("geneA", "miR")
)

cat(sprintf("  n=%d, k=%d, geneA=%d, cpg=%d, miR=%d\n", n, k, pA, pB, pC))
cat(sprintf("  Connections: %s\n",
            paste(sapply(connect_list, paste, collapse = "->"), collapse = ", ")))
cat(sprintf("  Signal features per cluster: geneA=30, cpg=15, miR=5\n"))
cat(sprintf("  Signal strength: geneA=2.5, cpg=2.0, miR=2.2\n\n"))


## ============================================================================
## Test 1: Full forest + get_multi_weights (baseline)
## ============================================================================
cat("--- Test 1: Full forest + IMD (baseline) ---\n")
t0 <- proc.time()

full_mods <- fit_multi_forest(
  dat.list = dat_list,
  connect_list = connect_list,
  ntree = 100,
  forest.wt = "all",
  yprob = 0.5,
  seed = 42
)

cat(sprintf("  fit_multi_forest done [%s]\n", elapsed(t0)))
cat("  Model structures:\n")
for (cn in names(full_mods)) {
  report_model_structure(full_mods[[cn]], cn)
}

t1 <- proc.time()
imd_full <- tryCatch(
  get_multi_weights(
    mod_list = full_mods,
    dat.list = dat_list,
    yprob = 0.5,
    parallel = FALSE,
    normalized = TRUE,
    seed = 42
  ),
  error = function(e) { fail(paste("get_multi_weights error:", e$message)); NULL }
)
cat(sprintf("  get_multi_weights done [%s]\n", elapsed(t1)))

if (!is.null(imd_full)) {
  check(!is.null(imd_full$weight_list), "weight_list exists")
  check(length(imd_full$weight_list) == length(dat_list),
        sprintf("weight_list has %d blocks", length(dat_list)))
  cat("  IMD weight distributions:\n")
  for (bn in names(dat_list)) {
    w <- imd_full$weight_list[[bn]]
    check(length(w) == ncol(dat_list[[bn]]),
          sprintf("  %s: weight length=%d matches ncol=%d", bn, length(w), ncol(dat_list[[bn]])))
    check(any(w > 0), sprintf("  %s: has non-zero weights", bn))
    report_weights(w, bn)
    report_top(w, bn)
  }
  cat(sprintf("  weight_list_init: %s\n",
              if (!is.null(imd_full$weight_list_init)) "present (per-tree)" else "NULL"))
}
cat(sprintf("  Total Test 1: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 2: Sub-MRF + compute_imd = TRUE
## ============================================================================
cat("--- Test 2: Sub-MRF + compute_imd = TRUE ---\n")
t0 <- proc.time()

sub_mods_imd <- fit_sub_multi_rfsrc(
  dat.list = dat_list,
  connect_list = connect_list,
  n_sub = 5,
  frac_response = 0.3,
  frac_predictor = 0.5,
  ntree_per_sub = 20,
  compute_imd = TRUE,
  yprob = 0.5,
  seed = 42,
  verbose = FALSE
)

cat(sprintf("  fit_sub_multi_rfsrc done [%s]\n", elapsed(t0)))

for (conn_name in names(sub_mods_imd)) {
  mod <- sub_mods_imd[[conn_name]]
  cat(sprintf("  Connection: %s\n", conn_name))
  report_model_structure(mod, conn_name)
  check(!is.null(mod$sub_mrf_info), sprintf("  %s: sub_mrf_info present", conn_name))
  check(!is.null(mod$imd_weights), sprintf("  %s: imd_weights present", conn_name))
  if (!is.null(mod$imd_weights)) {
    check(!is.null(mod$imd_weights$X), sprintf("  %s: imd_weights$X present", conn_name))
    check(!is.null(mod$imd_weights$Y), sprintf("  %s: imd_weights$Y present", conn_name))
    check(length(mod$imd_weights$X) == ncol(mod$xvar),
          sprintf("  %s: X weight length=%d matches xvar ncol=%d",
                  conn_name, length(mod$imd_weights$X), ncol(mod$xvar)))
    check(length(mod$imd_weights$Y) == ncol(mod$yvar),
          sprintf("  %s: Y weight length=%d matches yvar ncol=%d",
                  conn_name, length(mod$imd_weights$Y), ncol(mod$yvar)))
    check(any(mod$imd_weights$X > 0), sprintf("  %s: X has non-zero IMD", conn_name))
    check(any(mod$imd_weights$Y > 0), sprintf("  %s: Y has non-zero IMD", conn_name))
    cat("  IMD weight distributions (pre-computed):\n")
    report_weights(mod$imd_weights$X, paste0(conn_name, " X"))
    report_weights(mod$imd_weights$Y, paste0(conn_name, " Y"))
    report_top(mod$imd_weights$X, paste0(conn_name, " X"))
    report_top(mod$imd_weights$Y, paste0(conn_name, " Y"))
  }
}
cat(sprintf("  Total Test 2: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 3: Sub-MRF + compute_imd = FALSE + workflow guard
## ============================================================================
cat("--- Test 3: Sub-MRF without pre-computed IMD (should warn) ---\n")
t0 <- proc.time()

sub_mods_no_imd <- fit_sub_multi_rfsrc(
  dat.list = dat_list,
  connect_list = connect_list,
  n_sub = 3,
  frac_response = 0.3,
  frac_predictor = 0.5,
  ntree_per_sub = 10,
  compute_imd = FALSE,
  yprob = 0.5,
  seed = 42,
  verbose = FALSE
)

cat(sprintf("  fit_sub_multi_rfsrc done [%s]\n", elapsed(t0)))

for (conn_name in names(sub_mods_no_imd)) {
  mod <- sub_mods_no_imd[[conn_name]]
  check(is.null(mod$imd_weights), sprintf("  %s: imd_weights is NULL (as expected)", conn_name))
  cat(sprintf("    [%s] fields: %s\n", conn_name, paste(names(mod), collapse = ", ")))
}

## Simulate what the workflow does: detect sub-MRF without IMD and warn
warned <- FALSE
tryCatch(
  withCallingHandlers(
    {
      is_sub_mrf <- vapply(sub_mods_no_imd, function(m) !is.null(m$sub_mrf_info), logical(1))
      has_precomputed <- vapply(sub_mods_no_imd, function(m) !is.null(m$imd_weights), logical(1))
      if (any(is_sub_mrf & !has_precomputed)) {
        warning("Sub-MRF models detected without pre-computed IMD weights.", call. = FALSE)
      }
    },
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) fail(paste("Unexpected error:", e$message))
)
check(warned, "Warning issued for sub-MRF without pre-computed IMD")
cat(sprintf("  Total Test 3: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 4: Per-dimension threshold (mixed full/sub)
## ============================================================================
cat("--- Test 4: Per-dimension threshold (mixed full/sub per connection) ---\n")
t0 <- proc.time()

## miR has 50 features as response => below threshold 100 => use frac=1.0
## geneA has 500 as response => above threshold => use sub
mixed_mods <- fit_sub_multi_rfsrc(
  dat.list = dat_list,
  connect_list = connect_list,
  n_sub = 5,
  frac_response = 0.3,
  frac_predictor = 0.5,
  ntree_per_sub = 20,
  min_response_for_sub = 100,
  min_predictor_for_sub = 80,
  ntree_full = 100,
  compute_imd = TRUE,
  yprob = 0.5,
  seed = 42,
  verbose = TRUE
)

cat(sprintf("  fit_sub_multi_rfsrc done [%s]\n", elapsed(t0)))

for (conn_name in names(mixed_mods)) {
  mod <- mixed_mods[[conn_name]]
  has_sub <- !is.null(mod$sub_mrf_info)
  has_full_tree <- !is.null(mod$membership)
  cat(sprintf("  %s: sub_mrf=%s, full_forest=%s\n",
              conn_name, has_sub, has_full_tree))
  report_model_structure(mod, conn_name)

  ## For sub-MRF connections, check imd_weights
  if (has_sub) {
    check(!is.null(mod$imd_weights), sprintf("  %s: sub-MRF has imd_weights", conn_name))
    if (!is.null(mod$imd_weights)) {
      report_weights(mod$imd_weights$X, paste0(conn_name, " X"))
      report_weights(mod$imd_weights$Y, paste0(conn_name, " Y"))
    }
  }
  ## For full forest connections, check tree structure exists
  if (has_full_tree) {
    check(!is.null(mod$forest), sprintf("  %s: full forest has $forest", conn_name))
    cat(sprintf("    [%s] ntree=%d\n", conn_name, mod$ntree))
  }
}
cat(sprintf("  Total Test 4: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 5: mrf3_fit end-to-end with full forest + run_imd
## ============================================================================
cat("--- Test 5: mrf3_fit full forest + run_imd ---\n")
t0 <- proc.time()

fit_full <- tryCatch(
  mrf3_fit(
    dat.list = dat_list,
    ntree = 100,
    connect_list = connect_list,
    run_imd = TRUE,
    seed = 42
  ),
  error = function(e) { fail(paste("mrf3_fit error:", e$message)); NULL }
)

cat(sprintf("  mrf3_fit done [%s]\n", elapsed(t0)))

if (!is.null(fit_full)) {
  report_fit_structure(fit_full, "full forest")
  check(!is.null(fit_full$weights), "weights field present in output")
  if (!is.null(fit_full$weights)) {
    for (bn in names(dat_list)) {
      w <- fit_full$weights[[bn]]
      check(!is.null(w) && length(w) > 0,
            sprintf("  %s: weight present (len=%d)", bn, length(w)))
      report_weights(w, bn)
      report_top(w, bn)
    }
  }
}
cat(sprintf("  Total Test 5: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 6: mrf3_fit end-to-end with sub_mrf + run_imd (auto compute_imd)
## ============================================================================
cat("--- Test 6: mrf3_fit sub_mrf + run_imd (auto compute_imd injection) ---\n")
t0 <- proc.time()

fit_sub <- tryCatch(
  mrf3_fit(
    dat.list = dat_list,
    ntree = 100,
    connect_list = connect_list,
    run_imd = TRUE,
    sub_mrf = TRUE,
    sub_mrf_args = list(
      n_sub = 5,
      frac_response = 0.3,
      frac_predictor = 0.5,
      ntree_per_sub = 20
    ),
    seed = 42
  ),
  error = function(e) { fail(paste("mrf3_fit sub_mrf error:", e$message)); NULL }
)

cat(sprintf("  mrf3_fit done [%s]\n", elapsed(t0)))

if (!is.null(fit_sub)) {
  report_fit_structure(fit_sub, "sub-MRF")
  check(!is.null(fit_sub$weights), "weights field present in sub-MRF output")
  if (!is.null(fit_sub$weights)) {
    for (bn in names(dat_list)) {
      w <- fit_sub$weights[[bn]]
      check(!is.null(w) && length(w) > 0,
            sprintf("  %s: weight present (len=%d)", bn, length(w)))
      if (!is.null(w)) {
        check(any(w > 0), sprintf("  %s: has non-zero weights", bn))
        report_weights(w, bn)
        report_top(w, bn)
      }
    }
  }

  ## Check models are actually sub-MRF
  n_sub_models <- sum(vapply(fit_sub$models, function(m) !is.null(m$sub_mrf_info), logical(1)))
  cat(sprintf("  Sub-MRF models: %d / %d\n", n_sub_models, length(fit_sub$models)))
  check(n_sub_models > 0, "At least one model is sub-MRF")
}
cat(sprintf("  Total Test 6: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 7: IMD weight quality check (sub-MRF vs full forest)
## ============================================================================
cat("--- Test 7: IMD weight quality (sub-MRF vs full forest correlation) ---\n")
t0 <- proc.time()

if (!is.null(fit_full) && !is.null(fit_sub) &&
    !is.null(fit_full$weights) && !is.null(fit_sub$weights)) {

  for (bn in names(dat_list)) {
    w_full <- fit_full$weights[[bn]]
    w_sub  <- fit_sub$weights[[bn]]
    if (!is.null(w_full) && !is.null(w_sub)) {
      ## Align names
      common <- intersect(names(w_full), names(w_sub))
      if (length(common) > 10) {
        r <- cor(w_full[common], w_sub[common], method = "spearman")
        cat(sprintf("  %s: Spearman cor = %.3f (n=%d features)\n", bn, r, length(common)))
        check(r > 0.3, sprintf("  %s: reasonable IMD agreement (rho=%.3f > 0.3)", bn, r))

        ## Additional diagnostics: compare top features
        top_full <- names(sort(w_full[common], decreasing = TRUE))[1:min(20, length(common))]
        top_sub  <- names(sort(w_sub[common], decreasing = TRUE))[1:min(20, length(common))]
        overlap_top <- length(intersect(top_full, top_sub))
        cat(sprintf("  %s: top-20 overlap = %d/20 (%.0f%%)\n", bn, overlap_top, 100 * overlap_top / 20))

        ## Signal vs noise check: known signal features should rank higher
        ## Signal features for each block:
        ## geneA: geneA_1..geneA_90, cpg: cpg_1..cpg_45, miR: miR_1..miR_15
        n_signal <- if (bn == "geneA") 90 else if (bn == "cpg") 45 else 15
        signal_names <- paste0(sub("^(\\w+)$", "\\1", colnames(dat_list[[bn]])[1:n_signal]))
        signal_names <- colnames(dat_list[[bn]])[1:n_signal]
        if (all(signal_names %in% common)) {
          w_signal_full <- mean(w_full[signal_names])
          w_noise_full  <- mean(w_full[setdiff(common, signal_names)])
          w_signal_sub  <- mean(w_sub[signal_names])
          w_noise_sub   <- mean(w_sub[setdiff(common, signal_names)])
          cat(sprintf("  %s (full):  signal_mean=%.4f, noise_mean=%.4f, ratio=%.2f\n",
                      bn, w_signal_full, w_noise_full,
                      if (w_noise_full > 0) w_signal_full / w_noise_full else Inf))
          cat(sprintf("  %s (sub):   signal_mean=%.4f, noise_mean=%.4f, ratio=%.2f\n",
                      bn, w_signal_sub, w_noise_sub,
                      if (w_noise_sub > 0) w_signal_sub / w_noise_sub else Inf))
        }
      }
    }
  }
} else {
  cat("  Skipped (missing fit_full or fit_sub)\n")
}
cat(sprintf("  Total Test 7: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 8: Sub-MRF + compute_imd=TRUE with parallel=FALSE inside fit_one
## ============================================================================
cat("--- Test 8: Single connection fit_sub_mrf with compute_imd ---\n")
t0 <- proc.time()

single_mod <- tryCatch(
  fit_sub_mrf(
    X = dat_list$cpg,
    Y = dat_list$geneA,
    n_sub = 3,
    frac_response = 0.3,
    frac_predictor = 0.5,
    ntree_per_sub = 10,
    compute_imd = TRUE,
    seed = 42,
    verbose = FALSE
  ),
  error = function(e) { fail(paste("fit_sub_mrf error:", e$message)); NULL }
)

cat(sprintf("  fit_sub_mrf done [%s]\n", elapsed(t0)))

if (!is.null(single_mod)) {
  cat(sprintf("  Model fields: %s\n", paste(names(single_mod), collapse = ", ")))
  report_model_structure(single_mod, "cpg->geneA")

  check(!is.null(single_mod$imd_weights), "imd_weights present in single connection")
  check(!is.null(single_mod$imd_weights$X), "X weights present")
  check(!is.null(single_mod$imd_weights$Y), "Y weights present")
  check(length(single_mod$imd_weights$X) == ncol(dat_list$cpg),
        sprintf("X weights length=%d matches predictor ncol=%d",
                length(single_mod$imd_weights$X), ncol(dat_list$cpg)))
  check(length(single_mod$imd_weights$Y) == ncol(dat_list$geneA),
        sprintf("Y weights length=%d matches response ncol=%d",
                length(single_mod$imd_weights$Y), ncol(dat_list$geneA)))

  if (!is.null(single_mod$imd_weights)) {
    report_weights(single_mod$imd_weights$X, "X (predictor=cpg)")
    report_weights(single_mod$imd_weights$Y, "Y (response=geneA)")
    report_top(single_mod$imd_weights$X, "X top")
    report_top(single_mod$imd_weights$Y, "Y top")
  }

  ## Check coverage: features not included in any sub-model should have 0 weight
  cov_x <- single_mod$sub_mrf_info$col_coverage_predictor
  cov_y <- single_mod$sub_mrf_info$col_coverage_response
  uncov_x <- names(cov_x[cov_x == 0])
  uncov_y <- names(cov_y[cov_y == 0])
  cat(sprintf("  Coverage X: %d/%d covered | Coverage Y: %d/%d covered\n",
              sum(cov_x > 0), length(cov_x), sum(cov_y > 0), length(cov_y)))
  if (length(uncov_x) > 0) {
    cat(sprintf("  Uncovered X features: %d\n", length(uncov_x)))
    check(all(single_mod$imd_weights$X[uncov_x] == 0),
          "Uncovered X features have zero weight")
  } else {
    pass("All X features covered by at least one sub-model")
  }
  if (length(uncov_y) > 0) {
    cat(sprintf("  Uncovered Y features: %d\n", length(uncov_y)))
    check(all(single_mod$imd_weights$Y[uncov_y] == 0),
          "Uncovered Y features have zero weight")
  } else {
    pass("All Y features covered by at least one sub-model")
  }
}
cat(sprintf("  Total Test 8: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 9: mrf3_fit full forest + variable selection
## ============================================================================
cat("--- Test 9: mrf3_fit full forest + run_variable_selection ---\n")
t0 <- proc.time()

fit_full_vs <- tryCatch(
  mrf3_fit(
    dat.list = dat_list,
    ntree = 100,
    connect_list = connect_list,
    run_imd = TRUE,
    run_variable_selection = TRUE,
    seed = 42
  ),
  error = function(e) { fail(paste("mrf3_fit VS error:", e$message)); NULL }
)

cat(sprintf("  mrf3_fit done [%s]\n", elapsed(t0)))

if (!is.null(fit_full_vs)) {
  report_fit_structure(fit_full_vs, "full forest + VS")
  check(!is.null(fit_full_vs$weights), "weights present")
  check(!is.null(fit_full_vs$selected_vars), "selected_vars present")
  check(!is.null(fit_full_vs$vs_summary), "vs_summary present")
  if (!is.null(fit_full_vs$vs_summary)) {
    cat("  Variable selection summary:\n")
    for (i in seq_len(nrow(fit_full_vs$vs_summary))) {
      row <- fit_full_vs$vs_summary[i, ]
      cat(sprintf("    %s: %d / %d selected (%.1f%%)\n",
                  row$block, row$p_selected, row$p_total, row$selected_ratio * 100))
    }
    check(all(fit_full_vs$vs_summary$p_selected > 0),
          "All blocks have selected variables")
    check(all(fit_full_vs$vs_summary$p_selected <= fit_full_vs$vs_summary$p_total),
          "Selected <= total for all blocks")
  }
  check(!is.null(fit_full_vs$selected_data), "selected_data present")
  if (!is.null(fit_full_vs$selected_data)) {
    for (bn in names(fit_full_vs$selected_data)) {
      d <- fit_full_vs$selected_data[[bn]]
      if (!is.null(d)) {
        check(nrow(d) == n, sprintf("  %s: selected data has correct n=%d", bn, n))
        check(ncol(d) <= ncol(dat_list[[bn]]),
              sprintf("  %s: selected ncol=%d <= original=%d", bn, ncol(d), ncol(dat_list[[bn]])))
      }
    }
  }
  ## Report selected weights distribution (post-VS)
  if (!is.null(fit_full_vs$weights)) {
    cat("  Post-VS weight distributions:\n")
    for (bn in names(fit_full_vs$weights)) {
      report_weights(fit_full_vs$weights[[bn]], paste0(bn, " post-VS"))
    }
  }
}
cat(sprintf("  Total Test 9: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 10: mrf3_fit sub-MRF + variable selection
## ============================================================================
cat("--- Test 10: mrf3_fit sub_mrf + run_variable_selection ---\n")
t0 <- proc.time()

fit_sub_vs <- tryCatch(
  mrf3_fit(
    dat.list = dat_list,
    ntree = 100,
    connect_list = connect_list,
    run_imd = TRUE,
    run_variable_selection = TRUE,
    sub_mrf = TRUE,
    sub_mrf_args = list(
      n_sub = 5,
      frac_response = 0.3,
      frac_predictor = 0.5,
      ntree_per_sub = 20
    ),
    seed = 42
  ),
  error = function(e) { fail(paste("mrf3_fit sub_mrf VS error:", e$message)); NULL }
)

cat(sprintf("  mrf3_fit done [%s]\n", elapsed(t0)))

if (!is.null(fit_sub_vs)) {
  report_fit_structure(fit_sub_vs, "sub-MRF + VS")
  check(!is.null(fit_sub_vs$weights), "weights present (sub-MRF + VS)")
  check(!is.null(fit_sub_vs$selected_vars), "selected_vars present (sub-MRF + VS)")
  check(!is.null(fit_sub_vs$vs_summary), "vs_summary present (sub-MRF + VS)")
  if (!is.null(fit_sub_vs$vs_summary)) {
    cat("  Variable selection summary (sub-MRF):\n")
    for (i in seq_len(nrow(fit_sub_vs$vs_summary))) {
      row <- fit_sub_vs$vs_summary[i, ]
      cat(sprintf("    %s: %d / %d selected (%.1f%%)\n",
                  row$block, row$p_selected, row$p_total, row$selected_ratio * 100))
    }
    check(all(fit_sub_vs$vs_summary$p_selected > 0),
          "All blocks have selected variables (sub-MRF)")
  }
  ## Report post-VS weights
  if (!is.null(fit_sub_vs$weights)) {
    cat("  Post-VS weight distributions (sub-MRF):\n")
    for (bn in names(fit_sub_vs$weights)) {
      report_weights(fit_sub_vs$weights[[bn]], paste0(bn, " post-VS"))
    }
  }
}
cat(sprintf("  Total Test 10: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 11: mrf3_fit full forest + robust clustering
## ============================================================================
cat("--- Test 11: mrf3_fit full forest + run_robust_clustering ---\n")
t0 <- proc.time()

fit_full_robust <- tryCatch(
  mrf3_fit(
    dat.list = dat_list,
    ntree = 100,
    connect_list = connect_list,
    run_imd = TRUE,
    run_variable_selection = TRUE,
    run_robust_clustering = TRUE,
    seed = 42
  ),
  error = function(e) { fail(paste("mrf3_fit robust error:", e$message)); NULL }
)

cat(sprintf("  mrf3_fit done [%s]\n", elapsed(t0)))

if (!is.null(fit_full_robust)) {
  report_fit_structure(fit_full_robust, "full forest + robust")
  check(!is.null(fit_full_robust$weights), "weights present (robust)")
  check(!is.null(fit_full_robust$selected_vars), "selected_vars present (robust)")
  check(!is.null(fit_full_robust$robust_clusters), "robust_clusters present")
  if (!is.null(fit_full_robust$robust_clusters)) {
    tbl <- table(fit_full_robust$robust_clusters)
    n_cl <- length(tbl)
    cat(sprintf("  Robust clusters: %d unique labels, sizes: %s\n",
                n_cl, paste(sprintf("%s=%d", names(tbl), tbl), collapse = ", ")))
    check(n_cl >= 2, sprintf("  At least 2 robust clusters (got %d)", n_cl))
    check(length(fit_full_robust$robust_clusters) == n,
          sprintf("  Cluster vector length=%d matches n=%d",
                  length(fit_full_robust$robust_clusters), n))
    ## Compare with true clusters
    if (requireNamespace("mclust", quietly = TRUE)) {
      ari <- mclust::adjustedRandIndex(fit_full_robust$robust_clusters, cluster_true)
      cat(sprintf("  ARI vs true clusters: %.3f\n", ari))
    }
  }
  check(!is.null(fit_full_robust$robust_detail), "robust_detail present")
  if (!is.null(fit_full_robust$robust_detail)) {
    check(!is.null(fit_full_robust$robust_detail$shared), "robust shared branch present")
    cat(sprintf("  robust_detail fields: %s\n",
                paste(names(fit_full_robust$robust_detail), collapse = ", ")))
  }
}
cat(sprintf("  Total Test 11: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 12: mrf3_fit sub-MRF + robust clustering
## ============================================================================
cat("--- Test 12: mrf3_fit sub_mrf + run_robust_clustering ---\n")
t0 <- proc.time()

fit_sub_robust <- tryCatch(
  mrf3_fit(
    dat.list = dat_list,
    ntree = 100,
    connect_list = connect_list,
    run_imd = TRUE,
    run_variable_selection = TRUE,
    run_robust_clustering = TRUE,
    sub_mrf = TRUE,
    sub_mrf_args = list(
      n_sub = 5,
      frac_response = 0.3,
      frac_predictor = 0.5,
      ntree_per_sub = 20
    ),
    seed = 42
  ),
  error = function(e) { fail(paste("mrf3_fit sub_mrf robust error:", e$message)); NULL }
)

cat(sprintf("  mrf3_fit done [%s]\n", elapsed(t0)))

if (!is.null(fit_sub_robust)) {
  report_fit_structure(fit_sub_robust, "sub-MRF + robust")
  check(!is.null(fit_sub_robust$weights), "weights present (sub-MRF + robust)")
  check(!is.null(fit_sub_robust$selected_vars), "selected_vars present (sub-MRF + robust)")
  check(!is.null(fit_sub_robust$robust_clusters), "robust_clusters present (sub-MRF)")
  if (!is.null(fit_sub_robust$robust_clusters)) {
    tbl <- table(fit_sub_robust$robust_clusters)
    n_cl <- length(tbl)
    cat(sprintf("  Robust clusters (sub-MRF): %d unique labels, sizes: %s\n",
                n_cl, paste(sprintf("%s=%d", names(tbl), tbl), collapse = ", ")))
    check(n_cl >= 2, sprintf("  At least 2 robust clusters (got %d)", n_cl))
    ## Compare with true clusters
    if (requireNamespace("mclust", quietly = TRUE)) {
      ari <- mclust::adjustedRandIndex(fit_sub_robust$robust_clusters, cluster_true)
      cat(sprintf("  ARI vs true clusters: %.3f\n", ari))
    }
  }
  check(!is.null(fit_sub_robust$robust_detail), "robust_detail present (sub-MRF)")
  if (!is.null(fit_sub_robust$robust_detail)) {
    cat(sprintf("  robust_detail fields: %s\n",
                paste(names(fit_sub_robust$robust_detail), collapse = ", ")))
  }
}
cat(sprintf("  Total Test 12: %s\n\n", elapsed(t0)))


## ============================================================================
## Test 13: Variable selection quality (sub-MRF vs full forest agreement)
## ============================================================================
cat("--- Test 13: Variable selection agreement (sub-MRF vs full forest) ---\n")
t0 <- proc.time()

if (!is.null(fit_full_vs) && !is.null(fit_sub_vs) &&
    !is.null(fit_full_vs$selected_vars) && !is.null(fit_sub_vs$selected_vars)) {
  for (bn in names(dat_list)) {
    vars_full <- fit_full_vs$selected_vars[[bn]]
    vars_sub  <- fit_sub_vs$selected_vars[[bn]]
    if (!is.null(vars_full) && !is.null(vars_sub)) {
      overlap <- length(intersect(vars_full, vars_sub))
      union_n <- length(union(vars_full, vars_sub))
      jaccard <- if (union_n > 0) overlap / union_n else 0
      cat(sprintf("  %s: full=%d vars, sub=%d vars, overlap=%d, Jaccard=%.3f\n",
                  bn, length(vars_full), length(vars_sub), overlap, jaccard))
      check(jaccard > 0.1,
            sprintf("  %s: reasonable variable selection overlap (J=%.3f > 0.1)", bn, jaccard))

      ## Check how many known signal features are captured
      n_signal <- if (bn == "geneA") 90 else if (bn == "cpg") 45 else 15
      signal_names <- colnames(dat_list[[bn]])[1:n_signal]
      sig_full <- length(intersect(vars_full, signal_names))
      sig_sub  <- length(intersect(vars_sub, signal_names))
      cat(sprintf("  %s: signal captured: full=%d/%d (%.0f%%), sub=%d/%d (%.0f%%)\n",
                  bn, sig_full, n_signal, 100 * sig_full / n_signal,
                  sig_sub, n_signal, 100 * sig_sub / n_signal))
    }
  }
} else {
  cat("  Skipped (missing fit_full_vs or fit_sub_vs)\n")
}
cat(sprintf("  Total Test 13: %s\n\n", elapsed(t0)))


## ============================================================================
## Summary
## ============================================================================
total_elapsed <- proc.time() - TOTAL_START
cat("=== Test Suite Complete ===\n")
cat(sprintf("Total time: %.1fs\n", total_elapsed["elapsed"]))
if (FAILURES == 0L) {
  cat("All tests PASSED!\n")
} else {
  cat(sprintf("%d test(s) FAILED.\n", FAILURES))
}
