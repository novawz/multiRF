#' Simulate Shared + Specific Multi-Omics Data via InterSIM
#'
#' Generate multi-omics data where shared signal is provided by
#' `InterSIM::InterSIM()` and additional omics-specific signal is injected by
#' per-omics latent labels `U`.
#'
#' This function wraps the notebook simulation workflow into a reusable package
#' API and optionally returns workflow-ready `dat.list` and truth labels.
#'
#' @param n Number of samples.
#' @param Z_prop Cluster proportion for InterSIM shared clusters.
#' @param delta_methyl Shared-signal strength for methylation in InterSIM.
#' @param delta_expr Shared-signal strength for expression in InterSIM.
#' @param delta_protein Shared-signal strength for protein in InterSIM.
#' @param p_DMP Proportion of differential methylation features in InterSIM.
#' Default matches `InterSIM::InterSIM()` (`0.2`).
#' @param p_DEG Proportion of differential expression features in InterSIM.
#' Default matches `InterSIM::InterSIM()` (`NULL`).
#' @param p_DEP Proportion of differential protein features in InterSIM.
#' Default matches `InterSIM::InterSIM()` (`NULL`).
#' @param p_spec_feat Proportion of features used as specific set `T` per omics.
#' @param delta_spec Named numeric vector with specific shift size for
#' `methyl`, `gene`, and `protein`.
#' @param U_nclass Number of categories for each specific latent label `U`.
#' @param U_within_Z Logical; whether each `U` is generated within each shared
#' cluster `Z` (recommended).
#' @param max_features_per_block Maximum number of features kept per omics block
#' in returned `dat.list` when `return_workflow_data = TRUE`.
#' @param sample_prefix Prefix for generated sample IDs in `dat.list`.
#' @param return_workflow_data Logical; whether to return workflow-ready
#' `dat.list`, `truth`, and `k_ref`.
#' @param keep_sim_raw Logical; whether to keep original InterSIM output as
#' `sim_raw`.
#' @param seed Optional random seed. If `NULL`, current RNG state is used.
#'
#' @return A list containing:
#' - `X_list`: simulated matrices (`methyl`, `gene`, `protein`)
#' - `Z`: shared labels from InterSIM
#' - `U`: specific labels per omics
#' - `S`: shared feature index sets per omics
#' - `T`: specific feature index sets per omics
#' - `sim_raw`: raw InterSIM output (optional)
#' - `dat.list`, `truth`, `k_ref`: workflow-ready outputs (optional)
#'
simulate_intersim_shared_specific <- function(
    n = 500,
    Z_prop = c(0.3, 0.3, 0.4),
    delta_methyl = 2,
    delta_expr = 2,
    delta_protein = 2,
    p_DMP = 0.2,
    p_DEG = NULL,
    p_DEP = NULL,
    p_spec_feat = 0.05,
    delta_spec = c(methyl = 1, gene = 1, protein = 1),
    U_nclass = 2,
    U_within_Z = TRUE,
    max_features_per_block = 400,
    sample_prefix = "S",
    return_workflow_data = TRUE,
    keep_sim_raw = TRUE,
    seed = NULL) {

  if (!requireNamespace("InterSIM", quietly = TRUE)) {
    stop("Package `InterSIM` is required. Please install it first.")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("`seed` must be NULL or a single finite numeric value.")
    }
    set.seed(as.integer(seed))
  }

  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2) {
    stop("`n` must be a single integer >= 2.")
  }
  n <- as.integer(n)

  if (!is.numeric(U_nclass) || length(U_nclass) != 1L || !is.finite(U_nclass) || U_nclass < 2) {
    stop("`U_nclass` must be a single integer >= 2.")
  }
  U_nclass <- as.integer(U_nclass)

  if (!is.numeric(max_features_per_block) || length(max_features_per_block) != 1L ||
      !is.finite(max_features_per_block) || max_features_per_block < 1) {
    stop("`max_features_per_block` must be a single integer >= 1.")
  }
  max_features_per_block <- as.integer(max_features_per_block)

  if (!is.numeric(p_spec_feat) || length(p_spec_feat) != 1L || p_spec_feat < 0 || p_spec_feat > 1) {
    stop("`p_spec_feat` must be in [0, 1].")
  }
  check_prob <- function(x, name) {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
      stop("`", name, "` must be NULL or a numeric scalar in [0, 1].")
    }
    as.numeric(x)
  }
  p_DMP <- check_prob(p_DMP, "p_DMP")
  p_DEG <- check_prob(p_DEG, "p_DEG")
  p_DEP <- check_prob(p_DEP, "p_DEP")

  parse_three_named <- function(x, arg_name) {
    target <- c("methyl", "gene", "protein")
    if (!is.numeric(x)) {
      stop("`", arg_name, "` must be numeric.")
    }
    if (length(x) == 1L) {
      out <- rep(as.numeric(x), 3L)
      names(out) <- target
      return(out)
    }
    if (is.null(names(x))) {
      if (length(x) != 3L) {
        stop("`", arg_name, "` must have length 1 or 3, or be named by methyl/gene/protein.")
      }
      names(x) <- target
      return(x[target])
    }
    miss <- setdiff(target, names(x))
    if (length(miss) > 0L) {
      stop("`", arg_name, "` missing names: ", paste(miss, collapse = ", "))
    }
    as.numeric(x[target]) |> stats::setNames(target)
  }

  delta_spec <- parse_three_named(delta_spec, "delta_spec")

  get_sim_mats <- function(sim) {
    cand <- list(
      methyl  = c("methylation", "methyl", "X.methyl", "X_methyl", "dat.methyl"),
      gene    = c("gene", "expr", "expression", "X.expr", "X_expression", "dat.expr"),
      protein = c("protein", "prot", "X.protein", "X_protein", "dat.protein")
    )
    pick <- function(keys) {
      for (k in keys) {
        if (!is.null(sim[[k]])) {
          return(sim[[k]])
        }
      }
      NULL
    }
    Xm <- pick(cand$methyl)
    Xg <- pick(cand$gene)
    Xp <- pick(cand$protein)
    if (is.null(Xm) || is.null(Xg) || is.null(Xp)) {
      stop("Cannot find methyl/gene/protein matrices in InterSIM output.")
    }
    list(methyl = as.matrix(Xm), gene = as.matrix(Xg), protein = as.matrix(Xp))
  }

  coerce_samples_in_rows <- function(mat, n_samples) {
    x <- as.matrix(mat)
    if (nrow(x) == n_samples) {
      return(x)
    }
    if (ncol(x) == n_samples) {
      return(t(x))
    }
    stop("Cannot align matrix orientation: neither nrow nor ncol equals n_samples.")
  }

  make_U_within_Z <- function(Z, nclass) {
    Z <- as.factor(Z)
    U <- integer(length(Z))
    for (lvl in levels(Z)) {
      idx <- which(Z == lvl)
      idx <- sample(idx)
      cut_id <- cut(seq_along(idx), breaks = nclass, labels = FALSE)
      U[idx] <- as.integer(cut_id)
    }
    U
  }

  add_shift_by_label <- function(X, U, feat_idx, delta = 1, type = "methyl", eps = 1e-6, seed = NULL) {
    if (!is.null(seed)) {
      set.seed(as.integer(seed))
    }

    X2 <- as.matrix(X)
    type <- match.arg(type, c("methyl", "gene", "protein"))
    n <- nrow(X2)
    p <- ncol(X2)

    if (length(U) != n) {
      stop("`U` must have the same length as `nrow(X)`.")
    }

    feat_idx <- as.integer(feat_idx)
    feat_idx <- feat_idx[is.finite(feat_idx)]
    feat_idx <- unique(feat_idx[feat_idx >= 1L & feat_idx <= p])
    if (length(feat_idx) == 0L) {
      return(X2)
    }

    if (type == "methyl") {
      logit <- function(pv) log(pv / (1 - pv))
      inv_logit <- function(z) 1 / (1 + exp(-z))
      latent <- logit(pmin(pmax(X2, eps), 1 - eps))
      inverse_transform <- function(M) inv_logit(M)
    } else {
      X2_center <- colMeans(X2, na.rm = TRUE)
      X2_scale <- apply(X2, 2, stats::sd, na.rm = TRUE)
      X2_scale[!is.finite(X2_scale) | X2_scale <= 0] <- 1
      latent <- sweep(X2, 2, X2_center, FUN = "-")
      latent <- sweep(latent, 2, X2_scale, FUN = "/")
      inverse_transform <- function(M) {
        out <- sweep(M, 2, X2_scale, FUN = "*")
        sweep(out, 2, X2_center, FUN = "+")
      }
    } 

    delta_norm <- rnorm(length(feat_idx), mean = 0, sd = delta)

    if (U_nclass == 2L) {
      row_coef <- as.numeric(U == 2L)
    } else {
      lv <- sort(unique(U))
      lv <- lv[!is.na(lv)]
      if (length(lv) <= 1L) {
        return(X2)
      }
      row_coef <- rep(0, n)
      offset <- seq(-0.5, 0.5, length.out = length(lv))
      for (ii in seq_along(lv)) {
        row_coef[U == lv[ii]] <- offset[ii]
      }
    }

    shift_mat <- outer(row_coef, delta_norm)
    latent[, feat_idx] <- latent[, feat_idx, drop = FALSE] + shift_mat
    inverse_transform(latent)
  }

  pick_truth_from_sim <- function(sim_obj) {
    if (!"clustering.assignment" %in% names(sim_obj)) {
      stop("InterSIM output missing `clustering.assignment`.")
    }
    tb <- sim_obj$clustering.assignment
    cand <- c("cluster.id", "cluster_id", "cluster")
    cand <- cand[cand %in% names(tb)]
    if (length(cand) == 0L) {
      stop("Cannot find cluster label column in `clustering.assignment`.")
    }
    factor(tb[[cand[1]]])
  }

  # Force exact sum-to-1 to satisfy InterSIM's strict check
  Z_prop <- Z_prop / sum(Z_prop)

  sim <- InterSIM::InterSIM(
    n.sample = n,
    cluster.sample.prop = Z_prop,
    delta.methyl = delta_methyl,
    delta.expr = delta_expr,
    delta.protein = delta_protein,
    p.DMP = p_DMP,
    p.DEG = p_DEG,
    p.DEP = p_DEP
  )

  X <- get_sim_mats(sim)
  n_samples <- nrow(X$methyl)
  Z <- sim$clustering.assignment$cluster.id

  make_U <- function() {
    if (isTRUE(U_within_Z)) {
      make_U_within_Z(Z, nclass = U_nclass)
    } else {
      sample(rep(seq_len(U_nclass), length.out = n_samples))
    }
  }
  U <- list(
    methyl = make_U(),
    gene = make_U(),
    protein = make_U()
  )

  pick_feats <- function(p, frac) {
    m <- ceiling(p * frac)
    if (m <= 0L) {
      return(integer(0))
    }
    sample.int(p, size = m, replace = FALSE)
  }

  shared_frac <- list(
    methyl = p_DMP,
    gene = if (is.null(p_DEG)) p_DMP else p_DEG,
    protein = if (is.null(p_DEP)) p_DMP else p_DEP
  )

  S <- list(
    methyl = pick_feats(ncol(X$methyl), shared_frac$methyl),
    gene = pick_feats(ncol(X$gene), shared_frac$gene),
    protein = pick_feats(ncol(X$protein), shared_frac$protein)
  )

  pick_T <- function(p, frac, avoid_idx) {
    m <- ceiling(p * frac)
    if (m <= 0L) {
      return(integer(0))
    }
    pool <- setdiff(seq_len(p), avoid_idx)
    if (length(pool) == 0L) {
      return(integer(0))
    }
    sample(pool, size = min(m, length(pool)), replace = FALSE)
  }

  T <- list(
    methyl = pick_T(ncol(X$methyl), p_spec_feat, S$methyl),
    gene = pick_T(ncol(X$gene), p_spec_feat, S$gene),
    protein = pick_T(ncol(X$protein), p_spec_feat, S$protein)
  )

  Xm2 <- add_shift_by_label(X$methyl, U$methyl, T$methyl, delta = delta_spec[["methyl"]], type = "methyl")
  Xg2 <- add_shift_by_label(X$gene, U$gene, T$gene, delta = delta_spec[["gene"]], type = "gene")
  Xp2 <- add_shift_by_label(X$protein, U$protein, T$protein, delta = delta_spec[["protein"]], type = "protein")

  out <- list(
    X_list = list(methyl = Xm2, gene = Xg2, protein = Xp2),
    Z = Z,
    U = U,
    S = S,
    T = T
  )
  if (isTRUE(keep_sim_raw)) {
    out$sim_raw <- sim
  }

  if (isTRUE(return_workflow_data)) {
    truth <- pick_truth_from_sim(sim)
    raw_blocks <- out$X_list[c("gene", "methyl", "protein")]
    dat_sim <- lapply(raw_blocks, function(x) {
      x2 <- coerce_samples_in_rows(x, n_samples = length(truth))
      x2 <- x2[, seq_len(min(ncol(x2), max_features_per_block)), drop = FALSE]
      as.data.frame(x2, check.names = FALSE)
    })
    names(dat_sim) <- c("gene", "methyl", "protein")

    sample_ids <- sprintf("%s%03d", sample_prefix, seq_len(length(truth)))
    for (nm in names(dat_sim)) {
      rownames(dat_sim[[nm]]) <- sample_ids
    }

    out$dat.list <- dat_sim
    out$truth <- truth
    out$k_ref <- nlevels(truth)
  }

  out
}


#' Backward-compatible alias of `simulate_intersim_shared_specific()`
#'
#' @param ... Arguments passed to `simulate_intersim_shared_specific()`.
#' @return Same as `simulate_intersim_shared_specific()`.
simulate_shared_specific_intersim <- function(...) {
  simulate_intersim_shared_specific(...)
}
