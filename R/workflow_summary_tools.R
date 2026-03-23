# is an mrf3 fit-like output
is_mrf3_fit_object <- function(x) {
  inherits(x, "mrf3_fit") ||
    (is.list(x) && !is.null(x$reconstruction) && !is.null(x$models))
}

mrf3_sample_names <- function(x) {
  if (!is.null(x$reconstruction$sample_names)) {
    return(as.character(x$reconstruction$sample_names))
  }
  if (!is.null(x$data) && is.list(x$data) && length(x$data) > 0L) {
    rn <- rownames(x$data[[1]])
    if (!is.null(rn)) {
      return(as.character(rn))
    }
  }
  if (!is.null(x$shared$clustering$similarity)) {
    rn <- rownames(x$shared$clustering$similarity)
    if (!is.null(rn)) {
      return(as.character(rn))
    }
  }
  stop("Cannot infer sample names from `x`.")
}

align_labels_to_samples <- function(labels, sample_names) {
  labels <- as.character(labels)
  if (!is.null(names(labels)) && all(sample_names %in% names(labels))) {
    return(labels[sample_names])
  }
  if (length(labels) != length(sample_names)) {
    stop("`labels` could not be aligned to `sample_names`.")
  }
  stats::setNames(labels, sample_names)
}

cluster_subset <- function(mat, idx, k, method = c("auto", "PAM", "Spectral"), kind = c("similarity", "feature")) {
  method <- match.arg(method)
  kind <- match.arg(kind)
  x_sub <- as.matrix(mat)[idx, , drop = FALSE]
  if (identical(kind, "similarity")) {
    x_sub <- x_sub[, idx, drop = FALSE]
  }
  method_use <- if (identical(method, "auto")) "PAM" else method
  if (identical(kind, "similarity")) {
    fit <- cluster_similarity_matrix(S = x_sub, k = k, method = method_use)
  } else {
    s <- stats::cor(t(x_sub), use = "pairwise.complete.obs")
    s[!is.finite(s)] <- 0
    diag(s) <- 1
    fit <- cluster_similarity_matrix(S = s, k = k, method = method_use)
  }
  list(cl = fit$cl, method = fit$method)
}

update_coassign_mats <- function(co_num, co_den, idx, cl_sub) {
  idx <- as.integer(idx)
  same <- outer(cl_sub, cl_sub, "==") * 1
  co_num[idx, idx] <- co_num[idx, idx] + same
  co_den[idx, idx] <- co_den[idx, idx] + 1
  list(co_num = co_num, co_den = co_den)
}

cluster_metric_pair <- function(pred, ref) {
  list(
    ari = cluster_ari(pred, ref),
    jaccard = cluster_jaccard(pred, ref),
    nmi = cluster_nmi(pred, ref)
  )
}

cluster_consensus_matrix <- function(co_mat, k, method = c("PAM", "Spectral")) {
  method <- match.arg(method)
  fit <- cluster_similarity_matrix(S = co_mat, k = k, method = method)
  list(cl = fit$cl, method = fit$method)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  mean(x)
}


#' mrf3 Stability Evaluation (No Refit)
#'
#' Evaluate clustering stability by repeated subsampling/bootstrapping on
#' existing `mrf3_fit` outputs, without refitting forests.
#'
#' @param x An object from `mrf3_fit()`.
#' @param branches Branches to evaluate. Supported:
#' `"specific_shared"` and `"robust_clustering"`.
#' @param n_rep Number of repeats.
#' @param sample_frac Sampling fraction per repeat.
#' @param sample_mode Sampling mode: `"subsample"` or `"bootstrap"`.
#' @param specific_k Optional `k` override for `"specific_shared"` branch.
#' @param robust_k Optional `k` override for `"robust_clustering"` branch.
#' @param specific_method Clustering backend for `"specific_shared"` repeats:
#' `"auto"`, `"Spectral"`, or `"PAM"`.
#' @param robust_method Clustering backend for `"robust_clustering"` repeats:
#' `"auto"`, `"Spectral"`, or `"PAM"`.
#' @param run_consensus Logical; whether to derive consensus clustering.
#' @param consensus_method Consensus clustering backend: `"Spectral"` or `"PAM"`.
#' @param truth Optional truth labels for external metrics.
#' If names are provided, they are aligned to sample names.
#' @param seed Random seed.
#' @param verbose Logical; whether to print progress.
#'
#' @return A list with per-branch repeat metrics and optional consensus results.
#' @export
mrf3_stability <- function(x,
                               branches = c("specific_shared", "robust_clustering"),
                               n_rep = 50L,
                               sample_frac = 0.8,
                               sample_mode = c("subsample", "bootstrap"),
                               specific_k = NULL,
                               robust_k = NULL,
                               specific_method = c("auto", "PAM", "Spectral"),
                               robust_method = c("auto", "PAM", "Spectral"),
                               run_consensus = TRUE,
                               consensus_method = c("PAM", "Spectral"),
                               truth = NULL,
                               seed = 529,
                               verbose = TRUE) {
  if (!is_mrf3_fit_object(x)) {
    stop("`x` must be an `mrf3_fit` object.")
  }
  if (!is.numeric(n_rep) || length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 1) {
    stop("`n_rep` must be a single integer >= 1.")
  }
  n_rep <- as.integer(n_rep)
  if (!is.numeric(sample_frac) || length(sample_frac) != 1L || !is.finite(sample_frac) ||
      sample_frac <= 0 || sample_frac > 1) {
    stop("`sample_frac` must be in (0, 1].")
  }
  sample_mode <- match.arg(sample_mode)
  specific_method <- match.arg(specific_method)
  robust_method <- match.arg(robust_method)
  consensus_method <- match.arg(consensus_method)
  branches <- unique(as.character(branches))
  valid_branches <- c("specific_shared", "robust_clustering")
  if (length(setdiff(branches, valid_branches)) > 0L) {
    stop("Unsupported `branches`: ", paste(setdiff(branches, valid_branches), collapse = ", "))
  }

  sample_names <- mrf3_sample_names(x)
  n <- length(sample_names)
  n_sub <- max(2L, as.integer(floor(n * sample_frac)))

  truth_use <- NULL
  if (!is.null(truth)) {
    truth_use <- align_labels_to_samples(truth, sample_names = sample_names)
  }

  build_branch_input <- function(branch) {
    if (identical(branch, "specific_shared")) {
      S <- x$shared$clustering$similarity
      cl0 <- x$clusters
      if (is.null(S) || is.null(cl0)) {
        stop("`specific_shared` branch is unavailable in `mrf3_fit` output.")
      }
      S <- as.matrix(S)
      list(
        branch = branch,
        kind = "similarity",
        mat = S,
        base_cl = as.character(cl0),
        k = if (is.null(specific_k)) length(unique(cl0)) else as.integer(specific_k),
        method = specific_method
      )
    } else {
      dat <- x$robust_detail$shared$clustering$similarity
      cl0 <- x$robust_clusters
      if (is.null(dat) || is.null(cl0)) {
        stop("`robust_clustering` branch is unavailable in `mrf3_fit` output.")
      }
      dat <- as.matrix(dat)
      kind <- if (nrow(dat) == ncol(dat)) "similarity" else "feature"
      list(
        branch = branch,
        kind = kind,
        mat = dat,
        base_cl = as.character(cl0),
        k = if (is.null(robust_k)) length(unique(cl0)) else as.integer(robust_k),
        method = robust_method
      )
    }
  }

  branch_inputs <- lapply(branches, build_branch_input)
  names(branch_inputs) <- branches

  if (verbose) {
    message("Running stability repeats: n_rep = ", n_rep, ", sample_size = ", n_sub, "/", n)
  }
  set.seed(seed)

  out_branches <- vector("list", length(branch_inputs))
  names(out_branches) <- names(branch_inputs)

  for (b in names(branch_inputs)) {
    bi <- branch_inputs[[b]]
    if (verbose) {
      message("Stability branch: ", b)
    }
    n_sample_branch <- nrow(bi$mat)
    if (n_sample_branch != n) {
      stop("Branch `", b, "` sample size mismatch with mrf3 samples.")
    }

    co_num <- matrix(0, nrow = n, ncol = n, dimnames = list(sample_names, sample_names))
    co_den <- matrix(0, nrow = n, ncol = n, dimnames = list(sample_names, sample_names))
    metrics_rows <- vector("list", n_rep)
    cl_by_rep <- vector("list", n_rep)

    for (r in seq_len(n_rep)) {
      idx_raw <- if (identical(sample_mode, "bootstrap")) {
        sample.int(n, size = n_sub, replace = TRUE)
      } else {
        sample.int(n, size = n_sub, replace = FALSE)
      }
      idx <- sort(unique(idx_raw))
      if (length(idx) < 2L) {
        metrics_rows[[r]] <- data.frame(
          branch = b,
          rep = r,
          n_sub = length(idx),
          status = "skipped",
          method = NA_character_,
          ari_vs_base = NA_real_,
          jaccard_vs_base = NA_real_,
          nmi_vs_base = NA_real_,
          ari_vs_truth = NA_real_,
          jaccard_vs_truth = NA_real_,
          nmi_vs_truth = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }

      fit <- tryCatch(
        cluster_subset(
          mat = bi$mat,
          idx = idx,
          k = bi$k,
          method = bi$method,
          kind = bi$kind
        ),
        error = function(e) e
      )

      if (inherits(fit, "error")) {
        metrics_rows[[r]] <- data.frame(
          branch = b,
          rep = r,
          n_sub = length(idx),
          status = "failed",
          method = NA_character_,
          ari_vs_base = NA_real_,
          jaccard_vs_base = NA_real_,
          nmi_vs_base = NA_real_,
          ari_vs_truth = NA_real_,
          jaccard_vs_truth = NA_real_,
          nmi_vs_truth = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }

      cl_sub <- as.character(fit$cl)
      cl_full <- rep(NA_character_, n)
      cl_full[idx] <- cl_sub
      names(cl_full) <- sample_names
      cl_by_rep[[r]] <- cl_full

      co_res <- update_coassign_mats(co_num = co_num, co_den = co_den, idx = idx, cl_sub = cl_sub)
      co_num <- co_res$co_num
      co_den <- co_res$co_den

      base_res <- cluster_metric_pair(pred = cl_sub, ref = bi$base_cl[idx])
      truth_res <- if (is.null(truth_use)) {
        list(ari = NA_real_, jaccard = NA_real_, nmi = NA_real_)
      } else {
        cluster_metric_pair(pred = cl_sub, ref = truth_use[idx])
      }

      metrics_rows[[r]] <- data.frame(
        branch = b,
        rep = r,
        n_sub = length(idx),
        status = "ok",
        method = fit$method,
        ari_vs_base = base_res$ari,
        jaccard_vs_base = base_res$jaccard,
        nmi_vs_base = base_res$nmi,
        ari_vs_truth = truth_res$ari,
        jaccard_vs_truth = truth_res$jaccard,
        nmi_vs_truth = truth_res$nmi,
        stringsAsFactors = FALSE
      )
    }

    metrics_tb <- do.call(rbind, metrics_rows)
    co_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(sample_names, sample_names))
    nz <- co_den > 0
    co_mat[nz] <- co_num[nz] / co_den[nz]
    diag(co_mat) <- 1

    consensus <- NULL
    if (isTRUE(run_consensus)) {
      cons_fit <- tryCatch(
        cluster_consensus_matrix(
          co_mat = co_mat,
          k = bi$k,
          method = consensus_method
        ),
        error = function(e) e
      )
      if (!inherits(cons_fit, "error")) {
        c_base <- cluster_metric_pair(pred = cons_fit$cl, ref = bi$base_cl)
        c_truth <- if (is.null(truth_use)) {
          list(ari = NA_real_, jaccard = NA_real_, nmi = NA_real_)
        } else {
          cluster_metric_pair(pred = cons_fit$cl, ref = truth_use)
        }
        consensus <- list(
          cl = cons_fit$cl,
          k = bi$k,
          method = cons_fit$method,
          metrics = data.frame(
            branch = b,
            ari_vs_base = c_base$ari,
            jaccard_vs_base = c_base$jaccard,
            nmi_vs_base = c_base$nmi,
            ari_vs_truth = c_truth$ari,
            jaccard_vs_truth = c_truth$jaccard,
            nmi_vs_truth = c_truth$nmi,
            stringsAsFactors = FALSE
          )
        )
      }
    }

    ok <- metrics_tb$status == "ok"
    summary_tb <- data.frame(
      branch = b,
      n_rep = n_rep,
      n_ok = sum(ok),
      ari_vs_base_mean = safe_mean(metrics_tb$ari_vs_base[ok]),
      jaccard_vs_base_mean = safe_mean(metrics_tb$jaccard_vs_base[ok]),
      nmi_vs_base_mean = safe_mean(metrics_tb$nmi_vs_base[ok]),
      ari_vs_truth_mean = safe_mean(metrics_tb$ari_vs_truth[ok]),
      jaccard_vs_truth_mean = safe_mean(metrics_tb$jaccard_vs_truth[ok]),
      nmi_vs_truth_mean = safe_mean(metrics_tb$nmi_vs_truth[ok]),
      stringsAsFactors = FALSE
    )

    out_branches[[b]] <- list(
      metrics = metrics_tb,
      summary = summary_tb,
      coassign = co_mat,
      coassign_coverage = ifelse(co_den > 0, 1, 0),
      consensus = consensus,
      cluster_by_rep = cl_by_rep,
      settings = list(k = bi$k, method = bi$method, kind = bi$kind)
    )
  }

  out <- list(
    summary = do.call(rbind, lapply(out_branches, `[[`, "summary")),
    by_branch = out_branches,
    params = list(
      branches = branches,
      n_rep = n_rep,
      sample_frac = sample_frac,
      sample_mode = sample_mode,
      run_consensus = run_consensus,
      consensus_method = consensus_method,
      seed = seed
    )
  )
  class(out) <- c("mrf3_stability", "list")
  out
}


#' Print `mrf3_stability`
#'
#' @method print mrf3_stability
#' @param x Output from `mrf3_stability()`.
#' @param ... Unused.
#' @return Input object invisibly.
#' @export
print.mrf3_stability <- function(x, ...) {
  cat("mrf3_stability summary\n")
  if (is.data.frame(x$summary)) {
    print(x$summary, row.names = FALSE)
  } else {
    cat("No summary available.\n")
  }
  invisible(x)
}


#' Tidy cluster labels from mrf3_fit
#'
#' @param x An object from `mrf3_fit()`.
#'
#' @return A tibble with columns:
#' `sample`, `branch`, `omics`, `cluster`.
mrf3_tidy_clusters <- function(x) {
  if (!is_mrf3_fit_object(x)) {
    stop("`x` must be an `mrf3_fit` object.")
  }
  sample_names <- mrf3_sample_names(x)
  rows <- list()
  rid <- 1L

  push_rows <- function(branch, cl, omics = NA_character_) {
    if (is.null(cl)) {
      return(invisible(NULL))
    }
    cl <- as.character(cl)
    if (length(cl) != length(sample_names)) {
      return(invisible(NULL))
    }
    rows[[rid]] <<- tibble::tibble(
      sample = sample_names,
      branch = branch,
      omics = omics,
      cluster = cl
    )
    rid <<- rid + 1L
    invisible(NULL)
  }

  push_rows("specific_shared", x$clusters)
  push_rows("robust_clustering", x$robust_clusters)

  by_omics <- x$specific$clustering$by_omics
  if (is.list(by_omics) && length(by_omics) > 0L) {
    for (nm in names(by_omics)) {
      push_rows("specific_specific", by_omics[[nm]]$cl, omics = nm)
    }
  }

  if (length(rows) == 0L) {
    return(tibble::tibble(
      sample = character(0),
      branch = character(0),
      omics = character(0),
      cluster = character(0)
    ))
  }
  do.call(rbind, rows)
}


#' Tidy IMD weights from mrf3_fit
#'
#' @param x An object from `mrf3_fit()`.
#' @param top_n Optional top-n variables per block.
#' @param abs_weight Logical; whether to rank by absolute weight.
#' @param cluster Optional cluster id. If provided, export IMD from
#' `x$cluster_imd$by_cluster[[cluster]]`.
#'
#' @return A tibble with columns:
#' `source`, `cluster`, `block`, `variable`, `weight`, `abs_weight`, `rank`.
mrf3_tidy_imd <- function(x, top_n = NULL, abs_weight = TRUE, cluster = NULL) {
  if (!is_mrf3_fit_object(x)) {
    stop("`x` must be an `mrf3_fit` object.")
  }
  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1) {
      stop("`top_n` must be NULL or a single integer >= 1.")
    }
    top_n <- as.integer(top_n)
  }

  source <- "imd"
  w_list <- x$weights
  cl_name <- NA_character_
  if (!is.null(cluster)) {
    source <- "cluster_imd"
    cl_name <- as.character(cluster)[1]
    w_list <- x$cluster_imd$by_cluster[[cl_name]]$imd$weight_list
  }
  if (!is.list(w_list) || length(w_list) == 0L) {
    return(tibble::tibble(
      source = character(0),
      cluster = character(0),
      block = character(0),
      variable = character(0),
      weight = numeric(0),
      abs_weight = numeric(0),
      rank = integer(0)
    ))
  }

  rows <- lapply(names(w_list), function(b) {
    w <- as.numeric(w_list[[b]])
    nm <- names(w_list[[b]])
    if (is.null(nm) || length(nm) != length(w)) {
      nm <- paste0("V", seq_along(w))
    }
    aw <- if (isTRUE(abs_weight)) abs(w) else w
    ord <- order(aw, decreasing = TRUE)
    w <- w[ord]
    aw <- aw[ord]
    nm <- nm[ord]
    rk <- seq_along(w)
    if (!is.null(top_n)) {
      keep <- seq_len(min(top_n, length(w)))
      w <- w[keep]
      aw <- aw[keep]
      nm <- nm[keep]
      rk <- rk[keep]
    }
    tibble::tibble(
      source = source,
      cluster = cl_name,
      block = b,
      variable = nm,
      weight = w,
      abs_weight = aw,
      rank = as.integer(rk)
    )
  })
  do.call(rbind, rows)
}


#' Tidy shared-fraction summary from mrf3_fit
#'
#' @param x An object from `mrf3_fit()`.
#'
#' @return A tibble. Empty tibble when unavailable.
mrf3_tidy_shared <- function(x) {
  if (!is_mrf3_fit_object(x)) {
    stop("`x` must be an `mrf3_fit` object.")
  }
  tb <- x$shared$frac
  if (!is.data.frame(tb) || nrow(tb) == 0L) {
    return(tibble::tibble())
  }
  tibble::as_tibble(tb)
}
