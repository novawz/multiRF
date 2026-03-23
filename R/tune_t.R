#' Tune model-level top-v
#' @param dat.list A list of omics matrices used for clustering.
#' @param mod A fitted `mrf3` model object.
#' @param tmin Minimum model-level top-v cutoff to evaluate.
#' The effective upper bound is fixed to `ceiling(0.2 * n)` where `n` is
#' sample size.
#' @param by Base step size for model-level top-v grid construction.
#' @param k Optional fixed number of clusters.
#' @param sample_n Optional integer. If set, tune on a random subset of samples.
#' @param sample_frac Optional fraction in (0, 1]. Used when `sample_n` is `NULL`.
#' @param auto_sample_n Logical; when `sample_n` and `sample_frac` are both `NULL`,
#' automatically infer a tuning sample size from input data size.
#' @param max_candidates Optional maximum number of grid points to evaluate.
#' If exceeded, candidate spacing is expanded automatically. Default is `20`.
#' @param reuse_tuned_k Logical; when `k` is `NULL`, tune `k` once on baseline and
#' reuse it across all candidates for speed.
#' @param parallel Logical; whether to evaluate candidate grid values in parallel
#' (POSIX systems only). Default is `TRUE`.
#' @param cores Number of cores used when `parallel = TRUE`.
#' @param seed Random seed used for optional sample subsampling.
#' @param object Objective used to choose `model_top_v`
#' (`"entropy_elbow"` (default), `"diss"`, `"silhouette"`, or `"eigen"`).
#' @rdname tune_model_top_v
tune_model_top_v <- function(dat.list, mod, tmin = 10, by = 1, k = NULL,
                             sample_n = NULL, sample_frac = NULL,
                             auto_sample_n = TRUE,
                             max_candidates = 20,
                             reuse_tuned_k = TRUE,
                             parallel = TRUE,
                             cores = NULL,
                             seed = 529,
                             object = "entropy_elbow"){
  entropy_only <- identical(as.character(object)[1], "entropy_elbow")

  tune_prep <- prepare_tune_inputs(
    dat.list = dat.list,
    mod = mod,
    sample_n = sample_n,
    sample_frac = sample_frac,
    auto_sample_n = auto_sample_n,
    seed = seed
  )
  dat.list <- tune_prep$dat.list
  mod <- tune_prep$mod

  n_samp <- nrow(dat.list[[1]])
  tmax <- infer_fused_tune_vmax(n_samp)
  message("Auto model_top_v `tmax` = ", tmax, " (n = ", n_samp, ").")
  if (tmin > tmax) {
    message("Adjusting `tmin` from ", tmin, " to auto `tmax` = ", tmax, ".")
    tmin <- tmax
  }

  t_grid <- make_tune_grid(
    lower = tmin,
    upper = tmax,
    by = by,
    max_candidates = max_candidates
  )
  if (length(t_grid) == 0L) {
    stop("No valid `model_top_v` candidates were generated.")
  }
  baseline_v <- max(t_grid)
  tune_grid <- setdiff(t_grid, baseline_v)

  tune_ctx <- resolve_tune_mod_inputs(mod)
  rfit <- tune_ctx$rfit
  model_names <- names(rfit)
  alpha <- compute_tune_model_alpha(
    model_names = model_names,
    connection_score = tune_ctx$connection_score
  )
  cache_list <- build_model_weight_cache_list(
    rfit = rfit,
    vmax = max(c(t_grid, baseline_v))
  )

  eval_one <- function(tv, k_use = NULL) {
    W_all_raw <- build_fused_weight_from_cache(
      cache_list = cache_list,
      top_v = tv,
      alpha = alpha,
      keep_ties = TRUE
    )
    W_all <- postprocess_fused_weight(
      W = W_all_raw,
      top_v = NULL,
      row_normalize = TRUE,
      keep_ties = TRUE
    )
    ent <- calc_fused_weight_entropy(W_all)
    if (isTRUE(entropy_only)) {
      k_out <- if (!is.null(k_use) && is.finite(k_use)) as.integer(k_use)[1] else NA_integer_
      return(list(
        obj = NA_real_,
        sil = NA_real_,
        diffe = NA_real_,
        entropy = ent,
        model_top_v = as.integer(tv),
        k = k_out
      ))
    }
    S <- W_all %*% t(W_all)
    diag(S) <- 0
    stat <- evaluate_similarity_for_tuning(S, k_use = k_use)
    list(
      obj = stat$obj,
      sil = stat$sil,
      diffe = stat$diffe,
      entropy = ent,
      model_top_v = as.integer(tv),
      k = as.integer(stat$k)
    )
  }

  init <- eval_one(baseline_v, k_use = k)
  k_base <- init$k
  k_eval <- k
  if (!isTRUE(entropy_only) && is.null(k_eval) && isTRUE(reuse_tuned_k) && is.finite(k_base)) {
    k_eval <- k_base
  }

  rows <- list()
  if (length(tune_grid) > 0L) {
    rows <- eval_grid(
      grid = tune_grid,
      eval_fun = function(tv) eval_one(tv, k_use = k_eval),
      parallel = parallel,
      cores = cores
    )
  }

  df <- bind_tune_rows(
    rows = rows,
    value_col = "model_top_v"
  )
  df <- rbind(
    df,
    data.frame(
      obj = init$obj,
      sil = init$sil,
      diffe = init$diffe,
      entropy = init$entropy,
      model_top_v = baseline_v,
      k = k_base
    )
  )

  dfsumm <- select_tune_summary(df, object = object)

  list(
    tmax_tb = dfsumm,
    object = df
  )
}

#' Tune fused weight top-v truncation
#' @param dat.list A list of omics matrices used for clustering.
#' @param mod A fitted `mrf3` model object.
#' @param vmin Minimum `fused_top_v` cutoff to evaluate.
#' @param by Base step size for the `fused_top_v` grid.
#' @param vmax Maximum `fused_top_v` cutoff to evaluate.
#' If `NULL`, an adaptive small-v upper bound is used.
#' @param model_top_v Fixed model-level top-v cutoff used while tuning
#' `fused_top_v`.
#' @param k Optional fixed number of clusters.
#' @param sample_n Optional integer. If set, tune on a random subset of samples.
#' @param sample_frac Optional fraction in (0, 1]. Used when `sample_n` is `NULL`.
#' @param auto_sample_n Logical; when `sample_n` and `sample_frac` are both `NULL`,
#' automatically infer a tuning sample size from input data size.
#' @param max_candidates Optional maximum number of grid points to evaluate.
#' If exceeded, candidate spacing is expanded automatically. Default is `20`.
#' @param reuse_tuned_k Logical; when `k` is `NULL`, tune `k` once on baseline and
#' reuse it across all candidates for speed.
#' @param parallel Logical; whether to evaluate candidate grid values in parallel
#' (POSIX systems only). Default is `TRUE`.
#' @param cores Number of cores used when `parallel = TRUE`.
#' @param seed Random seed used for optional sample subsampling.
#' @param object Objective used to choose `fused_top_v`
#' (`"entropy_elbow"` (default), `"diss"`, `"silhouette"`, or `"eigen"`).
#' @param early_stop Logical; when `TRUE` and `object = "entropy_elbow"`,
#' stop tuning once a stable small-gain elbow is reached. Default is `FALSE`
#' so elbow selection is based on the full evaluated grid.
#' @param elbow_rel_tol Relative elbow threshold multiplier on max observed gain.
#' @param elbow_abs_tol Absolute lower bound for elbow threshold.
#' @param elbow_min_points Minimum evaluated points before early-stop is allowed.
#' @param elbow_min_frac Minimum evaluated fraction in `v` grid before early-stop.
#' @param elbow_patience Number of consecutive small-gain steps required to stop.
#' @param elbow_smooth_window Integer running-mean window used to smooth entropy
#' gains before elbow detection.
#' @rdname tune_fused_top_v
tune_fused_top_v <- function(dat.list, mod, vmin = 10, by = 1, vmax = NULL,
                             model_top_v = 10,
                             k = NULL,
                             sample_n = NULL, sample_frac = NULL,
                             auto_sample_n = TRUE,
                             max_candidates = 20,
                             reuse_tuned_k = TRUE,
                             parallel = TRUE,
                             cores = NULL,
                             seed = 529,
                             object = "entropy_elbow",
                             early_stop = FALSE,
                             elbow_rel_tol = 0.25,
                             elbow_abs_tol = 1e-4,
                             elbow_min_points = 4L,
                             elbow_min_frac = 0.2,
                             elbow_patience = 2L,
                             elbow_smooth_window = 3L){
  entropy_only <- identical(as.character(object)[1], "entropy_elbow")

  tune_prep <- prepare_tune_inputs(
    dat.list = dat.list,
    mod = mod,
    sample_n = sample_n,
    sample_frac = sample_frac,
    auto_sample_n = auto_sample_n,
    seed = seed
  )
  dat.list <- tune_prep$dat.list
  mod <- tune_prep$mod

  if (is.null(vmax)) {
    n_samp <- nrow(dat.list[[1]])
    vmax <- infer_fused_tune_vmax(n_samp)
    message("Auto fused_top_v `vmax` = ", vmax, " (n = ", n_samp, ").")
  }
  if (!is.numeric(vmin) || !is.numeric(vmax) || !is.numeric(by)) {
    stop("`vmin`, `vmax`, and `by` must be numeric.")
  }
  if (vmin <= 0 || vmax <= 0 || by <= 0) {
    stop("`vmin`, `vmax`, and `by` must be positive.")
  }
  if (vmin > vmax) {
    message("Adjusting `vmin` from ", vmin, " to auto `vmax` = ", vmax, ".")
    vmin <- vmax
  }
  if (!is.numeric(model_top_v) || length(model_top_v) != 1L ||
      !is.finite(model_top_v) || model_top_v <= 0) {
    stop("`model_top_v` must be a single positive numeric value.")
  }
  if (!is.logical(early_stop) || length(early_stop) != 1L || is.na(early_stop)) {
    stop("`early_stop` must be TRUE or FALSE.")
  }
  if (!is.numeric(elbow_rel_tol) || length(elbow_rel_tol) != 1L ||
      !is.finite(elbow_rel_tol) || elbow_rel_tol < 0) {
    stop("`elbow_rel_tol` must be a single non-negative numeric value.")
  }
  if (!is.numeric(elbow_abs_tol) || length(elbow_abs_tol) != 1L ||
      !is.finite(elbow_abs_tol) || elbow_abs_tol < 0) {
    stop("`elbow_abs_tol` must be a single non-negative numeric value.")
  }
  if (!is.numeric(elbow_min_points) || length(elbow_min_points) != 1L ||
      !is.finite(elbow_min_points) || elbow_min_points < 2) {
    stop("`elbow_min_points` must be a single integer >= 2.")
  }
  if (!is.numeric(elbow_min_frac) || length(elbow_min_frac) != 1L ||
      !is.finite(elbow_min_frac) || elbow_min_frac < 0 || elbow_min_frac > 1) {
    stop("`elbow_min_frac` must be a single numeric value in [0, 1].")
  }
  if (!is.numeric(elbow_patience) || length(elbow_patience) != 1L ||
      !is.finite(elbow_patience) || elbow_patience < 1) {
    stop("`elbow_patience` must be a single integer >= 1.")
  }
  if (!is.numeric(elbow_smooth_window) || length(elbow_smooth_window) != 1L ||
      !is.finite(elbow_smooth_window) || elbow_smooth_window < 1) {
    stop("`elbow_smooth_window` must be a single integer >= 1.")
  }
  elbow_min_points <- as.integer(elbow_min_points)
  elbow_patience <- as.integer(elbow_patience)
  elbow_smooth_window <- as.integer(elbow_smooth_window)
  model_top_v <- as.integer(model_top_v)

  tune_ctx <- resolve_tune_mod_inputs(mod)
  rfit <- tune_ctx$rfit
  model_names <- names(rfit)
  alpha <- compute_tune_model_alpha(
    model_names = model_names,
    connection_score = tune_ctx$connection_score
  )
  model_cache <- build_model_weight_cache_list(
    rfit = rfit,
    vmax = model_top_v
  )
  W_all_raw <- build_fused_weight_from_cache(
    cache_list = model_cache,
    top_v = model_top_v,
    alpha = alpha,
    keep_ties = TRUE
  )

  eval_one <- function(v, k_use = k) {
    if (is.null(v)) {
      W_eval <- postprocess_fused_weight(
        W = W_all_raw,
        top_v = NULL,
        row_normalize = TRUE,
        keep_ties = TRUE
      )
    } else {
      W_eval <- postprocess_fused_weight(
        W = W_all_raw,
        top_v = as.integer(v),
        row_normalize = TRUE,
        keep_ties = TRUE
      )
    }
    ent <- calc_fused_weight_entropy(W_eval)
    if (isTRUE(entropy_only)) {
      k_out <- if (!is.null(k_use) && is.finite(k_use)) as.integer(k_use)[1] else NA_integer_
      return(list(
        obj = NA_real_,
        sil = NA_real_,
        diffe = NA_real_,
        entropy = ent,
        k = k_out
      ))
    }
    S <- W_eval %*% t(W_eval)
    diag(S) <- 0
    stat <- evaluate_similarity_for_tuning(S, k_use = k_use)

    list(
      obj = stat$obj,
      sil = stat$sil,
      diffe = stat$diffe,
      entropy = ent,
      k = as.integer(stat$k)
    )
  }

  # Baseline: no fused top-v truncation.
  init <- eval_one(NULL, k_use = k)
  if (!isTRUE(entropy_only) && is.null(k) && isTRUE(reuse_tuned_k) && is.finite(init$k)) {
    k <- init$k
  }

  v_grid <- make_tune_grid(
    lower = vmin,
    upper = vmax,
    by = by,
    max_candidates = max_candidates
  )
  v_grid <- v_grid[v_grid > 0]
  v_grid <- v_grid[v_grid < ncol(W_all_raw)]
  if (length(v_grid) == 0L) {
    warning(
      "No valid finite `fused_top_v` candidates remain after filtering; returning no-trunc baseline only.",
      call. = FALSE
    )
  }

  rows <- list()
  early_stop_used <- FALSE
  if (length(v_grid) > 0L) {
    run_elbow_early_stop <- identical(object, "entropy_elbow") && isTRUE(early_stop)
    if (run_elbow_early_stop) {
      early_stop_used <- TRUE
      rows <- eval_grid_until_entropy_elbow(
        grid = v_grid,
        eval_fun = function(v) eval_one(v, k_use = k),
        rel_tol = elbow_rel_tol,
        abs_tol = elbow_abs_tol,
        min_points = elbow_min_points,
        min_frac = elbow_min_frac,
        patience = elbow_patience,
        smooth_window = elbow_smooth_window
      )
    } else {
      rows <- eval_grid(
        grid = v_grid,
        eval_fun = function(v) eval_one(v, k_use = k),
        parallel = parallel,
        cores = cores
      )
    }
  }
  rows <- lapply(seq_along(rows), function(i) {
    x <- rows[[i]]
    data.frame(
      obj = x$obj,
      sil = x$sil,
      diffe = x$diffe,
      entropy = x$entropy,
      fused_top_v = as.integer(v_grid[i]),
      k = x$k,
      is_no_trunc = FALSE
    )
  })

  df <- if (length(rows) > 0L) {
    do.call(rbind, rows)
  } else {
    data.frame(
      obj = numeric(0),
      sil = numeric(0),
      diffe = numeric(0),
      entropy = numeric(0),
      fused_top_v = numeric(0),
      k = integer(0),
      is_no_trunc = logical(0)
    )
  }
  df <- rbind(
    df,
    data.frame(
      obj = init$obj,
      sil = init$sil,
      diffe = init$diffe,
      entropy = init$entropy,
      fused_top_v = Inf,
      k = init$k,
      is_no_trunc = TRUE
    )
  )

  if (isTRUE(early_stop_used) && nrow(df) > 1L && length(v_grid) > 0L) {
    n_eval <- nrow(df) - 1L
    if (n_eval < length(v_grid)) {
      message("Entropy elbow early-stop evaluated ", n_eval, "/", length(v_grid), " fused_top_v candidates.")
    }
  }

  dfsumm <- select_tune_summary(
    df,
    object = object,
    elbow_rel_tol = elbow_rel_tol,
    elbow_abs_tol = elbow_abs_tol,
    elbow_min_points = elbow_min_points,
    elbow_min_frac = elbow_min_frac,
    elbow_patience = elbow_patience,
    elbow_smooth_window = elbow_smooth_window
  )

  list(
    vtop_tb = dfsumm,
    object = df
  )
}

resolve_tune_mod_inputs <- function(mod) {
  if (inherits(mod, "mrf3")) {
    rfit <- mod$mod
    connection_score <- mod$connection_score
  } else if (is.list(mod) && !is.null(mod$mod)) {
    rfit <- mod$mod
    connection_score <- mod$connection_score
  } else {
    rfit <- mod
    connection_score <- NULL
  }

  if (!is.list(rfit) || length(rfit) == 0L) {
    stop("`mod` must contain a non-empty model list.")
  }
  if (is.null(names(rfit)) || any(names(rfit) == "")) {
    stop("Model list in `mod` must be named as `response_predictor`.")
  }

  list(
    rfit = rfit,
    connection_score = connection_score
  )
}

compute_tune_model_alpha <- function(model_names, connection_score = NULL) {
  score <- match_model_scores(
    model_names = model_names,
    connection_score = connection_score
  )
  score <- pmax(score, 0)
  alpha <- normalize_fusion_weights(score, fallback_uniform = TRUE)
  names(alpha) <- model_names
  alpha
}

build_model_weight_cache_list <- function(rfit, vmax = NULL) {
  out <- lapply(
    rfit,
    function(mod) {
      fw <- mod$forest.wt
      if (is.null(fw)) {
        stop("Model is missing `forest.wt`.")
      }
      W_adj <- prepare_weight_matrix(
        W = fw,
        adjust = TRUE,
        top_v = NULL,
        row_normalize = TRUE,
        zero_diag = TRUE,
        keep_ties = TRUE
      )
      build_weight_topv_cache(W_adj, vmax = vmax)
    }
  )
  names(out) <- names(rfit)
  out
}

build_weight_topv_cache <- function(W, vmax = NULL) {
  W <- validate_weight_matrix(W)
  p <- ncol(W)
  if (is.null(vmax)) {
    vmax_use <- p
  } else {
    vmax_use <- as.integer(vmax)[1]
    if (!is.finite(vmax_use) || vmax_use <= 0L) {
      stop("`vmax` must be NULL or a positive finite integer.")
    }
    vmax_use <- min(vmax_use, p)
  }

  ord <- t(apply(W, 1, order, decreasing = TRUE))
  ord <- ord[, seq_len(vmax_use), drop = FALSE]
  row_idx <- rep(seq_len(nrow(W)), times = ncol(ord))
  sorted_values <- matrix(
    W[cbind(row_idx, as.vector(ord))],
    nrow = nrow(W),
    ncol = ncol(ord)
  )

  list(
    W = W,
    ord = ord,
    sorted_values = sorted_values,
    vmax = vmax_use
  )
}

materialize_weight_from_cache <- function(cache, v, keep_ties = TRUE) {
  if (is.null(v)) {
    return(cache$W)
  }
  v <- as.integer(v)[1]
  if (!is.finite(v)) {
    stop("`v` must be NULL or a finite integer.")
  }
  if (v <= 0L) {
    return(matrix(0, nrow = nrow(cache$W), ncol = ncol(cache$W), dimnames = dimnames(cache$W)))
  }
  if (v >= ncol(cache$W)) {
    return(cache$W)
  }
  v_use <- min(v, cache$vmax)

  out <- matrix(0, nrow = nrow(cache$W), ncol = ncol(cache$W), dimnames = dimnames(cache$W))
  if (!isTRUE(keep_ties)) {
    idx <- cache$ord[, seq_len(v_use), drop = FALSE]
    row_idx <- rep(seq_len(nrow(cache$W)), times = ncol(idx))
    out[cbind(row_idx, as.vector(idx))] <- cache$W[cbind(row_idx, as.vector(idx))]
    return(out)
  }

  for (i in seq_len(nrow(cache$W))) {
    thr <- cache$sorted_values[i, v_use]
    keep <- which(cache$W[i, ] >= thr)
    if (length(keep) > 0L) {
      out[i, keep] <- cache$W[i, keep]
    }
  }
  out
}

build_fused_weight_from_cache <- function(cache_list, top_v, alpha, keep_ties = TRUE) {
  model_names <- names(cache_list)
  alpha_use <- alpha[model_names]
  alpha_use[!is.finite(alpha_use)] <- 0
  if (sum(alpha_use) <= 0) {
    alpha_use <- rep(1 / length(alpha_use), length(alpha_use))
  } else {
    alpha_use <- alpha_use / sum(alpha_use)
  }
  names(alpha_use) <- model_names

  W_models <- lapply(
    model_names,
    function(m) materialize_weight_from_cache(cache_list[[m]], v = top_v, keep_ties = keep_ties)
  )
  names(W_models) <- model_names
  fuse_matrix_list(W_models, alpha_use)
}

evaluate_similarity_for_tuning <- function(S, k_use = NULL) {
  S <- as.matrix(S)
  S[!is.finite(S)] <- 0
  S <- (S + t(S)) / 2
  diag(S) <- 0
  k_fallback <- if (!is.null(k_use) && is.finite(k_use)) as.integer(k_use)[1] else NA_integer_
  if (!all(is.finite(S)) || nrow(S) < 2L || ncol(S) < 2L || sum(abs(S)) <= .Machine$double.eps) {
    return(list(
      obj = NA_real_,
      sil = NA_real_,
      diffe = NA_real_,
      k = k_fallback
    ))
  }
  clm <- tryCatch(
    {
      if (is.null(k_use)) {
        tune_k_clusters(S, return_cluster = TRUE, method = "Spectral")
      } else {
        spectral_cl(S, k_tune = as.integer(k_use))
      }
    },
    error = function(e) NULL
  )
  if (is.null(clm) || is.null(clm$cl)) {
    return(list(
      obj = NA_real_,
      sil = NA_real_,
      diffe = NA_real_,
      k = k_fallback
    ))
  }

  obj <- if (!is.null(clm$obj) && length(clm$obj) > 0L) as.numeric(clm$obj)[1] else NA_real_
  sil <- if (!is.null(clm$sil) && length(clm$sil) > 0L) as.numeric(clm$sil)[1] else NA_real_
  diffe <- if (!is.null(clm$diff_e) && length(clm$diff_e) > 0L) {
    max(as.numeric(clm$diff_e), na.rm = TRUE)
  } else {
    NA_real_
  }
  if (!is.finite(diffe)) {
    diffe <- NA_real_
  }

  list(
    obj = obj,
    sil = sil,
    diffe = diffe,
    k = as.integer(length(unique(clm$cl)))
  )
}

#' Compute effective neighbourhood size for each row of a weight matrix
#'
#' For each row \eqn{w_i}, the effective neighbourhood size is defined as
#' \eqn{n_{\mathrm{eff},i} = \exp(H(w_i))}, where
#' \eqn{H(w_i) = -\sum_j w_{ij} \log w_{ij}} is the Shannon entropy.
#' This equals the number of equally-weighted neighbours that would
#' produce the same entropy.
#'
#' @param W A square row-stochastic weight matrix (or will be row-normalised).
#' @param eps Small positive value to avoid log(0).
#' @return A numeric vector of length \code{nrow(W)} containing \eqn{n_{\mathrm{eff},i}}.
#' @keywords internal
effective_neighbourhood_size <- function(W, eps = 1e-12) {
  W <- as.matrix(W)
  W[!is.finite(W)] <- 0
  W <- pmax(W, 0)

  rs <- rowSums(W)
  prob <- W
  ok <- rs > eps
  if (any(ok)) {
    prob[ok, ] <- W[ok, , drop = FALSE] / rs[ok]
  }
  if (any(!ok)) {
    prob[!ok, ] <- 0
  }

  apply(prob, 1, function(v) {
    v <- v[v > eps]
    if (length(v) <= 1L) return(1)
    exp(-sum(v * log(v)))
  })
}

#' Select top-v from effective neighbourhood size
#'
#' Sets \eqn{v = \lceil \mathrm{quantile}_q(n_{\mathrm{eff},i}) \rceil}
#' where \eqn{n_{\mathrm{eff},i} = \exp(H(w_i))}.
#'
#' @param W A square weight matrix.
#' @param quantile_prob Quantile probability. Default 0.5 (median).
#' @param min_v Floor on the returned value.
#' @param eps Small positive value for entropy computation.
#' @return A single integer: the selected top-v.
#' @export
select_top_v_neff <- function(W, quantile_prob = 0.5, min_v = 10L, eps = 1e-12) {
  neff <- effective_neighbourhood_size(W, eps = eps)
  v <- as.integer(ceiling(stats::quantile(neff, probs = quantile_prob, na.rm = TRUE)))
  v <- max(v, as.integer(min_v))
  v <- min(v, ncol(W))
  v
}

calc_fused_weight_entropy <- function(W, eps = 1e-12) {
  W <- as.matrix(W)
  W[!is.finite(W)] <- 0
  W <- pmax(W, 0)

  p <- ncol(W)
  if (p <= 1L) {
    return(0)
  }

  rs <- rowSums(W)
  prob <- W
  ok <- rs > eps
  if (any(ok)) {
    prob[ok, ] <- W[ok, , drop = FALSE] / rs[ok]
  }
  if (any(!ok)) {
    prob[!ok, ] <- 0
  }

  row_entropy <- apply(
    prob,
    1,
    function(v) {
      v <- v[v > eps]
      if (length(v) <= 1L) {
        return(0)
      }
      -sum(v * log(v)) / log(p)
    }
  )
  e <- mean(row_entropy)
  if (!is.finite(e)) {
    return(NA_real_)
  }
  max(0, min(1, e))
}

bind_tune_rows <- function(rows, value_col) {
  if (length(rows) == 0L) {
    out <- data.frame(
      obj = numeric(0),
      sil = numeric(0),
      diffe = numeric(0),
      entropy = numeric(0),
      k = integer(0)
    )
    out[[value_col]] <- integer(0)
    return(out[, c("obj", "sil", "diffe", "entropy", value_col, "k"), drop = FALSE])
  }
  out <- do.call(
    rbind,
    lapply(rows, function(x) {
      data.frame(
        obj = x$obj,
        sil = x$sil,
        diffe = x$diffe,
        entropy = if (!is.null(x$entropy)) x$entropy else NA_real_,
        k = x$k,
        value = x[[value_col]],
        stringsAsFactors = FALSE
      )
    })
  )
  names(out)[names(out) == "value"] <- value_col
  out[[value_col]] <- as.integer(out[[value_col]])
  out$k <- as.integer(out$k)
  out
}

select_tune_summary <- function(df, object = "diss",
                                elbow_rel_tol = 0.25,
                                elbow_abs_tol = 1e-4,
                                elbow_min_points = 4L,
                                elbow_min_frac = 0.2,
                                elbow_patience = 2L,
                                elbow_smooth_window = 3L) {
  if (!is.data.frame(df) || nrow(df) == 0L) {
    stop("No tuning results available to summarize.")
  }
  object <- match.arg(object, c("diss", "silhouette", "eigen", "entropy_elbow"))

  if (object == "diss") {
    return(df %>% dplyr::group_by(k) %>% dplyr::slice_min(obj, with_ties = FALSE))
  }
  if (object == "silhouette") {
    return(df %>% dplyr::group_by(k) %>% dplyr::slice_max(sil, with_ties = FALSE))
  }
  if (object == "eigen") {
    return(df %>% dplyr::group_by(k) %>% dplyr::slice_max(diffe, with_ties = FALSE))
  }
  v_col <- if ("fused_top_v" %in% names(df)) "fused_top_v" else if ("model_top_v" %in% names(df)) "model_top_v" else NA_character_
  if (!is.character(v_col) || !nzchar(v_col)) {
    stop("`object = 'entropy_elbow'` is only available when tuning fused_top_v or model_top_v.")
  }
  if (!("entropy" %in% names(df))) {
    stop("`object = 'entropy_elbow'` requires an `entropy` column in tuning results.")
  }

  finite_rows <- is.finite(df[[v_col]]) & is.finite(df$entropy)
  if (!any(finite_rows)) {
    idx <- if ("is_no_trunc" %in% names(df)) which(df$is_no_trunc %in% TRUE)[1] else NA_integer_
    if (!is.finite(idx)) {
      ord0 <- order(df[[v_col]], na.last = NA)
      idx <- if (length(ord0) > 0L) ord0[1] else 1L
    }
    out <- df[idx, , drop = FALSE]
    out$entropy_delta <- NA_real_
    out$entropy_gain_smooth <- NA_real_
    out$elbow_threshold <- NA_real_
    out$elbow_rule <- "no_finite_entropy"
    out$elbow_selected <- TRUE
    return(out)
  }

  x <- df[finite_rows, , drop = FALSE]
  ord <- order(x[[v_col]], na.last = NA)
  x <- x[ord, , drop = FALSE]

  x$entropy_delta <- c(NA_real_, diff(x$entropy))
  if (nrow(x) <= 1L) {
    x$entropy_gain_smooth <- NA_real_
    x$elbow_threshold <- NA_real_
    x$elbow_rule <- "single_point"
    x$elbow_selected <- TRUE
    return(x[1, , drop = FALSE])
  }

  elbow_info <- detect_entropy_elbow(
    v = x[[v_col]],
    entropy = x$entropy,
    rel_tol = elbow_rel_tol,
    abs_tol = elbow_abs_tol,
    min_points = elbow_min_points,
    min_frac = elbow_min_frac,
    patience = elbow_patience,
    smooth_window = elbow_smooth_window,
    fallback = TRUE
  )
  elbow_idx <- elbow_info$elbow_idx
  if (!is.finite(elbow_idx) || elbow_idx < 1L || elbow_idx > nrow(x)) {
    elbow_idx <- min(max(1L, as.integer(elbow_min_points)), nrow(x))
  }

  x$entropy_delta <- elbow_info$entropy_delta
  x$entropy_gain_smooth <- elbow_info$entropy_gain_smooth
  x$elbow_threshold <- elbow_info$threshold
  x$elbow_rule <- elbow_info$rule
  x$elbow_selected <- FALSE
  x$elbow_selected[elbow_idx] <- TRUE
  x[elbow_idx, , drop = FALSE]
}

make_tune_grid <- function(lower, upper, by, max_candidates = NULL) {
  if (!is.numeric(lower) || !is.numeric(upper) || !is.numeric(by) ||
      length(lower) != 1L || length(upper) != 1L || length(by) != 1L ||
      !is.finite(lower) || !is.finite(upper) || !is.finite(by) ||
      by <= 0 || upper < lower) {
    return(integer(0))
  }

  step <- by
  if (!is.null(max_candidates) &&
      is.numeric(max_candidates) &&
      length(max_candidates) == 1L &&
      is.finite(max_candidates) &&
      max_candidates >= 2) {
    n_raw <- floor((upper - lower) / by) + 1L
    if (n_raw > max_candidates) {
      step <- ceiling((upper - lower) / (max_candidates - 1L))
    }
  }

  grid <- unique(as.integer(seq(lower, upper, by = step)))
  if (length(grid) == 0L || tail(grid, 1L) != as.integer(upper)) {
    grid <- unique(c(grid, as.integer(upper)))
  }
  grid
}

subset_rfit_to_samples <- function(rfit, sample_names = NULL, sample_idx = NULL) {
  subset_square <- function(mat) {
    if (is.null(mat)) return(NULL)
    if (!is.null(sample_names) && !is.null(rownames(mat))) {
      keep <- rownames(mat) %in% sample_names
      return(mat[keep, keep, drop = FALSE])
    }
    if (!is.null(sample_idx)) {
      return(mat[sample_idx, sample_idx, drop = FALSE])
    }
    mat
  }
  subset_rect <- function(mat) {
    if (is.null(mat)) return(NULL)
    if (!is.null(sample_names) && !is.null(rownames(mat))) {
      keep <- rownames(mat) %in% sample_names
      return(mat[keep, , drop = FALSE])
    }
    if (!is.null(sample_idx)) {
      return(mat[sample_idx, , drop = FALSE])
    }
    mat
  }

  purrr::map(
    rfit,
    function(m) {
      m2 <- m
      m2$xvar <- subset_rect(m2$xvar)
      m2$yvar <- subset_rect(m2$yvar)
      m2$forest.wt <- subset_square(m2$forest.wt)
      m2$proximity <- subset_square(m2$proximity)
      m2$membership <- subset_rect(m2$membership)
      m2
    }
  )
}

infer_tune_sample_n <- function(n) {
  if (!is.finite(n) || n <= 0) return(NA_integer_)
  n <- as.integer(n)
  if (n <= 600L) return(n)
  if (n <= 2000L) return(600L)
  800L
}

prepare_tune_inputs <- function(dat.list, mod, sample_n = NULL, sample_frac = NULL,
                                auto_sample_n = TRUE, seed = 529) {

  n <- nrow(dat.list[[1]])
  if (!all(vapply(dat.list, nrow, integer(1)) == n)) {
    warning("`dat.list` has inconsistent sample size across blocks. Skip subsampling.", call. = FALSE)
    return(list(dat.list = dat.list, mod = mod))
  }

  n_target <- sample_n
  if (is.null(n_target) && !is.null(sample_frac)) {
    n_target <- floor(n * sample_frac)
  }
  if (is.null(n_target) && isTRUE(auto_sample_n) && is.null(sample_frac)) {
    n_target <- infer_tune_sample_n(n)
    if (is.finite(n_target)) {
      message("Tuning speed-up: auto sample_n selected as ", n_target, " (n = ", n, ").")
    }
  }
  if (is.null(n_target)) {
    return(list(dat.list = dat.list, mod = mod))
  }
  if (!is.numeric(n_target) || length(n_target) != 1L || !is.finite(n_target)) {
    warning("Invalid `sample_n`/`sample_frac`. Skip subsampling.", call. = FALSE)
    return(list(dat.list = dat.list, mod = mod))
  }
  n_target <- as.integer(max(2L, min(n, n_target)))
  if (n_target >= n) {
    return(list(dat.list = dat.list, mod = mod))
  }

  set.seed(seed)
  idx <- sort(sample.int(n, n_target))
  dat_sub <- purrr::map(dat.list, ~.[idx, , drop = FALSE])
  sample_names <- rownames(dat_sub[[1]])
  if (is.null(sample_names)) {
    sample_names <- NULL
  }

  mod_sub <- mod
  if (!is.null(mod_sub$mod) && is.list(mod_sub$mod)) {
    mod_sub$mod <- subset_rfit_to_samples(
      mod_sub$mod,
      sample_names = sample_names,
      sample_idx = idx
    )
  }
  if (!is.null(mod_sub$recon)) {
    mod_sub$recon <- NULL
  }

  message("Tuning speed-up: subsampled ", n_target, "/", n, " samples.")
  list(dat.list = dat_sub, mod = mod_sub)
}

eval_grid <- function(grid, eval_fun, parallel = FALSE, cores = NULL) {
  if (!isTRUE(parallel) || length(grid) <= 1L || .Platform$OS.type == "windows") {
    return(lapply(grid, eval_fun))
  }
  cores <- sanitize_mc_cores(cores = cores, fallback = 1L)
  if (cores <= 1L) {
    return(lapply(grid, eval_fun))
  }
  parallel::mclapply(grid, eval_fun, mc.cores = cores)
}

smooth_running_mean <- function(x, window = 3L) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0L) {
    return(numeric(0))
  }
  window <- as.integer(max(1L, window))
  if (window <= 1L) {
    return(x)
  }

  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - window + 1L)
    seg <- x[lo:i]
    seg <- seg[is.finite(seg)]
    if (length(seg) > 0L) {
      out[i] <- mean(seg)
    }
  }
  out
}

pick_kneedle_fallback <- function(v, entropy, min_eval = 2L, smooth_window = 3L) {
  n <- length(entropy)
  min_eval <- min(max(2L, as.integer(min_eval)), n)
  if (n <= 2L) {
    return(list(idx = min_eval, rule = "min_eval_fallback"))
  }

  x <- as.numeric(v)
  y <- as.numeric(entropy)
  y_smooth <- smooth_running_mean(y, window = smooth_window)
  y_use <- ifelse(is.finite(y_smooth), y_smooth, y)

  xr <- range(x, na.rm = TRUE)
  yr <- range(y_use, na.rm = TRUE)
  if (!all(is.finite(xr)) || !all(is.finite(yr)) || xr[1] == xr[2] || yr[1] == yr[2]) {
    return(list(idx = min_eval, rule = "min_eval_fallback"))
  }

  x_norm <- (x - xr[1]) / (xr[2] - xr[1])
  y_norm <- (y_use - yr[1]) / (yr[2] - yr[1])
  dist <- y_norm - x_norm

  cand <- seq.int(max(2L, min_eval), max(2L, n - 1L))
  if (length(cand) == 0L) {
    return(list(idx = min_eval, rule = "min_eval_fallback"))
  }

  dist_cand <- dist[cand]
  dist_cand[!is.finite(dist_cand)] <- -Inf
  j <- which.max(dist_cand)
  if (length(j) == 0L || !is.finite(dist_cand[j]) || dist_cand[j] <= 0) {
    return(list(idx = min_eval, rule = "min_eval_fallback"))
  }

  list(idx = as.integer(cand[j]), rule = "kneedle_fallback")
}

detect_entropy_elbow <- function(v, entropy,
                                 rel_tol = 0.25,
                                 abs_tol = 1e-4,
                                 min_points = 4L,
                                 min_frac = 0.2,
                                 patience = 2L,
                                 smooth_window = 3L,
                                 fallback = TRUE) {
  v <- as.numeric(v)
  entropy <- as.numeric(entropy)
  n <- length(entropy)

  if (n == 0L) {
    return(list(
      elbow_idx = NA_integer_,
      threshold = NA_real_,
      entropy_delta = numeric(0),
      entropy_gain_smooth = numeric(0),
      rule = "none"
    ))
  }

  entropy_delta <- c(NA_real_, diff(entropy))
  gain <- pmax(entropy_delta, 0)
  gain_smooth <- c(NA_real_, smooth_running_mean(gain[-1], window = smooth_window))

  ref_gain <- suppressWarnings(max(gain_smooth[-1], na.rm = TRUE))
  if (!is.finite(ref_gain)) {
    ref_gain <- suppressWarnings(max(gain[-1], na.rm = TRUE))
  }
  if (!is.finite(ref_gain)) {
    ref_gain <- 0
  }
  threshold <- max(as.numeric(abs_tol), as.numeric(rel_tol) * ref_gain)

  min_eval <- max(
    as.integer(min_points),
    as.integer(ceiling(as.numeric(min_frac) * n)),
    2L
  )
  min_eval <- min(max(min_eval, 2L), n)
  patience <- max(1L, as.integer(patience))

  low_flag <- gain_smooth <= threshold
  low_flag[!is.finite(low_flag)] <- FALSE
  low_flag[1] <- FALSE

  elbow_idx <- NA_integer_
  max_start <- n - patience + 1L
  if (max_start >= min_eval) {
    for (i in seq.int(min_eval, max_start)) {
      idx <- seq.int(i, i + patience - 1L)
      if (all(low_flag[idx])) {
        elbow_idx <- i
        break
      }
    }
  }

  rule <- if (is.finite(elbow_idx)) "small_gain" else "none"
  if (!is.finite(elbow_idx) && isTRUE(fallback)) {
    fb <- pick_kneedle_fallback(
      v = v,
      entropy = entropy,
      min_eval = min_eval,
      smooth_window = smooth_window
    )
    elbow_idx <- fb$idx
    rule <- fb$rule
  }

  list(
    elbow_idx = as.integer(elbow_idx),
    threshold = as.numeric(threshold),
    entropy_delta = entropy_delta,
    entropy_gain_smooth = gain_smooth,
    rule = as.character(rule)[1]
  )
}

infer_fused_tune_vmax <- function(n,
                                  frac = 0.2) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0) {
    stop("`n` must be a single positive numeric value.")
  }
  if (!is.numeric(frac) || length(frac) != 1L ||
      !is.finite(frac) || frac <= 0) {
    stop("`frac` must be a single positive numeric value.")
  }

  n <- as.integer(n)
  target <- as.integer(ceiling(frac * n))
  vmax <- min(target, n)
  as.integer(max(1L, vmax))
}

eval_grid_until_entropy_elbow <- function(grid, eval_fun,
                                          rel_tol = 0.25,
                                          abs_tol = 1e-4,
                                          min_points = 4L,
                                          min_frac = 0.2,
                                          patience = 2L,
                                          smooth_window = 3L) {
  if (length(grid) == 0L) {
    return(list())
  }

  rows <- vector("list", length(grid))
  entropy_vals <- rep(NA_real_, length(grid))
  n_eval <- 0L

  for (i in seq_along(grid)) {
    rows[[i]] <- eval_fun(grid[[i]])
    n_eval <- i
    entropy_vals[i] <- as.numeric(rows[[i]]$entropy)[1]

    if (i < 2L) {
      next
    }

    elbow_info <- detect_entropy_elbow(
      v = grid[seq_len(i)],
      entropy = entropy_vals[seq_len(i)],
      rel_tol = rel_tol,
      abs_tol = abs_tol,
      min_points = min_points,
      min_frac = min_frac,
      patience = patience,
      smooth_window = smooth_window,
      fallback = FALSE
    )
    if (identical(elbow_info$rule, "small_gain") && is.finite(elbow_info$elbow_idx)) {
      break
    }
  }

  rows[seq_len(n_eval)]
}
