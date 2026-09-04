#' Tune model-level top-v
#' @param dat.list A list of omics matrices used for clustering.
#' @param mod A fitted `mrf3` model object.
#' @param tmin Minimum model-level top-v cutoff to evaluate.
#' The effective upper bound is the sample size, allowing the
#' `v >= 0.8 * n` no-truncation rule to be evaluated.
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
#' @param parallel Logical; whether to allow candidate-grid evaluation in fresh
#'   PSOCK worker processes. Parallel evaluation also requires an explicit
#'   `cores` value.
#' @param cores Number of cores used when `parallel = TRUE`. Default `NULL`
#'   keeps candidate-grid evaluation serial; supply a positive integer to use
#'   process-level workers.
#' @param seed Random seed used for optional sample subsampling.
#' @param object Objective used to choose `model_top_v`:
#' `"saturation"` (default), `"entropy_elbow"`, `"diss"`, `"silhouette"`, or
#' `"eigen"`. `"saturation"` selects the smallest `v` whose fused-weight row
#' entropy reaches a fraction `tau` of the no-truncation entropy (linearly
#' interpolated between grid points), so the choice does not depend on the
#' grid resolution. `"entropy_elbow"` keeps the previous small-gain elbow
#' heuristic; its elbow is selected among interior grid points, so the
#' smallest grid candidate cannot be selected directly.
#' @param tau Saturation fraction in (0, 1) used by `object = "saturation"`.
#' Default `0.9`.
#' @rdname tune_model_top_v
#' @export
tune_model_top_v <- function(dat.list, mod, tmin = 10, by = 1, k = NULL,
                             sample_n = NULL, sample_frac = NULL,
                             auto_sample_n = FALSE,
                             max_candidates = 20,
                             reuse_tuned_k = TRUE,
                             parallel = TRUE,
                             cores = NULL,
                             seed = 529,
                             object = "saturation",
                             tau = 0.9){
  object <- match.arg(
    as.character(object)[1L],
    c("saturation", "entropy_elbow", "diss", "silhouette", "eigen")
  )
  entropy_only <- object %in% c("saturation", "entropy_elbow")
  tau <- .check_saturation_tau(tau)

  if (!is.numeric(tmin) || !is.numeric(by)) {
    stop("`tmin` and `by` must be numeric.")
  }
  if (tmin <= 0 || by <= 0) {
    stop("`tmin` and `by` must be positive.")
  }

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
  t_grid <- t_grid[t_grid > 0]
  if (length(t_grid) == 0L) {
    stop("No valid `model_top_v` candidates were generated.")
  }
  baseline_v <- max(t_grid)
  tune_grid <- setdiff(t_grid, baseline_v)

  tune_ctx <- resolve_tune_mod_inputs(mod)
  rfit <- tune_ctx$rfit
  model_names <- names(rfit)
  cache_list <- build_model_weight_cache_list(
    rfit = rfit,
    vmax = max(c(t_grid, baseline_v))
  )

  eval_one <- function(tv, k_use = NULL) {
    W_all_raw <- build_response_fused_weight_from_cache(
      cache_list = cache_list,
      top_v = tv,
      connection_score = tune_ctx$connection_score,
      model_list = rfit,
      response_blocks = names(dat.list),
      recon_fusion = tune_ctx$recon_fusion,
      score_power = tune_ctx$score_power,
      score_floor = tune_ctx$score_floor,
      fallback_uniform = tune_ctx$fallback_uniform,
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

  df$entropy_frac <- entropy_fraction(df$entropy, reference = init$entropy)

  dfsumm <- select_tune_summary(df, object = object, tau = tau)

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
#' If `NULL`, the full sample size is used so the
#' `v >= 0.8 * n` no-truncation rule is represented in the grid.
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
#' @param parallel Logical; whether to allow candidate-grid evaluation in fresh
#'   PSOCK worker processes. Parallel evaluation also requires an explicit
#'   `cores` value.
#' @param cores Number of cores used when `parallel = TRUE`. Default `NULL`
#'   keeps candidate-grid evaluation serial; supply a positive integer to use
#'   process-level workers.
#' @param seed Random seed used for optional sample subsampling.
#' @param object Objective used to choose `fused_top_v`:
#' `"saturation"` (default), `"entropy_elbow"`, `"diss"`, `"silhouette"`, or
#' `"eigen"`. For the two entropy objectives the entropy of every truncation
#' level is obtained in closed form from the sorted fused weights (see
#' `fused_entropy_curve()`), so no candidate matrices are rebuilt.
#' `"saturation"` evaluates every integer `v` in `[vmin, vmax]` and selects
#' the smallest one whose entropy reaches `tau` times the no-truncation
#' entropy. `"entropy_elbow"` keeps the previous small-gain elbow heuristic on
#' the candidate grid; its elbow is selected among interior grid points, so
#' the smallest grid candidate and the no-truncation baseline cannot be
#' selected directly (no truncation is recovered separately via the
#' `v >= 0.8 * n` rule in the workflow).
#' @param tau Saturation fraction in (0, 1) used by `object = "saturation"`.
#' Default `0.9`.
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
#' @export
tune_fused_top_v <- function(dat.list, mod, vmin = 10, by = 1, vmax = NULL,
                             model_top_v = 10,
                             k = NULL,
                             sample_n = NULL, sample_frac = NULL,
                             auto_sample_n = FALSE,
                             max_candidates = 20,
                             reuse_tuned_k = TRUE,
                             parallel = TRUE,
                             cores = NULL,
                             seed = 529,
                             object = "saturation",
                             tau = 0.9,
                             early_stop = FALSE,
                             elbow_rel_tol = 0.25,
                             elbow_abs_tol = 1e-4,
                             elbow_min_points = 4L,
                             elbow_min_frac = 0.2,
                             elbow_patience = 2L,
                             elbow_smooth_window = 3L){
  object <- match.arg(
    as.character(object)[1L],
    c("saturation", "entropy_elbow", "diss", "silhouette", "eigen")
  )
  entropy_only <- object %in% c("saturation", "entropy_elbow")
  tau <- .check_saturation_tau(tau)

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
  model_cache <- build_model_weight_cache_list(
    rfit = rfit,
    vmax = model_top_v
  )
  W_all_raw <- build_response_fused_weight_from_cache(
    cache_list = model_cache,
    top_v = model_top_v,
    connection_score = tune_ctx$connection_score,
    model_list = rfit,
    response_blocks = names(dat.list),
    recon_fusion = tune_ctx$recon_fusion,
    score_power = tune_ctx$score_power,
    score_floor = tune_ctx$score_floor,
    fallback_uniform = tune_ctx$fallback_uniform,
    keep_ties = TRUE
  )

  # For the entropy objectives the whole entropy-versus-v curve follows in
  # closed form from the sorted fused weights, so candidates are looked up
  # instead of rebuilding a truncated matrix for each one.
  entropy_curve <- if (isTRUE(entropy_only)) {
    fused_entropy_curve(W_all_raw, keep_ties = TRUE)
  } else {
    NULL
  }

  eval_one <- function(v, k_use = k) {
    if (isTRUE(entropy_only)) {
      ent <- if (is.null(v)) {
        entropy_curve$entropy_inf
      } else {
        entropy_curve$entropy[as.integer(v)]
      }
      k_out <- if (!is.null(k_use) && is.finite(k_use)) as.integer(k_use)[1] else NA_integer_
      return(list(
        obj = NA_real_,
        sil = NA_real_,
        diffe = NA_real_,
        entropy = ent,
        k = k_out
      ))
    }
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

  v_grid <- if (identical(object, "saturation")) {
    # Closed-form lookups are free, so use the exact integer resolution.
    seq.int(as.integer(ceiling(vmin)), as.integer(floor(vmax)))
  } else {
    make_tune_grid(
      lower = vmin,
      upper = vmax,
      by = by,
      max_candidates = max_candidates
    )
  }
  v_grid <- v_grid[v_grid > 0]
  v_grid <- v_grid[v_grid <= ncol(W_all_raw)]
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
        parallel = isTRUE(parallel) && !isTRUE(entropy_only),
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
  rownames(df) <- NULL

  if (isTRUE(early_stop_used) && nrow(df) > 1L && length(v_grid) > 0L) {
    n_eval <- nrow(df) - 1L
    if (n_eval < length(v_grid)) {
      message("Entropy elbow early-stop evaluated ", n_eval, "/", length(v_grid), " fused_top_v candidates.")
    }
  }

  df$entropy_frac <- entropy_fraction(df$entropy, reference = init$entropy)

  dfsumm <- select_tune_summary(
    df,
    object = object,
    tau = tau,
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
  recon_fusion <- "weighted"
  score_power <- 1
  score_floor <- 0
  fallback_uniform <- TRUE
  response_blocks <- NULL
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

  if (is.list(mod)) {
    if (!is.null(mod$recon_fusion)) recon_fusion <- mod$recon_fusion
    if (!is.null(mod$score_power)) score_power <- mod$score_power
    if (!is.null(mod$score_floor)) score_floor <- mod$score_floor
    if (!is.null(mod$fallback_uniform)) fallback_uniform <- mod$fallback_uniform
    if (!is.null(mod$response_blocks)) response_blocks <- mod$response_blocks
  }

  if (!is.list(rfit) || length(rfit) == 0L) {
    stop("`mod` must contain a non-empty model list.")
  }
  if (is.null(names(rfit)) || any(names(rfit) == "")) {
    stop("Model list in `mod` must be named as `response_predictor`.")
  }
  if (!is.numeric(score_power) || length(score_power) != 1L ||
      !is.finite(score_power) || score_power <= 0) {
    stop("`score_power` must be a single positive numeric.")
  }
  if (!is.numeric(score_floor) || length(score_floor) != 1L ||
      !is.finite(score_floor) || score_floor < 0) {
    stop("`score_floor` must be a single non-negative numeric.")
  }

  list(
    rfit = rfit,
    connection_score = connection_score,
    recon_fusion = match.arg(as.character(recon_fusion)[1L], c("weighted", "uniform")),
    score_power = as.numeric(score_power)[1L],
    score_floor = as.numeric(score_floor)[1L],
    fallback_uniform = isTRUE(fallback_uniform),
    response_blocks = response_blocks
  )
}

# Build the Eq. 6-8 matrix used by response-stratified reconstruction. Scores
# are normalized within response blocks, and the response-level matrices are
# then averaged uniformly. A single global normalization over all directed
# connections would overweight response blocks having more models or larger
# raw modularity values.
build_response_fused_weight_from_cache <- function(cache_list, top_v,
                                                connection_score = NULL,
                                                model_list = NULL,
                                                response_blocks = NULL,
                                                recon_fusion = c("weighted", "uniform"),
                                                score_power = 1,
                                                score_floor = 0,
                                                fallback_uniform = TRUE,
                                                keep_ties = TRUE) {
  recon_fusion <- match.arg(recon_fusion)
  model_names <- names(cache_list)
  if (is.null(model_list)) model_list <- cache_list
  if (is.null(names(model_list))) names(model_list) <- model_names

  scores <- match_model_scores(
    model_names = model_names,
    connection_score = connection_score,
    model_list = model_list
  )
  if (identical(recon_fusion, "uniform")) {
    scores[] <- 1
  } else {
    scores <- pmax(scores, score_floor)^score_power
  }
  names(scores) <- model_names

  pairs <- lapply(seq_along(model_names), function(i) {
    parse_model_pair(model_names[[i]], model = model_list[[model_names[[i]]]])
  })
  responses <- vapply(pairs, function(pair) pair[[1L]], character(1))
  if (is.null(response_blocks)) response_blocks <- unique(responses)
  missing_response <- setdiff(response_blocks, responses)
  if (length(missing_response)) {
    stop(
      "Cannot construct Eq. 6-8 fusion: no fitted connection has response block(s): ",
      paste(missing_response, collapse = ", "), "."
    )
  }

  W_by_response <- lapply(response_blocks, function(response) {
    response_models <- model_names[responses == response]
    alpha <- normalize_fusion_weights(
      scores[response_models], fallback_uniform = fallback_uniform
    )
    names(alpha) <- response_models
    matrices <- lapply(response_models, function(model_name) {
      materialize_weight_from_cache(
        cache_list[[model_name]], v = top_v, keep_ties = keep_ties
      )
    })
    names(matrices) <- response_models
    fuse_matrix_list(matrices, alpha)
  })
  Reduce("+", W_by_response) / length(W_by_response)
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
  if (v > cache$vmax) {
    stop(
      "Requested top-v (", v, ") exceeds the cache `vmax` (", cache$vmax,
      "); rebuild the top-v cache with a larger `vmax`."
    )
  }
  v_use <- min(v, cache$vmax)

  out <- matrix(0, nrow = nrow(cache$W), ncol = ncol(cache$W), dimnames = dimnames(cache$W))
  if (!isTRUE(keep_ties)) {
    idx <- cache$ord[, seq_len(v_use), drop = FALSE]
    row_idx <- rep(seq_len(nrow(cache$W)), times = ncol(idx))
    out[cbind(row_idx, as.vector(idx))] <- cache$W[cbind(row_idx, as.vector(idx))]
    return(row_normalize_weights(out))
  }

  for (i in seq_len(nrow(cache$W))) {
    thr <- cache$sorted_values[i, v_use]
    keep <- which(cache$W[i, ] >= thr)
    if (length(keep) > 0L) {
      out[i, keep] <- cache$W[i, keep]
    }
  }
  row_normalize_weights(out)
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
                                tau = 0.9,
                                elbow_rel_tol = 0.25,
                                elbow_abs_tol = 1e-4,
                                elbow_min_points = 4L,
                                elbow_min_frac = 0.2,
                                elbow_patience = 2L,
                                elbow_smooth_window = 3L) {
  if (!is.data.frame(df) || nrow(df) == 0L) {
    stop("No tuning results available to summarize.")
  }
  object <- match.arg(
    object,
    c("diss", "silhouette", "eigen", "entropy_elbow", "saturation")
  )

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
    stop("`object = '", object, "'` is only available when tuning fused_top_v or model_top_v.")
  }
  if (!("entropy" %in% names(df))) {
    stop("`object = '", object, "'` requires an `entropy` column in tuning results.")
  }
  if (identical(object, "saturation")) {
    return(select_saturation_summary(df, v_col = v_col, tau = tau))
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

.check_saturation_tau <- function(tau) {
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
      tau <= 0 || tau >= 1) {
    stop("`tau` must be a single numeric value in (0, 1).", call. = FALSE)
  }
  as.numeric(tau)
}

entropy_fraction <- function(entropy, reference, eps = 1e-10) {
  reference <- as.numeric(reference)[1L]
  if (!is.finite(reference) || reference <= eps) {
    return(rep(NA_real_, length(entropy)))
  }
  as.numeric(entropy) / reference
}

# Saturation rule: the smallest v whose entropy reaches `tau` times the
# no-truncation entropy. Between two grid points the crossing is located by
# linear interpolation and rounded up, so the selected value does not depend
# on the grid resolution beyond interpolation error.
select_saturation_summary <- function(df, v_col, tau = 0.9, eps = 1e-10) {
  tau <- .check_saturation_tau(tau)
  v_all <- as.numeric(df[[v_col]])
  no_trunc <- if ("is_no_trunc" %in% names(df)) df$is_no_trunc %in% TRUE else rep(FALSE, nrow(df))
  no_trunc <- no_trunc | !is.finite(v_all)

  finite_rows <- which(is.finite(v_all) & is.finite(df$entropy))
  ref_idx <- if (any(no_trunc & is.finite(df$entropy))) {
    which(no_trunc & is.finite(df$entropy))[1L]
  } else if (length(finite_rows) > 0L) {
    finite_rows[which.max(v_all[finite_rows])]
  } else {
    NA_integer_
  }

  finish <- function(row, v_selected, frac, rule, interpolated) {
    out <- row
    out[[v_col]] <- v_selected
    out$entropy_frac <- frac
    out$saturation_tau <- tau
    out$saturation_target <- if (is.finite(ref_idx)) tau * df$entropy[ref_idx] else NA_real_
    out$saturation_rule <- rule
    out$saturation_interpolated <- interpolated
    rownames(out) <- NULL
    out
  }

  if (!is.finite(ref_idx)) {
    stop("`object = 'saturation'` requires at least one finite entropy value.")
  }
  h_inf <- df$entropy[ref_idx]
  if (!is.finite(h_inf) || h_inf <= eps) {
    # Flat curve: truncation changes nothing, so keep every neighbour.
    return(finish(df[ref_idx, , drop = FALSE], v_all[ref_idx], NA_real_,
                  "flat_entropy", FALSE))
  }

  x <- df[finite_rows, , drop = FALSE]
  vx <- v_all[finite_rows]
  ord <- order(vx)
  x <- x[ord, , drop = FALSE]
  vx <- vx[ord]
  frac <- x$entropy / h_inf
  target <- tau * h_inf

  hit <- which(x$entropy >= target - eps)
  if (length(hit) == 0L) {
    # Only the no-truncation baseline reaches the target.
    return(finish(df[ref_idx, , drop = FALSE], v_all[ref_idx], 1,
                  "not_reached", FALSE))
  }
  idx <- hit[1L]
  if (idx == 1L) {
    return(finish(x[idx, , drop = FALSE], vx[idx], frac[idx],
                  "first_candidate", FALSE))
  }
  if (vx[idx] - vx[idx - 1L] <= 1) {
    return(finish(x[idx, , drop = FALSE], vx[idx], frac[idx], "exact", FALSE))
  }
  h_lo <- x$entropy[idx - 1L]
  h_hi <- x$entropy[idx]
  v_lo <- vx[idx - 1L]
  v_hi <- vx[idx]
  v_star <- if (h_hi > h_lo) {
    v_lo + (target - h_lo) / (h_hi - h_lo) * (v_hi - v_lo)
  } else {
    v_hi
  }
  v_star <- min(max(ceiling(v_star), v_lo + 1), v_hi)
  selected_row <- x[idx, , drop = FALSE]
  h_selected <- if (v_hi > v_lo) {
    h_lo + (v_star - v_lo) / (v_hi - v_lo) * (h_hi - h_lo)
  } else {
    h_hi
  }
  selected_row$entropy <- h_selected
  finish(selected_row, as.numeric(v_star), h_selected / h_inf,
         "interpolated", TRUE)
}

#' Entropy of row-wise top-v truncated weights for every truncation level
#'
#' Computes, in closed form, the mean normalized row entropy that
#' `calc_fused_weight_entropy(postprocess_fused_weight(W, top_v = v))` would
#' return for every `v` from 1 to `ncol(W)`, plus the no-truncation value.
#' Each row is sorted once; the entropy of its top-`v` renormalized weights is
#' `log(S_v) - T_v / S_v` with `S_v` the cumulative weight and `T_v` the
#' cumulative `w * log(w)`, so the whole curve costs one sort per row.
#'
#' @param W Square numeric weight matrix (rows are samples).
#' @param keep_ties Logical; whether truncation keeps ties at the cutoff, as
#'   `truncate_top_v_rows()` does by default.
#' @param eps Weights at or below this value are treated as zero.
#' @return A list with `v` (integer vector `1:ncol(W)`), `entropy` (mean
#'   normalized row entropy at each `v`), and `entropy_inf` (no truncation).
#' @keywords internal
fused_entropy_curve <- function(W, keep_ties = TRUE, eps = 1e-12) {
  W <- as.matrix(W)
  if (!is.numeric(W) || length(dim(W)) != 2L) {
    stop("`W` must be a numeric matrix.", call. = FALSE)
  }
  n <- nrow(W)
  p <- ncol(W)
  if (n == 0L || p <= 1L) {
    return(list(v = seq_len(p), entropy = rep(0, p), entropy_inf = 0))
  }
  log_p <- log(p)
  total <- numeric(p)

  for (i in seq_len(n)) {
    w <- W[i, ]
    w[!is.finite(w)] <- 0
    w <- w[w > eps]
    m <- length(w)
    if (m <= 1L) {
      next
    }
    w <- sort(w, decreasing = TRUE)
    s <- cumsum(w)
    t <- cumsum(w * log(w))
    h <- (log(s) - t / s) / log_p
    h[1L] <- 0
    if (isTRUE(keep_ties)) {
      # Truncating at v keeps every entry >= w[v]; map v to the last index of
      # its tie group so tied neighbours at the cutoff are retained.
      last_in_group <- m + 1L - match(w, rev(w))
      h <- h[last_in_group]
    }
    if (m < p) {
      h <- c(h, rep(h[m], p - m))
    }
    total <- total + h
  }

  entropy <- unname(pmin(pmax(total / n, 0), 1))
  list(
    v = seq_len(p),
    entropy = entropy,
    entropy_inf = entropy[p]
  )
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

subset_rfit_to_samples <- function(rfit, sample_names = NULL, sample_idx = NULL,
                                   n_total = NULL) {
  # The positional path is exact when the matrix rows are still in the
  # original data order; name-based matching breaks silently on duplicated
  # rownames, so it is kept only as a checked fallback.
  check_name_match <- function(keep) {
    if (sum(keep) != length(sample_names)) {
      stop(
        "Cannot subset model matrices by sample names: matched ", sum(keep),
        " rows for ", length(sample_names),
        " subsampled names (duplicated or mismatched rownames)."
      )
    }
    keep
  }
  subset_square <- function(mat) {
    if (is.null(mat)) return(NULL)
    if (!is.null(sample_idx) && !is.null(n_total) && nrow(mat) == n_total) {
      return(mat[sample_idx, sample_idx, drop = FALSE])
    }
    if (!is.null(sample_names) && !is.null(rownames(mat))) {
      keep <- check_name_match(rownames(mat) %in% sample_names)
      return(mat[keep, keep, drop = FALSE])
    }
    if (!is.null(sample_idx)) {
      return(mat[sample_idx, sample_idx, drop = FALSE])
    }
    mat
  }
  subset_rect <- function(mat) {
    if (is.null(mat)) return(NULL)
    if (!is.null(sample_idx) && !is.null(n_total) && nrow(mat) == n_total) {
      return(mat[sample_idx, , drop = FALSE])
    }
    if (!is.null(sample_names) && !is.null(rownames(mat))) {
      keep <- check_name_match(rownames(mat) %in% sample_names)
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

  if (!is.null(sample_frac)) {
    if (!is.numeric(sample_frac) || length(sample_frac) != 1L ||
        !is.finite(sample_frac) || sample_frac <= 0 || sample_frac > 1) {
      stop("`sample_frac` must be a single numeric value in (0, 1].", call. = FALSE)
    }
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
  n_target <- as.integer(min(n, n_target))
  if (n_target >= n) {
    return(list(dat.list = dat.list, mod = mod))
  }
  if (n_target < 2L) {
    stop(
      "`sample_n`/`sample_frac` implies a degenerate tuning subset (",
      n_target, " < 2 samples).",
      call. = FALSE
    )
  }

  # Subsample under `seed` without clobbering the caller's global RNG state.
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(list = ".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
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
      sample_idx = idx,
      n_total = n
    )
  }
  if (!is.null(mod_sub$recon)) {
    mod_sub$recon <- NULL
  }

  message("Tuning speed-up: subsampled ", n_target, "/", n, " samples.")
  list(dat.list = dat_sub, mod = mod_sub)
}

eval_grid <- function(grid, eval_fun, parallel = FALSE, cores = NULL) {
  eval_fun <- match.fun(eval_fun)
  ## One parallel layer at a time: each candidate evaluation already uses
  ## thread-level forest parallelism, and the PSOCK workers
  ## must serialize the full closure (model + data) — measured slower
  ## than serial evaluation on realistic sizes. Outer process-level
  ## parallelism therefore requires an explicit `cores`.
  if (!isTRUE(parallel) || length(grid) <= 1L || is.null(cores)) {
    return(lapply(grid, eval_fun))
  }
  evaluated <- parallel_lapply_psock(
    grid,
    function(candidate) {
      tryCatch(
        list(value = eval_fun(candidate), error = NULL),
        error = function(e) list(value = NULL, error = conditionMessage(e))
      )
    },
    cores = cores
  )
  failed <- vapply(
    evaluated,
    function(result) !is.null(result$error),
    logical(1)
  )
  if (any(failed)) {
    first <- which(failed)[1L]
    stop(
      "Parallel grid evaluation failed for candidate value(s) ",
      paste(unlist(grid[failed]), collapse = ", "),
      ": ", evaluated[[first]]$error,
      call. = FALSE
    )
  }
  lapply(evaluated, `[[`, "value")
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

  if (min_eval > n - 1L) {
    return(list(idx = min_eval, rule = "min_eval_fallback"))
  }
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
                                  frac = 1.0) {
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
  # Enforce `min_frac` against the FULL candidate grid, not the evaluated
  # prefix, so early-stop cannot trigger before the documented fraction of
  # the grid has been evaluated.
  min_eval_full <- min(
    max(
      as.integer(min_points),
      as.integer(ceiling(as.numeric(min_frac) * length(grid))),
      2L
    ),
    length(grid)
  )

  for (i in seq_along(grid)) {
    rows[[i]] <- eval_fun(grid[[i]])
    n_eval <- i
    entropy_vals[i] <- as.numeric(rows[[i]]$entropy)[1]

    if (i < 2L || i < min_eval_full) {
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
