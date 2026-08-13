#' MRF variable selection
#'
#' @param mod A fitted `mrf3_fit` object, or an `mrf3`-like object that already
#'   contains raw forest IMD in `mod$imd`.
#' @param dat.list A named list of omics matrices used for feature selection
#'   and optional refitting.
#' @param method Feature-selection rule. `"filter"` adaptively tunes the
#'   cutoff `tau * sd(IMD)` using OOB error; `"mixture"` fits a
#'   point-mass/two-component model; and `"test"` (alias
#'   `"transformation"`) implements the Eq. 16 Student-t transformation.
#'   `"thres"` applies a fixed `se * sd(IMD)` cutoff.
#' @param signal Which signal component to select on: `"shared"` (cross-modal),
#'   `"specific"` (per-block residual), or `"all"` (the union of both).
#' @param se The fixed value of tau used by `method = "thres"`.
#' @param c1 Distribution family for the lower-IMD (noise) component in
#'   mixture mode: `"normal"` or `"truncn"`.
#' @param c2 Distribution family for the higher-IMD component in mixture mode:
#'   `"normal"` or `"gamma"`.
#' @param level One-sided significance level for the transformation, or the
#'   maximum posterior noise probability for mixture selection. A single
#'   numeric value is applied to every block; a named vector/list supplies
#'   block-specific values. `"auto"` retains the optional Bayesian-FDR rule in
#'   mixture mode.
#' @param tscore Deprecated compatibility switch. Use
#'   `method = "transformation"` instead.
#' @param use_distribution Deprecated compatibility argument.
#' @param re_weights Logical; whether selected IMD values are supplied as
#'   variable-sampling weights during the final refit.
#' @param re_fit Logical; whether to refit forests after feature selection.
#' @param ntree Number of trees used in the optional final refit. `NULL` reuses
#'   the fitted object's tree count, falling back to 300.
#' @param scale Logical; whether to standardize selected data before refitting.
#' @param k Number of repeated forest fits at each candidate filtering cutoff.
#' @param tol Tolerable adjacent change in mean OOB normalized MSE used to
#'   choose the filtering cutoff.
#' @param iter Maximum EM iterations in mixture mode.
#' @param eps EM convergence tolerance in mixture mode.
#' @param normalized Logical; whether to L2-normalize the retained weights
#'   after selection. These methods are always evaluated on raw IMD;
#'   therefore the default is `FALSE`.
#' @param select Data block(s) to which selection is applied; default `"ALL"`.
#' @param ... Additional arguments passed to forest fitting.
#'
#' @return An object inheriting from `mrf3` and `vs`. For `signal = "all"`,
#'   `dat.list` and `imd` contain the union selection and the two component
#'   results are retained in `signal_results`.
#' @export
mrf3_vs <- function(mod,
                    dat.list = NULL,
                    method = "filter",
                    signal = c("shared", "specific", "all"),
                    se = 1,
                    c1 = "normal",
                    c2 = "normal",
                    level = 0.05,
                    tscore = FALSE,
                    use_distribution = TRUE,
                    re_weights = FALSE,
                    re_fit = TRUE,
                    ntree = NULL,
                    scale = FALSE,
                    k = 3,
                    tol = 0.01,
                    iter = 1000,
                    eps = 1e-05,
                    normalized = FALSE,
                    select = "ALL",
                    ...) {

  signal <- match.arg(signal)
  method_requested <- match.arg(
    as.character(method)[1L],
    c("filter", "mixture", "test", "transformation", "thres")
  )
  method <- if (identical(method_requested, "transformation")) "test" else method_requested

  if (isTRUE(tscore)) {
    warning(
      "`tscore` is deprecated; using the Eq. 16 Student-t transformation.",
      call. = FALSE
    )
    method <- "test"
    method_requested <- "transformation"
  }
  select_all <- length(select) == 1L && identical(as.character(select), "ALL")
  if (!select_all) {
    select <- unique(as.character(select))
    if (!length(select) || anyNA(select) || any(!nzchar(select))) {
      stop("`select` must be \"ALL\" or one or more valid block names.")
    }
  }
  block_is_selected <- function(block_name) select_all || block_name %in% select

  # `signal = "all"` is a usable union result, rather than a bare pair of
  # unrelated objects. Fit at most once, after the two selections are merged.
  if (identical(signal, "all")) {
    if (identical(method, "filter")) {
      warning(
        "Adaptive OOB filtering is not defined for residual-specific IMD; ",
        "the shared component will use adaptive OOB filtering and the ",
        "specific component will use the fixed cutoff `se * sd(IMD)`.",
        call. = FALSE
      )
      shared_method <- "filter"
      specific_method <- "thres"
    } else {
      shared_method <- specific_method <- method_requested
    }
    source_dat <- dat.list
    if (is.null(source_dat) && inherits(mod, "mrf3_fit")) source_dat <- mod$data
    if (is.null(source_dat)) {
      stop("`signal = 'all'` requires `dat.list` or data retained in `mod$data`.")
    }

    cl <- match.call(expand.dots = TRUE)
    cl$method <- shared_method
    cl$signal <- "shared"
    cl$re_fit <- FALSE
    shared_res <- eval(cl, parent.frame())
    cl$method <- specific_method
    cl$signal <- "specific"
    specific_res <- eval(cl, parent.frame())

    dat_names <- names(source_dat)
    union_vars <- lapply(dat_names, function(block) {
      unique(c(
        colnames(shared_res$dat.list[[block]]),
        colnames(specific_res$dat.list[[block]])
      ))
    })
    names(union_vars) <- dat_names
    union_vars <- lapply(dat_names, function(block) {
      cols <- colnames(source_dat[[block]])
      cols[cols %in% union_vars[[block]]]
    }) |>
      stats::setNames(dat_names)

    combined_dat <- lapply(dat_names, function(block) {
      source_dat[[block]][, union_vars[[block]], drop = FALSE]
    }) |>
      stats::setNames(dat_names)

    combine_weight <- function(block) {
      cols <- colnames(source_dat[[block]])
      a <- .vs_align_vector(shared_res$imd[[block]], cols, fill = 0)
      b <- .vs_align_vector(specific_res$imd[[block]], cols, fill = 0)
      out <- pmax(a, b)
      out[!cols %in% union_vars[[block]]] <- 0
      out
    }
    combined_weights <- lapply(dat_names, combine_weight) |>
      stats::setNames(dat_names)

    out <- shared_res
    out$imd <- combined_weights
    out$dat.list <- if (isTRUE(scale)) {
      lapply(combined_dat, function(x) base::scale(x))
    } else {
      combined_dat
    }
    out$selected_vars <- union_vars
    out$selected_weights <- lapply(dat_names, function(block) {
      combined_weights[[block]][union_vars[[block]]]
    }) |>
      stats::setNames(dat_names)
    out$signal <- "all"
    out$selection_method <- method_requested
    out$signal_results <- list(shared = shared_res, specific = specific_res)
    out$thres <- list(shared = shared_res$thres, specific = specific_res$thres)
    out$refit_performed <- FALSE

    if (isTRUE(re_fit)) {
      if (any(vapply(out$dat.list, ncol, integer(1)) == 0L)) {
        warning(
          "The union selection contains an empty block; the final refit was skipped.",
          call. = FALSE
        )
      } else {
        refit_weights <- NULL
        if (isTRUE(re_weights)) {
          refit_weights <- lapply(dat_names, function(block) {
            combined_weights[[block]][colnames(out$dat.list[[block]])]
          }) |>
            stats::setNames(dat_names)
        }
        refit <- fit_multi_forest(
          out$dat.list,
          connect_list = out$connection,
          ntree = out$ntree,
          type = out$type,
          var.wt = refit_weights,
          ytry = out$ytry,
          ...
        )
        out$mod <- refit
        out$oob_err <- .vs_mean_model_oob_nmse(refit)
        out$refit_performed <- TRUE
      }
    }
    class(out) <- unique(c("mrf3_vs_all", class(out), "vs", "mrf3"))
    return(out)
  }

  if (inherits(mod, "mrf3_fit")) {
    wf <- mod
    if (is.null(dat.list)) dat.list <- wf$data
    if (is.null(dat.list)) {
      stop(
        "`mrf3_vs()` with `mrf3_fit` requires data in `mod$data`. ",
        "Run `mrf3_fit(..., return_data = TRUE)` or provide `dat.list`."
      )
    }

    if (identical(signal, "specific")) {
      spec <- wf$specific
      if (is.null(spec) || is.null(spec$imd)) {
        stop(
          "`signal = 'specific'` requires residual-forest IMD in ",
          "`mod$specific$imd`."
        )
      }
      wf_weights <- spec$imd
      wf_weights_ls <- NULL
      # Selecting on residual-specific IMD still refits the original
      # cross-omics forest structure; self-connections are used only to map the
      # per-block transformation records below.
      refit_connection_vs <- wf$connection
      spec_pt <- spec$imd_per_tree
      if (is.list(spec_pt) && length(spec_pt)) {
        block_names <- names(spec_pt)
        wf_weights_ls <- lapply(block_names, function(block) {
          pt_mat <- spec_pt[[block]]
          if (is.null(pt_mat)) return(NULL)
          pt_mat <- as.matrix(pt_mat)
          lapply(seq_len(ncol(pt_mat)), function(tree) {
            ans <- list(pt_mat[, tree])
            names(ans) <- block
            ans
          })
        })
        names(wf_weights_ls) <- block_names
        wf_weights_ls <- Filter(Negate(is.null), wf_weights_ls)
      }
      connect_list_vs <- lapply(names(wf_weights), function(block) c(block, block))
    } else {
      wf_weights <- wf$imd
      wf_weights_ls <- wf$imd_init
      connect_list_vs <- wf$connection
    }

    if (is.null(wf_weights)) {
      stop(
        "`mrf3_vs()` requires raw IMD weights. Run ",
        "`mrf3_fit(..., run_imd = TRUE)` first."
      )
    }

    mod <- list(
      imd = wf_weights,
      imd_ls = wf_weights_ls,
      connection = connect_list_vs,
      ytry = wf$config$ytry,
      ntree = wf$config$ntree,
      type = wf$type,
      oob_err = wf$oob_err,
      mod = wf$models,
      refit_connection = if (exists("refit_connection_vs", inherits = FALSE)) {
        refit_connection_vs
      } else {
        NULL
      }
    )
    class(mod) <- "mrf3"
  }

  if (is.null(dat.list) || !is.list(dat.list) || !length(dat.list)) {
    stop("`dat.list` must be a non-empty named list.")
  }
  if (is.null(names(dat.list)) || any(!nzchar(names(dat.list)))) {
    stop("`dat.list` must have non-empty block names.")
  }
  dat_names <- names(dat.list)
  if (!select_all) {
    unknown_select <- setdiff(select, dat_names)
    if (length(unknown_select)) {
      stop("Unknown block(s) in `select`: ", paste(unknown_select, collapse = ", "))
    }
  }

  weights <- mod$imd
  if (is.null(weights)) weights <- mod$weights
  if (is.null(weights)) stop("`mod` does not contain IMD weights in `$imd`.")

  # Align once by feature name. The selection methods below operate on raw
  # Eq. 8 forest IMD, whose domain is [0, 1].
  weights <- lapply(dat_names, function(block) {
    if (is.null(weights[[block]])) {
      stop("Missing IMD weights for block `", block, "`.")
    }
    out <- .vs_align_vector(weights[[block]], colnames(dat.list[[block]]), fill = 0)
    if (any(!is.finite(out))) stop("Non-finite IMD in block `", block, "`.")
    if (any(out < -sqrt(.Machine$double.eps)) ||
        any(out > 1 + sqrt(.Machine$double.eps))) {
      stop("Raw forest IMD must lie in [0, 1]; invalid values in block `", block, "`.")
    }
    out[out < 0] <- 0
    out[out > 1] <- 1
    out
  }) |>
    stats::setNames(dat_names)

  ntree_use <- ntree
  if (is.null(ntree_use)) ntree_use <- mod$ntree
  if (is.null(ntree_use) || !length(ntree_use) || !is.finite(ntree_use[1L])) {
    ntree_use <- 300L
  }
  ntree_use <- as.integer(ntree_use[1L])
  if (ntree_use < 1L) stop("`ntree` must be positive.")

  message("Variable selection..")
  thres <- NULL

  if (identical(method, "thres")) {
    thres <- chooss_thres3(weights, se = se)
  }

  if (identical(method, "test")) {
    if (is.null(mod$imd_ls)) {
      stop(
        "The IMD transformation requires per-tree IMD in `mod$imd_ls`; ",
        "it cannot be reconstructed from forest means."
      )
    }
    thres <- test_fn(
      wl = mod$imd_ls,
      connection = mod$connection,
      dat_names = dat_names,
      sig.thres = level,
      feature_names = lapply(dat.list, colnames)
    )
  }

  if (identical(method, "mixture")) {
    c1 <- match.arg(c1, c("normal", "truncn"))
    c2 <- match.arg(c2, c("normal", "gamma"))
    thres <- lapply(dat_names, function(block) {
      .fit_imd_mixture(
        weights[[block]],
        iter = iter,
        eps = eps,
        c1 = c1,
        c2 = c2
      )
    }) |>
      stats::setNames(dat_names)
  }

  if (identical(method, "filter")) {
    if (identical(signal, "specific")) {
      warning(
        "Adaptive OOB filtering is defined for cross-modal forests; ",
        "using the fixed cutoff `se * sd(IMD)` for specific IMD.",
        call. = FALSE
      )
      method <- "thres"
      method_requested <- "thres"
      thres <- chooss_thres3(weights, se = se)
    } else {
      thres <- choose_thres2(
        weights = weights,
        connection = mod$connection,
        new_dat = dat.list,
        ytry = mod$ytry,
        ntree = ntree_use,
        type = mod$type,
        oob_init = mod$oob_err,
        models_init = mod$mod,
        k = k,
        select = select,
        tol = tol,
        ...
      )
    }
  }

  weights_new <- lapply(dat_names, function(block) {
    w <- weights[[block]]
    if (!block_is_selected(block)) return(w)

    if (method %in% c("filter", "thres")) {
      cutoff <- unname(thres[[block]])
      keep <- is.finite(w) & w > cutoff
      w[!keep] <- 0
    } else if (identical(method, "test")) {
      keep <- thres$keep_idx[[block]]
      keep <- .vs_align_vector(keep, names(w), fill = 1) > 0
      w[!keep] <- 0
    } else if (identical(method, "mixture")) {
      post <- thres[[block]]
      # The explicit zero-mass component is part of the null/noise class.
      prob_noise <- post[, "noise"] + post[, "zero"]
      if (identical(level, "auto")) {
        target_fdr <- 0.05
        positive <- which(w > 0)
        keep <- rep(FALSE, length(w))
        if (length(positive)) {
          ord <- positive[order(prob_noise[positive])]
          cumulative_fdr <- cumsum(prob_noise[ord]) / seq_along(ord)
          n_keep <- max(c(0L, which(cumulative_fdr <= target_fdr)))
          if (n_keep > 0L) keep[ord[seq_len(n_keep)]] <- TRUE
        }
      } else {
        pr <- .vs_level_for_block(level, block)
        keep <- w > 0 & prob_noise < pr
      }
      w[!keep] <- 0
    }

    if (isTRUE(normalized)) {
      denom <- sqrt(sum(w^2))
      if (is.finite(denom) && denom > 0) w <- w / denom
    }
    w
  }) |>
    stats::setNames(dat_names)

  new_dat <- lapply(dat_names, function(block) {
    if (!block_is_selected(block)) return(dat.list[[block]])
    keep <- names(weights_new[[block]])[weights_new[[block]] > 0]
    keep <- colnames(dat.list[[block]])[colnames(dat.list[[block]]) %in% keep]
    dat.list[[block]][, keep, drop = FALSE]
  }) |>
    stats::setNames(dat_names)

  new_dat_fit <- if (isTRUE(scale)) {
    lapply(new_dat, function(x) base::scale(x))
  } else {
    new_dat
  }

  mod$imd <- weights_new
  mod$weights <- NULL
  mod$dat.list <- new_dat_fit
  mod$selected_vars <- lapply(new_dat, colnames)
  mod$selected_weights <- lapply(dat_names, function(block) {
    weights_new[[block]][colnames(new_dat[[block]])]
  }) |>
    stats::setNames(dat_names)
  mod$thres <- thres
  mod$ntree <- ntree_use
  mod$signal <- signal
  mod$selection_method <- method_requested
  mod$refit_performed <- FALSE

  if (isTRUE(re_fit)) {
    refit_connection <- if (!is.null(mod$refit_connection)) {
      mod$refit_connection
    } else {
      mod$connection
    }
    if (is.null(refit_connection) || !length(refit_connection)) {
      warning("Model refitting was skipped because no connection is available.", call. = FALSE)
    } else if (any(vapply(new_dat_fit, ncol, integer(1)) == 0L)) {
      warning("At least one selected block is empty; model refitting was skipped.", call. = FALSE)
    } else {
      message("Refit model..")
      refit_weights <- NULL
      if (isTRUE(re_weights)) {
        refit_weights <- lapply(dat_names, function(block) {
          if (!block_is_selected(block)) return(NULL)
          weights_new[[block]][colnames(new_dat_fit[[block]])]
        }) |>
          stats::setNames(dat_names)
      }
      refit <- fit_multi_forest(
        new_dat_fit,
        connect_list = refit_connection,
        ntree = ntree_use,
        type = mod$type,
        var.wt = refit_weights,
        ytry = mod$ytry,
        ...
      )
      mod$mod <- refit
      mod$oob_err <- .vs_mean_model_oob_nmse(refit)
      mod$refit_performed <- TRUE
    }
  }

  class(mod) <- unique(c(class(mod), "vs", "mrf3"))
  mod
}

# Align feature vectors by name without recycling. If names are unavailable,
# positional alignment is accepted only when lengths agree.
.vs_align_vector <- function(x, target_names, fill = NA_real_) {
  target_names <- as.character(target_names)
  out <- rep(fill, length(target_names))
  names(out) <- target_names
  if (is.null(x)) return(out)
  x_names <- names(x)
  x <- as.vector(x)
  if (is.null(x_names) || !any(nzchar(x_names))) {
    if (length(x) != length(out)) {
      stop("Unnamed feature vector has length ", length(x),
           ", expected ", length(out), ".")
    }
    out[] <- x
    return(out)
  }
  idx <- match(target_names, x_names)
  present <- !is.na(idx)
  out[present] <- x[idx[present]]
  out
}

.vs_level_for_block <- function(level, block, default = 0.05) {
  if (identical(level, "auto")) return(default)
  if (is.list(level)) level <- unlist(level, use.names = TRUE)
  if (!is.numeric(level) || !length(level)) {
    stop("`level` must be numeric, a named numeric vector/list, or \"auto\".")
  }
  value <- if (!is.null(names(level)) && block %in% names(level)) {
    level[[block]]
  } else {
    level[[1L]]
  }
  if (length(value) != 1L || !is.finite(value) || value <= 0 || value >= 1) {
    stop("Selection level for block `", block, "` must lie strictly between 0 and 1.")
  }
  as.numeric(value)
}

# Mean feature-level normalized OOB MSE, averaged again across directed
# connections. `get_oob_nmse()` already averages all p + q coordinates.
.vs_mean_model_oob_nmse <- function(models) {
  if (is.null(models) || !is.list(models) || !length(models)) return(NA_real_)
  score <- vapply(models, function(model) {
    value <- tryCatch(get_oob_nmse(model), error = function(e) NA_real_)
    as.numeric(value)
  }, numeric(1))
  if (!any(is.finite(score))) return(NA_real_)
  mean(score[is.finite(score)])
}

#' @keywords internal
choose_thres2 <- function(weights, connection, new_dat, ytry, ntree, type,
                          oob_init = NULL, models_init = NULL, k = 3,
                          tol = 0.01, select = "ALL",
                          tau_grid = seq(0.8, 3.1, by = 0.1),
                          .fit_fun = fit_multi_forest,
                          .score_fun = .vs_mean_model_oob_nmse, ...) {
  k <- as.integer(k)
  if (k < 1L) stop("`k` must be positive.")
  tau_grid <- as.numeric(tau_grid)
  if (!length(tau_grid) || any(!is.finite(tau_grid)) || any(tau_grid < 0)) {
    stop("`tau_grid` must contain non-negative finite values.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("`tol` must be one non-negative finite number.")
  }
  select_all <- length(select) == 1L && identical(as.character(select), "ALL")
  selected_block <- function(block) select_all || block %in% as.character(select)

  # Kept as an unused formal for source compatibility with pre-alignment code.
  invisible(oob_init)

  dots <- list(...)
  base_seed <- if (!is.null(dots$seed)) as.integer(dots$seed)[1L] else 529L
  dots$seed <- NULL
  models_have_oob_inputs <- is.list(models_init) && length(models_init) &&
    all(vapply(models_init, function(model) {
      !is.null(model$membership) && !is.null(model$inbag)
    }, logical(1)))
  baseline <- if (models_have_oob_inputs) {
    .score_fun(models_init)
  } else {
    NA_real_
  }
  if (!is.finite(baseline)) {
    # Historical `oob_init` values were often computed from all-tree weights,
    # so do not mix them with true OOB scores. Refit the unfiltered data using
    # exactly the same OOB protocol as the candidate cutoffs.
    baseline_repeat <- rep(NA_real_, k)
    for (replicate_id in seq_len(k)) {
      fit_args <- utils::modifyList(
        list(
          dat.list = new_dat,
          connect_list = connection,
          ntree = ntree,
          type = type,
          ytry = ytry,
          var.wt = NULL,
          forest.wt = "oob",
          seed = base_seed + replicate_id - 1L
        ),
        dots
      )
      # OOB scoring is part of the method and cannot be overridden by `...`.
      fit_args$dat.list <- new_dat
      fit_args$connect_list <- connection
      fit_args$forest.wt <- "oob"
      fit_args$seed <- base_seed + replicate_id - 1L
      baseline_fit <- do.call(.fit_fun, fit_args)
      baseline_repeat[[replicate_id]] <- .score_fun(baseline_fit)
    }
    if (any(is.finite(baseline_repeat))) {
      baseline <- mean(baseline_repeat[is.finite(baseline_repeat)])
    }
  }
  candidate_error <- rep(Inf, length(tau_grid))

  for (i in seq_along(tau_grid)) {
    tau <- tau_grid[[i]]
    keep_by_block <- lapply(names(weights), function(block) {
      if (!selected_block(block)) return(rep(TRUE, length(weights[[block]])))
      cutoff <- tau * stats::sd(weights[[block]])
      if (!is.finite(cutoff)) cutoff <- 0
      weights[[block]] > cutoff
    }) |>
      stats::setNames(names(weights))

    dat_reduced <- lapply(names(new_dat), function(block) {
      if (!selected_block(block)) return(new_dat[[block]])
      keep_names <- names(weights[[block]])[keep_by_block[[block]]]
      keep_names <- colnames(new_dat[[block]])[colnames(new_dat[[block]]) %in% keep_names]
      new_dat[[block]][, keep_names, drop = FALSE]
    }) |>
      stats::setNames(names(new_dat))

    required_blocks <- unique(unlist(connection, use.names = FALSE))
    if (any(vapply(dat_reduced[required_blocks], ncol, integer(1)) == 0L)) next

    repeat_error <- rep(NA_real_, k)
    for (replicate_id in seq_len(k)) {
      fit_args <- utils::modifyList(
        list(
          dat.list = dat_reduced,
          connect_list = connection,
          ntree = ntree,
          type = type,
          ytry = ytry,
          var.wt = NULL,
          forest.wt = "oob",
          seed = base_seed + replicate_id - 1L
        ),
        dots
      )
      fit_args$dat.list <- dat_reduced
      fit_args$connect_list <- connection
      fit_args$forest.wt <- "oob"
      fit_args$seed <- base_seed + replicate_id - 1L
      refit <- do.call(.fit_fun, fit_args)
      repeat_error[[replicate_id]] <- .score_fun(refit)
    }
    if (any(is.finite(repeat_error))) {
      candidate_error[[i]] <- mean(repeat_error[is.finite(repeat_error)])
    }
  }

  # Adaptive filtering looks for the first stable plateau between
  # adjacent candidate cutoffs. The unfiltered baseline is diagnostic only: it
  # may be a single existing fit and must not be compared with a k-fit mean.
  # The first grid point is the companion algorithm's anchor (`oob[-1]`), so a
  # plateau can start at the second candidate once two subsequent candidates
  # have tolerably similar errors.
  plateau_pool <- seq_along(candidate_error)[-1L]
  stable_start <- integer()
  if (length(plateau_pool) >= 2L) {
    left <- plateau_pool[-length(plateau_pool)]
    right <- plateau_pool[-1L]
    stable_start <- left[
      is.finite(candidate_error[left]) &
        is.finite(candidate_error[right]) &
        abs(candidate_error[right] - candidate_error[left]) <= tol
    ]
  }
  if (length(stable_start)) {
    chosen_index <- stable_start[[1L]]
  } else if (any(is.finite(candidate_error))) {
    chosen_index <- which.min(candidate_error)
  } else {
    stop("No filtering cutoff produced a fit with a finite OOB error.")
  }
  tau <- tau_grid[[chosen_index]]
  message("Choose ", format(tau), " times sd")

  thresholds <- vapply(weights, function(w) {
    value <- tau * stats::sd(w)
    if (is.finite(value)) value else 0
  }, numeric(1))
  attr(thresholds, "tau") <- tau
  attr(thresholds, "oob_trace") <- data.frame(
    tau = tau_grid,
    mean_oob_nmse = candidate_error
  )
  attr(thresholds, "baseline_oob_nmse") <- baseline
  thresholds
}

# Fixed threshold tau * sd(IMD). Kept under the historical helper name
# for compatibility with code that calls it via `multiRF:::`.
chooss_thres3 <- function(weights, se) {
  if (!is.numeric(se) || length(se) != 1L || !is.finite(se) || se < 0) {
    stop("`se` (tau) must be one non-negative finite number.")
  }
  out <- vapply(weights, function(w) {
    value <- se * stats::sd(w)
    if (is.finite(value)) value else 0
  }, numeric(1))
  names(out) <- names(weights)
  out
}

# Fit an explicit point mass at zero plus a two-component model on
# positive IMD. Components are relabelled after fitting so column `noise`
# always refers to the lower-mean component.
.fit_imd_mixture <- function(x, iter = 1000, eps = 1e-5,
                             c1 = c("normal", "truncn"),
                             c2 = c("normal", "gamma")) {
  c1 <- match.arg(c1)
  c2 <- match.arg(c2)
  x_names <- names(x)
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(x < 0) || any(x > 1)) {
    stop("Mixture input must be finite raw IMD in [0, 1].")
  }
  n <- length(x)
  post <- matrix(0, nrow = n, ncol = 3L,
                 dimnames = list(x_names, c("noise", "signal", "zero")))
  zero <- x == 0
  post[zero, "zero"] <- 1
  xp <- x[!zero]
  if (!length(xp)) return(post)
  if (length(xp) == 1L || stats::sd(xp) < sqrt(.Machine$double.eps)) {
    warning(
      "Positive IMD values do not contain enough variation to identify a ",
      "two-component mixture; conservatively treating them as noise.",
      call. = FALSE
    )
    post[!zero, "noise"] <- 1
    attr(post, "degenerate") <- TRUE
    attr(post, "component_means") <- rep(mean(xp), 2L)
    attr(post, "loglik_trace") <- numeric()
    return(post)
  }

  min_sd <- max(stats::sd(xp) * 1e-3, 1e-6)
  q <- as.numeric(stats::quantile(xp, c(0.3, 0.7), names = FALSE, type = 8))
  mu <- q
  sigma <- rep(max(stats::sd(xp) / 2, min_sd), 2L)
  gamma_shape <- max((mu[2L] / sigma[2L])^2, 1e-3)
  gamma_scale <- max(sigma[2L]^2 / mu[2L], 1e-6)
  mix <- 0.5

  log_density_one <- function(values, mean, sd, family) {
    sd <- max(sd, min_sd)
    if (identical(family, "truncn")) {
      # Eq. 13 is left-truncated at zero (not doubly truncated at one).
      log_survival_zero <- stats::pnorm(
        0, mean = mean, sd = sd, lower.tail = FALSE, log.p = TRUE
      )
      stats::dnorm(values, mean = mean, sd = sd, log = TRUE) -
        log_survival_zero
    } else {
      stats::dnorm(values, mean = mean, sd = sd, log = TRUE)
    }
  }
  log_density_two <- function(values, mean, sd, family, shape, scale) {
    if (identical(family, "gamma")) {
      stats::dgamma(
        values,
        shape = max(shape, 1e-6),
        scale = max(scale, .Machine$double.xmin),
        log = TRUE
      )
    } else {
      stats::dnorm(values, mean = mean, sd = max(sd, min_sd), log = TRUE)
    }
  }

  weighted_normal_mle <- function(values, weights) {
    total <- sum(weights)
    mean <- sum(weights * values) / total
    sd <- sqrt(sum(weights * (values - mean)^2) / total)
    c(mean = mean, sd = max(sd, min_sd))
  }

  weighted_truncn_mle <- function(values, weights, start_mean, start_sd) {
    objective <- function(par) {
      mean <- par[[1L]]
      sd <- exp(par[[2L]])
      log_density <- log_density_one(values, mean, sd, "truncn")
      value <- -sum(weights * log_density)
      if (is.finite(value)) value else .Machine$double.xmax / 100
    }
    initial <- c(start_mean, log(max(start_sd, min_sd)))
    fit <- tryCatch(
      stats::optim(
        initial,
        objective,
        method = "L-BFGS-B",
        lower = c(-20, log(min_sd)),
        upper = c(1, log(10)),
        control = list(maxit = 200L)
      ),
      error = function(e) NULL
    )
    if (is.null(fit) || !is.finite(fit$value)) {
      return(c(mean = start_mean, sd = max(start_sd, min_sd)))
    }
    c(mean = fit$par[[1L]], sd = max(exp(fit$par[[2L]]), min_sd))
  }

  weighted_gamma_mle <- function(values, weights, start_shape, start_scale) {
    total <- sum(weights)
    weighted_mean <- sum(weights * values) / total
    weighted_log_mean <- sum(weights * log(values)) / total
    target <- max(log(weighted_mean) - weighted_log_mean, 0)

    # Weighted gamma MLE: log(k) - digamma(k) = log(E_w X) - E_w log(X).
    if (target <= 1e-12) {
      shape <- 1e6
    } else {
      score <- function(shape) log(shape) - digamma(shape) - target
      lower <- 1e-6
      upper <- max(1, start_shape)
      while (score(upper) > 0 && upper < 1e8) upper <- upper * 2
      shape <- tryCatch(
        stats::uniroot(score, lower = lower, upper = upper, tol = 1e-10)$root,
        error = function(e) start_shape
      )
    }
    shape <- min(max(shape, 1e-6), 1e8)
    # Preserve the weighted MLE mean exactly, including very small raw IMD.
    scale <- max(weighted_mean / shape, .Machine$double.xmin)
    c(shape = shape, scale = scale)
  }

  log_sum_exp_two <- function(left, right) {
    maximum <- pmax(left, right)
    maximum + log(exp(left - maximum) + exp(right - maximum))
  }

  component_mean_one <- function(mean, sd, family) {
    if (!identical(family, "truncn")) return(mean)
    alpha <- -mean / sd
    mean + sd * exp(stats::dnorm(alpha, log = TRUE) -
                      stats::pnorm(alpha, lower.tail = FALSE, log.p = TRUE))
  }
  component_mean_two <- function(mean, family, shape, scale) {
    if (identical(family, "gamma")) shape * scale else mean
  }

  loglik_trace <- numeric()
  rolled_back <- FALSE
  responsibility <- rep(0.5, length(xp))
  for (step in seq_len(as.integer(iter))) {
    log_a <- log(mix) + log_density_one(xp, mu[1L], sigma[1L], c1)
    log_b <- log1p(-mix) + log_density_two(
      xp, mu[2L], sigma[2L], c2, gamma_shape, gamma_scale
    )
    log_denominator <- log_sum_exp_two(log_a, log_b)
    responsibility <- exp(log_a - log_denominator)
    loglik <- sum(log_denominator)
    if (!length(loglik_trace) ||
        abs(loglik - loglik_trace[[length(loglik_trace)]]) > 1e-12) {
      loglik_trace <- c(loglik_trace, loglik)
    }

    weight_one <- sum(responsibility)
    weight_two <- sum(1 - responsibility)
    if (weight_one <= 1e-8 || weight_two <= 1e-8) break

    old_parameters <- list(
      mix = mix,
      mu = mu,
      sigma = sigma,
      gamma_shape = gamma_shape,
      gamma_scale = gamma_scale
    )
    mix <- min(1 - 1e-6, max(1e-6, mean(responsibility)))
    first_fit <- if (identical(c1, "truncn")) {
      weighted_truncn_mle(xp, responsibility, mu[1L], sigma[1L])
    } else {
      weighted_normal_mle(xp, responsibility)
    }
    mu[1L] <- first_fit[["mean"]]
    sigma[1L] <- first_fit[["sd"]]

    if (identical(c2, "gamma")) {
      second_fit <- weighted_gamma_mle(
        xp, 1 - responsibility, gamma_shape, gamma_scale
      )
      gamma_shape <- second_fit[["shape"]]
      gamma_scale <- second_fit[["scale"]]
      mu[2L] <- gamma_shape * gamma_scale
      sigma[2L] <- sqrt(gamma_shape) * gamma_scale
    } else {
      second_fit <- weighted_normal_mle(xp, 1 - responsibility)
      mu[2L] <- second_fit[["mean"]]
      sigma[2L] <- second_fit[["sd"]]
    }

    proposed_a <- log(mix) + log_density_one(xp, mu[1L], sigma[1L], c1)
    proposed_b <- log1p(-mix) + log_density_two(
      xp, mu[2L], sigma[2L], c2, gamma_shape, gamma_scale
    )
    proposed_loglik <- sum(log_sum_exp_two(proposed_a, proposed_b))

    # A numerically capped gamma MLE can cease to be a true maximizer for
    # boundary-scale data. Never accept an M-step that lowers the observed-data
    # likelihood; retain the last valid parameters instead.
    if (!is.finite(proposed_loglik) || proposed_loglik < loglik) {
      mix <- old_parameters$mix
      mu <- old_parameters$mu
      sigma <- old_parameters$sigma
      gamma_shape <- old_parameters$gamma_shape
      gamma_scale <- old_parameters$gamma_scale
      rolled_back <- TRUE
      break
    }
    if (abs(proposed_loglik - loglik) <= eps) {
      if (abs(proposed_loglik - loglik_trace[[length(loglik_trace)]]) > 1e-12) {
        loglik_trace <- c(loglik_trace, proposed_loglik)
      }
      break
    }
  }

  log_a <- log(mix) + log_density_one(xp, mu[1L], sigma[1L], c1)
  log_b <- log1p(-mix) + log_density_two(
    xp, mu[2L], sigma[2L], c2, gamma_shape, gamma_scale
  )
  log_denominator <- log_sum_exp_two(log_a, log_b)
  positive_post <- cbind(exp(log_a - log_denominator),
                         exp(log_b - log_denominator))
  final_loglik <- sum(log_denominator)
  if (!length(loglik_trace) ||
      abs(final_loglik - loglik_trace[[length(loglik_trace)]]) > 1e-12) {
    loglik_trace <- c(loglik_trace, final_loglik)
  }
  component_means <- c(
    component_mean_one(mu[1L], sigma[1L], c1),
    component_mean_two(mu[2L], c2, gamma_shape, gamma_scale)
  )
  if (component_means[[1L]] > component_means[[2L]]) {
    positive_post <- positive_post[, 2:1, drop = FALSE]
  }
  post[!zero, c("noise", "signal")] <- positive_post
  attr(post, "loglik_trace") <- loglik_trace
  attr(post, "component_means") <- sort(component_means)
  attr(post, "em_rollback") <- rolled_back
  post
}

# Eq. 16 transformation for one per-tree IMD matrix (features x trees).
.imd_transformation_one <- function(mat, alpha) {
  mat <- as.matrix(mat)
  if (ncol(mat) < 2L) stop("IMD transformation requires at least two trees.")
  forest_mean <- rowMeans(mat)
  global_mean <- mean(mat)
  standard_error <- apply(mat, 1L, stats::sd) / sqrt(ncol(mat))
  difference <- forest_mean - global_mean
  t_score <- difference / standard_error
  zero_se <- !is.finite(t_score)
  if (any(zero_se)) {
    t_score[zero_se & difference > 0] <- Inf
    t_score[zero_se & difference < 0] <- -Inf
    t_score[zero_se & difference == 0] <- 0
  }
  p_value <- stats::pt(t_score, df = ncol(mat) - 1L, lower.tail = FALSE)
  keep <- as.integer(p_value < alpha)
  names(keep) <- names(p_value) <- names(t_score) <- rownames(mat)
  list(keep_idx = keep, pval = p_value, ts = t_score, ntree = ncol(mat))
}

#' @keywords internal
test_fn <- function(wl, connection, dat_names, sig.thres = 0.05,
                    feature_names = NULL) {
  if (is.null(feature_names)) {
    feature_names <- stats::setNames(vector("list", length(dat_names)), dat_names)
  }
  if (is.null(names(feature_names))) names(feature_names) <- dat_names
  records <- stats::setNames(vector("list", length(dat_names)), dat_names)

  for (i in seq_along(wl)) {
    trees <- wl[[i]]
    if (!is.list(trees) || !length(trees)) next
    valid_trees <- trees[vapply(trees, is.list, logical(1))]
    if (!length(valid_trees)) next
    side_keys <- unique(unlist(lapply(valid_trees, names), use.names = FALSE))
    if (!length(side_keys)) next

    conn <- if (length(connection) >= i) as.character(connection[[i]]) else character()
    if (all(side_keys %in% dat_names)) {
      side_to_block <- stats::setNames(side_keys, side_keys)
    } else {
      mapped <- rev(conn)
      mapped <- mapped[seq_len(min(length(mapped), length(side_keys)))]
      side_to_block <- stats::setNames(mapped, side_keys[seq_along(mapped)])
    }

    for (side in names(side_to_block)) {
      block <- unname(side_to_block[[side]])
      if (!block %in% dat_names) next
      raw <- lapply(valid_trees, function(tree) tree[[side]])
      raw <- raw[vapply(raw, function(x) !is.null(x) && length(x) > 0L, logical(1))]
      if (!length(raw)) next

      target <- feature_names[[block]]
      if (is.null(target) || !length(target)) {
        target <- unique(unlist(lapply(raw, names), use.names = FALSE))
        if (!length(target)) target <- paste0("V", seq_len(length(raw[[1L]])))
      }
      aligned <- lapply(raw, .vs_align_vector, target_names = target, fill = 0)
      mat <- do.call(cbind, aligned)
      rownames(mat) <- target
      alpha <- .vs_level_for_block(sig.thres, block)
      records[[block]][[length(records[[block]]) + 1L]] <-
        .imd_transformation_one(mat, alpha)
    }
  }

  aggregate_block <- function(block) {
    target <- feature_names[[block]]
    block_records <- records[[block]]
    if (!length(block_records)) {
      if (is.null(target)) target <- character()
      warning(
        "No connected forest supplies per-tree IMD for block `", block,
        "`; retaining its variables because the transformation was not evaluated.",
        call. = FALSE
      )
      return(list(
        keep_idx = stats::setNames(rep(1L, length(target)), target),
        pval = stats::setNames(rep(NA_real_, length(target)), target),
        ts = stats::setNames(rep(NA_real_, length(target)), target),
        n_models = 0L
      ))
    }
    if (is.null(target) || !length(target)) target <- names(block_records[[1L]]$keep_idx)
    keep_mat <- do.call(cbind, lapply(block_records, function(x) {
      .vs_align_vector(x$keep_idx, target, fill = 0)
    }))
    p_mat <- do.call(cbind, lapply(block_records, function(x) {
      .vs_align_vector(x$pval, target, fill = NA_real_)
    }))
    t_mat <- do.call(cbind, lapply(block_records, function(x) {
      .vs_align_vector(x$ts, target, fill = NA_real_)
    }))
    mean_na <- function(mat) {
      out <- rowMeans(mat, na.rm = TRUE)
      out[!is.finite(out)] <- NA_real_
      stats::setNames(out, target)
    }
    list(
      # Apply a strict majority across connected forests, row-wise and
      # separately for every feature.
      keep_idx = stats::setNames(as.integer(rowMeans(keep_mat) > 0.5), target),
      pval = mean_na(p_mat),
      ts = mean_na(t_mat),
      n_models = length(block_records)
    )
  }

  aggregated <- lapply(dat_names, aggregate_block) |>
    stats::setNames(dat_names)
  list(
    keep_idx = lapply(aggregated, `[[`, "keep_idx"),
    pval = lapply(aggregated, `[[`, "pval"),
    ts = lapply(aggregated, `[[`, "ts"),
    n_models = vapply(aggregated, `[[`, integer(1), "n_models")
  )
}
