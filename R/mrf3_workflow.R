#' Stage-based multiRF fitting pipeline
#'
#' @param dat.list A named list of omics matrices (samples in rows, features in columns).
#' @param ntree Number of trees for RF fitting.
#' @param scale Logical; whether to z-standardize each feature before fitting.
#' @param ytry Number of response variables sampled per split. `NULL` delegates
#'   to the engine default.
#' @param samptype Sampling scheme passed to forest fitting: `"swor"` or `"swr"`.
#' @param connect_list Optional predefined connections (`list(c(response, predictor), ...)`).
#' If `NULL`, directional connections are selected automatically by `find_connection()`
#' from fitted forests.
#' @param filter_mode Feature filtering mode passed to `filter_omics()`.
#' @param filter_method Feature dispersion metric passed to `filter_omics()`.
#' @param top_n_by_type Optional auto-filter overrides.
#' @param top_n_manual Optional manual top-n filtering configuration.
#' @param filter_verbose Logical; whether to print filtering diagnostics.
#' @param main_clustering Global clustering strategy applied consistently to
#' shared, specific, and robust branches: `"similarity"` (default),
#' `"proximity"`, or `"enhanced_proximity"`.
#' @param shared_specific_args A named list of additional arguments passed to
#' `get_shared_specific_weights()`. When this branch is created, per-omics
#' shared fraction (`1 - ||R||_F^2 / ||X||_F^2`) is also computed by
#' `get_shared_frac()` and attached as `shared_frac`.
#' By default, this function passes `specific_top_v = selected_fused_top_v`
#' into this branch.
#' @param clustering_args A named list of additional arguments passed to
#' the shared/specific clustering stage.
#' @param run_imd Logical; whether to run `get_multi_weights()` as a pipeline stage.
#' @param run_cluster_imd Logical; whether to run `cluster_imd()` after global
#' IMD when `run_imd = TRUE`. If `NULL` (default), inherits from `run_imd`.
#' Set to `FALSE` to compute global IMD weights only without the per-cluster
#' refit, which can be very memory-intensive for sub-MRF models.
#' @param imd_args A named list of additional arguments passed to `get_multi_weights()`.
#' By default, this function sets `parallel = TRUE` for IMD unless overridden here.
#' @param run_variable_selection Logical; whether to run variable selection
#' using `mrf3_vs()` after IMD weights are available.
#' @param variable_selection_args A named list of additional arguments passed to
#' `mrf3_vs()`. `re_fit` is managed by workflow and ignored if provided:
#' it is forced to `TRUE` when `run_robust_clustering = TRUE`, otherwise `FALSE`.
#' @param run_robust_clustering Logical; whether to run a robust clustering
#' branch that first selects variables from IMD weights and then re-clusters.
#' @param robust_clustering_args A named list of clustering arguments
#' (`shared_k`, `specific_k`, etc.) specific to the robust clustering branch.
#' By default, robust clustering inherits the actual `k` chosen in Stage 4;
#' values in this list override those inherited defaults.
#' @param cluster_imd_args A named list of additional arguments passed to
#' `cluster_imd()` when `run_imd = TRUE`.
#' By default, this function reuses cluster labels from
#' `clusters`.
#' @param top_v Optional unified top-v cutoff applied to both `model_top_v`
#' and `fused_top_v`.
#' @param model_top_v Model-level top-v cutoff used on each single-model
#' forest weight matrix before fusion. Default `Inf` = no truncation (use all
#' weights). Set to `NULL` to auto-tune via `tune_model_top_v()`.
#' @param recon_fusion Reconstruction fusion mode passed to `get_reconstr_matrix()`:
#' `"weighted"` (default) or `"uniform"`.
#' @param score_power Exponent applied to connection scores for weighted reconstruction.
#' @param score_floor Non-negative floor applied to connection scores before weighting.
#' @param fallback_uniform Logical; whether weighted reconstruction falls back to uniform averaging when scores are unavailable.
#' @param fused_top_v Row-wise top-v truncation for fused weights. Default `Inf`
#' = no truncation. Set to `NULL` to auto-tune, `FALSE` to skip entirely.
#' @param fused_row_normalize Logical; whether to row-normalize fused weights after optional truncation.
#' @param fused_keep_ties Logical; whether fused top-v truncation keeps ties at cutoff.
#' @param top_v_method Strategy used when auto-selecting `top_v`.
#' @param neff_quantile Quantile of effective neighborhood size used by the
#'   `"neff"` top-v rule.
#' @param model_top_v_tune_args A named list of additional arguments passed to
#' `tune_model_top_v()` (e.g., `tmin`, `by`, `k`).
#' Workflow always uses `object = "entropy_elbow"`.
#' @param fused_top_v_tune_args A named list of additional arguments passed to
#' `tune_fused_top_v()` (e.g., `vmin`, `by`, `vmax`, `k`).
#' Workflow always uses `object = "entropy_elbow"`.
#' @param return_data Logical; whether to include filtered/scaled data in output objects.
#' @param compact_output Logical; whether to drop heavy duplicated objects from
#' output for lower memory usage (for example fitted model copies in robust
#' clustering and full `vs_fit` object). Useful when keeping multiple
#' `mrf3_fit` objects in memory.
#' @param verbose Logical; whether to print stage-level progress messages.
#' @param seed Random seed.
#' @param ... Additional arguments passed to `mrf3_init()`.
#'
#' @return A compact `mrf3_fit` object with top-level components:
#' `config`, `init`, `tuning`, `reconstruction`, `clusters`,
#' `shared`, `specific`, `imd`, `cluster_imd`,
#' `variable_selection`, `robust_clustering`, and optional `data`.
#' @export
mrf3_fit <- function(dat.list,
                          ntree = 500,
                          scale = TRUE,
                          ytry = NULL,
                          samptype = c("swor", "swr"),
                          connect_list = NULL,
                          filter_mode = c("auto", "none", "manual"),
                          filter_method = c("mad", "variance"),
                          top_n_by_type = NULL,
                          top_n_manual = NULL,
                          filter_verbose = TRUE,
                          main_clustering = c("similarity", "proximity", "enhanced_proximity"),
                          shared_specific_args = list(),
                          clustering_args = list(),
                          run_imd = FALSE,
                          run_cluster_imd = FALSE,
                          imd_args = list(),
                          run_variable_selection = FALSE,
                          variable_selection_args = list(),
                          run_robust_clustering = FALSE,
                          robust_clustering_args = list(),
                          cluster_imd_args = list(),
                          top_v = NULL,
                          model_top_v = NULL,
                          recon_fusion = c("weighted", "uniform"),
                          score_power = 1,
                          score_floor = 0,
                          fallback_uniform = TRUE,
                          fused_top_v = NULL,
                          fused_row_normalize = TRUE,
                          fused_keep_ties = TRUE,
                          top_v_method = c("entropy_elbow", "neff"),
                          neff_quantile = 0.5,
                          model_top_v_tune_args = list(),
                          fused_top_v_tune_args = list(),
                          return_data = FALSE,
                          compact_output = FALSE,
                          verbose = TRUE,
                          seed = 529,
                          ...){


  dots <- list(...)
  if (!is.null(dots$auto_tune_model_top_v) || !is.null(dots$auto_tune_fused_top_v)) {
    stop(
      "`auto_tune_model_top_v` and `auto_tune_fused_top_v` have been removed. ",
      "Set `model_top_v`/`fused_top_v` to numeric values to skip tuning; leave them NULL to auto-tune."
    )
  }
  if (!is.null(dots$connection)) dots$connection <- NULL
  if (!is.null(dots$connection_strategy)) dots$connection_strategy <- NULL
  if (!is.null(dots$run_clustering)) dots$run_clustering <- NULL
  if (!is.null(dots$clustering_strategy)) dots$clustering_strategy <- NULL
  if (!is.null(dots$run_shared_specific)) dots$run_shared_specific <- NULL
  if (!is.null(dots$run_specific_clustering)) dots$run_specific_clustering <- NULL
  if (!is.null(variable_selection_args$re_fit)) variable_selection_args$re_fit <- NULL
  ## robust_clustering_args is now a formal parameter; ignore stale dots entry
  if (!is.null(dots$robust_clustering_args)) dots$robust_clustering_args <- NULL
  if (!is.null(dots$clustering_args)) dots$clustering_args <- NULL
  if (is.character(main_clustering) && length(main_clustering) == 1L &&
      identical(main_clustering, "off")) {
    stop(
      "`main_clustering = 'off'` has been removed. ",
      "Use one of: 'similarity', 'proximity', 'enhanced_proximity'."
    )
  }
  main_clustering <- match.arg(main_clustering)
  samptype <- match.arg(samptype)
  filter_mode <- match.arg(filter_mode)
  filter_method <- match.arg(filter_method)
  recon_fusion <- match.arg(recon_fusion)
  if (!is.null(top_v)) {
    if (!is.numeric(top_v) || length(top_v) != 1L || !is.finite(top_v) || top_v <= 0) {
      stop("`top_v` must be NULL or a single positive integer.")
    }
    top_v <- as.integer(top_v)
    model_top_v <- top_v
    fused_top_v <- top_v
  }
  connect_list <- normalize_connect_list(
    connect_list = connect_list,
    n_blocks = length(dat.list),
    valid_names = names(dat.list)
  )
  # model_top_v: NULL = auto-tune, Inf = no truncation, integer = fixed cutoff
  if (isTRUE(model_top_v)) {
    model_top_v <- NULL
  }
  if (!is.null(model_top_v) && !identical(model_top_v, Inf)) {
    if (!is.numeric(model_top_v) || length(model_top_v) != 1L ||
        !is.finite(model_top_v) || model_top_v <= 0) {
      stop("`model_top_v` must be NULL, Inf, or a single positive integer.")
    }
    model_top_v <- as.integer(model_top_v)
  }
  # fused_top_v: NULL = auto-tune, Inf = no truncation, FALSE = skip, integer = fixed
  if (isTRUE(fused_top_v)) {
    fused_top_v <- NULL
  }
  if (!is.null(fused_top_v) && !isFALSE(fused_top_v) && !identical(fused_top_v, Inf)) {
    if (!is.numeric(fused_top_v) || length(fused_top_v) != 1L ||
        !is.finite(fused_top_v) || fused_top_v <= 0) {
      stop("`fused_top_v` must be NULL, Inf, FALSE, or a single positive integer.")
    }
    fused_top_v <- as.integer(fused_top_v)
  }

  # Stage 1: initialization through mrf3_init
  if (verbose) message("Fitting forests..")

  ## When IMD is needed and sub_mrf is active, inject compute_imd = TRUE
  ## into the sub_mrf_args so IMD is pre-computed during fitting.
  need_imd_early <- isTRUE(run_imd) || isTRUE(run_variable_selection) || isTRUE(run_robust_clustering)
  if (need_imd_early && isTRUE(dots$sub_mrf)) {
    sma <- if (!is.null(dots$sub_mrf_args)) dots$sub_mrf_args else list()
    if (is.null(sma$compute_imd)) sma$compute_imd <- TRUE
    dots$sub_mrf_args <- sma
  }

  # Skip proximity computation when clustering doesn't need it
  prox_arg <- if (main_clustering == "similarity") "none" else "all"

  # When enhanced_proximity is requested and the native engine is active,
  # compute enhanced proximity inside C++ during tree building.
  # This avoids the slow R-level cl_forest() foreach loop.
  if (identical(main_clustering, "enhanced_proximity") &&
      identical(getOption("multiRF.engine", "native"), "native")) {
    if (is.null(dots$enhanced_prox)) dots$enhanced_prox <- TRUE
    # Forward sibling_gamma from clustering_args if set there
    if (is.null(dots$sibling_gamma) && !is.null(clustering_args$sibling_gamma)) {
      dots$sibling_gamma <- clustering_args$sibling_gamma
    }
    if (is.null(dots$leaf_embed_dim) && !is.null(clustering_args$leaf_embed_dim)) {
      dots$leaf_embed_dim <- clustering_args$leaf_embed_dim
    }
    # Also inject into shared_specific_args so the residual unsupervised
    # forests (for specific clustering) also get enhanced_prox computed in C++.
    if (is.null(shared_specific_args$enhanced_prox)) {
      shared_specific_args$enhanced_prox <- TRUE
    }
    if (!is.null(dots$sibling_gamma) && is.null(shared_specific_args$sibling_gamma)) {
      shared_specific_args$sibling_gamma <- dots$sibling_gamma
    }
    if (!is.null(dots$leaf_embed_dim) && is.null(shared_specific_args$leaf_embed_dim)) {
      shared_specific_args$leaf_embed_dim <- dots$leaf_embed_dim
    }
  }

  init_args <- c(
    list(
      dat.list = dat.list,
      ntree = ntree,
      scale = scale,
      ytry = ytry,
      samptype = samptype,
      proximity = prox_arg,
      connect_list = connect_list,
      filter_mode = filter_mode,
      filter_method = filter_method,
      top_n_by_type = top_n_by_type,
      top_n_manual = top_n_manual,
      filter_verbose = filter_verbose,
      return_data = TRUE,
      verbose = verbose,
      seed = seed
    ),
    dots
  )
  init_mod <- do.call(mrf3_init, init_args)

  dat_filtered <- init_mod$dat.list
  if (is.null(dat_filtered)) {
    stop("`mrf3_fit()` expects `mrf3_init(return_data = TRUE)` to provide filtered data.")
  }
  if (isTRUE(scale)) {
    dat_fit <- purrr::map(dat_filtered, ~as.data.frame(base::scale(.), check.names = FALSE))
  } else {
    dat_fit <- dat_filtered
  }

  mod_list <- init_mod$mod
  oob_err <- init_mod$oob_err
  type <- init_mod$type

  conn_strategy <- if (length(dat_fit) == 1L) {
    "single_block"
  } else if (!is.null(connect_list)) {
    "provided"
  } else {
    "auto"
  }
  conn_note <- switch(
    conn_strategy,
    single_block = "Single-block mode: connection fixed to itself.",
    provided = "Using user-provided connect_list.",
    auto = "Auto-selected by forest-weight quality ranking."
  )
  stage_connection <- list(
    strategy = conn_strategy,
    ready = TRUE,
    model_connection = names(mod_list),
    connect_list = init_mod$connection,
    score = init_mod$connection_score,
    top_v_used = init_mod$connection_top_v_used,
    note = conn_note
  )
  if (is.null(stage_connection$connect_list)) {
    stage_connection$connect_list <- enumerate_connections(names(dat_fit))
    stage_connection$note <- paste(
      "Connection list missing from init output; using enumerated directed pairs."
    )
  }

  connect_for_downstream <- stage_connection$connect_list

  shared_k_for_tune <- if (!is.null(clustering_args$shared_k)) {
    clustering_args$shared_k
  } else {
    NULL
  }

  disable_fused_top_v <- isFALSE(fused_top_v)

  top_v_method <- match.arg(top_v_method)
  top_v_main <- resolve_top_v_values(
    dat_input = dat_fit,
    mod_input = mod_list,
    connection_input = connect_for_downstream,
    connection_score = stage_connection$score,
    model_top_v_input = model_top_v,
    fused_top_v_input = fused_top_v,
    disable_fused_top_v = disable_fused_top_v,
    shared_k_for_tune = shared_k_for_tune,
    top_v_method = top_v_method,
    neff_quantile = neff_quantile,
    model_top_v_tune_args = model_top_v_tune_args,
    fused_top_v_tune_args = fused_top_v_tune_args,
    stage_prefix = "[Stage 3/5]",
    verbose = verbose
  )
  stage_reconstruction_tuning <- top_v_main$tuning
  final_model_top_v <- top_v_main$model_top_v
  final_fused_top_v <- top_v_main$fused_top_v

  if (verbose) message("Reconstructing and clustering..")
  main_run <- run_branch_pipeline(
    dat_input = dat_fit,
    mod_input = mod_list,
    model_top_v_use = final_model_top_v,
    fused_top_v_use = final_fused_top_v,
    connection_score = stage_connection$score,
    recon_args = list(
      recon_fusion = recon_fusion,
      score_power = score_power,
      score_floor = score_floor,
      fallback_uniform = fallback_uniform,
      fused_row_normalize = fused_row_normalize,
      fused_keep_ties = fused_keep_ties
    ),
    shared_specific_args = shared_specific_args,
    clustering_args = clustering_args,
    main_clustering = main_clustering
  )
  stage_reconstruction <- main_run$reconstruction
  main_branch <- main_run$branch
  main_clusters <- main_branch$clusters
  main_shared <- main_branch$shared
  main_specific <- main_branch$specific

  # Stage 5: IMD / variable selection
  stage_imd <- NULL
  stage_cluster_imd <- NULL
  stage_variable_selection <- NULL
  branch_robust_clustering <- NULL
  do_cluster_imd <- FALSE
  run_variable_selection_exec <- isTRUE(run_variable_selection) || isTRUE(run_robust_clustering)
  need_imd <- isTRUE(run_imd) || isTRUE(run_variable_selection_exec)
  if (need_imd) {
    if (isTRUE(run_imd)) {
      if (verbose) message("Computing IMD..")
    } else {
      if (verbose) message("Computing IMD..")
    }

    ## Classify each model: sub-MRF (no tree structure) vs full forest
    is_sub_mrf <- vapply(mod_list, function(m) !is.null(m$sub_mrf_info), logical(1))
    has_precomputed <- vapply(mod_list, function(m) !is.null(m$imd_weights), logical(1))

    if (any(is_sub_mrf)) {
      ## Mixed or all sub-MRF case
      ## Sub-MRF models MUST have pre-computed IMD; full forests can be computed normally
      sub_without_imd <- is_sub_mrf & !has_precomputed
      if (any(sub_without_imd)) {
        warning(
          "Sub-MRF models detected without pre-computed IMD weights. ",
          "Set `compute_imd = TRUE` in `sub_mrf_args` to enable. Skipping IMD.",
          call. = FALSE
        )
        stage_imd <- NULL
      } else {
        if (verbose) message("  Using pre-computed IMD.")

        ## Compute IMD per-connection
        weight_l <- lapply(names(mod_list), function(m_name) {
          mod <- mod_list[[m_name]]
          if (!is.null(mod$imd_weights)) {
            ## Pre-computed from sub-MRF
            iw <- mod$imd_weights
          } else {
            ## Full forest: compute via standard tree traversal
            imp_out <- get_imp_forest(mod, parallel = TRUE, robust = FALSE,
                                       calc = "Both", normalized = FALSE,
                                       seed = seed)
            iw <- imp_out$imp_ls
          }
          m_name_sep <- rev(unlist(stringr::str_split(m_name, "_")))
          names(iw) <- m_name_sep
          iw
        })

        ## Merge across connections: average per block
        block_names <- names(dat_fit)
        weight_list <- lapply(block_names, function(bn) {
          ww <- purrr::compact(purrr::map(weight_l, bn))
          if (length(ww) == 0L) return(setNames(numeric(ncol(dat_fit[[bn]])), colnames(dat_fit[[bn]])))
          w <- Reduce("+", ww) / length(ww)
          denom <- sqrt(sum(w^2))
          if (is.finite(denom) && denom > 0) w <- w / denom
          w
        })
        names(weight_list) <- block_names
        stage_imd <- list(weight_list = weight_list, weight_list_init = NULL, net = NULL)
      }
    } else {
      ## All full forests: standard path
      default_imd_args <- list(
        mod_list = mod_list,
        dat.list = dat_fit,
        ytry = ytry,
        parallel = TRUE,
        seed = seed
      )
      final_imd_args <- utils::modifyList(default_imd_args, imd_args)
      stage_imd <- do.call(get_multi_weights, final_imd_args)
    }
    ## Resolve run_cluster_imd: default is FALSE unless explicitly requested
    do_cluster_imd <- isTRUE(run_cluster_imd)

    if (isTRUE(run_imd) && do_cluster_imd) {
      cluster_labels <- NULL
      if (!is.null(main_clusters)) {
        cluster_labels <- main_clusters
      }
      if (!is.null(cluster_imd_args$cluster)) {
        cluster_labels <- cluster_imd_args$cluster
      }

      if (is.null(cluster_labels)) {
        warning(
          "`run_imd = TRUE`: skipped `cluster_imd` because no cluster labels were found. ",
          "Pass `cluster` via `cluster_imd_args`, or ensure shared clustering returns labels.",
          call. = FALSE
        )
      } else {
        if (verbose) message("  Cluster-specific IMD..")
        wf_for_cluster_imd <- list(
          models = mod_list,
          connection = connect_for_downstream,
          config = list(
            ntree = ntree,
            ytry = ytry
          ),
          clusters = main_clusters,
          shared = main_shared,
          specific = main_specific,
          data = dat_fit
        )
        class(wf_for_cluster_imd) <- c("mrf3_fit", "list")

        default_cluster_imd_args <- list(
          x = wf_for_cluster_imd,
          cluster = cluster_labels,
          dat.list = dat_fit,
          parallel = TRUE,
          seed = seed
        )
        final_cluster_imd_args <- utils::modifyList(default_cluster_imd_args, cluster_imd_args)
        stage_cluster_imd <- do.call(cluster_imd, final_cluster_imd_args)
      }
    }
  }

  if (isTRUE(run_variable_selection_exec)) {
    if (is.null(stage_imd) || is.null(stage_imd$weight_list)) {
      warning(
        "Variable-selection branch requires valid IMD weights. Skipping variable selection branch.",
        call. = FALSE
      )
    } else {
      if (verbose) message("Variable selection..")
      force_refit <- isTRUE(run_robust_clustering)
      vs_defaults <- list(
        method = "test",
        re_fit = force_refit,
        re_weights = FALSE,
        normalized = TRUE
      )
      final_vs_args <- utils::modifyList(vs_defaults, variable_selection_args)
      final_vs_args$re_fit <- force_refit
      if (identical(as.character(final_vs_args$method)[1], "filter") && !isTRUE(final_vs_args$re_fit)) {
        warning(
          "`mrf3_vs(method = 'filter')` performs iterative refits. ",
          "Switching to `method = 'thres'` for no-refit pipeline mode.",
          call. = FALSE
        )
        final_vs_args$method <- "thres"
      }

      ## `method = 'test'` requires per-tree weight distributions (weight_list_init).
      ## Sub-MRF pre-computed IMD does not produce these, so fall back to 'mixture'.
      if (identical(as.character(final_vs_args$method)[1], "test") &&
          is.null(stage_imd$weight_list_init)) {
        if (verbose) message("  Falling back to method = 'mixture'.")
        final_vs_args$method <- "mixture"
      }

      mod_vs <- list(
        weights = stage_imd$weight_list,
        weights_ls = stage_imd$weight_list_init,
        connection = connect_for_downstream,
        ytry = ytry,
        ntree = ntree,
        type = type,
        oob_err = oob_err,
        mod = mod_list
      )
      class(mod_vs) <- "mrf3"
      vs_fit <- do.call(mrf3_vs, c(list(mod = mod_vs, dat.list = dat_fit), final_vs_args))
      dat_selected <- lapply(vs_fit$dat.list, function(x) {
        if (is.null(x)) {
          return(NULL)
        }
        as.data.frame(x, check.names = FALSE)
      })
      dat_names <- names(dat_fit)
      summary_rows <- lapply(dat_names, function(nm) {
        p_total <- ncol(dat_fit[[nm]])
        p_selected <- if (is.null(dat_selected[[nm]])) 0L else as.integer(ncol(dat_selected[[nm]]))
        data.frame(
          block = nm,
          p_total = p_total,
          p_selected = p_selected,
          selected_ratio = ifelse(p_total > 0, p_selected / p_total, NA_real_),
          stringsAsFactors = FALSE
        )
      })
      stage_variable_selection <- list(
        method = as.character(final_vs_args$method)[1],
        re_fit = isTRUE(final_vs_args$re_fit),
        args = final_vs_args,
        summary = do.call(rbind, summary_rows),
        selected_vars = lapply(dat_selected, colnames),
        dat_selected = dat_selected,
        vs_fit = if (isTRUE(compact_output)) NULL else vs_fit,
        mod_list = if (isTRUE(final_vs_args$re_fit)) vs_fit$mod else mod_list,
        connection = if (!is.null(vs_fit$connection)) vs_fit$connection else connect_for_downstream
      )
    }
  }

  if (isTRUE(run_robust_clustering)) {
    if (is.null(stage_imd) || is.null(stage_imd$weight_list)) {
      warning(
        "`run_robust_clustering = TRUE` requires valid IMD weights. Skipping robust clustering.",
        call. = FALSE
      )
    } else if (is.null(stage_variable_selection$dat_selected)) {
      warning(
        "`run_robust_clustering = TRUE` requires variable-selection output from `mrf3_vs`. Skipping robust clustering.",
        call. = FALSE
      )
    } else {
      if (verbose) message("Robust clustering..")
      robust_dat <- stage_variable_selection$dat_selected

      robust_mod_list <- mod_list
      robust_connection <- connect_for_downstream
      if (isTRUE(stage_variable_selection$re_fit) && is.list(stage_variable_selection$mod_list)) {
        robust_mod_list <- stage_variable_selection$mod_list
        robust_connection <- stage_variable_selection$connection
      }

      top_v_robust <- list(
        model_top_v = final_model_top_v,
        fused_top_v = final_fused_top_v
      )

      ## Inherit main clustering k, then let robust_clustering_args override
      robust_cl_args <- clustering_args
      ## Extract actual k from Stage 4
      main_shared_k <- main_shared$clustering$k
      if (!is.null(main_shared_k) && is.null(robust_cl_args$shared_k)) {
        robust_cl_args$shared_k <- main_shared_k
      }
      main_specific_k <- lapply(main_specific$clustering$by_omics, `[[`, "k")
      main_specific_k <- Filter(Negate(is.null), main_specific_k)
      if (length(main_specific_k) > 0L && is.null(robust_cl_args$specific_k)) {
        robust_cl_args$specific_k <- main_specific_k
      }
      ## User overrides via robust_clustering_args take precedence
      robust_cl_args <- utils::modifyList(robust_cl_args, robust_clustering_args)

      robust_run <- run_branch_pipeline(
        dat_input = robust_dat,
        mod_input = robust_mod_list,
        model_top_v_use = top_v_robust$model_top_v,
        fused_top_v_use = top_v_robust$fused_top_v,
        connection_score = stage_connection$score,
        recon_args = list(
          recon_fusion = recon_fusion,
          score_power = score_power,
          score_floor = score_floor,
          fallback_uniform = fallback_uniform,
          fused_row_normalize = fused_row_normalize,
          fused_keep_ties = fused_keep_ties
        ),
        shared_specific_args = shared_specific_args,
        clustering_args = robust_cl_args,
        main_clustering = main_clustering
      )
      robust_branch <- robust_run$branch
      branch_robust_clustering <- list(
        clusters = robust_branch$clusters,
        shared = robust_branch$shared,
        specific = robust_branch$specific
      )
    }
  }

  if (isTRUE(compact_output)) {
    if (is.list(stage_variable_selection)) {
      stage_variable_selection$mod_list <- NULL
      stage_variable_selection$vs_fit <- NULL
    }
    mod_list_out <- NULL
  } else {
    mod_list_out <- mod_list
  }

  out <- list(
    call = match.call(),

    ## ---------- Config (parameters) ----------
    config = list(
      ntree = ntree,
      scale = scale,
      ytry = ytry,
      filter_mode = filter_mode,
      filter_method = filter_method,
      main_clustering = main_clustering,
      run_imd = run_imd,
      run_variable_selection = run_variable_selection_exec,
      run_robust_clustering = run_robust_clustering,
      run_cluster_imd = do_cluster_imd,
      imd_computed = !is.null(stage_imd),
      compact_output = isTRUE(compact_output),
      recon_fusion = recon_fusion,
      score_power = score_power,
      score_floor = score_floor,
      fallback_uniform = fallback_uniform,
      fused_row_normalize = fused_row_normalize,
      fused_keep_ties = fused_keep_ties,
      seed = seed
    ),

    ## ---------- Init (flat) ----------
    type = type,
    oob_err = oob_err,
    models = mod_list_out,
    filter_summary = init_mod$filter_summary,
    connection = connect_for_downstream,
    connection_score = init_mod$connection_score,

    ## ---------- Tuning ----------
    model_top_v = final_model_top_v,
    fused_top_v = final_fused_top_v,
    tuning_detail = list(
      model_top_v_candidates = stage_reconstruction_tuning$model_top_v,
      fused_top_v_candidates = stage_reconstruction_tuning$fused_top_v
    ),

    ## ---------- Reconstruction ----------
    reconstruction = stage_reconstruction,

    ## ---------- Clustering ----------
    clusters = main_clusters,
    shared = main_shared,
    specific = main_specific,

    ## ---------- IMD ----------
    weights = if (!is.null(stage_imd)) stage_imd$weight_list else NULL,
    weights_init = if (!is.null(stage_imd)) stage_imd$weight_list_init else NULL,
    imd_net = if (!is.null(stage_imd)) stage_imd$net else NULL,

    ## ---------- Cluster IMD ----------
    cluster_imd = NULL,

    ## ---------- Variable Selection ----------
    selected_vars = if (!is.null(stage_variable_selection)) stage_variable_selection$selected_vars else NULL,
    selected_data = if (!is.null(stage_variable_selection)) stage_variable_selection$dat_selected else NULL,
    vs_summary = if (!is.null(stage_variable_selection)) stage_variable_selection$summary else NULL,
    vs_detail = if (!is.null(stage_variable_selection)) {
      list(
        method = stage_variable_selection$method,
        re_fit = stage_variable_selection$re_fit,
        args = stage_variable_selection$args,
        vs_fit = stage_variable_selection$vs_fit,
        mod_list = stage_variable_selection$mod_list,
        connection = stage_variable_selection$connection
      )
    } else {
      NULL
    },

    ## ---------- Robust Clustering ----------
    robust_clusters = if (!is.null(branch_robust_clustering)) branch_robust_clustering$clusters else NULL,
    robust_detail = if (!is.null(branch_robust_clustering)) {
      list(shared = branch_robust_clustering$shared, specific = branch_robust_clustering$specific)
    } else {
      NULL
    },

    ## ---------- Data ----------
    data = if (return_data) dat_fit else NULL
  )
  if (!is.null(stage_cluster_imd)) {
    out$cluster_imd <- stage_cluster_imd
  }
  class(out) <- c("mrf3_fit", "list")
  out
}

  
