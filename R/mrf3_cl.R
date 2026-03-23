pick_tuned_value <- function(tb, column, allow_infinite = FALSE) {
  if (!is.data.frame(tb) || nrow(tb) == 0L || !column %in% names(tb)) {
    stop("Tuning did not return any valid candidate.")
  }
  value <- tb[[column]][1]
  if (!is.finite(value)) {
    if (isTRUE(allow_infinite)) {
      return(NULL)
    }
    stop("Selected tuning value in `", column, "` is not finite.")
  }
  as.integer(value)
}

#' Run Shared/Specific Clustering Pipeline
#'
#' Main clustering entrypoint for the current `mrf3_fit` workflow.
#'
#' @param recon Reconstruction object from `get_reconstr_matrix()`.
#' @param shared_specific Shared/specific weight object from
#'   `get_shared_specific_weights()`.
#' @param mod_for_shared Random-forest model list used by proximity-based shared
#'   clustering.
#' @param clustering_args Named list of clustering arguments (for example
#'   `shared_k`, `specific_k`, `specific_prox_method_cl`, `tune_method`,
#'   `gap_w`).
#' @param cluster_method One of `"similarity"`, `"proximity"`, or
#'   `"enhanced_proximity"`.
#'
#' @return A list with `shared` and `specific` clustering outputs.
#' @export
run_cluster_pipeline <- function(recon,
                                 shared_specific,
                                 mod_for_shared,
                                 clustering_args,
                                 cluster_method) {
  cluster_method <- match.arg(
    as.character(cluster_method)[1],
    c("similarity", "proximity", "enhanced_proximity")
  )

  args <- utils::modifyList(
    list(recon = recon, shared_specific = shared_specific),
    clustering_args
  )

  # Remove forest-fitting params that should never reach clustering functions
  args$mtry <- NULL
  args$ytry <- NULL

  if (identical(cluster_method, "similarity")) {
    if (is.null(args$shared_method)) args$shared_method <- "PAM"
    if (is.null(args$specific_method)) args$specific_method <- "PAM"
    return(do.call(mrf3_specific_clustering, args))
  }

  enhanced <- identical(cluster_method, "enhanced_proximity")
  prox_label <- if (enhanced) "Enhanced_Proximity" else "Proximity"

  shared_k <- args$shared_k
  specific_k <- args$specific_k
  prox_method_cl <- if (!is.null(args$specific_prox_method_cl)) {
    as.character(args$specific_prox_method_cl)[1]
  } else {
    "PAM"
  }

  prox_passthrough <- args
  drop_names <- c(
    "recon", "shared_specific",
    "k", "method", "enhanced", "method_cl", "prox_method_cl", "rfit",
    "shared_mode", "shared_dat_use", "shared_k", "shared_method", "shared_similarity_type",
    "specific_k", "specific_method", "specific_prox_method_cl", "specific_similarity_type",
    "tune_method", "gap_w", "gamma", "alpha_init", "ao_max_iter", "ao_tol",
    "knn_q", "hollow", "ao_symm", "ao_verbose",
    "mtry", "ytry"
  )
  prox_passthrough[intersect(names(prox_passthrough), drop_names)] <- NULL

  shared_prox <- do.call(
    mrf3_cl_prox,
    c(
      list(
        rfit = mod_for_shared,
        k = shared_k,
        enhanced = enhanced,
        method_cl = prox_method_cl
      ),
      prox_passthrough
    )
  )
  shared_out <- list(
    mode = "proximity",
    dat_used = "mod_list",
    similarity = shared_prox$dat,
    similarity_type = if (enhanced) "enhanced_proximity" else "proximity",
    cl = shared_prox$cl,
    cl_mod = shared_prox$cl_mod,
    k = shared_prox$k,
    method = prox_label,
    embed = NULL,
    ao_fit = NULL
  )

  specific_out <- do.call(
    cluster_specific_similarity,
    c(
      list(
        shared_specific = shared_specific,
        k = specific_k,
        method = prox_label,
        prox_method_cl = prox_method_cl
      ),
      prox_passthrough
    )
  )

  list(
    shared = shared_out,
    specific = specific_out
  )
}

build_shared_specific <- function(dat_input,
                                  recon_input,
                                  specific_top_v,
                                  shared_specific_args = list()) {
  default_shared_specific_args <- list(
    dat.list = dat_input,
    recon = recon_input,
    specific_top_v = specific_top_v
  )
  final_shared_specific_args <- utils::modifyList(default_shared_specific_args, shared_specific_args)
  out <- do.call(get_shared_specific_weights, final_shared_specific_args)
  out$shared_frac <- get_shared_frac(
    dat.list = dat_input,
    shared_specific = out
  )
  out
}

merge_cluster_outputs <- function(shared_specific, specific_clustering) {
  shared_out <- list(
    weights = shared_specific$shared,
    frac = shared_specific$shared_frac,
    clustering = specific_clustering$shared
  )
  clusters <- shared_out$clustering$cl
  shared_out$clustering$cl <- NULL
  specific_out <- list(
    weights = shared_specific$specific,
    clustering = specific_clustering$specific
  )
  list(
    clusters = clusters,
    shared = shared_out,
    specific = specific_out
  )
}

run_branch_pipeline <- function(dat_input,
                                mod_input,
                                model_top_v_use,
                                fused_top_v_use,
                                connection_score,
                                recon_args = list(),
                                shared_specific_args = list(),
                                clustering_args = list(),
                                main_clustering = c("similarity", "proximity", "enhanced_proximity")) {
  main_clustering <- match.arg(main_clustering)
  recon_defaults <- list(
    rfit = mod_input,
    model_top_v = model_top_v_use,
    connection_score = connection_score
  )
  final_recon_args <- utils::modifyList(recon_defaults, recon_args)
  final_recon_args$fused_top_v <- fused_top_v_use

  recon <- do.call(get_reconstr_matrix, final_recon_args)
  shared_specific <- build_shared_specific(
    dat_input = dat_input,
    recon_input = recon,
    specific_top_v = fused_top_v_use,
    shared_specific_args = shared_specific_args
  )
  specific_clustering <- run_cluster_pipeline(
    recon = recon,
    shared_specific = shared_specific,
    mod_for_shared = mod_input,
    clustering_args = clustering_args,
    cluster_method = main_clustering
  )
  list(
    reconstruction = recon,
    branch = merge_cluster_outputs(
      shared_specific = shared_specific,
      specific_clustering = specific_clustering
    )
  )
}
