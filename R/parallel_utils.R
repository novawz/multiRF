sanitize_mc_cores <- function(cores = NULL, fallback = 1L) {
  detect <- suppressWarnings(parallel::detectCores())
  if (length(detect) != 1L || is.na(detect) || !is.finite(detect) || detect < 1L) {
    detect <- fallback
  }
  if (length(cores) != 1L || is.na(cores) || !is.finite(cores) || cores < 1L) {
    cores <- fallback
  }
  as.integer(max(1L, min(as.integer(cores), as.integer(detect))))
}

initialize_psock_workers <- function(cluster) {
  lib_paths <- .libPaths()
  option_names <- intersect(
    c("multiRF.engine", "multiRF.nthread"),
    names(options())
  )
  package_options <- options()[option_names]

  parallel::clusterCall(
    cluster,
    function(paths, opts) {
      .libPaths(paths)
      if (length(opts)) {
        options(opts)
      }
      invisible(NULL)
    },
    lib_paths,
    package_options
  )
  invisible(cluster)
}

# Apply a function in fresh R processes.  PSOCK is used explicitly so this
# helper is safe in front ends, such as Positron, that do not allow the active
# R session to be forked.  parLapply() preserves the order of X.
parallel_lapply_psock <- function(X, FUN, cores = NULL, ...) {
  FUN <- match.fun(FUN)
  n_tasks <- length(X)
  if (n_tasks == 0L) {
    return(lapply(X, FUN, ...))
  }

  if (is.null(cores)) {
    cores <- max(1L, parallel::detectCores() - 1L)
  }
  cores <- sanitize_mc_cores(cores = cores, fallback = 1L)
  cores <- min(cores, n_tasks)
  if (cores <= 1L) {
    return(lapply(X, FUN, ...))
  }

  cluster <- parallel::makePSOCKcluster(cores)
  on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
  initialize_psock_workers(cluster)
  parallel::parLapply(cluster, X, FUN, ...)
}
