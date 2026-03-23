#' Filter multi-omics feature blocks before RF fitting
#'
#' @param dat.list A named list containing multi-omics datasets with samples in
#' rows and features in columns.
#' @param filter_mode Feature filtering mode. `"auto"` applies omics-aware
#' defaults, `"none"` keeps all features, and `"manual"` uses `top_n_manual`.
#' @param filter_method Feature dispersion metric used for ranking. Choose from
#' `"mad"` or `"variance"`.
#' @param top_n_by_type Optional named numeric/list object used in `"auto"` mode.
#' Supported names are `rna`, `mirna`, `methylation`, `protein`, and `unknown`.
#' Defaults are RNA = 5000, miRNA = adaptive 500-1000, methylation = adaptive
#' 5000-20000, protein = all, unknown = all.
#' @param top_n_manual Optional named or positional numeric/list object used in
#' `"manual"` mode. Names should match `dat.list` block names.
#' @param return_summary Logical; whether to return filtering summary metadata.
#' @param verbose Logical; whether to print per-block filtering messages.
#' @param ... Additional arguments (currently ignored).
#'
#' @return If `return_summary = FALSE`, returns filtered `dat.list`.
#' Otherwise returns a list with `dat_filtered` and `filter_summary`.
#' @export filter_omics
filter_omics <- function(dat.list,
                         filter_mode = c("auto", "none", "manual"),
                         filter_method = c("mad", "variance"),
                         top_n_by_type = NULL,
                         top_n_manual = NULL,
                         return_summary = TRUE,
                         verbose = TRUE,
                         ...){

  if (!is.list(dat.list) || length(dat.list) == 0) {
    stop("`dat.list` must be a non-empty list of omics blocks.")
  }

  if (is.null(names(dat.list)) || any(names(dat.list) == "")) {
    names(dat.list) <- paste0("omics", seq_along(dat.list))
  }

  filter_mode <- match.arg(filter_mode)
  filter_method <- match.arg(filter_method)
  top_n_map <- normalize_top_n_by_type(top_n_by_type)

  dat_filtered <- vector("list", length(dat.list))
  names(dat_filtered) <- names(dat.list)
  filter_rows <- vector("list", length(dat.list))

  for (i in seq_along(dat.list)) {
    block_name <- names(dat.list)[i]
    dat_block <- sanitize_omics_block(dat.list[[i]], block_name = block_name)

    omics_type <- infer_omics_type(block_name)
    n_samples <- nrow(dat_block)
    n_features <- ncol(dat_block)
    top_n <- resolve_top_n(
      block_name = block_name,
      block_idx = i,
      omics_type = omics_type,
      n_samples = n_samples,
      n_features = n_features,
      filter_mode = filter_mode,
      top_n_map = top_n_map,
      top_n_manual = top_n_manual
    )

    filter_out <- filter_omics_block(
      dat_block = dat_block,
      top_n = top_n,
      filter_method = filter_method
    )

    dat_filtered[[block_name]] <- filter_out$dat

    if (verbose) {
      message(
        sprintf(
          "Filtering %s [%s]: %d -> %d features (%s)",
          block_name, omics_type, n_features, ncol(filter_out$dat),
          ifelse(filter_out$applied, filter_method, "none")
        )
      )
    }

    filter_rows[[i]] <- data.frame(
      block = block_name,
      omics_type = omics_type,
      n_samples = n_samples,
      n_features_in = n_features,
      n_features_out = ncol(filter_out$dat),
      top_n = top_n,
      filter_applied = filter_out$applied,
      filter_method = ifelse(filter_out$applied, filter_method, "none"),
      stringsAsFactors = FALSE
    )
  }

  if (!return_summary) {
    return(dat_filtered)
  }

  filter_summary <- do.call(rbind, filter_rows)
  list(dat_filtered = dat_filtered, filter_summary = filter_summary)
}


sanitize_omics_block <- function(dat_block, block_name) {
  dat_block <- as.data.frame(dat_block, check.names = FALSE)
  if (ncol(dat_block) == 0) {
    stop("Block `", block_name, "` has zero features.")
  }

  is_num <- vapply(dat_block, is.numeric, logical(1))
  if (!all(is_num)) {
    stop("All columns in block `", block_name, "` must be numeric.")
  }

  dat_block
}


infer_omics_type <- function(block_name) {
  n <- tolower(block_name)
  if (grepl("mirna|micro.?rna", n)) {
    return("mirna")
  }
  if (grepl("meth|methyl|dnam|cpg|450k|850k", n)) {
    return("methylation")
  }
  if (grepl("protein|rppa|prot\\b", n)) {
    return("protein")
  }
  if (grepl("rna|mrna|gene|expr|transcript", n)) {
    return("rna")
  }
  "unknown"
}


default_mirna_top_n <- function(n_samples) {
  if (n_samples < 150) 500L else 1000L
}


default_methylation_top_n <- function(n_samples) {
  if (n_samples <= 100) {
    return(5000L)
  }
  if (n_samples <= 250) {
    return(10000L)
  }
  if (n_samples <= 500) {
    return(15000L)
  }
  20000L
}


normalize_top_n_by_type <- function(top_n_by_type) {
  top_n_map <- c(
    rna = 5000,
    mirna = NA_real_,
    methylation = NA_real_,
    protein = Inf,
    unknown = Inf
  )

  if (is.null(top_n_by_type)) {
    return(top_n_map)
  }

  top_n_vec <- unlist(top_n_by_type, use.names = TRUE)
  if (is.null(names(top_n_vec))) {
    stop("`top_n_by_type` must be named.")
  }

  nm <- tolower(names(top_n_vec))
  alias <- c(
    mrna = "rna",
    expression = "rna",
    transcript = "rna",
    gene = "rna",
    micro_rna = "mirna",
    mi_rna = "mirna",
    methyl = "methylation",
    methy = "methylation",
    dnam = "methylation",
    prot = "protein",
    rppa = "protein"
  )
  idx_alias <- nm %in% names(alias)
  nm[idx_alias] <- alias[nm[idx_alias]]

  valid <- nm %in% names(top_n_map)
  if (any(!valid)) {
    warning(
      "Ignoring unknown `top_n_by_type` names: ",
      paste(unique(names(top_n_vec)[!valid]), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(valid)) {
    top_n_map[nm[valid]] <- as.numeric(top_n_vec[valid])
  }

  top_n_map
}


resolve_manual_top_n <- function(top_n_manual, block_name, block_idx) {
  if (is.null(top_n_manual)) {
    return(Inf)
  }

  manual_vec <- unlist(top_n_manual, use.names = TRUE)
  if (is.null(names(manual_vec))) {
    if (length(manual_vec) >= block_idx) {
      return(as.numeric(manual_vec[[block_idx]]))
    }
    return(Inf)
  }

  if (block_name %in% names(manual_vec)) {
    return(as.numeric(manual_vec[[block_name]]))
  }
  idx_name <- as.character(block_idx)
  if (idx_name %in% names(manual_vec)) {
    return(as.numeric(manual_vec[[idx_name]]))
  }
  Inf
}


resolve_top_n <- function(block_name, block_idx, omics_type, n_samples, n_features,
                          filter_mode, top_n_map, top_n_manual) {
  if (filter_mode == "none") {
    return(n_features)
  }

  if (filter_mode == "manual") {
    top_n <- resolve_manual_top_n(top_n_manual, block_name, block_idx)
  } else {
    top_n <- top_n_map[[omics_type]]
    if (omics_type == "mirna" && (is.na(top_n) || !is.finite(top_n))) {
      top_n <- default_mirna_top_n(n_samples)
    }
    if (omics_type == "methylation" && (is.na(top_n) || !is.finite(top_n))) {
      top_n <- default_methylation_top_n(n_samples)
    }
  }

  if (is.null(top_n) || is.na(top_n) || !is.finite(top_n) || top_n <= 0) {
    return(n_features)
  }

  as.integer(min(n_features, max(1, floor(top_n))))
}


filter_omics_block <- function(dat_block, top_n, filter_method) {
  n_features <- ncol(dat_block)
  if (top_n >= n_features) {
    return(list(dat = dat_block, applied = FALSE))
  }

  m <- as.matrix(dat_block)
  if (filter_method == "mad") {
    scores <- matrixStats::colMads(m, na.rm = TRUE, constant = 1)
  } else {
    scores <- matrixStats::colVars(m, na.rm = TRUE)
  }
  scores[!is.finite(scores)] <- -Inf
  keep <- order(scores, decreasing = TRUE)[seq_len(top_n)]

  list(dat = dat_block[, keep, drop = FALSE], applied = TRUE)
}


enumerate_connections <- function(dat_names) {
  if (length(dat_names) <= 1) {
    return(list(c(dat_names)))
  }

  connect_list <- list()
  for (response_name in dat_names) {
    predictor_names <- dat_names[dat_names != response_name]
    for (predictor_name in predictor_names) {
      connect_list[[length(connect_list) + 1L]] <- c(response_name, predictor_name)
    }
  }

  connect_list
}
