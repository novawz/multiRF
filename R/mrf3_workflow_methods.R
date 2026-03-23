#' Summarize `mrf3_fit` output
#'
#' @method summary mrf3_fit
#' @param object Output from `mrf3_fit()`.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A structured summary object with key metadata and branch status.
#' @export
summary.mrf3_fit <- function(object, ...) {
  if (!is.list(object) || is.null(object$reconstruction)) {
    stop("`object` must be an output from `mrf3_fit()`.")
  }

  recon <- object$reconstruction

  fit_mod <- object$models
  if (is.list(fit_mod)) {
    model_names <- names(fit_mod)
    n_models <- length(fit_mod)
  } else {
    model_names <- character(0)
    n_models <- NA_integer_
  }

  sample_n <- NA_integer_
  if (!is.null(recon$sample_names)) {
    sample_n <- length(recon$sample_names)
  } else if (!is.null(recon$W$W_all)) {
    sample_n <- nrow(recon$W$W_all)
  }

  block_names <- character(0)
  if (length(model_names) > 0L) {
    pairs <- strsplit(model_names, "_", fixed = TRUE)
    block_names <- unique(unlist(lapply(pairs, function(x) x[seq_len(min(2L, length(x)))])))
  }
  if (length(block_names) == 0L && !is.null(recon$fused_mat)) {
    block_names <- names(recon$fused_mat)
  }
  if (length(block_names) == 0L && !is.null(object$specific$weights$W)) {
    block_names <- names(object$specific$weights$W)
  }

  n_connections <- if (is.list(object$connection)) length(object$connection) else NA_integer_

  branch_tb <- data.frame(
    branch = c("shared", "specific", "imd", "cluster_imd", "variable_selection", "robust_clustering"),
    ran = c(
      !is.null(object$shared),
      !is.null(object$specific),
      !is.null(object$weights),
      !is.null(object$cluster_imd),
      !is.null(object$selected_vars),
      !is.null(object$robust_clusters)
    ),
    detail = c(
      if (!is.null(object$shared)) {
        "computed shared weights and shared clustering"
      } else {
        "skipped"
      },
      if (!is.null(object$specific)) {
        sp_method <- object$specific$clustering$method
        if (is.null(sp_method)) "specific=NA" else paste0("specific=", sp_method)
      } else {
        "skipped"
      },
      if (!is.null(object$weights)) {
        paste0("net=", length(object$imd_net))
      } else {
        "skipped"
      },
      if (!is.null(object$cluster_imd)) {
        ok_n <- NA_integer_
        if (!is.null(object$cluster_imd$summary$status)) {
          ok_n <- sum(object$cluster_imd$summary$status == "ok")
        }
        vs_tag <- ""
        if (isTRUE(object$cluster_imd$params$run_vs)) {
          n_vs <- sum(vapply(object$cluster_imd$by_cluster, function(cl) !is.null(cl$vs), logical(1)))
          vs_tag <- paste0(", vs=", n_vs, "/", ok_n)
        }
        paste0("clusters=", ifelse(is.na(ok_n), "NA", as.character(ok_n)), vs_tag)
      } else {
        "skipped"
      },
      if (!is.null(object$selected_vars)) {
        sel_tb <- object$vs_summary
        ptxt <- if (is.data.frame(sel_tb) && nrow(sel_tb) > 0L) {
          paste0("p_selected=", paste(sel_tb$p_selected, collapse = "/"))
        } else {
          "p_selected=NA"
        }
        vs_method <- if (!is.null(object$vs_detail$method)) object$vs_detail$method else "NA"
        paste0("method=", vs_method, ", ", ptxt)
      } else {
        "skipped"
      },
      if (!is.null(object$robust_clusters)) {
        k_use <- length(unique(object$robust_clusters))
        paste0("k=", k_use)
      } else {
        "skipped"
      }
    ),
    stringsAsFactors = FALSE
  )

  shared_frac <- NULL
  if (!is.null(object$shared$frac) &&
      is.data.frame(object$shared$frac)) {
    shared_frac <- object$shared$frac
    if ("shared_frac" %in% names(shared_frac)) {
      o <- order(shared_frac$shared_frac, decreasing = TRUE, na.last = TRUE)
      shared_frac <- shared_frac[o, , drop = FALSE]
    }
  }

  out <- list(
    metadata = list(
      n_blocks = length(block_names),
      block_names = block_names,
      n_samples = sample_n,
      n_models = n_models,
      model_names = model_names
    ),
    connection = list(
      n_connections = n_connections
    ),
    selection = list(
      model_top_v = object$model_top_v,
      fused_top_v = object$fused_top_v
    ),
    branches = branch_tb,
    shared_frac = shared_frac
  )
  class(out) <- "summary.mrf3_fit"
  out
}


#' Print mrf3 summary
#'
#' @method print summary.mrf3_fit
#' @param x Summary object from `summary.mrf3_fit()`.
#' @param digits Digits used for numeric formatting.
#' @param max_shared_frac_rows Maximum number of shared-fraction rows to print.
#' @param ... Additional arguments (unused).
#'
#' @return The input object, invisibly.
#' @export
print.summary.mrf3_fit <- function(x, digits = 3, max_shared_frac_rows = 8L, ...) {
  meta <- x$metadata
  conn <- x$connection
  sel <- x$selection
  branch_tb <- x$branches

  cat("mrf3_fit Summary\n")
  cat(strrep("-", 16), "\n")
  cat("  Blocks      : ", meta$n_blocks, "\n", sep = "")
  if (length(meta$block_names) > 0L) {
    cat("  Block names : ", paste(meta$block_names, collapse = ", "), "\n", sep = "")
  }
  cat("  Samples     : ", ifelse(is.na(meta$n_samples), "NA", as.character(meta$n_samples)), "\n", sep = "")
  cat("  Models      : ", ifelse(is.na(meta$n_models), "NA", as.character(meta$n_models)), "\n", sep = "")
  cat("  Connections : ", ifelse(is.na(conn$n_connections), "NA", as.character(conn$n_connections)), "\n", sep = "")

  cat("\n  Top-v: model=", sel$model_top_v,
      ", fused=", ifelse(is.null(sel$fused_top_v), "NULL", as.character(sel$fused_top_v)),
      "\n", sep = "")

  cat("\nBranches\n")
  for (i in seq_len(nrow(branch_tb))) {
    cat(
      "  - ", branch_tb$branch[i], ": ",
      ifelse(branch_tb$ran[i], "ran", "skipped"),
      " (", branch_tb$detail[i], ")\n",
      sep = ""
    )
  }

  if (!is.null(x$shared_frac) && nrow(x$shared_frac) > 0L) {
    cat("\nShared Fraction (top blocks)\n")
    tb <- x$shared_frac
    keep <- seq_len(min(nrow(tb), as.integer(max_shared_frac_rows)))
    show_cols <- intersect(c("data", "shared_frac", "specific_ratio", "method"), names(tb))
    print(utils::head(tb[keep, show_cols, drop = FALSE], n = length(keep)), digits = digits, row.names = FALSE)
  }

  invisible(x)
}


#' Print `mrf3_fit` object in compact form
#'
#' @method print mrf3_fit
#' @param x Output from `mrf3_fit()`.
#' @param ... Additional arguments passed to `print.summary.mrf3_fit()`.
#'
#' @return The input object, invisibly.
#' @export
print.mrf3_fit <- function(x, ...) {
  sx <- summary.mrf3_fit(x)
  print.summary.mrf3_fit(sx, ...)
  cat("\nKey fields: $clusters, $weights, $selected_vars, $robust_clusters,\n")
  cat("  $models, $connection, $shared, $specific, $cluster_imd, $reconstruction, $data\n")
  cat("Accessors: get_clusters(), get_weights(), get_selected_vars(), get_top_vars(), ...\n")
  invisible(x)
}


# ---------- Backward compatibility -------------------------------------------

#' @export
`$.mrf3_fit` <- function(x, name) {
  # Map old nested paths to new flat fields
  compat <- list(
    init = function() {
      .Deprecated(msg = paste0(
        "Access via `$init` is deprecated.\n",
        "Use `$type`, `$oob_err`, `$models`, `$connection`, `$filter_summary` directly."
      ))
      list(
        type = .subset2(x, "type"),
        oob_err = .subset2(x, "oob_err"),
        mod_list = .subset2(x, "models"),
        connection = list(
          connect_list = .subset2(x, "connection"),
          score = .subset2(x, "connection_score")
        ),
        filter_summary = .subset2(x, "filter_summary")
      )
    },
    tuning = function() {
      .Deprecated(msg = paste0(
        "Access via `$tuning` is deprecated.\n",
        "Use `$model_top_v`, `$fused_top_v`, `$tuning_detail` directly."
      ))
      list(
        selected_model_top_v = .subset2(x, "model_top_v"),
        selected_fused_top_v = .subset2(x, "fused_top_v"),
        model_top_v_candidates = .subset2(x, "tuning_detail")$model_top_v_candidates,
        fused_top_v_candidates = .subset2(x, "tuning_detail")$fused_top_v_candidates
      )
    },
    imd = function() {
      .Deprecated(msg = paste0(
        "Access via `$imd` is deprecated.\n",
        "Use `$weights`, `$weights_init`, `$imd_net` directly."
      ))
      list(
        weight_list = .subset2(x, "weights"),
        weight_list_init = .subset2(x, "weights_init"),
        net = .subset2(x, "imd_net")
      )
    },
    variable_selection = function() {
      .Deprecated(msg = paste0(
        "Access via `$variable_selection` is deprecated.\n",
        "Use `$selected_vars`, `$selected_data`, `$vs_summary`, `$vs_detail` directly."
      ))
      list(
        selected_vars = .subset2(x, "selected_vars"),
        dat_selected = .subset2(x, "selected_data"),
        summary = .subset2(x, "vs_summary"),
        method = .subset2(x, "vs_detail")$method,
        vs_fit = .subset2(x, "vs_detail")$vs_fit,
        mod_list = .subset2(x, "vs_detail")$mod_list
      )
    },
    robust_clustering = function() {
      .Deprecated(msg = paste0(
        "Access via `$robust_clustering` is deprecated.\n",
        "Use `$robust_clusters`, `$robust_detail` directly."
      ))
      detail <- .subset2(x, "robust_detail")
      c(list(clusters = .subset2(x, "robust_clusters")), if (is.list(detail)) detail else list())
    }
  )

  if (name %in% names(compat)) {
    return(compat[[name]]())
  }

  .subset2(x, name)
}
