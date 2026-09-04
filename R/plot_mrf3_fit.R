# ============================================================================
# Unified S3 plot entry point for `mrf3_fit` objects
#
# `plot(fit, type = "...")` extracts the right components from the fit and
# forwards them (plus `...`) to the existing exported plot functions, applying
# one consistent cluster palette across all plot types.
# ============================================================================


#' Fixed cluster palette for `mrf3_fit` plots (internal)
#'
#' Returns a fixed, colour-blind-safe colour vector used as the default
#' cluster colouring by `plot.mrf3_fit()`. The first eight colours are
#' identical to the internal palette used by the individual plot functions
#' (`plot_network()`, `plot_circos()`, `plot_km()`, ...), so the same cluster
#' id receives the same colour in every plot type. Requests beyond the palette
#' length are recycled with a warning.
#'
#' @param n Number of colours requested. `NULL` returns the full palette.
#'
#' @return A character vector of hex colours of length `n`.
#' @noRd
mrf3_palette <- function(n = NULL) {
  base_palette <- c(
    "#3182BD", # signal blue
    "#33B5A5", # signal teal
    "#E28E2C", # accent orange
    "#D24B40", # accent red
    "#756BB1", # violet
    "#56B4E9", # sky blue
    "#CC79A7", # rose
    "#8C8C42"  # olive
  )
  if (is.null(n)) {
    return(base_palette)
  }
  n <- as.integer(n)[1]
  if (!is.finite(n) || n < 1L) {
    return(character(0))
  }
  if (n <= length(base_palette)) {
    return(base_palette[seq_len(n)])
  }
  warning(
    "`mrf3_palette()` provides ", length(base_palette),
    " distinct colours but ", n, " were requested; colours are recycled.",
    call. = FALSE
  )
  rep_len(base_palette, n)
}


# Default sample-level group labels for a plot type: shared clusters, robust
# clusters when the robust similarity is being plotted, or block-specific
# clusters when a block-specific source is requested.
.mrf3_plot_default_group <- function(x, source = NULL, omics = NULL) {
  if (identical(source, "robust_clustering")) {
    return(get_clusters(x, which = "robust"))
  }
  if (identical(source, "specific_specific") && !is.null(omics)) {
    cl <- get_clusters(x, block = omics)
    if (!is.null(cl)) {
      return(cl)
    }
  }
  get_clusters(x)
}


# Append the fixed cluster palette to a ggplot embedding as a named manual
# scale so that a cluster id always maps to the same colour, whichever plot
# type displays it. The appended scale replaces the (identical for <= 8
# groups) scale the plot function set; users can still add their own
# `scale_fill_manual()` / `scale_colour_manual()` afterwards, which replaces
# this one.
.mrf3_apply_cluster_scale <- function(p, group) {
  if (is.null(group) || !inherits(p, "ggplot")) {
    return(p)
  }
  gf <- droplevels(factor(group))
  if (nlevels(gf) < 1L) {
    return(p)
  }
  values <- stats::setNames(mrf3_palette(nlevels(gf)), levels(gf))
  mapped <- names(p$mapping)
  scale_use <- if ("fill" %in% mapped) {
    ggplot2::scale_fill_manual(values = values, drop = FALSE)
  } else if (any(c("colour", "color") %in% mapped)) {
    ggplot2::scale_colour_manual(values = values, drop = FALSE)
  } else {
    NULL
  }
  if (is.null(scale_use)) {
    return(p)
  }
  suppressMessages(p + scale_use)
}


#' Plot an `mrf3_fit` object
#'
#' Unified plotting entry point for fitted multiRF pipelines. The method
#' extracts the components each display needs from the fit (the learned
#' sample similarity, cluster labels via [get_clusters()], IMD weights via
#' [get_weights()], pairwise feature adjacency via [pairwise_imd()]) and
#' forwards them, together with `...`, to the corresponding exported plot
#' function, so every option of the underlying function stays reachable.
#'
#' Supported `type` values and their targets:
#'
#' * `"tsne"` (default): [plot_tsne()] on the fitted similarity matrix
#'   (t-SNE), coloured by shared cluster.
#' * `"umap"`: [plot_umap()] on the fitted similarity matrix.
#' * `"embed"`: deprecated pairs display via [plot_embed()], drawn on the top
#'   principal components of the fitted similarity matrix (control the number
#'   of components with `embed_dims`, default 3).
#' * `"network"`: [plot_network()] sample network from the fitted similarity
#'   matrix, vertices coloured by shared cluster.
#' * `"circos"`: [pairwise_imd()] on the fit, then [plot_circos()] of the
#'   feature-level adjacency (`source = "adj_var"` by default). Requires IMD
#'   (`run_imd = TRUE`) and retained models (`compact_output = FALSE`).
#' * `"composition"`: [plot_cluster_composition()] of shared clusters against
#'   a user-supplied `annotation` vector.
#' * `"km"`: [plot_km()] Kaplan--Meier curves of the shared clusters; the fit
#'   holds no survival data, so `time_var`, `event_var`, and `pheno_mat` must
#'   be supplied.
#' * `"weights"`: [plot_weights()] of the per-block IMD weights. Requires
#'   `run_imd = TRUE` (or cluster-level IMD via `run_cluster_imd = TRUE`).
#'
#' For every type the default sample grouping is `get_clusters(x)`; pass
#' `group` through `...` to override it (for `"tsne"`, `"umap"`, `"embed"`,
#' `"network"`, `"composition"`, and `"km"`). When `source =
#' "robust_clustering"` is passed through `...`, robust cluster labels are
#' used instead, and `source = "specific_specific"` with `omics` uses that
#' block's specific cluster labels.
#'
#' Cluster colours are consistent across all types: cluster ids are mapped in
#' sorted order onto a fixed colour-blind-safe palette (the same palette the
#' individual plot functions use), so cluster 2 has the same colour in the
#' t-SNE, network, and Kaplan--Meier displays. A user-supplied colour setting
#' passed through `...` (for example `palette` for `"km"` or `colours` for
#' `"composition"`) always wins over this default, and the returned ggplot
#' objects accept a replacement `scale_fill_manual()` /
#' `scale_colour_manual()`.
#'
#' @param x An `mrf3_fit` object returned by [mrf3_fit()] (or [mrf3()]).
#' @param type Plot type; one of `"tsne"`, `"umap"`, `"embed"`, `"network"`,
#'   `"circos"`, `"composition"`, `"km"`, `"weights"`. Defaults to `"tsne"`.
#' @param ... Additional arguments forwarded to the target plot function
#'   (e.g. `seed`, `perplexity`, `cutoff`, `top`, `source`, `omics`,
#'   `group`, `cut.off`, `risk.table`). For `type = "circos"`,
#'   `feature_source` and `normalized` are forwarded to [pairwise_imd()].
#'   For `type = "embed"`, `embed_dims` sets the number of principal
#'   components displayed.
#' @param annotation Required for `type = "composition"`: a vector or factor
#'   of external sample annotations (e.g. tumour subtypes), one value per
#'   sample and in the same sample order as the fitted data. The fit itself
#'   stores no annotations.
#' @param time_var Required for `type = "km"`: name of the survival/follow-up
#'   time column in `pheno_mat`.
#' @param event_var Required for `type = "km"`: name of the binary 0/1 event
#'   column in `pheno_mat`.
#' @param pheno_mat Required for `type = "km"`: data frame of sample-level
#'   phenotype data containing `time_var` and `event_var`, with one row per
#'   sample in the same sample order as the fitted data.
#'
#' @return Whatever the target plot function returns: a `ggplot` object for
#'   `"tsne"`, `"umap"`, `"composition"`, and `"weights"`; a `ggsurvplot`
#'   object for `"km"`; invisibly, an `igraph` object for `"network"`, the
#'   chord-diagram data for `"circos"`, and `NULL` for `"embed"`.
#'
#' @seealso [plot_tsne()], [plot_umap()], [plot_network()], [plot_circos()],
#'   [plot_cluster_composition()], [plot_km()], [plot_weights()],
#'   [get_clusters()], [pairwise_imd()]
#'
#' @examples
#' \donttest{
#' data("tcga_brca_data", package = "multiRF")
#' ids <- rownames(tcga_brca[[1]])[seq_len(80)]
#' dat <- lapply(tcga_brca, function(x) {
#'   x[ids, seq_len(min(30L, ncol(x))), drop = FALSE]
#' })
#'
#' fit <- mrf3_fit(
#'   dat,
#'   ntree = 30,
#'   filter_mode = "none",
#'   clustering_args = list(shared_k = 3),
#'   run_imd = TRUE,
#'   return_data = TRUE,
#'   filter_verbose = FALSE,
#'   verbose = FALSE,
#'   seed = 529
#' )
#'
#' ## Embeddings and networks of the learned similarity
#' plot(fit, type = "tsne", seed = 529)
#' plot(fit, type = "network", seed = 529)
#'
#' ## IMD feature weights
#' plot(fit, type = "weights", top = 8)
#'
#' ## Cluster composition against an external annotation
#' annotation <- factor(rep(c("Subtype A", "Subtype B"), length.out = 80))
#' plot(fit, type = "composition", annotation = annotation)
#'
#' ## Kaplan-Meier curves need user-supplied survival data
#' if (requireNamespace("survival", quietly = TRUE) &&
#'     requireNamespace("survminer", quietly = TRUE)) {
#'   pheno <- data.frame(
#'     os_time = stats::rexp(80, rate = 1 / 50),
#'     os_event = stats::rbinom(80, 1, 0.6)
#'   )
#'   plot(fit, type = "km", time_var = "os_time", event_var = "os_event",
#'        pheno_mat = pheno, risk.table = FALSE)
#' }
#' }
#'
#' @method plot mrf3_fit
#' @export
plot.mrf3_fit <- function(x,
                          type = c("tsne", "umap", "embed", "network",
                                   "circos", "composition", "km", "weights"),
                          ...,
                          annotation = NULL,
                          time_var = NULL,
                          event_var = NULL,
                          pheno_mat = NULL) {
  type <- match.arg(type)
  dots <- list(...)

  ## Types that colour/group samples by cluster labels
  grouped_types <- c("tsne", "umap", "embed", "network", "composition", "km")
  group_supplied <- "group" %in% names(dots)
  group <- NULL
  if (type %in% grouped_types) {
    if (group_supplied) {
      group <- dots$group
    } else {
      group <- .mrf3_plot_default_group(
        x, source = dots$source, omics = dots$omics
      )
    }
    dots$group <- NULL
  }

  switch(
    type,

    tsne = {
      p <- do.call(plot_tsne, c(list(mod = x, group = group), dots))
      .mrf3_apply_cluster_scale(p, group)
    },

    umap = {
      p <- do.call(plot_umap, c(list(dat = x, group = group), dots))
      .mrf3_apply_cluster_scale(p, group)
    },

    embed = {
      embed_dims <- dots$embed_dims
      if (is.null(embed_dims)) embed_dims <- 3L
      if (length(embed_dims) != 1L || !is.numeric(embed_dims) ||
          !is.finite(embed_dims) || embed_dims < 2L) {
        stop("`embed_dims` must be a single integer >= 2.", call. = FALSE)
      }
      embed_dims <- as.integer(embed_dims)
      mat <- extract_mrf3_plot_matrix(
        x,
        source = if (is.null(dots$source)) "auto" else dots$source,
        omics = dots$omics,
        cluster = dots$cluster
      )
      scores <- stats::prcomp(mat)$x
      scores <- scores[, seq_len(min(embed_dims, ncol(scores))), drop = FALSE]
      colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
      dots <- dots[setdiff(
        names(dots), c("embed_dims", "source", "omics", "cluster")
      )]
      invisible(do.call(plot_embed, c(list(dat = scores, group = group), dots)))
    },

    network = {
      res <- do.call(plot_network, c(list(dat = x, group = group), dots))
      invisible(res)
    },

    circos = {
      if (is.null(x$imd)) {
        stop(
          "`type = \"circos\"` requires IMD weights, but this fit has none. ",
          "Refit with `run_imd = TRUE` (or `run_variable_selection = TRUE`).",
          call. = FALSE
        )
      }
      has_pairwise <- is.list(x$models) &&
        any(vapply(x$models, function(m) !is.null(m$pairwise_xy), logical(1)))
      if (!has_pairwise && is.null(x$imd_net)) {
        stop(
          "`type = \"circos\"` requires pairwise IMD input: fitted models ",
          "with pre-computed `pairwise_xy` or an IMD network (`imd_net`). ",
          "Refit with `run_imd = TRUE` and `compact_output = FALSE`.",
          call. = FALSE
        )
      }
      pairwise_args <- dots[intersect(
        names(dots), c("feature_source", "normalized")
      )]
      dots <- dots[setdiff(names(dots), c("feature_source", "normalized"))]
      pia <- do.call(pairwise_imd, c(list(x = x), pairwise_args))
      res <- do.call(plot_circos, c(list(mat = pia), dots))
      invisible(res)
    },

    composition = {
      if (is.null(annotation)) {
        stop(
          "`type = \"composition\"` needs an `annotation` vector: external ",
          "sample labels (e.g. tumour subtypes), one per sample, in the same ",
          "sample order as the fitted data. The fit does not store ",
          "annotations. Example:\n",
          "  plot(fit, type = \"composition\", annotation = subtype)",
          call. = FALSE
        )
      }
      if (is.null(group)) {
        stop(
          "`type = \"composition\"` needs cluster labels, but ",
          "`get_clusters()` found none in this fit. ",
          "Pass labels explicitly via `group`.",
          call. = FALSE
        )
      }
      if (length(annotation) != length(group)) {
        stop(
          "`annotation` has length ", length(annotation),
          " but the fit provides ", length(group), " cluster labels. ",
          "Supply one annotation per fitted sample, in the same order.",
          call. = FALSE
        )
      }
      dots <- dots[setdiff(names(dots), c("source", "omics", "cluster"))]
      do.call(
        plot_cluster_composition,
        c(list(cluster = group, annotation = annotation), dots)
      )
    },

    km = {
      missing_args <- c(
        if (is.null(time_var)) "`time_var`",
        if (is.null(event_var)) "`event_var`",
        if (is.null(pheno_mat)) "`pheno_mat`"
      )
      if (length(missing_args) > 0L) {
        stop(
          "`type = \"km\"` needs survival data that the fit does not store; ",
          "missing: ", paste(missing_args, collapse = ", "), ". Example:\n",
          "  plot(fit, type = \"km\", time_var = \"os_time\", ",
          "event_var = \"os_event\", pheno_mat = clinical)",
          call. = FALSE
        )
      }
      if (is.null(group)) {
        stop(
          "`type = \"km\"` needs sample groups, but `get_clusters()` found ",
          "no cluster labels in this fit. Pass groups explicitly via `group`.",
          call. = FALSE
        )
      }
      if (!is.data.frame(pheno_mat)) pheno_mat <- as.data.frame(pheno_mat)
      for (nm in c(time_var, event_var)) {
        if (length(nm) != 1L || !is.character(nm) || !nm %in% names(pheno_mat)) {
          stop(
            "`time_var` and `event_var` must each name one column of ",
            "`pheno_mat`. Column `", nm, "` was not found.",
            call. = FALSE
          )
        }
      }
      if (length(group) != nrow(pheno_mat)) {
        stop(
          "`pheno_mat` has ", nrow(pheno_mat), " rows but the fit provides ",
          length(group), " cluster labels. Supply one phenotype row per ",
          "fitted sample, in the same order.",
          call. = FALSE
        )
      }
      ## Cluster labels from the fit are group labels, never a continuous
      ## score; make that explicit for plot_km().
      test_var_use <- if (group_supplied) group else factor(group)
      dots <- dots[setdiff(names(dots), c("source", "omics", "cluster"))]
      if (!"palette" %in% names(dots)) {
        complete <- !is.na(test_var_use) &
          !is.na(pheno_mat[[time_var]]) & !is.na(pheno_mat[[event_var]])
        is_continuous <- is.numeric(test_var_use) &&
          length(unique(test_var_use[complete])) > 5L
        if (is_continuous) {
          ## plot_km() will dichotomize into Low/High
          dots$palette <- mrf3_palette(2L)
        } else {
          f_all <- factor(test_var_use)
          pal <- stats::setNames(mrf3_palette(nlevels(f_all)), levels(f_all))
          keep_levels <- levels(droplevels(f_all[complete]))
          dots$palette <- unname(pal[keep_levels])
        }
      }
      do.call(
        plot_km,
        c(
          list(
            test_var = test_var_use,
            time_var = time_var,
            event_var = event_var,
            pheno_mat = pheno_mat
          ),
          dots
        )
      )
    },

    weights = {
      if (is.null(x$imd) && is.null(x$cluster_imd)) {
        stop(
          "`type = \"weights\"` requires IMD weights, but this fit has none. ",
          "Refit with `run_imd = TRUE` (or `run_variable_selection = TRUE`, ",
          "or `run_cluster_imd = TRUE` for cluster-level IMD).",
          call. = FALSE
        )
      }
      if (!is.null(dots$source) && is.null(dots$weight_source)) {
        dots$weight_source <- dots$source
        dots$source <- NULL
      }
      do.call(plot_weights, c(list(weights = x), dots))
    }
  )
}
