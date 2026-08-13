#' Multiple plot functions
#' @param dat A data frame or matrix
#' @param weights Weight vector or weight-carrying object for `plot_weights()`.
#' @param mat Adjacency-like matrix for `plot_circos()`.
#' @param mod A fitted model object used to compute embeddings internally.
#' Can also be an `mrf3_fit` object.
#' @param group Class group
#' @param label_group Logical; whether to draw group labels on embeddings.
#' @param position Legend position.
#' @param perplexity tSNE perplexity used by `Rtsne::Rtsne()`.
#' @param pca Logical; whether to run PCA before embedding.
#' @param ncomp Number of principal components retained when `pca = TRUE`.
#' @param size Point size in scatter plots.
#' @param shape Point shape in embedding plots.
#' @param base_size Base font size used by the multiRF plotting theme.
#' @param base_family Font family used by the multiRF plotting theme.
#' @param seed Optional non-negative integer seed for stochastic embedding or
#' network-layout initialization. The caller's RNG state is restored.
#' @param source Source used when `dat`/`mod`/`weights` is an `mrf3_fit`.
#' Supported values for `mrf3_fit` inputs:
#' `"auto"`, `"specific_shared"`, `"specific_specific"`, `"clustering"`,
#' `"robust_clustering"`, `"reconstruction_weight_all"`,
#' `"reconstruction_weight_block"`, `"reconstruction_fused_block"`.
#' `plot_network()` and `plot_circos()` additionally accept a
#' `pairwise_imd_analysis` object with `source = "adj_var"` or `"adj_dat"`.
#' The native t-SNE model path requires a square sample-similarity source;
#' use the `dat` argument for a raw rectangular feature matrix.
#' For `plot_weights()`, supported values are:
#' `"auto"`, `"imd"`, `"cluster_imd"`, `"mrf"`.
#' @param omics Optional omics block name when `source` is block-specific.
#' @param cluster Optional cluster label when `source` uses cluster-level IMD output.
#' @param layout Network layout function for graph plotting.
#' @param vertex.size Vertex size in network plots.
#' @param label.dist Label distance in network plots.
#' @param edge.width Edge width in network plots.
#' @param vertex.frame.color Vertex frame color in network plots.
#' @param vertex.label.color Vertex label color in network plots.
#' @param vertex.label.cex Vertex label font size in network plots.
#' @param vertex.label.dist Vertex label distance in network plots.
#' @param edge.curved Edge curvature in network plots.
#' @param vertex.label.degree Vertex label angle in network plots.
#' @param plot.which Which omics block(s) to plot in `plot_weights()`.
#' @param ncol Number of facet columns in `plot_weights()`. By default, plots
#' with more than two blocks use two columns to preserve readable feature
#' labels.
#' @param labels Optional panel labels in `plot_weights()`.
#' @param weight_source Alias of `source` used by `plot_weights()`.
#' @param names.list Named list of variable groups for `plot_circos()`.
#' @param cut.off Cutoff applied before plotting a circos graph.
#' @param highlight Vector of sectors to highlight in `plot_circos()`.
#' @param pch Point symbol for base plotting functions.
#' @param config UMAP configuration object passed to `umap::umap()`.
#' @param method UMAP backend method passed to `umap::umap()`.
#' @inheritParams graphics::plot
#' @return `plot_tsne()`, `plot_umap()`, and `plot_weights()` return a
#' `ggplot` object. `plot_network()` invisibly returns an `igraph` object,
#' `plot_circos()` invisibly returns the chord-diagram result, and the
#' deprecated `plot_embed()` invisibly returns `NULL` after drawing.
#'
#' @name plot_tsne
NULL


# ---- mrf3-aware helpers ----
mrf3_plot_palette <- function(n = NULL) {
  palette <- c(
    signal_blue = "#3182BD",
    signal_teal = "#33B5A5",
    accent_orange = "#E28E2C",
    accent_red = "#D24B40",
    violet = "#756BB1",
    sky_blue = "#56B4E9",
    rose = "#CC79A7",
    olive = "#8C8C42",
    neutral_mid = "#767676",
    neutral_dark = "#272727",
    neutral_light = "#D8D8D8"
  )
  if (is.null(n)) {
    return(palette)
  }
  n <- as.integer(n)[1]
  if (!is.finite(n) || n < 1L) {
    return(character(0))
  }
  cluster_palette <- palette[seq_len(8L)]
  if (n <= length(cluster_palette)) {
    return(unname(cluster_palette[seq_len(n)]))
  }
  grDevices::hcl.colors(n, palette = "Dark 3")
}

theme_mrf3 <- function(base_size = 8, base_family = "sans") {
  if (length(base_size) != 1L || !is.numeric(base_size) ||
      !is.finite(base_size) || base_size <= 0) {
    stop("`base_size` must be a finite positive number.", call. = FALSE)
  }
  if (length(base_family) != 1L || !is.character(base_family) ||
      is.na(base_family)) {
    stop("`base_family` must be a single character string.", call. = FALSE)
  }
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#272727"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "#272727"),
      axis.title = ggplot2::element_text(size = base_size, colour = "#272727"),
      axis.text = ggplot2::element_text(size = max(base_size - 0.5, 1), colour = "#272727"),
      legend.title = ggplot2::element_text(size = max(base_size - 0.3, 1)),
      legend.text = ggplot2::element_text(size = max(base_size - 0.7, 1)),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = max(base_size - 0.3, 1), face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      panel.grid = ggplot2::element_blank()
    )
}

plot_embedding_gg <- function(embed, group = NULL, axes = c("Dim 1", "Dim 2"),
                              main = NULL, label_group = TRUE, position = "right",
                              size = 1.2, shape = 16, base_size = 8,
                              base_family = "sans") {
  embed <- as_numeric_matrix(embed, "embedding")
  if (ncol(embed) < 2L || nrow(embed) < 1L || any(!is.finite(embed[, 1:2, drop = FALSE]))) {
    stop("`embedding` must contain at least two finite columns and one row.")
  }
  n <- nrow(embed)
  if (is.null(group)) {
    group <- factor(rep("Samples", n))
  } else {
    if (length(group) != n) {
      stop("Length of `group` (", length(group), ") does not match embedding rows (", n, ").")
    }
    group <- droplevels(factor(group))
    if (anyNA(group)) {
      stop("`group` must not contain missing values.")
    }
  }
  df <- data.frame(dim1 = embed[, 1], dim2 = embed[, 2], group = group)
  palette <- stats::setNames(mrf3_plot_palette(nlevels(group)), levels(group))
  use_fill <- identical(as.integer(shape), 21L)
  if (use_fill) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, fill = group)) +
      ggplot2::geom_point(
        size = size, shape = shape, alpha = 0.92,
        colour = "#272727", stroke = 0.18
      ) +
      ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
      ggplot2::labs(x = axes[1], y = axes[2], title = main, fill = NULL)
    legend_aesthetic <- "fill"
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, colour = group)) +
      ggplot2::geom_point(size = size, shape = shape, alpha = 0.9) +
      ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
      ggplot2::labs(x = axes[1], y = axes[2], title = main, colour = NULL)
    legend_aesthetic <- "colour"
  }
  p <- p +
    ggplot2::coord_equal() +
    theme_mrf3(base_size = base_size, base_family = base_family) +
    ggplot2::theme(legend.position = position)
  if (nlevels(group) == 1L || isTRUE(label_group)) {
    p <- p + if (identical(legend_aesthetic, "fill")) {
      ggplot2::guides(fill = "none")
    } else {
      ggplot2::guides(colour = "none")
    }
  } else {
    guide <- ggplot2::guide_legend(
      override.aes = list(size = max(size, 2.2), alpha = 1)
    )
    p <- p + if (identical(legend_aesthetic, "fill")) {
      ggplot2::guides(fill = guide)
    } else {
      ggplot2::guides(colour = guide)
    }
  }
  if (isTRUE(label_group) && nlevels(group) > 1L) {
    centroids <- stats::aggregate(cbind(dim1, dim2) ~ group, data = df, FUN = mean)
    centers <- do.call(rbind, lapply(seq_len(nrow(centroids)), function(i) {
      candidates <- df[df$group == centroids$group[i], , drop = FALSE]
      nearest <- which.min(
        (candidates$dim1 - centroids$dim1[i])^2 +
          (candidates$dim2 - centroids$dim2[i])^2
      )
      candidates[nearest, c("dim1", "dim2", "group"), drop = FALSE]
    }))
    p <- p + ggrepel::geom_text_repel(
      data = centers,
      ggplot2::aes(label = group),
      size = 2.4,
      colour = "#272727",
      min.segment.length = 0,
      seed = 529,
      show.legend = FALSE,
      max.overlaps = Inf
    )
  }
  p
}

as_numeric_matrix <- function(x, label = "input") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`", label, "` must be a numeric matrix/data.frame.")
  }
  if (length(x) == 0L || any(!is.finite(x))) {
    stop("`", label, "` must be non-empty and contain only finite values.")
  }
  x
}

resolve_choice <- function(value, choices, what = "value") {
  if (length(choices) == 0L) {
    stop("No available ", what, ".")
  }
  if (is.null(value)) {
    return(choices[[1]])
  }
  value <- as.character(value)[1]
  if (!value %in% choices) {
    stop("Unknown ", what, ": `", value, "`. Available: ", paste(choices, collapse = ", "))
  }
  value
}

extract_mrf3_plot_matrix <- function(x,
                                         source = c(
                                           "auto",
                                           "specific_shared",
                                           "specific_specific",
                                           "clustering",
                                           "robust_clustering",
                                           "reconstruction_weight_all",
                                           "reconstruction_weight_block",
                                           "reconstruction_fused_block"
                                         ),
                                         omics = NULL,
                                         cluster = NULL) {
  source <- match.arg(source)

  get_matrix <- function(src) {
    if (identical(src, "specific_shared")) {
      S <- x$shared$clustering$similarity
      if (is.null(S)) stop("`shared$clustering$similarity` is unavailable.")
      return(as_numeric_matrix(S, "specific_shared"))
    }

    if (identical(src, "specific_specific")) {
      by_omics <- x$specific$clustering$by_omics
      if (!is.list(by_omics) || length(by_omics) == 0L) {
        stop("`specific$clustering$by_omics` is unavailable.")
      }
      block <- resolve_choice(omics, names(by_omics), what = "omics")
      S <- by_omics[[block]]$similarity
      if (is.null(S)) stop("`specific_specific` similarity is unavailable for omics `", block, "`.")
      return(as_numeric_matrix(S, paste0("specific_specific:", block)))
    }

    if (identical(src, "clustering")) {
      warning("`source = 'clustering'` is an alias of `specific_shared`; use `specific_shared`.", call. = FALSE)
      S <- x$shared$clustering$similarity
      if (is.null(S)) stop("`shared$clustering$similarity` is unavailable.")
      return(as_numeric_matrix(S, "shared$clustering$similarity"))
    }

    if (identical(src, "robust_clustering")) {
      S <- x$robust_detail$shared$clustering$similarity
      if (is.null(S)) stop("`robust_detail$shared$clustering$similarity` is unavailable.")
      return(as_numeric_matrix(S, "robust_detail$shared$clustering$similarity"))
    }

    if (identical(src, "reconstruction_weight_all")) {
      S <- x$reconstruction$W$W_all
      if (is.null(S)) stop("`reconstruction$W$W_all` is unavailable.")
      return(as_numeric_matrix(S, "reconstruction$W$W_all"))
    }

    if (identical(src, "reconstruction_weight_block")) {
      by_block <- x$reconstruction$W$W_per_response
      if (!is.list(by_block) || length(by_block) == 0L) {
        stop("`reconstruction$W$W_per_response` is unavailable.")
      }
      block <- resolve_choice(omics, names(by_block), what = "omics")
      return(as_numeric_matrix(by_block[[block]], paste0("reconstruction$W$W_per_response:", block)))
    }

    if (identical(src, "reconstruction_fused_block")) {
      fused <- x$reconstruction$fused_mat
      if (!is.list(fused) || length(fused) == 0L) {
        stop("`reconstruction$fused_mat` is unavailable.")
      }
      block <- resolve_choice(omics, names(fused), what = "omics")
      return(as_numeric_matrix(fused[[block]], paste0("reconstruction$fused_mat:", block)))
    }

    stop("Unsupported mrf3 source: `", src, "`.")
  }

  if (!identical(source, "auto")) {
    return(get_matrix(source))
  }

  auto_candidates <- c(
    "specific_shared",
    "clustering",
    "robust_clustering",
    "reconstruction_weight_all",
    "reconstruction_fused_block"
  )
  for (src in auto_candidates) {
    candidate <- tryCatch(get_matrix(src), error = function(e) NULL)
    if (!is.null(candidate)) {
      return(candidate)
    }
  }
  stop("Cannot auto-resolve a plotting matrix from this mrf3_fit object. Set `source` explicitly.")
}

extract_plot_matrix <- function(x, source = "auto", omics = NULL, cluster = NULL) {
  if (inherits(x, "pairwise_imd_analysis")) {
    source <- match.arg(source, c("auto", "adj_var", "adj_dat"))
    if (identical(source, "auto") || identical(source, "adj_var")) {
      return(as_numeric_matrix(x$adj_var_mat, "pairwise_imd_analysis$adj_var_mat"))
    }
    return(as_numeric_matrix(x$adj_dat_mat, "pairwise_imd_analysis$adj_dat_mat"))
  }
  if (is_mrf3_fit_object(x)) {
    return(extract_mrf3_plot_matrix(
      x = x,
      source = source,
      omics = omics,
      cluster = cluster
    ))
  }
  as_numeric_matrix(x, "dat")
}

as_weight_list <- function(x) {
  if (is.numeric(x) && !is.list(x)) {
    x <- list(weights = x)
  }
  if (!is.list(x) || length(x) == 0L) {
    stop("`weights` must be a non-empty list (or numeric vector).")
  }
  out <- lapply(seq_along(x), function(i) {
    w <- x[[i]]
    if (is.null(w)) {
      return(NULL)
    }
    if (!is.numeric(w)) {
      stop("Every weight block must be numeric.")
    }
    w0 <- as.numeric(w)
    if (any(!is.finite(w0))) {
      stop("Weight blocks must contain only finite values.")
    }
    nm <- names(w)
    if (is.null(nm) || any(!nzchar(nm))) {
      nm <- paste0("V", seq_along(w0))
    }
    names(w0) <- make.unique(nm)
    w0
  })
  names(out) <- names(x)
  out <- purrr::compact(out)
  if (length(out) == 0L) {
    stop("No valid numeric weights found.")
  }
  if (is.null(names(out)) || any(!nzchar(names(out)))) {
    names(out) <- paste0("omics", seq_along(out))
  }
  out
}

extract_plot_weights <- function(x, source = c("auto", "imd", "cluster_imd", "mrf"), cluster = NULL) {
  source <- match.arg(source)

  if (identical(source, "mrf")) {
    wts <- x$imd
    if (!is.list(x) || is.null(wts)) {
      stop("`source = 'mrf'` requires an mrf3-like object with `$imd`.")
    }
    return(as_weight_list(wts))
  }

  if (is_mrf3_fit_object(x)) {
    if (source %in% c("auto", "imd")) {
      wts <- x$imd
      if (!is.null(wts)) {
        return(as_weight_list(wts))
      }
      if (identical(source, "imd")) {
        stop("`weights` is unavailable in mrf3_fit object.")
      }
    }

    if (source %in% c("auto", "cluster_imd")) {
      by_cluster <- x$cluster_imd$by_cluster
      if (is.list(by_cluster) && length(by_cluster) > 0L) {
        cl <- resolve_choice(cluster, names(by_cluster), what = "cluster")
        w <- by_cluster[[cl]]$imd$weight_list
        if (!is.null(w)) {
          return(as_weight_list(w))
        }
      }
      if (identical(source, "cluster_imd")) {
        stop("`cluster_imd` weights are unavailable in mrf3_fit object.")
      }
    }

    stop("Cannot auto-resolve weights from mrf3_fit object. Use `source = 'imd'` or `source = 'cluster_imd'`.")
  }

  if (is.list(x) && (inherits(x, "mrf3") || identical(source, "auto"))) {
    wts <- x$imd
    if (!is.null(wts)) return(as_weight_list(wts))
  }

  as_weight_list(x)
}


# plot tSNE
#' @rdname plot_tsne
#' @export
plot_tsne <- function(dat = NULL, mod = NULL, group = NULL, label_group = TRUE,
                      position = "right", perplexity = 30, pca = FALSE,
                      ncomp = 70, main = "t-SNE", size = 1.4, shape = 21,
                      source = "auto", omics = NULL, cluster = NULL,
                      seed = NULL, base_size = 8, base_family = "sans", ...){

  perplexity_was_supplied <- !missing(perplexity)
  if (!is.logical(pca) || length(pca) != 1L || is.na(pca)) {
    stop("`pca` must be TRUE or FALSE.", call. = FALSE)
  }
  ncomp <- .tsne_integer_scalar(ncomp, "ncomp")
  if (length(perplexity) != 1L || !is.numeric(perplexity) ||
      !is.finite(perplexity) || perplexity <= 0) {
    stop("`perplexity` must be a finite positive number.", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .tsne_integer_scalar(seed, "seed", minimum = 0)
  }

  if (!is.null(mod)) {
    if (is_mrf3_fit_object(mod)) {
      mod <- list(dat = extract_mrf3_plot_matrix(mod, source = source, omics = omics, cluster = cluster))
    } else if (is.matrix(mod) || is.data.frame(mod)) {
      mod <- list(dat = as_numeric_matrix(mod, "mod"))
    }
    if (is.list(mod) && !is.null(mod$dat)) {
      mod$dat <- as_numeric_matrix(mod$dat, "mod$dat")
    } else {
      stop("`mod` must contain numeric `mod$dat`, or be a matrix, or be an `mrf3_fit` object.")
    }
    embed <- mrf3_tsne(mod, seed = seed, ...)
  } else {
    if (is_mrf3_fit_object(dat)) {
      dat_wf <- extract_mrf3_plot_matrix(dat, source = source, omics = omics, cluster = cluster)
      embed <- mrf3_tsne(list(dat = dat_wf), seed = seed, ...)
    } else {
      dat <- as_numeric_matrix(dat, "dat")
      if (!requireNamespace("Rtsne", quietly = TRUE)) {
        stop("Package `Rtsne` is required when `plot_tsne()` receives a raw data matrix.",
             call. = FALSE)
      }
      if (nrow(dat) < 5L) {
        stop("Raw-data t-SNE requires at least five observations.", call. = FALSE)
      }
      max_perplexity <- max(1L, floor((nrow(dat) - 2L) / 3L))
      if (perplexity > max_perplexity) {
        if (perplexity_was_supplied) {
          warning(
            "`perplexity` exceeds the value supported by this sample size; using ",
            max_perplexity, ".", call. = FALSE
          )
        }
        perplexity <- max_perplexity
      }
      if (!is.null(seed)) {
        old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit({
          if (old_seed_exists) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
          } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        }, add = TRUE)
        set.seed(seed)
      }
      t <- Rtsne::Rtsne(dat, perplexity = perplexity, pca = pca,
                        initial_dims = min(ncomp, ncol(dat)), ...)
      embed <- t$Y
    }
  }
  plot_embedding_gg(
    embed = embed,
    group = group,
    axes = c("t-SNE 1", "t-SNE 2"),
    main = main,
    label_group = label_group,
    position = position,
    size = size,
    shape = shape,
    base_size = base_size,
    base_family = base_family
  )
}

#' @param dat A data frame or matrix
#' @param group Class group
#' @param label A logical parameter that determine whether to show labels or not.
#' @param cutoff Edge-weight cutoff. For a regular sample matrix the default
#' is the mean positive upper-triangle weight; for pairwise IMD it is the 75th
#' percentile of positive upper-triangle weights.
#' @param ... Additional arguments passed to the underlying plotting routine.
#'
#' @export
#' @rdname plot_tsne

# plot sample or feature network
plot_network <- function(dat, group = NULL, label = FALSE, cutoff = NULL,
                         layout = igraph::layout_with_fr, vertex.size = 5, label.dist = NULL,
                         edge.width = NULL, vertex.frame.color = "white", vertex.label.color = "black",
                         vertex.label.cex = 0.5, vertex.label.dist = 0,
                         edge.curved = 0.5,vertex.label.degree =  -pi/2,
                         position = "bottomright", source = "auto", omics = NULL,
                         cluster = NULL, seed = NULL, ...){
  if (!is.null(label.dist)) {
    warning("`label.dist` is deprecated; use `vertex.label.dist`.", call. = FALSE)
    vertex.label.dist <- label.dist
  }
  if (inherits(dat, "pairwise_imd_analysis")) {
    source <- match.arg(source, c("auto", "adj_var", "adj_dat"))
    mat <- if (identical(source, "adj_dat")) dat$adj_dat_mat else dat$adj_var_mat
    var_use <- dat$var_use

    if (is.null(group) && !is.null(var_use)) {
      feat_group <- character(nrow(mat))
      names(feat_group) <- rownames(mat)
      for (block in names(var_use)) {
        feat_group[intersect(var_use[[block]], rownames(mat))] <- block
      }
      feat_group[feat_group == ""] <- "other"
      group <- feat_group[rownames(mat)]
    }

    mat <- as_numeric_matrix(mat, "pairwise adjacency")
    if (is.null(cutoff)) {
      positive <- mat[upper.tri(mat) & mat > 0]
      cutoff <- if (length(positive)) stats::quantile(positive, 0.75, na.rm = TRUE) else Inf
    }
  } else {
    mat <- extract_plot_matrix(dat, source = source, omics = omics, cluster = cluster)
    if (is.null(cutoff)) {
      positive <- mat[upper.tri(mat) & mat > 0]
      cutoff <- if (length(positive)) mean(positive) else Inf
    }
  }
  if (nrow(mat) != ncol(mat)) {
    stop("`plot_network()` requires a square adjacency/similarity matrix.")
  }
  if (length(cutoff) != 1L || !is.numeric(cutoff) || is.na(cutoff)) {
    stop("`cutoff` must be NULL or a single non-missing numeric value.", call. = FALSE)
  }
  if (!is.null(group)) {
    if (length(group) != nrow(mat)) stop("`group` must have one value per matrix row.")
    group <- as.character(group)
    if (anyNA(group)) stop("`group` must not contain missing values.")
  }
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 0
  mat[mat < cutoff] <- 0
  network <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  keep <- igraph::degree(network) > 0
  if (!any(keep)) {
    warning("No edges remain after applying `cutoff`; returning an empty graph.", call. = FALSE)
    network <- igraph::delete_vertices(network, igraph::V(network))
    return(invisible(network))
  }
  if (!is.null(group)) {
    group <- group[keep]
  }
  network <- igraph::delete_vertices(network, which(!keep))
  if (is.null(group)) {
    vertex_col <- mrf3_plot_palette()[["signal_blue"]]
    group_palette <- NULL
  } else {
    levels_group <- sort(unique(group))
    group_palette <- stats::setNames(mrf3_plot_palette(length(levels_group)), levels_group)
    vertex_col <- unname(group_palette[group])
    igraph::V(network)$group <- group
  }
  igraph::V(network)$color <- vertex_col
  edge_weights <- igraph::E(network)$weight
  edge_width_use <- edge.width
  if (is.null(edge_width_use)) {
    if (length(edge_weights) && diff(range(edge_weights)) > 0) {
      edge_width_use <- 0.25 + 1.25 *
        (edge_weights - min(edge_weights)) / diff(range(edge_weights))
    } else {
      edge_width_use <- rep(0.6, length(edge_weights))
    }
  }
  vlab <- if (isTRUE(label)) igraph::V(network)$name else NA
  if (!is.null(seed)) {
    seed <- .tsne_integer_scalar(seed, "seed", minimum = 0)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  graphics::plot(
    network,
    vertex.color = vertex_col,
    vertex.size = vertex.size,
    vertex.label = vlab,
    vertex.label.cex = vertex.label.cex,
    vertex.label.color = vertex.label.color,
    vertex.label.dist = vertex.label.dist,
    vertex.frame.color = vertex.frame.color,
    vertex.label.degree = vertex.label.degree,
    edge.width = edge_width_use,
    edge.color = grDevices::adjustcolor(mrf3_plot_palette()[["neutral_mid"]], alpha.f = 0.55),
    edge.curved = edge.curved,
    layout = layout,
    ...
  )
  if (!is.null(group_palette)) {
    graphics::legend(position, fill = group_palette, legend = names(group_palette),
                     border = NA, bty = "n", cex = 0.75)
  }
  invisible(network)
}


#' @param dat A data frame or matrix
#' @param group Class group
#' @param top The number of top weighted variables to show in the plot. The default is 20. Can be chosen from numerical values or all for setting the parameter = NULL.
#' @import ggplot2
#'
#' @export
#' @rdname plot_tsne

plot_weights <- function(weights, plot.which = "all", top = 20, labels = NULL,
                         weight_source = "auto", cluster = NULL,
                         ncol = NULL, base_size = 8, base_family = "sans"){

  weights <- extract_plot_weights(weights, source = weight_source, cluster = cluster)
  if (length(plot.which) == 1L && identical(toupper(as.character(plot.which)), "ALL")) {
    imp <- weights
  } else {
    plot.which <- as.character(plot.which)
    missing_blocks <- setdiff(plot.which, names(weights))
    if (length(missing_blocks)) {
      stop("Unknown `plot.which`: ", paste(missing_blocks, collapse = ", "),
           ". Available: ", paste(names(weights), collapse = ", "), ".")
    }
    imp <- weights[plot.which]
  }
  if (!is.null(top)) {
    if (length(top) != 1L || !is.numeric(top) || !is.finite(top) || top < 1 || top != floor(top)) {
      stop("`top` must be NULL or a positive integer.")
    }
    top <- as.integer(top)
  }
  panel_names <- names(imp)
  if (!is.null(labels)) {
    if (length(labels) != length(panel_names)) stop("`labels` must match the selected number of blocks.")
    panel_names <- as.character(labels)
  }
  if (is.null(ncol)) {
    ncol <- if (length(imp) <= 2L) length(imp) else 2L
  } else {
    ncol <- .tsne_integer_scalar(ncol, "ncol")
  }
  dfs <- lapply(seq_along(imp), function(i) {
    w <- imp[[i]]
    ord <- order(w, decreasing = TRUE, na.last = NA)
    if (!is.null(top)) ord <- utils::head(ord, top)
    data.frame(
      block = panel_names[i],
      variable = names(w)[ord],
      weight = unname(w[ord]),
      rank = seq_along(ord),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, dfs)
  if (!nrow(df)) stop("No weights are available to plot.")
  df$variable_panel <- interaction(df$block, df$variable, drop = TRUE, lex.order = TRUE)
  level_order <- rev(df$variable_panel[order(df$block, df$rank)])
  df$variable_panel <- factor(df$variable_panel, levels = unique(level_order))
  ggplot2::ggplot(df, ggplot2::aes(x = variable_panel, y = weight)) +
    ggplot2::geom_col(width = 0.72, fill = mrf3_plot_palette()[["signal_blue"]]) +
    ggplot2::facet_wrap(~block, scales = "free_y", ncol = ncol) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(labels = stats::setNames(df$variable, df$variable_panel)) +
    ggplot2::labs(x = NULL, y = "IMD weight") +
    theme_mrf3(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.spacing = grid::unit(6, "pt")
    )
}



#' @export
#' @rdname plot_tsne
plot_circos <- function(mat, names.list = NULL, group = NULL, cut.off = NULL, highlight = NULL, source = c("adj_var", "adj_dat"), ...){
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package `circlize` is required for `plot_circos()`. Install it with install.packages('circlize').")
  }
  source <- match.arg(source)
  if (inherits(mat, "pairwise_imd_analysis")) {
    obj <- mat
    mat <- if (identical(source, "adj_var")) obj$adj_var_mat else obj$adj_dat_mat
    if (is.null(names.list) && identical(source, "adj_var")) {
      names.list <- obj$var_use
    }
  }
  if (is.null(names.list) || !is.list(names.list) || length(names.list) == 0L) {
    stop("`names.list` must be a non-empty named list. For `pairwise_imd_analysis`, omit it only when `source = 'adj_var'`.")
  }

  vimp <- as_numeric_matrix(mat, "mat")
  if (nrow(vimp) != ncol(vimp) || is.null(rownames(vimp)) || is.null(colnames(vimp))) {
    stop("`mat` must be a square named matrix.")
  }
  if (!identical(rownames(vimp), colnames(vimp))) {
    stop("`mat` must use the same row and column names in the same order.", call. = FALSE)
  }
  if (is.null(names(names.list)) || any(!nzchar(names(names.list)))) {
    stop("`names.list` must have non-empty group names.", call. = FALSE)
  }
  listed <- unlist(names.list, use.names = FALSE)
  if (anyDuplicated(listed) || !setequal(listed, rownames(vimp))) {
    stop("`names.list` must assign every matrix name to exactly one group.", call. = FALSE)
  }
  if (is.null(group)) group <- names(names.list)
  if (length(group) != length(names.list)) stop("`group` must have one label per entry in `names.list`.")
  group <- as.character(group)

  diag(vimp) <- 0
  if(!is.null(cut.off)){
    if (length(cut.off) != 1L || !is.numeric(cut.off) || is.na(cut.off)) {
      stop("`cut.off` must be NULL or a single non-missing numeric value.", call. = FALSE)
    }
    vimp[vimp <= cut.off] <- 0
  }
  keep <- rowSums(abs(vimp)) != 0 | colSums(abs(vimp)) != 0
  vimp <- vimp[keep, keep, drop = FALSE]
  names.list <- lapply(names.list, function(l) l[l %in% colnames(vimp)])
  keep_groups <- lengths(names.list) > 0L
  names.list <- names.list[keep_groups]
  group <- group[keep_groups]
  if (!length(names.list) || !length(vimp) || !any(vimp != 0)) {
    stop("No non-zero links remain after applying `cut.off`.")
  }

  # Circular Network Diagram Plot

  group.list <- lapply(seq_along(group), function(g){
    structure(rep(group[g], length(names.list[[g]])), names =  names.list[[g]])
  }
  )


  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  sector_group <- unlist(group.list)
  group_palette <- stats::setNames(mrf3_plot_palette(length(unique(group))), unique(group))
  sector_colors <- unname(group_palette[sector_group])
  names(sector_colors) <- names(sector_group)
  out <- circlize::chordDiagram(
    vimp,
    annotationTrack = "grid",
    group = sector_group,
    grid.col = sector_colors,
    preAllocateTracks = list(
      track.height = circlize::mm_h(4),
      track.margin = c(circlize::mm_h(4), 0)
    ),
    ...
  )


  circlize::circos.track(track.index = 2, panel.fun = function(x, y) {
    sector.index.all = circlize::get.cell.meta.data("sector.index")
    for(i in seq_along(names.list)){

      sector.index = sector.index.all[sector.index.all %in% names.list[[i]]]
      xlim = circlize::get.cell.meta.data("xlim")
      ylim = circlize::get.cell.meta.data("ylim")

      circlize::circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.55,
                            facing = "clockwise", niceFacing = TRUE)

    }}, bg.border = NA)



  for(i in seq_along(names.list)){
    circlize::highlight.sector(names.list[[i]], track.index = 1,
                     col = grDevices::adjustcolor(mrf3_plot_palette(length(names.list))[i], alpha.f = 0.22),
                     text = group[i], cex = 0.8, text.col = "#272727", niceFacing = TRUE)
  }

  if(!is.null(highlight)){
    highlight <- unique(intersect(highlight, unlist(names.list, use.names = FALSE)))
    if (length(highlight) > 0L) {
      circlize::highlight.sector(
        highlight,
        col = grDevices::adjustcolor(
          mrf3_plot_palette()[["accent_orange"]], alpha.f = 0.28
        )
      )
    }
  }


  invisible(out)
}



#' @export
#' @rdname plot_tsne
plot_embed <- function(dat, group = NULL, position = "bottom", pch = 20, ...){
  .Deprecated("plot_tsne", package = "multiRF",
              msg = "`plot_embed()` is deprecated; use `plot_tsne()` or `plot_umap()`.")
  dat <- as_numeric_matrix(dat, "dat")
  if (!is.null(group) && length(group) != nrow(dat)) stop("`group` must have one value per row of `dat`.")
  old_par <- graphics::par(c("xpd", "oma"))
  on.exit(graphics::par(old_par), add = TRUE)
  if(!is.null(group)){
    group <- as.character(group)
    l <- sort(na.omit(unique(group)))
    get_p <- mrf3_plot_palette(length(l))
    names(get_p) <- l
    col <- get_p[group]
  } else {col = 1}

  graphics::pairs(dat, col = col, pch = pch, oma = c(6, 3, 3, 3), ...)
  graphics::par(xpd = TRUE)
  if(!is.null(group)){
    graphics::legend(position, fill = get_p, legend = l, horiz=TRUE,
           xpd=TRUE, bty="n", cex = .5)
  }
  invisible(NULL)
}


#' @export
#' @rdname plot_tsne
plot_umap <- function(dat, group = NULL, main = "UMAP", label_group = TRUE,
                      pca = TRUE, ncomp = 70, position = "right", pch = 21,
                      size = 1.4, config = umap::umap.defaults,
                      method = "naive", source = "auto", omics = NULL,
                      cluster = NULL, seed = NULL, base_size = 8,
                      base_family = "sans", ...){

  if (!requireNamespace("umap", quietly = TRUE)) {
    stop("Package `umap` is required for `plot_umap()`.", call. = FALSE)
  }
  dat <- extract_plot_matrix(dat, source = source, omics = omics, cluster = cluster)
  if (nrow(dat) < 3L) {
    stop("UMAP requires at least three observations.", call. = FALSE)
  }
  if (!is.logical(pca) || length(pca) != 1L || is.na(pca)) {
    stop("`pca` must be TRUE or FALSE.", call. = FALSE)
  }
  ncomp <- .tsne_integer_scalar(ncomp, "ncomp")
  if (!is.null(seed)) {
    seed <- .tsne_integer_scalar(seed, "seed", minimum = 0)
  }

  if(pca) {
    pca_fit <- stats::prcomp(dat)
    x <- pca_fit$x[, seq_len(min(ncomp, ncol(pca_fit$x))), drop = FALSE]
  } else {
    x <- dat
  }

  config_use <- config
  if (!is.list(config_use) || is.null(config_use$n_neighbors) ||
      length(config_use$n_neighbors) != 1L ||
      !is.numeric(config_use$n_neighbors) ||
      !is.finite(config_use$n_neighbors)) {
    stop("`config` must be a UMAP configuration with a finite `n_neighbors` value.",
         call. = FALSE)
  }
  config_use$n_neighbors <- min(
    max(2L, as.integer(config_use$n_neighbors)),
    nrow(x) - 1L
  )
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  t <- umap::umap(x, config = config_use, method = method, ...)
  plot_embedding_gg(
    embed = t$layout,
    group = group,
    axes = c("UMAP 1", "UMAP 2"),
    main = main,
    label_group = label_group,
    position = position,
    size = size,
    shape = pch,
    base_size = base_size,
    base_family = base_family
  )
}
