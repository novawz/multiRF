#' Multiple plot functions
#' @param dat A data frame or matrix
#' @param weights Weight vector or weight-carrying object for `plot_weights()`.
#' @param mat Adjacency-like matrix for `plot_circos()`.
#' @param mod A fitted model object used to compute embeddings internally.
#' Can also be an `mrf3_fit` object.
#' @param group Class group
#' @param label_group Logical; whether to draw group labels on embeddings.
#' @param position Legend position. The default is "bottomright"
#' @param perplexity tSNE perplexity used by `Rtsne::Rtsne()`.
#' @param pca Logical; whether to run PCA before embedding.
#' @param ncomp Number of principal components retained when `pca = TRUE`.
#' @param size Point size in scatter plots.
#' @param source Source used when `dat`/`mod`/`weights` is an `mrf3_fit`.
#' Supported values for matrix-based plots:
#' `"auto"`, `"specific_shared"`, `"specific_specific"`, `"clustering"`,
#' `"robust_clustering"`, `"reconstruction_weight_all"`,
#' `"reconstruction_weight_block"`, `"reconstruction_fused_block"`.
#' When input is a `pairwise_imd_analysis` object, supported values are
#' `"auto"`, `"adj_var"`, and `"adj_dat"`.
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
#' @param labels Optional panel labels in `plot_weights()`.
#' @param weight_source Alias of `source` used by `plot_weights()`.
#' @param names.list Named list of variable groups for `plot_circos()`.
#' @param cut.off Cutoff applied before plotting a circos graph.
#' @param highlight Vector of sectors to highlight in `plot_circos()`.
#' @param pch Point symbol for base plotting functions.
#' @param config UMAP configuration object passed to `umap::umap()`.
#' @param method UMAP backend method passed to `umap::umap()`.
#' @inheritParams graphics::plot
#' @return A plot object
#'
#' @name plot_tsne
NULL


# ---- mrf3-aware helpers ----
as_numeric_matrix <- function(x, label = "input") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`", label, "` must be a numeric matrix/data.frame.")
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
      S <- x$reconstruction$W$W_all
      if (is.null(S)) stop("`reconstruction$W$W_all` is unavailable.")
      return(as_numeric_matrix(S, "reconstruction$W$W_all"))
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
    w0 <- as.numeric(w)
    names(w0) <- names(w)
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
    if (!is.list(x) || is.null(x$weights)) {
      stop("`source = 'mrf'` requires an mrf3-like object with `$weights`.")
    }
    return(as_weight_list(x$weights))
  }

  if (is_mrf3_fit_object(x)) {
    if (source %in% c("auto", "imd")) {
      if (!is.null(x$weights)) {
        return(as_weight_list(x$weights))
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

  if (is.list(x) && !is.null(x$weights) && (inherits(x, "mrf3") || identical(source, "auto"))) {
    return(as_weight_list(x$weights))
  }

  as_weight_list(x)
}


# plot tSNE
#' @rdname plot_tsne
#' @export
plot_tsne <- function(dat = NULL, mod = NULL, group = NULL, label_group = TRUE, position = "right", perplexity = 30, pca = FALSE, ncomp = 70, main = "tSNE", size = .5, source = "auto", omics = NULL, cluster = NULL, ...){

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
    embed <- mrf3_tsne(mod, ...)
  } else {
    if (is_mrf3_fit_object(dat)) {
      dat_wf <- extract_mrf3_plot_matrix(dat, source = source, omics = omics, cluster = cluster)
      embed <- mrf3_tsne(list(dat = dat_wf), ...)
    } else {
      dat <- as_numeric_matrix(dat, "dat")
      t <- Rtsne::Rtsne(dat, perplexity = perplexity, pca = pca, initial_dims = ncomp, ...)
      embed <- t$Y
    }
  }

  if (is.null(group)) {
    group <- rep("black", nrow(embed))
  } else if (length(group) != nrow(embed)) {
    stop("Length of `group` (", length(group), ") does not match embedding rows (", nrow(embed), ").")
  }
  group <- as.factor(group)

  df <- data.frame(tSNE1 = embed[,1], tSNE2 = embed[,2], group = group)
  df <- na.omit(df)
  p1 <- ggplot(df, aes(tSNE1, tSNE2)) +
    geom_point(aes(color = group), size = size) +
    xlab("tSNE dim 1") +
    ylab("tSNE dim 2") +
    ggtitle(main) +
    theme_classic()

  if(length(unique(group)) == 1) {
    p1 <- p1 + guides(color = "none")
  } else {
    p1 <- p1 +
      ggsci::scale_color_igv() +
      guides(color = guide_legend(override.aes = list(size=3))) +
      theme(legend.position = position,
            legend.title = element_blank(),
            legend.text = element_text(size = 6))
  }

  if (label_group && length(unique(group)) != 1) {
    data2 <- df %>% group_by(group) %>% dplyr::select(tSNE1, tSNE2) %>% dplyr::summarise(dplyr::across(dplyr::everything(), mean))

    p1 <- p1 +
      ggrepel::geom_text_repel(data = data2, aes(label = group),
                               size = 2.5, color = "grey20",
                               max.overlaps = 15)
  }


  p1


}

#' @param dat A data frame or matrix
#' @param group Class group
#' @param label A logical parameter that determine whether to show labels or not.
#' @param cutoff Cutoff of edgeweights. The default is mean value of dat.
#' @param position Legend position. The default is "bottomright"
#' @param ... Additional arguments passed to the underlying plotting routine.
#'
#' @export
#' @rdname plot_tsne

# plot sample or feature network
plot_network <- function(dat, group = NULL, label = FALSE, cutoff = NULL,
                         layout = igraph::layout_with_fr, vertex.size = 5, label.dist =1,
                         edge.width = 0.3, vertex.frame.color = "white", vertex.label.color = "black",
                         vertex.label.cex = 0.5, vertex.label.dist = 0,
                         edge.curved = 0.5,vertex.label.degree =  -pi/2,
                         position = "bottomright", source = "auto", omics = NULL, cluster = NULL, ...){

  # Feature network mode: pairwise_imd_analysis input
  if (inherits(dat, "pairwise_imd_analysis")) {
    mat <- if (identical(source, "adj_dat")) dat$adj_dat_mat else dat$adj_var_mat
    var_use <- dat$var_use

    # Build group mapping: feature name -> omics block
    if (!is.null(var_use)) {
      feat_group <- character(nrow(mat))
      names(feat_group) <- rownames(mat)
      for (block in names(var_use)) {
        feat_group[intersect(var_use[[block]], rownames(mat))] <- block
      }
      feat_group[feat_group == ""] <- "other"
      group <- feat_group[rownames(mat)]
    }

    diag(mat) <- 0
    if (is.null(cutoff)) cutoff <- quantile(mat[mat > 0], 0.75, na.rm = TRUE)
    mat[mat < cutoff] <- 0

    network <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)

    if (!is.null(group)) {
      group <- as.character(group)
      n_class <- length(unique(group))
      get_p <- ggpubr::get_palette("jco", k = n_class)
      names(get_p) <- sort(unique(group))
      igraph::V(network)$class <- get_p[group]
    }

    # Remove isolated vertices
    network <- igraph::delete_vertices(network, igraph::degree(network) == 0)

    vlab <- if (isTRUE(label)) igraph::V(network)$name else NA

    plot(network,
         vertex.color = igraph::V(network)$class,
         vertex.size = vertex.size,
         vertex.label = vlab,
         vertex.label.cex = vertex.label.cex,
         vertex.label.color = vertex.label.color,
         vertex.label.dist = vertex.label.dist,
         vertex.frame.color = vertex.frame.color,
         edge.width = edge.width,
         edge.curved = edge.curved,
         layout = layout,
         ...)
    if (!is.null(group)) {
      legend(position, fill = get_p, legend = sort(unique(group)), bty = "n", cex = 0.8)
    }
    return(invisible(NULL))
  }

  # Sample network mode (original behavior)
  dat <- extract_plot_matrix(dat, source = source, omics = omics, cluster = cluster)
  if (nrow(dat) != ncol(dat)) {
    stop("`plot_network()` requires a square adjacency/similarity matrix.")
  }

  diag(dat) <- 0
  if (is.null(cutoff)) cutoff <- mean(dat)
  network <- igraph::graph_from_adjacency_matrix(dat,
                                                 mode = "undirected",
                                                 weighted = TRUE)

  if(!is.null(group)){
    group <- as.character(group)
    n_class <- length(unique(group))
    get_p <- ggpubr::get_palette("jco", k = n_class)
    names(get_p) <- sort(unique(group))
    col <- get_p[group]
    igraph::V(network)$class <- col
  }

  if(!label){
    label <- NA
  } else {label <- igraph::V(network)$name}

  network <- igraph::delete_edges(network, igraph::E(network)[which(igraph::E(network)$weight < cutoff)])
  network <- igraph::delete_vertices(network, igraph::degree(network) == 0)

  plot(network,
       vertex.color = igraph::V(network)$class,
       vertex.size = vertex.size,
       vertex.label = label,
       label.dist = label.dist,
       edge.width = edge.width,
       vertex.frame.color = vertex.frame.color,
       vertex.label.color = vertex.label.color,
       vertex.label.cex = vertex.label.cex,
       vertex.label.dist = vertex.label.dist,
       edge.curved = edge.curved,
       vertex.label.degree = -pi/2,
       layout = layout,
       ...)
  if(!is.null(group)){
    legend(position, fill = get_p, legend = unique(group))
  }

}


#' @param dat A data frame or matrix
#' @param group Class group
#' @param top The number of top weighted variables to show in the plot. The default is 20. Can be chosen from numerical values or all for setting the parameter = NULL.
#' @import ggplot2
#'
#' @export
#' @rdname plot_tsne

plot_weights <- function(weights, plot.which = "all", top = 20, labels = NULL, weight_source = "auto", cluster = NULL){

  weights <- extract_plot_weights(weights, source = weight_source, cluster = cluster)

  if(toupper(plot.which) == "ALL"){
    imp <- weights
  } else {
    imp <- weights[plot.which]
  }

  df <-  purrr::map(imp,
                    ~data.frame(Weights = sort(., decreasing = TRUE),
                                Var_names = names(.)[order(., decreasing = TRUE)],
                                Ranks = c(1:length(.))))

  if(!is.null(top)){

    df <- purrr::map(df, function(tb) {
      if (nrow(tb) == 0L) {
        return(tb)
      }
      tb[seq_len(min(as.integer(top), nrow(tb))), , drop = FALSE]
    })

    p <- plyr::llply(
      df,
      .fun = function(g){
        p <- ggplot(
          data = g,
          aes(x = reorder(Var_names, -order(Ranks, decreasing = FALSE)), y = Weights, fill = Weights)
        ) +
          geom_bar(stat = 'identity') +
          ggpubr::theme_pubr(base_size = 9) +
          theme(axis.text.x = element_text(size = .5)) +
          xlab("Variables") +
          coord_flip() +
          scale_fill_gradientn(colours = colorRampPalette(RColorBrewer::brewer.pal(8, "Oranges"))(10))

      }
    )

  } else {

    p <- plyr::llply(
      df,
      .fun = function(g){
        p <- ggplot(
          data = g,
          aes(x = Var_names, y = Weights)
        ) +
          geom_point() +
          ggpubr::theme_pubr(base_size = 9) +
          theme(axis.text.x = element_blank()) +
          xlab("Variables")

      }
    )

  }
  if(is.null(labels)) labels <- names(df)
  ggpubr::ggarrange(plotlist = p, nrow = 1, labels = labels)

}



#' @export
#' @rdname plot_tsne
plot_circos <- function(mat, names.list = NULL, group = NULL, cut.off = NULL, highlight = NULL, source = c("adj_var", "adj_dat"), ...){

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

  vimp <- mat

  group <- names(names.list)

  if(!is.null(cut.off)){
    vimp[vimp <= cut.off] <- 0
    vimp <- vimp[ rowSums(vimp) != 0 ,]
    vimp <- vimp[, colnames(vimp) %in% rownames(vimp) ]
    names.list <- lapply(names.list, function(l) l[l %in% colnames(vimp)])
  }

  # Circular Network Diagram Plot

  group.list <- lapply(1:length(group), function(g){
    structure(rep(group[g], length(names.list[[g]])), names =  names.list[[g]])
  }
  )


  circos <- circlize::chordDiagram(vimp, annotationTrack = "grid",
                         group = unlist(group.list),
                         preAllocateTracks = list(
                           track.height = circlize::mm_h(4),
                           track.margin = c(circlize::mm_h(4), 0),
                           ...)
  )


  circlize::circos.track(track.index = 2, panel.fun = function(x, y) {
    sector.index.all = circlize::get.cell.meta.data("sector.index")
    for(i in c(1:length(names.list))){

      sector.index = sector.index.all[sector.index.all %in% names.list[[i]]]
      xlim = circlize::get.cell.meta.data("xlim")
      ylim = circlize::get.cell.meta.data("ylim")

      circlize::circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.4 , facing = "clockwise", niceFacing = TRUE)

    }}, bg.border = NA)



  for(i in c(1:length(names.list))){
    circlize::highlight.sector(names.list[[i]], track.index = 1, col = i,
                     text = group[i], cex = 0.8, text.col = "white", niceFacing = TRUE)
  }

  if(!is.null(highlight)){
    highlight <- unique(intersect(highlight, unlist(names.list, use.names = FALSE)))
    if (length(highlight) > 0L) {
      circlize::highlight.sector(highlight, col = "#00FF0040")
    }
  }


  circlize::circos.clear()


}



#' @export
#' @rdname plot_tsne
plot_embed <- function(dat, group = NULL, position = "bottom", pch = 20, ...){


  if(!is.null(group)){
    group <- as.character(group)
    l <- sort(na.omit(unique(group)))
    get_p <- ggpubr::get_palette("jco", k = length(l))
    names(get_p) <- l
    col <- get_p[group]
  } else {col = 1}

  pairs(dat, col = col, pch = pch,  oma=c(10,4,6,4), ...)
  par(xpd=TRUE)
  if(!is.null(group)){
    legend(position, fill = get_p, legend = l, horiz=TRUE,
           xpd=TRUE, bty="n", cex = .5)
  }


}


#' @export
#' @rdname plot_tsne
plot_umap <- function(dat, group = NULL, main = "UMAP", label_group = TRUE,
                      pca = TRUE, ncomp = 70, position = "right", pch = 20, config = umap::umap.defaults,  method = "umap-learn", source = "auto", omics = NULL, cluster = NULL, ...){

  dat <- extract_plot_matrix(dat, source = source, omics = omics, cluster = cluster)

  if(pca) {
    x <- prcomp(dat)$x[, seq_len(min(ncomp, ncol(dat))), drop = FALSE]
  } else {
    x <- dat
  }

  t <- umap::umap(x, config = config, method = method)

  if(is.null(group)) group <- rep("black", nrow(dat))
  group <- as.factor(group)

  df <- data.frame(umap1 = t$layout[,1], umap2 = t$layout[,2], group = group)
  df <- na.omit(df)

  p1 <- ggplot(df, aes(umap1, umap2)) +
    geom_point(aes(color = group), size = .5) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    ggtitle(main) +
    theme_classic()

  if(length(unique(group)) == 1) {
    p1 <- p1 + guides(color = "none")
  } else {
    p1 <- p1 +
      ggsci::scale_color_igv() +
      guides(color = guide_legend(override.aes = list(size=3))) +
      theme(legend.position = position,
            legend.title = element_blank(),
            legend.text = element_text(size = 6))
  }

  if (label_group && length(unique(group)) != 1) {
    data2 <- df %>%
      group_by(group) %>%
      dplyr::select(umap1, umap2) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), mean))

    p1 <- p1 +
      ggrepel::geom_text_repel(data = data2, aes(label = group),
                               size = 2.5, color = "grey20",
                               max.overlaps = 15)
  }


  p1



}
