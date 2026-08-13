#' Plot cluster composition as a heatmap
#'
#' Build a composition heatmap from paired cluster labels and an external
#' annotation. The default shows the distribution of annotation classes
#' within each cluster. Sample labels are paired by position; align them by
#' sample identifier before calling this function when necessary.
#'
#' @param cluster Cluster labels. A factor preserves its declared level order;
#'   other atomic vectors use first-occurrence order.
#' @param annotation Annotation or reference labels paired positionally with
#'   `cluster`. Factor levels control the displayed column order.
#' @param normalize Normalization direction. `"cluster"` gives row fractions,
#'   `"annotation"` gives column fractions, `"total"` gives fractions of all
#'   retained samples, and `"none"` displays raw counts.
#' @param label Cell-label style: sample count and percentage, count only,
#'   percentage only, or no label. Percentage labels are unavailable when
#'   `normalize = "none"`.
#' @param na.rm Logical; if `TRUE`, remove pairs with a missing cluster or
#'   annotation. If `FALSE`, missing values produce an error.
#' @param drop Logical; whether to drop unused factor levels. With
#'   `drop = FALSE`, zero-count combinations remain visible.
#' @param percent_digits Integer from 0 to 6 giving the number of decimal
#'   places used in percentage labels.
#' @param main,xlab,ylab Plot title and axis labels. Use `NULL` to suppress a
#'   label.
#' @param colours Character vector of at least two colours defining the
#'   sequential fill gradient.
#' @param base_size Base font size.
#' @param base_family Font family.
#'
#' @return A `ggplot` object. Its `data` component contains `cluster`,
#'   `annotation`, `count`, `fraction`, `value`, and `label`, which can be used
#'   as figure source data. `fraction` is `NA` when `normalize = "none"`.
#' @export
#'
#' @examples
#' cluster <- factor(c("C1", "C1", "C1", "C2", "C2"))
#' subtype <- factor(c("A", "A", "B", "B", "C"), levels = c("A", "B", "C"))
#' plot_cluster_composition(cluster, subtype)
plot_cluster_composition <- function(
    cluster,
    annotation,
    normalize = c("cluster", "annotation", "total", "none"),
    label = c("count_percent", "count", "percent", "none"),
    na.rm = TRUE,
    drop = TRUE,
    percent_digits = 0L,
    main = "Cluster composition",
    xlab = NULL,
    ylab = "Cluster",
    colours = c("#F7FBFF", "#B4C0E4", "#484878"),
    base_size = 8,
    base_family = "sans") {
  normalize <- match.arg(normalize)
  label <- match.arg(label)

  validate_composition_vector(cluster, "cluster")
  validate_composition_vector(annotation, "annotation")
  if (length(cluster) != length(annotation)) {
    stop("`cluster` and `annotation` must have the same length.", call. = FALSE)
  }
  if (length(cluster) == 0L) {
    stop("`cluster` and `annotation` must not be empty.", call. = FALSE)
  }
  if (length(na.rm) != 1L || !is.logical(na.rm) || is.na(na.rm)) {
    stop("`na.rm` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(drop) != 1L || !is.logical(drop) || is.na(drop)) {
    stop("`drop` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(percent_digits) != 1L || !is.numeric(percent_digits) ||
      !is.finite(percent_digits) || percent_digits < 0 ||
      percent_digits > 6 || percent_digits != as.integer(percent_digits)) {
    stop("`percent_digits` must be an integer from 0 to 6.", call. = FALSE)
  }
  percent_digits <- as.integer(percent_digits)
  validate_composition_text(main, "main", allow_null = TRUE)
  validate_composition_text(xlab, "xlab", allow_null = TRUE)
  validate_composition_text(ylab, "ylab", allow_null = TRUE)
  composition_theme <- theme_mrf3(
    base_size = base_size,
    base_family = base_family
  )
  if (!is.character(colours) || length(colours) < 2L ||
      anyNA(colours) || any(!nzchar(colours))) {
    stop("`colours` must contain at least two valid colour strings.", call. = FALSE)
  }
  colour_ok <- tryCatch({
    grDevices::col2rgb(colours)
    TRUE
  }, error = function(e) FALSE)
  if (!colour_ok) {
    stop("`colours` contains an invalid colour.", call. = FALSE)
  }
  if (normalize == "none" && label %in% c("count_percent", "percent")) {
    stop(
      "Percentage labels require a normalized fill; use `label = \"count\"` when `normalize = \"none\"`.",
      call. = FALSE
    )
  }

  complete <- !is.na(cluster) & !is.na(annotation)
  if (!isTRUE(na.rm) && any(!complete)) {
    stop(
      "Missing labels are present; use `na.rm = TRUE` or recode them as an explicit category.",
      call. = FALSE
    )
  }
  cluster_use <- cluster[complete]
  annotation_use <- annotation[complete]
  if (length(cluster_use) == 0L) {
    stop("No complete cluster-annotation pairs remain to plot.", call. = FALSE)
  }

  cluster_levels <- composition_levels(cluster, cluster_use, drop = drop)
  annotation_levels <- composition_levels(annotation, annotation_use, drop = drop)
  cluster_chr <- as.character(cluster_use)
  annotation_chr <- as.character(annotation_use)
  if (any(!nzchar(trimws(cluster_levels))) ||
      any(!nzchar(trimws(annotation_levels)))) {
    stop("Cluster and annotation labels must not be blank.", call. = FALSE)
  }

  counts <- as.data.frame(
    table(
      cluster = factor(cluster_chr, levels = cluster_levels),
      annotation = factor(annotation_chr, levels = annotation_levels)
    ),
    stringsAsFactors = TRUE,
    responseName = "count"
  )
  counts$count <- as.integer(counts$count)

  denominator <- switch(
    normalize,
    cluster = ave(counts$count, counts$cluster, FUN = sum),
    annotation = ave(counts$count, counts$annotation, FUN = sum),
    total = rep(sum(counts$count), nrow(counts)),
    none = rep(NA_real_, nrow(counts))
  )
  counts$fraction <- ifelse(denominator > 0, counts$count / denominator, 0)
  counts$value <- if (normalize == "none") counts$count else counts$fraction

  percent_text <- formatC(
    100 * counts$fraction,
    digits = percent_digits,
    format = "f"
  )
  counts$label <- switch(
    label,
    count_percent = sprintf("%d\n%s%%", counts$count, percent_text),
    count = as.character(counts$count),
    percent = paste0(percent_text, "%"),
    none = rep("", nrow(counts))
  )
  counts$label[counts$count == 0L] <- ""

  contrast_value <- counts$value
  if (normalize == "none" && max(contrast_value) > 0) {
    contrast_value <- contrast_value / max(contrast_value)
  }
  counts$text_on_dark <- contrast_value > 0.45

  legend_title <- switch(
    normalize,
    cluster = "Within-cluster\nfraction",
    annotation = "Within-annotation\nfraction",
    total = "Overall fraction",
    none = "Sample count"
  )
  fill_limits <- if (normalize == "none") NULL else c(0, 1)

  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = annotation, y = cluster, fill = value)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, colour = text_on_dark),
      size = max(base_size / 3.4, 1.8),
      lineheight = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      colours = colours,
      limits = fill_limits,
      name = legend_title
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "#272727", "TRUE" = "white")
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = main, x = xlab, y = ylab) +
    composition_theme +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
    )
}


validate_composition_vector <- function(x, name) {
  if (!is.atomic(x) || !is.null(dim(x))) {
    stop(sprintf("`%s` must be an atomic vector or factor.", name), call. = FALSE)
  }
  if (is.numeric(x) && any(!is.finite(x) & !is.na(x))) {
    stop(sprintf("`%s` must not contain infinite values.", name), call. = FALSE)
  }
  invisible(TRUE)
}


validate_composition_text <- function(x, name, allow_null = FALSE) {
  if (isTRUE(allow_null) && is.null(x)) {
    return(invisible(TRUE))
  }
  if (length(x) != 1L || !is.character(x) || is.na(x)) {
    stop(sprintf("`%s` must be a single character string or NULL.", name), call. = FALSE)
  }
  invisible(TRUE)
}


composition_levels <- function(original, retained, drop) {
  retained_chr <- as.character(retained)
  if (is.factor(original)) {
    levels_out <- levels(original)
    levels_out <- levels_out[!is.na(levels_out)]
    if (isTRUE(drop)) {
      levels_out <- levels_out[levels_out %in% retained_chr]
    }
    return(levels_out)
  }
  unique(retained_chr)
}
