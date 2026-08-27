# Tests for the unified S3 plot entry point: plot(fit, type = "...")

dispatch_fit_env <- new.env(parent = emptyenv())

get_dispatch_tcga <- function() {
  tcga <- tryCatch(multiRF::tcga_brca, error = function(e) NULL)
  if (is.null(tcga)) {
    data_env <- new.env()
    utils::data("tcga_brca_data", package = "multiRF", envir = data_env)
    tcga <- data_env$tcga_brca
  }
  tcga
}

get_dispatch_fit <- function() {
  if (!is.null(dispatch_fit_env$fit)) {
    return(dispatch_fit_env$fit)
  }
  tcga <- get_dispatch_tcga()
  ids <- rownames(tcga[[1L]])[seq_len(150L)]
  dat <- lapply(tcga, function(x) {
    x[ids, seq_len(min(60L, ncol(x))), drop = FALSE]
  })
  fit <- suppressWarnings(suppressMessages(mrf3_fit(
    dat.list = dat,
    ntree = 30,
    filter_mode = "none",
    clustering_args = list(shared_k = 3L, specific_k = 2L),
    run_imd = TRUE,
    run_variable_selection = TRUE,
    variable_selection_args = list(
      method = "mixture", signal = "shared", level = 0.05, re_fit = FALSE
    ),
    return_data = TRUE,
    nthread = 1,
    filter_verbose = FALSE,
    verbose = FALSE,
    seed = 529
  )))
  dispatch_fit_env$fit <- fit
  fit
}

with_null_pdf <- function(code) {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  force(code)
}


test_that("mrf3_palette is fixed, colour-blind safe, and recycles with warning", {
  pal <- multiRF:::mrf3_palette()
  expect_length(pal, 8L)
  expect_identical(pal[1:4], c("#3182BD", "#33B5A5", "#E28E2C", "#D24B40"))
  ## Same first colours as the package-wide internal plotting palette
  expect_identical(pal, unname(multiRF:::mrf3_plot_palette())[1:8])
  expect_identical(multiRF:::mrf3_palette(3), pal[1:3])
  expect_identical(multiRF:::mrf3_palette(NULL), pal)
  expect_length(multiRF:::mrf3_palette(0), 0L)
  expect_warning(long <- multiRF:::mrf3_palette(10), "recycled")
  expect_identical(long, rep_len(pal, 10L))
})

test_that("plot.mrf3_fit is registered and validates `type`", {
  fit <- get_dispatch_fit()
  expect_true(inherits(fit, "mrf3_fit"))
  expect_error(plot(fit, type = "nope"), "'arg' should be one of")
})

test_that("type = 'tsne' (the default) dispatches to plot_tsne with clusters", {
  fit <- get_dispatch_fit()
  p <- plot(fit, type = "tsne", seed = 7, max_iter = 60, learning_rate = 50)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "t-SNE 1")
  expect_equal(nrow(p$data), 150L)
  expect_identical(
    sort(levels(p$data$group)),
    sort(as.character(unique(get_clusters(fit))))
  )

  ## Default type is tsne
  p_default <- plot(fit, seed = 7, max_iter = 30, learning_rate = 50)
  expect_s3_class(p_default, "ggplot")
  expect_identical(p_default$labels$x, "t-SNE 1")
})

test_that("type = 'umap' dispatches to plot_umap with consistent colours", {
  skip_if_not_installed("umap")
  fit <- get_dispatch_fit()
  cl <- get_clusters(fit)
  k <- length(unique(cl))
  p <- plot(fit, type = "umap", seed = 7)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "UMAP 1")
  built <- ggplot2::ggplot_build(p)
  fill_by_cluster <- tapply(
    built$data[[1]]$fill, p$data$group, function(z) unique(z)[1]
  )
  expect_identical(
    as.character(fill_by_cluster[as.character(sort(unique(cl)))]),
    multiRF:::mrf3_palette(k)
  )
})

test_that("type = 'embed' dispatches to plot_embed on similarity PCs", {
  fit <- get_dispatch_fit()
  with_null_pdf({
    expect_warning(plot(fit, type = "embed"), "deprecated")
  })
  expect_error(plot(fit, type = "embed", embed_dims = 1), "embed_dims")
})

test_that("type = 'network' dispatches to plot_network and returns an igraph", {
  fit <- get_dispatch_fit()
  g <- with_null_pdf(plot(fit, type = "network", seed = 7))
  expect_s3_class(g, "igraph")
  expect_gt(igraph::vcount(g), 0L)
  expect_true(all(
    igraph::V(g)$group %in% as.character(unique(get_clusters(fit)))
  ))
})

test_that("the same cluster id has the same colour in tsne and network plots", {
  fit <- get_dispatch_fit()
  cl <- get_clusters(fit)
  k <- length(unique(cl))

  p <- plot(fit, type = "tsne", seed = 7, max_iter = 60, learning_rate = 50)
  built <- ggplot2::ggplot_build(p)
  fill_by_cluster <- tapply(
    built$data[[1]]$fill, p$data$group, function(z) unique(z)[1]
  )
  ## The dispatcher palette itself, keyed by sorted cluster id
  expect_identical(
    as.character(fill_by_cluster[as.character(sort(unique(cl)))]),
    multiRF:::mrf3_palette(k)
  )

  g <- with_null_pdf(plot(fit, type = "network", seed = 7))
  colour_by_cluster <- tapply(
    igraph::V(g)$color, igraph::V(g)$group, function(z) unique(z)[1]
  )
  common <- intersect(names(fill_by_cluster), names(colour_by_cluster))
  expect_gt(length(common), 0L)
  expect_identical(
    as.character(fill_by_cluster[common]),
    as.character(colour_by_cluster[common])
  )
})

test_that("type = 'circos' runs pairwise IMD and draws a chord diagram", {
  skip_if_not_installed("circlize")
  fit <- get_dispatch_fit()
  res <- with_null_pdf(
    suppressWarnings(suppressMessages(plot(fit, type = "circos")))
  )
  expect_s3_class(res, "data.frame")
  expect_length(circlize::get.all.sector.index(), 0L)
})

test_that("type = 'circos' gives an informative error without IMD", {
  fit <- get_dispatch_fit()
  fit_no_imd <- fit
  fit_no_imd$imd <- NULL
  expect_error(plot(fit_no_imd, type = "circos"), "run_imd = TRUE")
})

test_that("type = 'composition' needs an annotation and errors informatively", {
  fit <- get_dispatch_fit()
  expect_error(plot(fit, type = "composition"), "`annotation`")
  expect_error(
    plot(fit, type = "composition", annotation = factor(c("A", "B"))),
    "one annotation per fitted sample"
  )

  annotation <- factor(rep(c("Subtype A", "Subtype B", "Subtype C"),
                           length.out = 150L))
  p <- plot(fit, type = "composition", annotation = annotation)
  expect_s3_class(p, "ggplot")
  expect_identical(
    sort(levels(p$data$cluster)),
    sort(as.character(unique(get_clusters(fit))))
  )
  expect_identical(levels(p$data$annotation), levels(annotation))
})

test_that("type = 'km' needs survival inputs and errors informatively", {
  fit <- get_dispatch_fit()
  expect_error(plot(fit, type = "km"), "`time_var`")
  expect_error(
    plot(fit, type = "km", time_var = "t", event_var = "e"),
    "`pheno_mat`"
  )

  skip_if_not_installed("survival")
  skip_if_not_installed("survminer")
  set.seed(11)
  pheno <- data.frame(
    os_time = stats::rexp(150L, rate = 1 / 50),
    os_event = stats::rbinom(150L, 1L, 0.6)
  )
  expect_error(
    plot(fit, type = "km", time_var = "bad_col", event_var = "os_event",
         pheno_mat = pheno),
    "was not found"
  )
  expect_error(
    plot(fit, type = "km", time_var = "os_time", event_var = "os_event",
         pheno_mat = pheno[1:10, , drop = FALSE]),
    "one phenotype row per"
  )

  p <- suppressWarnings(
    plot(fit, type = "km", time_var = "os_time", event_var = "os_event",
         pheno_mat = pheno, risk.table = FALSE)
  )
  expect_s3_class(p, "ggsurvplot")

  ## Default curve colours come from the shared cluster palette
  k <- length(unique(get_clusters(fit)))
  curve_cols <- unique(stats::na.omit(
    ggplot2::ggplot_build(p$plot)$data[[1]]$colour
  ))
  expect_true(all(curve_cols %in% multiRF:::mrf3_palette(k)))

  ## A user-supplied palette wins over the default
  custom <- c("#101010", "#202020", "#303030")[seq_len(k)]
  p_custom <- suppressWarnings(
    plot(fit, type = "km", time_var = "os_time", event_var = "os_event",
         pheno_mat = pheno, risk.table = FALSE, palette = custom)
  )
  custom_cols <- unique(stats::na.omit(
    ggplot2::ggplot_build(p_custom$plot)$data[[1]]$colour
  ))
  expect_true(all(custom_cols %in% custom))
})

test_that("type = 'weights' dispatches to plot_weights and checks IMD", {
  fit <- get_dispatch_fit()
  p <- plot(fit, type = "weights", top = 5)
  expect_s3_class(p, "ggplot")
  expect_identical(
    sort(unique(p$data$block)),
    sort(names(get_weights(fit)))
  )

  fit_no_imd <- fit
  fit_no_imd$imd <- NULL
  fit_no_imd$cluster_imd <- NULL
  expect_error(plot(fit_no_imd, type = "weights"), "run_imd = TRUE")
})

test_that("a user-supplied `group` overrides the default cluster labels", {
  fit <- get_dispatch_fit()
  custom_group <- factor(rep(c("X", "Y"), length.out = 150L))
  p <- plot(fit, type = "tsne", group = custom_group,
            seed = 7, max_iter = 30, learning_rate = 50)
  expect_identical(levels(p$data$group), c("X", "Y"))
})
