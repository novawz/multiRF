test_that("plot theme and palette are stable and composable", {
  pal <- multiRF:::mrf3_plot_palette()
  expect_identical(
    unname(pal[c("signal_blue", "signal_teal", "accent_orange", "accent_red")]),
    c("#3182BD", "#33B5A5", "#E28E2C", "#D24B40")
  )
  expect_s3_class(multiRF:::theme_mrf3(), "theme")
})

test_that("plot_cluster_composition returns auditable row composition", {
  cluster <- factor(c("C1", "C1", "C1", "C2", "C2"), levels = c("C1", "C2", "C3"))
  annotation <- factor(
    c("A", "A", "B", "B", "C"),
    levels = c("A", "B", "C", "D")
  )
  p <- plot_cluster_composition(cluster, annotation, drop = FALSE)
  expect_s3_class(p, "ggplot")
  expect_identical(levels(p$data$cluster), c("C1", "C2", "C3"))
  expect_identical(levels(p$data$annotation), c("A", "B", "C", "D"))
  expect_equal(nrow(p$data), 12L)
  expect_equal(
    unname(as.numeric(tapply(p$data$fraction, p$data$cluster, sum))),
    c(1, 1, 0)
  )
  cell <- p$data[p$data$cluster == "C1" & p$data$annotation == "A", ]
  expect_identical(cell$count, 2L)
  expect_equal(cell$fraction, 2 / 3)
  expect_identical(cell$label, "2\n67%")
  expect_true(all(p$data$label[p$data$count == 0L] == ""))
  expect_equal(p$scales$get_scales("fill")$limits, c(0, 1))
})

test_that("plot_cluster_composition supports normalization choices", {
  cluster <- c("C1", "C1", "C2", "C2", "C2")
  annotation <- c("A", "B", "A", "A", "B")

  p_annotation <- plot_cluster_composition(
    cluster, annotation, normalize = "annotation", label = "percent"
  )
  expect_equal(
    unname(as.numeric(tapply(
      p_annotation$data$fraction, p_annotation$data$annotation, sum
    ))),
    c(1, 1)
  )

  p_total <- plot_cluster_composition(cluster, annotation, normalize = "total")
  expect_equal(sum(p_total$data$fraction), 1)

  p_count <- plot_cluster_composition(
    cluster, annotation, normalize = "none", label = "count"
  )
  expect_identical(p_count$data$value, p_count$data$count)
  expect_true(all(is.na(p_count$data$fraction)))
  expect_null(p_count$scales$get_scales("fill")$limits)
  expect_error(
    plot_cluster_composition(cluster, annotation, normalize = "none"),
    "Percentage labels require"
  )
})

test_that("plot_cluster_composition validates paired labels", {
  p <- plot_cluster_composition(
    factor(c("C1", NA, "C2"), levels = c("C1", "C2")),
    factor(c("A", "B", "A"), levels = c("A", "B"))
  )
  expect_equal(sum(p$data$count), 2L)
  expect_error(
    plot_cluster_composition(c("C1", NA), c("A", "B"), na.rm = FALSE),
    "Missing labels"
  )
  expect_error(plot_cluster_composition(c("C1"), c("A", "B")), "same length")
  expect_error(plot_cluster_composition(character(), character()), "must not be empty")
  expect_error(plot_cluster_composition(matrix("C1"), "A"), "atomic vector")
  expect_error(plot_cluster_composition("C1", "A", colours = "red"), "at least two")
  expect_error(plot_cluster_composition("", "A"), "must not be blank")
  expect_error(plot_cluster_composition("C1", "A", base_size = "large"), "base_size")
  expect_error(plot_cluster_composition("C1", "A", percent_digits = 7), "0 to 6")
})

test_that("plot_weights validates inputs and returns one faceted ggplot", {
  weights <- list(
    gene = stats::setNames(c(0.8, 0.4, 0.2), c("G1", "G2", "G3")),
    methy = c(0.7, 0.1)
  )
  p <- plot_weights(weights, plot.which = c("gene", "methy"), top = 2)
  expect_s3_class(p, "ggplot")
  expect_identical(sort(unique(p$data$block)), c("gene", "methy"))
  expect_true(all(c("V1", "V2") %in% p$data$variable))
  built <- ggplot2::ggplot_build(p)
  expect_length(built$layout$panel_scales_x, 2L)
  expect_error(plot_weights(weights, plot.which = "missing"), "Unknown `plot.which`")
  expect_error(plot_weights(weights, top = 0), "positive integer")
  bad <- weights
  bad$gene[1] <- NA_real_
  expect_error(plot_weights(bad), "finite")
})

test_that("mrf3 block plotting source uses per-response weights", {
  W_all <- diag(3)
  W_gene <- matrix(1 / 3, 3, 3)
  fit <- structure(
    list(
      models = list(dummy = list()),
      reconstruction = list(
        W = list(W_all = W_all, W_per_response = list(gene = W_gene)),
        fused_mat = list(gene = matrix(1, 3, 2))
      ),
      shared = list(clustering = list(similarity = W_all))
    ),
    class = "mrf3_fit"
  )
  expect_identical(
    multiRF:::extract_plot_matrix(fit, source = "reconstruction_weight_block", omics = "gene"),
    W_gene
  )
})

test_that("embedding plots use the shared visual contract", {
  S <- matrix(c(
    0, 0.9, 0.2, 0.1,
    0.9, 0, 0.1, 0.2,
    0.2, 0.1, 0, 0.8,
    0.1, 0.2, 0.8, 0
  ), 4, 4, byrow = TRUE)
  p <- plot_tsne(
    mod = list(dat = S), group = c("A", "A", "B", "B"),
    seed = 17, max_iter = 10, learning_rate = 10
  )
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "t-SNE 1")
  expect_equal(nrow(p$data), 4L)
  expect_error(
    plot_tsne(mod = list(dat = matrix(1, 4, 2)), max_iter = 2),
    "square matrix"
  )

  skip_if_not_installed("Rtsne")
  set.seed(91)
  raw <- matrix(stats::rnorm(100), nrow = 20)
  rng_before <- .Random.seed
  expect_warning(
    raw_plot <- plot_tsne(raw, seed = 17, max_iter = 250),
    NA
  )
  expect_s3_class(raw_plot, "ggplot")
  expect_identical(.Random.seed, rng_before)
  expect_error(plot_tsne(raw, seed = 1.5), "seed")
})

test_that("plot_umap defaults to the native R backend", {
  skip_if_not_installed("umap")
  set.seed(4)
  x <- matrix(stats::rnorm(60), nrow = 12, ncol = 5)
  rng_before <- .Random.seed
  p <- plot_umap(x, group = rep(c("A", "B"), each = 6), pca = TRUE,
                 seed = 17)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "UMAP 1")
  expect_identical(.Random.seed, rng_before)
  expect_error(plot_umap(x, group = 1:3), "does not match")
})

test_that("plot_network validates groups and returns an igraph", {
  S <- matrix(c(
    1, 0.8, 0.1, 0.1,
    0.8, 1, 0.1, 0.1,
    0.1, 0.1, 1, 0.7,
    0.1, 0.1, 0.7, 1
  ), 4, 4, byrow = TRUE)
  dimnames(S) <- list(paste0("S", 1:4), paste0("S", 1:4))
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  set.seed(81)
  rng_before <- .Random.seed
  graph <- plot_network(S, group = c("A", "A", "B", "B"), cutoff = 0.5,
                        seed = 17)
  expect_s3_class(graph, "igraph")
  expect_equal(igraph::vcount(graph), 4L)
  expect_identical(.Random.seed, rng_before)
  expect_error(plot_network(S, group = c("A", "B"), cutoff = 0.5), "one value per")
  expect_warning(empty <- plot_network(diag(4), cutoff = 1), "No edges remain")
  expect_equal(igraph::vcount(empty), 0L)
})

test_that("plot_embed is deprecated and restores graphics state", {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  old <- graphics::par(c("xpd", "oma"))
  expect_warning(
    plot_embed(matrix(stats::rnorm(30), 10, 3), group = rep(c("A", "B"), each = 5)),
    "deprecated"
  )
  expect_equal(graphics::par(c("xpd", "oma")), old)
})

test_that("plot_circos clears state and rejects empty links", {
  skip_if_not_installed("circlize")
  mat <- matrix(c(0, 0.8, 0.8, 0), 2, 2, byrow = TRUE)
  dimnames(mat) <- list(c("A1", "B1"), c("A1", "B1"))
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  expect_no_error(plot_circos(mat, names.list = list(A = "A1", B = "B1")))
  expect_length(circlize::get.all.sector.index(), 0L)
  expect_error(
    plot_circos(mat, names.list = list(A = "A1", B = "B1"), cut.off = 1),
    "No non-zero links"
  )
  expect_length(circlize::get.all.sector.index(), 0L)

  sparse <- matrix(0, 3, 3, dimnames = list(c("A1", "B1", "C1"), c("A1", "B1", "C1")))
  sparse[1, 2] <- sparse[2, 1] <- 0.8
  expect_no_error(
    plot_circos(
      sparse,
      names.list = list(A = "A1", B = "B1", C = "C1")
    )
  )
  expect_length(circlize::get.all.sector.index(), 0L)
})

test_that("plot_km validates survival inputs and uses restrained colours", {
  skip_if_not_installed("survival")
  skip_if_not_installed("survminer")
  pheno <- data.frame(time = seq_len(12), event = rep(c(0, 1), 6))
  p <- suppressWarnings(
    plot_km(rep(c("A", "B"), each = 6), "time", "event", pheno,
            risk.table = FALSE)
  )
  expect_s3_class(p, "ggsurvplot")
  expect_error(plot_km(1:3, "time", "event", pheno), "one value per row")
  bad <- pheno
  bad$event[1] <- 2
  expect_error(plot_km(rep(c("A", "B"), each = 6), "time", "event", bad), "0/1 events")
})
