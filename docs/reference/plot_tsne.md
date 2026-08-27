# Multiple plot functions

Multiple plot functions

## Usage

``` r
plot_tsne(
  dat = NULL,
  mod = NULL,
  group = NULL,
  label_group = TRUE,
  position = "right",
  perplexity = 30,
  pca = FALSE,
  ncomp = 70,
  main = "t-SNE",
  size = 1.4,
  shape = 21,
  source = "auto",
  omics = NULL,
  cluster = NULL,
  seed = NULL,
  base_size = 8,
  base_family = "sans",
  ...
)

plot_network(
  dat,
  group = NULL,
  label = FALSE,
  cutoff = NULL,
  layout = igraph::layout_with_fr,
  vertex.size = 5,
  label.dist = NULL,
  edge.width = NULL,
  vertex.frame.color = "white",
  vertex.label.color = "black",
  vertex.label.cex = 0.5,
  vertex.label.dist = 0,
  edge.curved = 0.5,
  vertex.label.degree = -pi/2,
  position = "bottomright",
  source = "auto",
  omics = NULL,
  cluster = NULL,
  seed = NULL,
  ...
)

plot_weights(
  weights,
  plot.which = "all",
  top = 20,
  labels = NULL,
  weight_source = "auto",
  cluster = NULL,
  ncol = NULL,
  base_size = 8,
  base_family = "sans"
)

plot_circos(
  mat,
  names.list = NULL,
  group = NULL,
  cut.off = NULL,
  highlight = NULL,
  source = c("adj_var", "adj_dat"),
  ...
)

plot_embed(dat, group = NULL, position = "bottom", pch = 20, ...)

plot_umap(
  dat,
  group = NULL,
  main = "UMAP",
  label_group = TRUE,
  pca = TRUE,
  ncomp = 70,
  position = "right",
  pch = 21,
  size = 1.4,
  config = umap::umap.defaults,
  method = "naive",
  source = "auto",
  omics = NULL,
  cluster = NULL,
  seed = NULL,
  base_size = 8,
  base_family = "sans",
  ...
)
```

## Arguments

- dat:

  A data frame or matrix

- mod:

  A fitted model object used to compute embeddings internally. Can also
  be an `mrf3_fit` object.

- group:

  Class group

- label_group:

  Logical; whether to draw group labels on embeddings.

- position:

  Legend position.

- perplexity:

  tSNE perplexity used by
  [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html).

- pca:

  Logical; whether to run PCA before embedding.

- ncomp:

  Number of principal components retained when `pca = TRUE`.

- main:

  a main title for the plot, see also
  [`title`](https://rdrr.io/r/graphics/title.html).

- size:

  Point size in scatter plots.

- shape:

  Point shape in embedding plots.

- source:

  Source used when `dat`/`mod`/`weights` is an `mrf3_fit`. Supported
  values for `mrf3_fit` inputs: `"auto"`, `"specific_shared"`,
  `"specific_specific"`, `"clustering"`, `"robust_clustering"`,
  `"reconstruction_weight_all"`, `"reconstruction_weight_block"`,
  `"reconstruction_fused_block"`. `plot_network()` and `plot_circos()`
  additionally accept a `pairwise_imd_analysis` object with
  `source = "adj_var"` or `"adj_dat"`. The t-SNE model path requires a
  square sample-similarity source; use the `dat` argument for a raw
  rectangular feature matrix. For `plot_weights()`, supported values
  are: `"auto"`, `"imd"`, `"cluster_imd"`, `"mrf"`.

- omics:

  Optional omics block name when `source` is block-specific.

- cluster:

  Optional cluster label when `source` uses cluster-level IMD output.

- seed:

  Optional non-negative integer seed for stochastic embedding or
  network-layout initialization. The caller's RNG state is restored.

- base_size:

  Base font size used by the multiRF plotting theme.

- base_family:

  Font family used by the multiRF plotting theme.

- ...:

  Additional arguments passed to the underlying plotting routine.

- label:

  A logical parameter that determine whether to show labels or not.

- cutoff:

  Edge-weight cutoff. For a regular sample matrix the default is the
  mean positive upper-triangle weight; for pairwise IMD it is the 75th
  percentile of positive upper-triangle weights.

- layout:

  Network layout function for graph plotting.

- vertex.size:

  Vertex size in network plots.

- label.dist:

  Label distance in network plots.

- edge.width:

  Edge width in network plots.

- vertex.frame.color:

  Vertex frame color in network plots.

- vertex.label.color:

  Vertex label color in network plots.

- vertex.label.cex:

  Vertex label font size in network plots.

- vertex.label.dist:

  Vertex label distance in network plots.

- edge.curved:

  Edge curvature in network plots.

- vertex.label.degree:

  Vertex label angle in network plots.

- weights:

  Weight vector or weight-carrying object for `plot_weights()`.

- plot.which:

  Which omics block(s) to plot in `plot_weights()`.

- top:

  The number of top weighted variables to show in the plot. The default
  is 20. Can be chosen from numerical values or all for setting the
  parameter = NULL.

- labels:

  Optional panel labels in `plot_weights()`.

- weight_source:

  Alias of `source` used by `plot_weights()`.

- ncol:

  Number of facet columns in `plot_weights()`. By default, plots with
  more than two blocks use two columns to preserve readable feature
  labels.

- mat:

  Adjacency-like matrix for `plot_circos()`.

- names.list:

  Named list of variable groups for `plot_circos()`.

- cut.off:

  Cutoff applied before plotting a circos graph.

- highlight:

  Vector of sectors to highlight in `plot_circos()`.

- pch:

  Point symbol for base plotting functions.

- config:

  UMAP configuration object passed to
  [`umap::umap()`](https://rdrr.io/pkg/umap/man/umap.html).

- method:

  UMAP backend method passed to
  [`umap::umap()`](https://rdrr.io/pkg/umap/man/umap.html).

## Value

`plot_tsne()`, `plot_umap()`, and `plot_weights()` return a `ggplot`
object. `plot_network()` invisibly returns an `igraph` object,
`plot_circos()` invisibly returns the chord-diagram result, and the
deprecated `plot_embed()` invisibly returns `NULL` after drawing.
