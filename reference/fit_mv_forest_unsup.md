# Fit an unsupervised random forest

Emulates rfsrc unsupervised mode: at each tree, randomly partition the
columns of `X` into pseudo-predictor and pseudo-response halves, then
fit a multivariate regression tree. Only `forest.wt` and `proximity` are
meaningful; there is no real Y.

## Usage

``` r
fit_mv_forest_unsup(
  X,
  ntree = 500L,
  ytry = NULL,
  nsplit = 10L,
  forest.wt = "all",
  proximity = c("all", "inbag", "oob", "none"),
  nodesize = 3L,
  max_depth = 0L,
  seed = 529L,
  samptype = c("swor", "swr"),
  nthread = getOption("multiRF.nthread", 0L),
  enhanced_prox = FALSE,
  sibling_gamma = 0.5,
  leaf_embed_dim = 10L
)
```

## Arguments

- X:

  Data frame or matrix of features (n x p).

- ntree:

  Number of trees.

- ytry:

  Number of candidate pseudo-Y columns per split. Default `NULL` = 15.

- nsplit:

  Number of candidate cutpoints per pseudo-predictor; `0` scans all
  distinct cutpoints.

- forest.wt:

  Forest-weight mode: `"all"` uses every sample in a terminal node,
  `"inbag"` uses bootstrap donors, and `"oob"` uses only trees in which
  the target sample is out of bag.

- proximity:

  Proximity output mode.

- nodesize:

  Minimum terminal node size.

- max_depth:

  Maximum tree depth (0 = unlimited).

- seed:

  Random seed.

- samptype:

  Sampling scheme: `"swor"` or `"swr"`.

- nthread:

  Number of threads used for fitting.

- enhanced_prox:

  Logical; whether to compute enhanced proximity.

- sibling_gamma:

  Strength of the sibling-leaf correction used by enhanced proximity.

- leaf_embed_dim:

  Embedding dimension used by enhanced proximity.

## Value

A list with `forest.wt`, `proximity`, `membership`, `xvar`,
`xvar.names`, `ntree`, and `engine = "multiRF"`.
