# Fit a multivariate regression forest

Drop-in replacement for
[`fit_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
when `type = "regression"`. Returns an object with the same interface
(`$forest.wt`, `$proximity`, `$membership`, `$xvar`, `$yvar`, `$ntree`,
`$xvar.names`) so downstream code in multiRF works without changes.

## Usage

``` r
fit_mv_forest(
  X,
  Y,
  ntree = 500L,
  mtry = NULL,
  ytry = NULL,
  nsplit = 10L,
  forest.wt = "all",
  proximity = c("all", "inbag", "oob", "none"),
  nodesize = 5L,
  max_depth = 0L,
  seed = 529L,
  samptype = c("swor", "swr"),
  nthread = getOption("multiRF.nthread", 0L),
  xvar.wt = NULL,
  yvar.wt = NULL,
  enhanced_prox = FALSE,
  sibling_gamma = 0.5,
  leaf_embed_dim = 10L
)
```

## Arguments

- X:

  Data frame or matrix of predictor features (n x px).

- Y:

  Data frame or matrix of response features (n x qy).

- ntree:

  Number of trees.

- mtry:

  Number of candidate X variables per split. Accepts an integer, a
  formula string (`"sqrt(p)"`, `"p/3"`, `"p/2"`), or `NULL` for the
  default (`ceiling(px/3)` for regression, `ceiling(sqrt(px))` for
  classification). In formulas, `p` is the number of predictor columns.

- ytry:

  Number of candidate Y variables per split. Accepts an integer, a
  formula string (`"sqrt(p)"`, `"p/3"`), or `NULL` for the default
  (`ceiling(qy/3)`). In formulas, `p` is the number of Y columns.

- nsplit:

  Number of candidate numeric cutpoints per split variable.

- forest.wt:

  Forest-weight mode: `"all"`, `"inbag"`, or `"oob"`.

- proximity:

  Proximity output mode.

- nodesize:

  Minimum terminal node size. Default: 5.

- max_depth:

  Maximum tree depth (0 = unlimited). Default: 0.

- seed:

  Random seed.

- samptype:

  Sampling scheme: `"swor"` or `"swr"`.

- nthread:

  Number of threads used for fitting.

- xvar.wt:

  Optional non-negative predictor sampling weights.

- yvar.wt:

  Optional non-negative response sampling weights.

- enhanced_prox:

  Logical; whether to compute enhanced proximity.

- sibling_gamma:

  Strength of the sibling-leaf correction used by enhanced proximity.

- leaf_embed_dim:

  Embedding dimension used by enhanced proximity.

## Value

A list compatible with rfsrc output, containing:

- forest.wt:

  n x n forest weight matrix

- proximity:

  n x n proximity matrix

- membership:

  n x ntree terminal node membership

- xvar:

  predictor data frame

- yvar:

  response data frame

- xvar.names:

  character vector of predictor names

- ntree:

  number of trees

- tree_info:

  per-tree structure for IMD
