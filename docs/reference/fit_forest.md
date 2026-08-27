# Fit random forest model and create model lists

Fit random forest model and create model lists

## Usage

``` r
fit_forest(
  X,
  Y = NULL,
  type = c("regression", "classification", "unsupervised"),
  nodedepth = NULL,
  max_depth = NULL,
  nodesize = NULL,
  ntree = 200,
  forest.wt = "all",
  proximity = "all",
  mtry = NULL,
  ytry = NULL,
  nsplit = 10,
  samptype = c("swor", "swr"),
  nthread = getOption("multiRF.nthread", 0L),
  xvar.wt = NULL,
  yvar.wt = NULL,
  split.wt = NULL,
  case.wt = NULL,
  seed = 529L,
  engine = getOption("multiRF.engine", "native"),
  enhanced_prox = FALSE,
  sibling_gamma = 0.5,
  leaf_embed_dim = 10L,
  ...
)

fit_multi_forest(
  dat.list,
  connect_list = NULL,
  var.wt = NULL,
  yprob = 1,
  ytry = NULL,
  seed = 529L,
  parallel_connections = FALSE,
  cores_connections = NULL,
  ...
)

fit_rfsrc(
  X,
  Y = NULL,
  type = "regression",
  nodedepth = NULL,
  max_depth = NULL,
  nodesize = NULL,
  ntree = 200,
  forest.wt = "all",
  proximity = "all",
  mtry = NULL,
  ytry = NULL,
  nsplit = 10,
  samptype = c("swor", "swr"),
  nthread = getOption("multiRF.nthread", 0L),
  xvar.wt = NULL,
  yvar.wt = NULL,
  split.wt = NULL,
  case.wt = NULL,
  seed = 529L,
  engine = getOption("multiRF.engine", "native"),
  enhanced_prox = FALSE,
  sibling_gamma = 0.5,
  leaf_embed_dim = 10L,
  ...
)

fit_multi_rfsrc(
  dat.list,
  connect_list = NULL,
  var.wt = NULL,
  yprob = 1,
  ytry = NULL,
  seed = 529L,
  ...
)
```

## Arguments

- X:

  A data frame that consider as predictor set

- Y:

  A data frame that consider as response set. If Y = NULL, an
  unsupervised RF is conducted

- type:

  Select the type of RF model. The default is regression. Can select
  from "regression", "classification", and "unsupervised"

- nodedepth:

  Backward-compatible alias for `max_depth`.

- max_depth:

  Maximum tree depth; `0` means unlimited.

- nodesize:

  Minimum terminal-node size. Defaults to 5 for regression, 1 for
  classification, and 3 for unsupervised forests.

- ntree:

  Number of trees.

- forest.wt:

  Forest-weight output mode: `"all"` gives all-sample terminal-node
  weights, `"inbag"` restricts donors to bootstrap samples, and `"oob"`
  evaluates each target only on trees where it is OOB.

- proximity:

  Proximity output mode: `"all"`, `"inbag"`, `"oob"`, or `"none"`.

- mtry:

  Number of candidate X variables per split. Default `NULL` =
  `ceiling(px/3)` for multivariate regression.

- ytry:

  Number of response variables randomly selected per split. Default
  `NULL` means the native engine uses `ceiling(qy/3)`. Set to a specific
  integer to override (e.g., `ytry = ncol(Y) / 2`).

- nsplit:

  Number of candidate numeric cutpoints evaluated per variable. Native
  and `randomForestSRC` both default to `10`; set `0` to scan all.

- samptype:

  Sampling scheme: `"swor"` (without replacement) or `"swr"` (with
  replacement).

- nthread:

  Number of threads used by the native engine.

- xvar.wt, yvar.wt:

  Optional non-negative predictor/response sampling weights used by the
  native multivariate forest.

- split.wt, case.wt:

  Optional split/case weights. These are forwarded to `randomForestSRC`;
  the native engine fails explicitly because it does not yet implement
  them.

- seed:

  Random seed passed to the selected engine.

- engine:

  Forest backend. Default is `getOption("multiRF.engine", "native")`.
  Native is the default and recommended engine. `randomForestSRC` is
  used only as a non-native fallback when explicitly requested.

- enhanced_prox:

  Logical; whether to compute enhanced proximity in the native engine.

- sibling_gamma:

  Strength of the sibling-leaf correction used by enhanced proximity.

- leaf_embed_dim:

  Embedding dimension used by the native enhanced proximity path.

- ...:

  Additional arguments passed to `randomForestSRC::rfsrc()` when
  `engine != "native"`.

- dat.list:

  A list that contains multi-omics datasets with samples in rows and
  features in columns. Samples should be matched in each dataset.

- connect_list:

  A pre-defined connection list between datasets. If `NULL`, all
  directed pairwise connections are enumerated.

- var.wt:

  Optional named list of variable weights; names must cover the block
  names in `names(dat.list)`.

- yprob:

  Deprecated. Use `ytry` directly instead.

- parallel_connections:

  Logical; whether directed connections are fitted in parallel.

- cores_connections:

  Optional core budget for connection-level parallelism.

## Value

A model list

## Details

`fit_forest()` now defaults to the package-native engine for
classification, multivariate regression, and unsupervised fitting.
`randomForestSRC` is optional and is only used when
`engine != "native"`. If `type` is omitted and `Y = NULL`, unsupervised
fitting is selected.
