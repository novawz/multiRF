# Fit an initial MRF model

Fit an initial MRF model

## Usage

``` r
mrf3_init(
  dat.list,
  ntree = 300,
  scale = TRUE,
  ytry = NULL,
  samptype = c("swor", "swr"),
  connect_list = NULL,
  filter_mode = c("auto", "none", "manual"),
  filter_method = c("mad", "variance"),
  top_n_by_type = NULL,
  top_n_manual = NULL,
  filter_verbose = TRUE,
  return_data = FALSE,
  sub_mrf = FALSE,
  sub_mrf_args = list(),
  select_connection = FALSE,
  compute_oob = FALSE,
  verbose = TRUE,
  seed = 529,
  ...
)
```

## Arguments

- dat.list:

  A list containing multi-omics datasets with samples in rows and
  features in columns.

- ntree:

  Number of trees for fitting MRF model. Default is 300.

- scale:

  Whether to z-standardize each feature. Default is TRUE.

- ytry:

  Number of response variables sampled at each split. `NULL` delegates
  to the engine default.

- samptype:

  Sampling scheme passed to forest fitting: `"swor"` or `"swr"`.

- connect_list:

  Optional pre-defined connection list. If `NULL`, all directed pairwise
  models are retained. Their modularity scores are used as soft fusion
  weights; set `select_connection = TRUE` only to request the legacy
  direction-selection step.

- filter_mode:

  Feature filtering mode passed to
  [`filter_omics()`](https://novawz.github.io/multiRF/reference/filter_omics.md):
  `"auto"`, `"none"`, or `"manual"`.

- filter_method:

  Feature dispersion metric passed to
  [`filter_omics()`](https://novawz.github.io/multiRF/reference/filter_omics.md):
  `"mad"` or `"variance"`.

- top_n_by_type:

  Optional auto-filter overrides passed to
  [`filter_omics()`](https://novawz.github.io/multiRF/reference/filter_omics.md).

- top_n_manual:

  Optional manual top-n configuration passed to
  [`filter_omics()`](https://novawz.github.io/multiRF/reference/filter_omics.md).

- filter_verbose:

  Logical; whether to print filtering diagnostics.

- return_data:

  Whether to return the data list. Default is FALSE.

- sub_mrf:

  Logical; if `TRUE`, uses
  [`fit_sub_multi_rfsrc()`](https://novawz.github.io/multiRF/reference/fit_sub_multi_rfsrc.md)
  instead of
  [`fit_multi_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
  for the final model fit. This sub-samples response and predictor
  features per connection, fits smaller MRFs, and averages the resulting
  n x n matrices. Useful when p is large. Default is `FALSE`.

- sub_mrf_args:

  Named list of arguments forwarded to
  [`fit_sub_multi_rfsrc()`](https://novawz.github.io/multiRF/reference/fit_sub_multi_rfsrc.md)
  when `sub_mrf = TRUE`. Common options: `n_sub` (number of replicates,
  default 15), `frac_response` (default 0.2), `frac_predictor` (default
  0.2), `ntree_per_sub` (by default `ceiling(ntree / n_sub)`),
  `min_response_for_sub` (default 500; blocks smaller than this use all
  features), `min_predictor_for_sub` (default 500), and `ytry` (the
  main-model value). Unless explicitly overridden, sub-MRF fits also
  inherit the main forest's `samptype`, `mtry`, engine, split count,
  thread count, node size, and depth.

- select_connection:

  Logical; whether to apply the legacy quality-based
  connection-selection step. The default is `FALSE`, retaining all
  connections for modularity-weighted fusion.

- compute_oob:

  Logical; whether to compute OOB normalized MSE during connection
  selection. When `TRUE`, quality score = modularity - oob_nmse; when
  `FALSE` (default), quality score = modularity only, which is faster.

- verbose:

  Logical; whether to print progress messages.

- seed:

  Random seed.

- ...:

  Additional arguments passed to forest fitting helpers.

## Value

mrf3 object

## Details

`mrf3_init()` now performs initialization only (filtering, forest
fitting, and optional connection selection). IMD weights are computed
downstream via
[`get_multi_weights()`](https://novawz.github.io/multiRF/reference/get_multi_weights.md)
(for example inside `mrf3_fit(..., run_imd = TRUE)`).

When `sub_mrf = TRUE`, the final model fit uses a sub-sampling ensemble
strategy: for each connection, features are randomly sub-sampled from
both response and predictor blocks, smaller MRFs are fitted, and the n x
n forest-weight and proximity matrices are averaged. This trades a small
amount of matrix approximation accuracy for substantial speed gains when
p is large (e.g., 6-8x faster at p = 2000). When `connect_list = NULL`,
all directed connections are retained by default and modularity is used
only as a fusion weight. Legacy direction selection is opt-in through
`select_connection = TRUE`.
