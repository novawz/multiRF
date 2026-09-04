# Stage-based multiRF fitting pipeline

Stage-based multiRF fitting pipeline

## Usage

``` r
mrf3_fit(
  dat.list,
  ntree = 500,
  scale = TRUE,
  ytry = NULL,
  samptype = c("swor", "swr"),
  connect_list = NULL,
  filter_mode = c("auto", "none", "manual"),
  filter_method = c("mad", "variance"),
  top_n_by_type = NULL,
  top_n_manual = NULL,
  filter_verbose = TRUE,
  main_clustering = c("similarity", "proximity", "enhanced_proximity"),
  shared_specific_args = list(),
  clustering_args = list(),
  run_imd = FALSE,
  run_cluster_imd = FALSE,
  imd_args = list(),
  run_variable_selection = FALSE,
  variable_selection_args = list(),
  run_robust_clustering = FALSE,
  robust_clustering_args = list(),
  cluster_imd_args = list(),
  top_v = NULL,
  model_top_v = NULL,
  recon_fusion = c("weighted", "uniform"),
  global_fusion = c("average", "pmin"),
  score_power = 1,
  score_floor = 0,
  fallback_uniform = TRUE,
  fused_top_v = NULL,
  fused_row_normalize = TRUE,
  fused_keep_ties = TRUE,
  top_v_method = c("saturation", "entropy_elbow", "neff"),
  top_v_tau = 0.9,
  neff_quantile = 0.5,
  model_top_v_tune_args = list(),
  fused_top_v_tune_args = list(),
  select_connection = FALSE,
  return_data = FALSE,
  compact_output = FALSE,
  verbose = TRUE,
  seed = 529,
  ...
)
```

## Arguments

- dat.list:

  A named list of omics matrices (samples in rows, features in columns).

- ntree:

  Number of trees for RF fitting.

- scale:

  Logical; whether to z-standardize each feature before fitting.

- ytry:

  Number of response variables sampled per split. `NULL` delegates to
  the engine default.

- samptype:

  Sampling scheme passed to forest fitting: `"swor"` or `"swr"`.

- connect_list:

  Optional predefined connections (`list(c(response, predictor), ...)`).
  If `NULL`, all directed pairwise connections are retained and their
  modularity scores are used as soft fusion weights.

- filter_mode:

  Feature filtering mode passed to
  [`filter_omics()`](https://novawz.github.io/multiRF/reference/filter_omics.md).

- filter_method:

  Feature dispersion metric passed to
  [`filter_omics()`](https://novawz.github.io/multiRF/reference/filter_omics.md).

- top_n_by_type:

  Optional auto-filter overrides.

- top_n_manual:

  Optional manual top-n filtering configuration.

- filter_verbose:

  Logical; whether to print filtering diagnostics.

- main_clustering:

  Global clustering strategy applied consistently to shared, specific,
  and robust branches: `"similarity"` (default), `"proximity"`, or
  `"enhanced_proximity"`.

- shared_specific_args:

  A named list of additional arguments passed to
  [`get_shared_specific_weights()`](https://novawz.github.io/multiRF/reference/get_shared_specific_weights.md).
  When this branch is created, per-omics shared fraction
  (`1 - ||R||_F^2 / ||X||_F^2`) is also computed by
  [`get_shared_frac()`](https://novawz.github.io/multiRF/reference/get_shared_frac.md)
  and attached as `shared_frac`. By default, this function passes
  `specific_top_v = selected_fused_top_v` into this branch. Residual
  unsupervised forests inherit the main `ntree`, `samptype`, `ytry`,
  `proximity`, `nsplit`, `nthread`, `nodesize`, and maximum depth
  settings unless a corresponding `specific_*` value is supplied here.

- clustering_args:

  A named list of additional arguments passed to the shared/specific
  clustering stage.

- run_imd:

  Logical; whether to run
  [`get_multi_weights()`](https://novawz.github.io/multiRF/reference/get_multi_weights.md)
  as a pipeline stage.

- run_cluster_imd:

  Logical; whether to run
  [`cluster_imd()`](https://novawz.github.io/multiRF/reference/cluster_imd.md)
  after global IMD when `run_imd = TRUE`. The default is `FALSE`; set it
  explicitly to `TRUE` to request the cluster-specific stage.

- imd_args:

  A named list of additional arguments passed to
  [`get_multi_weights()`](https://novawz.github.io/multiRF/reference/get_multi_weights.md).
  By default, this function sets `parallel = TRUE` for IMD unless
  overridden here. When any model is a sub-MRF, entries shared with
  [`get_imp_forest()`](https://novawz.github.io/multiRF/reference/get_imp_forest.md)
  (e.g. `parallel`, `normalized`) are also applied to the per-connection
  IMD computation of full-forest models.

- run_variable_selection:

  Logical; whether to run variable selection using
  [`mrf3_vs()`](https://novawz.github.io/multiRF/reference/mrf3_vs.md)
  after IMD weights are available.

- variable_selection_args:

  A named list of additional arguments passed to
  [`mrf3_vs()`](https://novawz.github.io/multiRF/reference/mrf3_vs.md).
  The workflow defaults to raw-IMD filtering. An explicit `re_fit` is
  honored unless `run_robust_clustering = TRUE`, which requires and
  therefore forces refitting.

- run_robust_clustering:

  Logical; whether to run a robust clustering branch that first selects
  variables from IMD weights and then re-clusters.

- robust_clustering_args:

  A named list of clustering arguments (`shared_k`, `specific_k`, etc.)
  specific to the robust clustering branch. By default, robust
  clustering inherits the actual `k` chosen in Stage 4; values in this
  list override those inherited defaults.

- cluster_imd_args:

  A named list of additional arguments passed to
  [`cluster_imd()`](https://novawz.github.io/multiRF/reference/cluster_imd.md)
  when `run_imd = TRUE`. By default, this function reuses cluster labels
  from `clusters`.

- top_v:

  Optional unified top-v cutoff applied to both `model_top_v` and
  `fused_top_v`.

- model_top_v:

  Model-level top-v cutoff used on each single-model forest weight
  matrix before fusion. The default `NULL` auto-tunes via
  [`tune_model_top_v()`](https://novawz.github.io/multiRF/reference/tune_model_top_v.md);
  use `Inf` for no truncation. A selected or fixed value at least
  `0.8 * n` is treated as no truncation.

- recon_fusion:

  Reconstruction fusion mode passed to
  [`get_reconstr_matrix()`](https://novawz.github.io/multiRF/reference/mrf3_reconstr.md):
  `"weighted"` (default) or `"uniform"`.

- global_fusion:

  Global fusion across response blocks. `"average"` implements Eq. 8;
  `"pmin"` is an optional intersection extension and disables fused
  top-v tuning.

- score_power:

  Exponent applied to connection scores for weighted reconstruction.

- score_floor:

  Non-negative floor applied to connection scores before weighting.

- fallback_uniform:

  Logical; whether weighted reconstruction falls back to uniform
  averaging when scores are unavailable.

- fused_top_v:

  Row-wise top-v truncation for fused weights. Default `NULL` =
  auto-tune via
  [`tune_fused_top_v()`](https://novawz.github.io/multiRF/reference/tune_fused_top_v.md).
  Set to `Inf` for no truncation or `FALSE` to skip it. Values at least
  `0.8 * n` are treated as no truncation.

- fused_row_normalize:

  Logical; whether to row-normalize fused weights after optional
  truncation.

- fused_keep_ties:

  Logical; whether fused top-v truncation keeps ties at cutoff.

- top_v_method:

  Strategy used when auto-selecting `top_v`: `"saturation"` (default),
  `"entropy_elbow"`, or `"neff"`. `"saturation"` keeps the smallest
  number of neighbours whose fused-weight row entropy reaches
  `top_v_tau` times the no-truncation entropy; it is independent of the
  candidate grid and has a single interpretable parameter.
  `"entropy_elbow"` is the previous small-gain elbow heuristic and
  `"neff"` uses the effective neighbourhood size without any grid
  search.

- top_v_tau:

  Saturation fraction in (0, 1) used by `top_v_method = "saturation"`.
  Default `0.9`.

- neff_quantile:

  Quantile of effective neighborhood size used by the `"neff"` top-v
  rule.

- model_top_v_tune_args:

  A named list of additional arguments passed to
  [`tune_model_top_v()`](https://novawz.github.io/multiRF/reference/tune_model_top_v.md)
  (e.g., `tmin`, `by`, `k`, `max_candidates`). The objective always
  follows `top_v_method`.

- fused_top_v_tune_args:

  A named list of additional arguments passed to
  [`tune_fused_top_v()`](https://novawz.github.io/multiRF/reference/tune_fused_top_v.md)
  (e.g., `vmin`, `by`, `vmax`, `k`). The objective always follows
  `top_v_method`.

- select_connection:

  Logical; whether to opt into the legacy quality-based
  direction-selection step. The default `FALSE` preserves all directed
  connections for modularity-weighted fusion.

- return_data:

  Logical; whether to include filtered/scaled data in output objects.

- compact_output:

  Logical; whether to drop heavy duplicated objects from output for
  lower memory usage (for example fitted model copies in robust
  clustering and full `vs_fit` object). Useful when keeping multiple
  `mrf3_fit` objects in memory.

- verbose:

  Logical; whether to print stage-level progress messages.

- seed:

  Random seed.

- ...:

  Additional arguments passed to
  [`mrf3_init()`](https://novawz.github.io/multiRF/reference/mrf3_init.md).

## Value

An `mrf3_fit` list with flat top-level components. Core fields are
`config`, `models`, `connection`, `connection_score`, `model_top_v`,
`fused_top_v`, `tuning_detail`, `reconstruction`, `clusters`, `shared`,
and `specific`. Optional stages populate `imd`, `cluster_imd`,
`selected_vars`, `selected_weights`, `selected_data`, `vs_detail`,
`robust_clusters`, and `robust_detail`; `data` is included when
`return_data = TRUE`.
