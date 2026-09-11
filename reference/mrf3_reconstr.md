# MRF unsupervised clustering – Reconstruction method

MRF unsupervised clustering – Reconstruction method

## Usage

``` r
mrf3_reconstr(
  recon = NULL,
  rfit = NULL,
  k = NULL,
  weights_cluster = TRUE,
  dat_use_to_cluster = "ALL",
  model_top_v = 10,
  recon_fusion = c("weighted", "uniform"),
  global_fusion = c("average", "pmin"),
  connection_score = NULL,
  response_blocks = NULL,
  score_power = 1,
  score_floor = 0,
  fallback_uniform = TRUE,
  fused_top_v = NULL,
  fused_row_normalize = TRUE,
  fused_keep_ties = TRUE,
  fusion_mode = "average",
  gamma = 0.1,
  alpha_init = NULL,
  ao_max_iter = 20,
  ao_tol = 1e-04,
  knn_q = NULL,
  hollow = TRUE,
  ao_symm = TRUE,
  ao_verbose = FALSE,
  ...
)

get_reconstr_matrix(
  rfit,
  model_top_v = 10,
  recon_fusion = c("weighted", "uniform"),
  global_fusion = c("average", "pmin"),
  connection_score = NULL,
  response_blocks = NULL,
  score_power = 1,
  score_floor = 0,
  fallback_uniform = TRUE,
  fused_top_v = NULL,
  fused_row_normalize = TRUE,
  fused_keep_ties = TRUE
)
```

## Arguments

- recon:

  Optional reconstruction object returned by `get_reconstr_matrix()`. If
  `NULL`, reconstruction is built from `rfit`.

- rfit:

  A model list of random forest models.

- k:

  Pre-defined number of clusters. The default is selecting the optimal k
  by tuning method.

- weights_cluster:

  Logical; whether to cluster based on reconstructed weight affinity.

- dat_use_to_cluster:

  If clustering method selected as Reconstr, the user can choose the
  data to conduct clustering. In weight-clustering mode, shared
  clustering uses `W_all`; this argument only applies when
  `weights_cluster = FALSE`.

- model_top_v:

  Model-level top-v cutoff applied to each single-model forest weight
  matrix before fusion.

- recon_fusion:

  Reconstruction fusion mode. `"weighted"` (default) uses
  `connection_score`; `"uniform"` uses equal weights.

- global_fusion:

  Fusion across block-specific matrices. `"average"` implements a
  uniform average; `"pmin"` is an optional element-wise intersection
  extension.

- connection_score:

  Optional directional score matrix from
  `find_connection(return_score = TRUE)`.

- response_blocks:

  Optional character vector containing the data blocks represented in
  the fused matrix. For each block, reconstruction first uses forests
  where it is the response, then forests where it is the predictor. A
  block absent from every fitted connection uses the average of the
  available block-specific matrices.

- score_power:

  Exponent applied to raw connection scores before normalization.

- score_floor:

  Non-negative floor applied to raw connection scores before
  normalization.

- fallback_uniform:

  Logical; whether to fallback to uniform averaging when weighted scores
  are unavailable.

- fused_top_v:

  Optional integer. If set, apply row-wise top-v truncation to fused
  weights (`W_all`) after fusion. If `FALSE`, no fused top-v truncation
  is applied.

- fused_row_normalize:

  Logical; whether to row-normalize fused weights after optional
  truncation.

- fused_keep_ties:

  Logical; whether fused top-v truncation keeps ties at cutoff.

- fusion_mode:

  Clustering-affinity fusion mode. `"average"` keeps current behavior;
  `"ao_reg"` uses alternating optimization with simplex-regularized
  weights.

- gamma:

  L2 regularization strength for the simplex weights in `"ao_reg"` mode.

- alpha_init:

  Optional initial fusion weights for `"ao_reg"` mode.

- ao_max_iter:

  Maximum alternating-optimization iterations.

- ao_tol:

  Convergence tolerance for `||alpha_t - alpha_{t-1}||_1`.

- knn_q:

  Optional kNN sparsification per base similarity (rows keep top-q).

- hollow:

  Logical; whether to zero out diagonal in base/fused similarities.

- ao_symm:

  Logical; whether to symmetrize kNN-sparsified similarities by
  `pmax(S, t(S))`.

- ao_verbose:

  Logical; whether to print AO iteration diagnostics.

- ...:

  Additional arguments passed to clustering backends.

## Value

mrf3 clustering object
