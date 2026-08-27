# Cluster Shared Similarity From Reconstruction Weights

Cluster Shared Similarity From Reconstruction Weights

## Usage

``` r
cluster_shared_similarity(
  recon,
  mode = c("average", "ao_reg"),
  dat_use = "ALL",
  k = NULL,
  method = c("Spectral", "PAM"),
  tune_method = "silhouette",
  gap_w = "uniform",
  similarity_type = c("second", "first"),
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
```

## Arguments

- recon:

  Reconstruction output from
  [`get_reconstr_matrix()`](mrf3_reconstr.md).

- mode:

  Shared clustering mode: `"average"` or `"ao_reg"`.

- dat_use:

  Shared weight source for `"average"` mode. Shared clustering now
  always uses `W_all`; this argument is kept for compatibility.

- k:

  Optional number of clusters. Required for `"ao_reg"` mode.

- method:

  Clustering backend for `"average"` mode.

- tune_method:

  Tuning criterion for PAM in `"average"` mode when `k` is `NULL`.

- gap_w:

  Weighting scheme for spectral eigengap in `"average"` mode when `k` is
  `NULL`.

- similarity_type:

  Similarity structure used in `"average"` mode.

- gamma:

  L2 regularization strength for `"ao_reg"` simplex weights.

- alpha_init:

  Optional initial AO fusion weights.

- ao_max_iter:

  Maximum AO iterations.

- ao_tol:

  AO convergence tolerance.

- knn_q:

  Optional kNN sparsification per base similarity in AO mode.

- hollow:

  Logical; whether to zero diagonal in AO mode.

- ao_symm:

  Logical; whether to symmetrize kNN sparsification in AO mode.

- ao_verbose:

  Logical; whether to print AO diagnostics.

- ...:

  Additional arguments passed to clustering backends.

## Value

A list with shared similarity, clustering result, and AO details (if
used).
