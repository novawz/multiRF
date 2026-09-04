# Joint Shared/Specific Similarity Clustering Wrapper

Joint Shared/Specific Similarity Clustering Wrapper

## Usage

``` r
mrf3_specific_clustering(
  recon,
  shared_specific,
  shared_mode = c("average", "ao_reg"),
  shared_dat_use = "ALL",
  shared_k = NULL,
  shared_method = c("Spectral", "PAM"),
  shared_similarity_type = c("second", "first"),
  specific_k = NULL,
  specific_method = c("Spectral", "PAM", "Proximity", "Enhanced_Proximity"),
  specific_prox_method_cl = c("PAM", "Spectral"),
  specific_similarity_type = c("second", "first"),
  tune_method = "silhouette",
  gap_w = "uniform",
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
  [`get_reconstr_matrix()`](https://novawz.github.io/multiRF/reference/mrf3_reconstr.md).

- shared_specific:

  Output from
  [`get_shared_specific_weights()`](https://novawz.github.io/multiRF/reference/get_shared_specific_weights.md).

- shared_mode:

  Shared clustering mode passed to
  [`cluster_shared_similarity()`](https://novawz.github.io/multiRF/reference/cluster_shared_similarity.md).

- shared_dat_use:

  Shared weight source for `"average"` mode.

- shared_k:

  Optional `k` for shared clustering. Required for `"ao_reg"` mode.

- shared_method:

  Clustering backend for shared `"average"` mode.

- shared_similarity_type:

  Similarity structure for shared `"average"` mode.

- specific_k:

  Optional `k` configuration passed to
  [`cluster_specific_similarity()`](https://novawz.github.io/multiRF/reference/cluster_specific_similarity.md).

- specific_method:

  Clustering backend for specific clustering. Choose from `"Spectral"`,
  `"PAM"`, `"Proximity"`, and `"Enhanced_Proximity"`.

- specific_prox_method_cl:

  Proximity clustering backend used when `specific_method` is
  `"Proximity"` or `"Enhanced_Proximity"`.

- specific_similarity_type:

  Similarity structure for specific clustering.

- tune_method:

  Tuning criterion for PAM when tuning `k`.

- gap_w:

  Weighting scheme for spectral eigengap when tuning `k`.

- gamma:

  L2 regularization for shared AO fusion.

- alpha_init:

  Optional initial AO fusion weights.

- ao_max_iter:

  Maximum AO iterations.

- ao_tol:

  AO convergence tolerance.

- knn_q:

  Optional kNN sparsification for shared AO fusion.

- hollow:

  Logical; whether to zero diagonal in shared AO fusion.

- ao_symm:

  Logical; whether to symmetrize kNN sparsification in AO fusion.

- ao_verbose:

  Logical; whether to print AO diagnostics.

- ...:

  Additional arguments passed to clustering backends.

## Value

A list with `shared` and `specific` clustering outputs.
