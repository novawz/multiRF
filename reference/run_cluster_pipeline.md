# Run Shared/Specific Clustering Pipeline

Main clustering entrypoint for the current `mrf3_fit` workflow.

## Usage

``` r
run_cluster_pipeline(
  recon,
  shared_specific,
  mod_for_shared,
  clustering_args,
  cluster_method
)
```

## Arguments

- recon:

  Reconstruction object from
  [`get_reconstr_matrix()`](https://novawz.github.io/multiRF/reference/mrf3_reconstr.md).

- shared_specific:

  Shared/specific weight object from
  [`get_shared_specific_weights()`](https://novawz.github.io/multiRF/reference/get_shared_specific_weights.md).

- mod_for_shared:

  Random-forest model list used by proximity-based shared clustering.

- clustering_args:

  Named list of clustering arguments (for example `shared_k`,
  `specific_k`, `specific_prox_method_cl`, `tune_method`, `gap_w`).

- cluster_method:

  One of `"similarity"`, `"proximity"`, or `"enhanced_proximity"`.

## Value

A list with `shared` and `specific` clustering outputs.
