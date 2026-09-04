# Fit sub-MRF ensemble for a multi-omics connection list

Higher-level wrapper that matches the interface of
[`fit_multi_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md).
For each directed connection (response_block -\> predictor_block),
features are sub-sampled from both blocks independently, and the
resulting n x n matrices are averaged across replicates.

## Usage

``` r
fit_sub_multi_rfsrc(
  dat.list,
  connect_list,
  n_sub = 15L,
  frac_response = 0.2,
  frac_predictor = 0.2,
  ntree_per_sub = 25L,
  mtry = NULL,
  ytry = NULL,
  min_response = 20L,
  min_predictor = 30L,
  min_response_for_sub = 500L,
  min_predictor_for_sub = 500L,
  ntree_full = 300L,
  enhanced = FALSE,
  compute_imd = FALSE,
  imd_args = list(),
  seed = 529L,
  parallel = FALSE,
  cores = 2L,
  parallel_connections = FALSE,
  cores_connections = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- dat.list:

  Named list of data frames (samples x features).

- connect_list:

  A list of connections, each a character vector of length 2:
  `c(response_name, predictor_name)`. Same format as used by
  [`fit_multi_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md).

- n_sub:

  Integer; number of sub-MRF replicates to fit.

- frac_response:

  Fraction of response columns to sample per replicate (default 0.2).
  Convergence analysis shows symmetric low fractions (0.1–0.3 for both
  response and predictor) yield the best approximation of the full
  forest-weight matrix.

- frac_predictor:

  Fraction of predictor columns to sample per replicate (default 0.2).
  Avoid setting this to 1.0 when `frac_response` is small, as the
  asymmetric case leads to poor convergence (each sub-model overfits to
  its response subset).

- ntree_per_sub:

  Number of trees per sub-MRF (default 25 for direct calls; the workflow
  inherits the main tree budget by using `ceiling(ntree / n_sub)` unless
  explicitly overridden).

- mtry:

  Number of candidate X variables per split. Passed through to
  [`fit_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md).

- ytry:

  Number of candidate Y variables per split (0 = `min(qy, ceiling(p/3))`
  default).

- min_response:

  Minimum number of response columns per sub-MRF.

- min_predictor:

  Minimum number of predictor columns per sub-MRF.

- min_response_for_sub:

  Minimum response-block size below which all response features are used
  instead of sub-sampling.

- min_predictor_for_sub:

  Minimum predictor-block size below which all predictor features are
  used instead of sub-sampling.

- ntree_full:

  Number of trees used when both blocks are below their sub-sampling
  thresholds and a full forest is fitted instead.

- enhanced:

  Logical; if `TRUE`, compute soft enhanced proximity using full-data
  sample embeddings within each sub-model. The result is stored as
  `$enhanced_prox` in the output so that
  [`mrf3_cl_prox()`](https://novawz.github.io/multiRF/reference/mrf3_cl_prox.md)
  can use it directly without re-traversing trees. Default `FALSE`
  because the extra computation is non-trivial.

- compute_imd:

  Logical; whether to pre-compute IMD weights across sub-models.

- imd_args:

  Named list of arguments passed to
  [`get_imp_forest()`](https://novawz.github.io/multiRF/reference/get_imp_forest.md)
  when `compute_imd = TRUE`.

- seed:

  Base random seed.

- parallel:

  Logical; if `TRUE`, fit replicates in fresh PSOCK worker processes.

- cores:

  Number of cores used for within-connection parallelism. The default
  for this wrapper is `2L`.

- parallel_connections:

  Logical; whether distinct directed connections may be fitted in fresh
  PSOCK worker processes.

- cores_connections:

  Optional total core budget for connection-level parallelism.

- verbose:

  Logical; print progress messages.

- ...:

  Passed to
  [`fit_sub_mrf()`](https://novawz.github.io/multiRF/reference/fit_sub_mrf.md)
  and then
  [`fit_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md).

## Value

A **named list** of rfsrc-compatible objects, one per connection. Names
follow the `"response_predictor"` convention used by
[`fit_multi_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md).
Each object contains `forest.wt`, `forest.wt.oob`, `proximity`, `xvar`,
`yvar`.

## Details

Sub-MRF uses multiRF for its internal fits. OOB forest-weight
reconstruction uses the returned `inbag` matrix.
