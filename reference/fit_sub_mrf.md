# Fit an ensemble of sub-sampled MRF models for a single connection

When the number of response and/or predictor features is large, fitting
a single MRF with all features is computationally expensive.
`fit_sub_mrf()` addresses this by repeatedly sub-sampling features from
both the response and predictor blocks of a **single connection**,
fitting a smaller MRF each time, and averaging the resulting
forest-weight and proximity matrices.

## Usage

``` r
fit_sub_mrf(
  X,
  Y,
  n_sub = 15L,
  frac_response = 0.2,
  frac_predictor = 0.2,
  ntree_per_sub = 25L,
  mtry = NULL,
  ytry = NULL,
  min_response = 20L,
  min_predictor = 30L,
  enhanced = FALSE,
  compute_imd = FALSE,
  imd_args = list(),
  seed = 529L,
  parallel = FALSE,
  cores = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- X:

  Predictor data frame (samples x features) — a single omics block.

- Y:

  Response data frame (samples x features) — a single omics block.

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
  [`fit_forest()`](fit_forest.md).

- ytry:

  Number of candidate Y variables per split. `NULL` delegates to the
  forest engine (the multiRF default is `ceiling(qy / 3)`).

- min_response:

  Minimum number of response columns per sub-MRF.

- min_predictor:

  Minimum number of predictor columns per sub-MRF.

- enhanced:

  Logical; if `TRUE`, compute soft enhanced proximity using full-data
  sample embeddings within each sub-model. The result is stored as
  `$enhanced_prox` in the output so that
  [`mrf3_cl_prox()`](mrf3_cl_prox.md) can use it directly without
  re-traversing trees. Default `FALSE` because the extra computation is
  non-trivial.

- compute_imd:

  Logical; whether to pre-compute IMD weights across sub-models.

- imd_args:

  Named list of arguments passed to
  [`get_imp_forest()`](get_imp_forest.md) when `compute_imd = TRUE`.

- seed:

  Base random seed.

- parallel:

  Logical; if `TRUE`, fit replicates in fresh PSOCK worker processes.

- cores:

  Number of cores when `parallel = TRUE`; `NULL` (default) uses
  `parallel::detectCores() - 1`.

- verbose:

  Logical; print progress messages.

- ...:

  Additional arguments forwarded to [`fit_forest()`](fit_forest.md).

## Value

An rfsrc-like list with:

- forest.wt:

  Averaged n x n forest-weight matrix (all samples).

- forest.wt.oob:

  Averaged n x n OOB forest-weight matrix.

- proximity:

  Averaged n x n proximity matrix.

- xvar:

  Full predictor data (all columns of `X`).

- yvar:

  Full response data (all columns of `Y`).

- ntree:

  Total tree count across all sub-MRFs.

- sub_mrf_info:

  List with coverage and timing metadata.

## Details

Both `forest.wt` and `proximity` are n x n matrices that do NOT depend
on feature dimension. Averaging across sub-models with different feature
subsets produces a smoothed estimate.

`forest.wt.oob` is computed from each sub-model's `membership` and
`inbag` matrices: for sample i, only trees where i is out-of-bag
contribute to its weight row. This provides a regularized version
suitable for connection selection via
[`find_connection()`](find_connection.md).
