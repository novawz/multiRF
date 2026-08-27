# Build Shared/Specific Weights From Reconstruction

Build Shared/Specific Weights From Reconstruction

## Usage

``` r
get_shared_specific_weights(
  dat.list,
  recon,
  per_response_recon = TRUE,
  specific_top_v = NULL,
  specific_keep_ties = TRUE,
  specific_row_normalize = TRUE,
  specific_seed = NULL,
  specific_n_consensus = 1L,
  specific_ntree = 500L,
  specific_samptype = c("swor", "swr"),
  specific_ytry = NULL,
  specific_proximity = "all",
  specific_nsplit = 10L,
  specific_nthread = getOption("multiRF.nthread", 0L),
  specific_nodesize = 3L,
  specific_max_depth = 0L,
  specific_forest_wt = "all",
  ...
)
```

## Arguments

- dat.list:

  A named list of omics matrices (samples in rows, features in columns).

- recon:

  Reconstruction output from
  [`get_reconstr_matrix()`](mrf3_reconstr.md).

- per_response_recon:

  Logical; when `TRUE` (default), the shared reconstruction for each
  block `k` uses only the forest weight matrices from connections where
  block `k` is the *response*. Specifically, `W^(k) = sum_m alpha_m W_m`
  (modularity-weighted) over models `m` whose response is block `k`, and
  then `X_hat^(k) = W^(k) X^(k)`. This ensures that each block's
  residual captures genuinely block-specific variation rather than
  reconstruction artefacts from connections where the block served as
  predictor. When `FALSE`, reverts to the legacy behaviour that uses the
  global fused reconstruction (`recon$fused_mat`) or `W_{all} %*% X`.

- specific_top_v:

  Optional integer. If set, each row of fused specific residual-forest
  weights keeps top-v entries.

- specific_keep_ties:

  Logical; whether specific top-v truncation keeps ties at cutoff.

- specific_row_normalize:

  Logical; whether to row-normalize fused specific weights after
  optional truncation.

- specific_seed:

  Integer seed for the residual unsupervised forests. When `NULL`
  (default), falls back to 529 with a warning.
  [`mrf3_fit()`](mrf3_fit.md) automatically injects the pipeline seed
  here, so users calling [`mrf3()`](mrf3.md) or
  [`mrf3_fit()`](mrf3_fit.md) need not set this. Direct callers of
  `get_shared_specific_weights()` should pass an explicit positive seed
  for reproducibility.

- specific_n_consensus:

  Integer; number of consensus runs for the residual unsupervised
  forests. When `> 1`, fits the residual RF `specific_n_consensus` times
  with seeds `specific_seed, specific_seed + 1, ...` and averages the
  forest weights and IMD across runs. This reduces variance in the
  specific signal at the cost of proportionally more computation.
  Default `1L` (single run, no consensus averaging).

- specific_ntree, specific_samptype, specific_ytry, specific_proximity:

  RF structure inherited from the main forest workflow for residual
  unsupervised forests.

- specific_nsplit, specific_nthread:

  Candidate cutpoint and thread counts inherited from the main workflow.

- specific_nodesize, specific_max_depth:

  Residual-tree stopping settings.

- specific_forest_wt:

  Forest-weight mode for residual forests. The bioRxiv clustering
  definition uses `"all"`.

- ...:

  Deprecated/ignored arguments kept for compatibility.

## Value

A list with `shared` and `specific` components:

- `shared$W_all`: shared fused weights (global, for shared clustering).

- `shared$W_per_response`: named list of per-response fused weight
  matrices (only present when `per_response_recon = TRUE`).

- `specific$residual`: residual omics matrices `R = X - X_pred`.

- `specific$predicted`: predicted omics matrices `X_pred` from
  per-response (or global) reconstruction.

- `specific$residual_mod`: unsupervised RF models fitted on residual
  matrices.

- `specific$W`: specific residual weights from residual RF models after
  adjusted/truncate/row-sum-normalize.

- `specific$imd`: named list of per-block variable-level IMD weights
  (named numeric vectors) from unsupervised RF on residuals.
