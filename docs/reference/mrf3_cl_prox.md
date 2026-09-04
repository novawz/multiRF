# MRF unsupervised clustering – Proximity method

MRF unsupervised clustering – Proximity method

## Usage

``` r
mrf3_cl_prox(
  rfit,
  k = NULL,
  enhanced = TRUE,
  size_min = 5,
  use = "X",
  symm = TRUE,
  leaf_embed_dim = 10,
  merge_quantile = 0.9,
  merge_mode = c("soft", "hard"),
  sibling_gamma = 0.5,
  sibling_fun = c("constant", "positive", "shift01", "abs", "identity"),
  sibling_cap = TRUE,
  hard_prox_mode = c("soft_enhanced", "binary"),
  parallel = TRUE,
  sparse = FALSE,
  method_cl = "PAM",
  tune_method = "silhouette",
  gap_w = "uniform",
  cores = NULL,
  ...
)
```

## Arguments

- rfit:

  A model list of random forest models.

- k:

  Pre-defined number of clusters. By default, the optimal value is
  selected by the requested tuning method.

- enhanced:

  Logical; whether to calculate enhanced proximity.

- size_min:

  Minimum terminal-node size when computing enhanced proximity.

- use:

  Which importance side to use in tree traversal, `"X"` or `"Y"`.

- symm:

  Logical; whether to symmetrize directional proximities.

- leaf_embed_dim:

  Dimension used for low-dimensional leaf embedding.

- merge_quantile:

  Quantile threshold for sibling-pair merges in enhanced mode. `0.9`
  keeps the top 10\\ iteration. Used only when `merge_mode = "hard"`.

- merge_mode:

  Enhanced-proximity merge mode: `"soft"` (default) or `"hard"`.
  `"hard"` performs iterative sibling merges; `"soft"` directly builds
  per-tree similarity with same-leaf = 1 and sibling-leaf =
  `sibling_gamma * f(corr)`. `"soft"` is generally more stable.

- sibling_gamma:

  Multiplicative weight `gamma` used in soft merge mode.

- sibling_fun:

  Correlation transform `f(corr)` used in soft merge mode. One of
  `"constant"`, `"positive"`, `"shift01"`, `"abs"`, or `"identity"`.

- sibling_cap:

  Logical; whether to cap soft sibling similarities at 1.

- hard_prox_mode:

  Proximity construction used after hard merge: `"soft_enhanced"` uses
  same-enhanced-leaf = 1 and sibling-enhanced-leaf =
  `sibling_gamma * f(corr)`, while `"binary"` uses same-enhanced-leaf
  only.

- parallel:

  Logical; whether to allow forest-level computation in fresh PSOCK
  worker processes. Parallel execution also requires an explicit `cores`
  value.

- sparse:

  Logical; whether to sparsify the enhanced proximity matrix.

- method_cl:

  Clustering backend (`"PAM"` or `"Spectral"`).

- tune_method:

  Tuning criterion for PAM when `k` is `NULL` (`"silhouette"` or
  `"ratio"`).

- gap_w:

  Weighting scheme for spectral eigengap when `k` is `NULL` (`"uniform"`
  or `"log"`).

- cores:

  Number of CPU cores used by parallel steps. Default `NULL` keeps
  matrix-heavy forest reductions serial to limit peak memory.

- ...:

  Additional arguments passed to downstream clustering helpers.

## Value

An `mrf3` clustering object.
