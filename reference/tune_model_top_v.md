# Tune model-level top-v

Tune model-level top-v

## Usage

``` r
tune_model_top_v(
  dat.list,
  mod,
  tmin = 10,
  by = 1,
  k = NULL,
  sample_n = NULL,
  sample_frac = NULL,
  auto_sample_n = FALSE,
  max_candidates = 20,
  reuse_tuned_k = TRUE,
  parallel = TRUE,
  cores = NULL,
  seed = 529,
  object = "saturation",
  tau = 0.9
)
```

## Arguments

- dat.list:

  A list of omics matrices used for clustering.

- mod:

  A fitted `mrf3` model object.

- tmin:

  Minimum model-level top-v cutoff to evaluate. The effective upper
  bound is the sample size, allowing the `v >= 0.8 * n` no-truncation
  rule to be evaluated.

- by:

  Base step size for model-level top-v grid construction.

- k:

  Optional fixed number of clusters.

- sample_n:

  Optional integer. If set, tune on a random subset of samples.

- sample_frac:

  Optional fraction in (0, 1\]. Used when `sample_n` is `NULL`.

- auto_sample_n:

  Logical; when `sample_n` and `sample_frac` are both `NULL`,
  automatically infer a tuning sample size from input data size.

- max_candidates:

  Optional maximum number of grid points to evaluate. If exceeded,
  candidate spacing is expanded automatically. Default is `20`.

- reuse_tuned_k:

  Logical; when `k` is `NULL`, tune `k` once on baseline and reuse it
  across all candidates for speed.

- parallel:

  Logical; whether to allow candidate-grid evaluation in fresh PSOCK
  worker processes. Parallel evaluation also requires an explicit
  `cores` value.

- cores:

  Number of cores used when `parallel = TRUE`. Default `NULL` keeps
  candidate-grid evaluation serial; supply a positive integer to use
  process-level workers.

- seed:

  Random seed used for optional sample subsampling.

- object:

  Objective used to choose `model_top_v`: `"saturation"` (default),
  `"entropy_elbow"`, `"diss"`, `"silhouette"`, or `"eigen"`.
  `"saturation"` selects the smallest `v` whose fused-weight row entropy
  reaches a fraction `tau` of the no-truncation entropy (linearly
  interpolated between grid points), so the choice does not depend on
  the grid resolution. `"entropy_elbow"` keeps the previous small-gain
  elbow heuristic; its elbow is selected among interior grid points, so
  the smallest grid candidate cannot be selected directly.

- tau:

  Saturation fraction in (0, 1) used by `object = "saturation"`. Default
  `0.9`.
