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
  object = "entropy_elbow"
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

  Logical; whether to evaluate candidate grid values in parallel (POSIX
  systems only). Default is `TRUE`.

- cores:

  Number of cores used when `parallel = TRUE`. Default `NULL` uses
  `max(1, parallel::detectCores() - 1)` cores.

- seed:

  Random seed used for optional sample subsampling.

- object:

  Objective used to choose `model_top_v` (`"entropy_elbow"` (default),
  `"diss"`, `"silhouette"`, or `"eigen"`). For `"entropy_elbow"`, the
  elbow is selected among interior grid points, so the smallest grid
  candidate cannot be selected directly.
