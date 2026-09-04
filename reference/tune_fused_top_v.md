# Tune fused weight top-v truncation

Tune fused weight top-v truncation

## Usage

``` r
tune_fused_top_v(
  dat.list,
  mod,
  vmin = 10,
  by = 1,
  vmax = NULL,
  model_top_v = 10,
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
  tau = 0.9,
  early_stop = FALSE,
  elbow_rel_tol = 0.25,
  elbow_abs_tol = 1e-04,
  elbow_min_points = 4L,
  elbow_min_frac = 0.2,
  elbow_patience = 2L,
  elbow_smooth_window = 3L
)
```

## Arguments

- dat.list:

  A list of omics matrices used for clustering.

- mod:

  A fitted `mrf3` model object.

- vmin:

  Minimum `fused_top_v` cutoff to evaluate.

- by:

  Base step size for the `fused_top_v` grid.

- vmax:

  Maximum `fused_top_v` cutoff to evaluate. If `NULL`, the full sample
  size is used so the `v >= 0.8 * n` no-truncation rule is represented
  in the grid.

- model_top_v:

  Fixed model-level top-v cutoff used while tuning `fused_top_v`.

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

  Objective used to choose `fused_top_v`: `"saturation"` (default),
  `"entropy_elbow"`, `"diss"`, `"silhouette"`, or `"eigen"`. For the two
  entropy objectives the entropy of every truncation level is obtained
  in closed form from the sorted fused weights (see
  [`fused_entropy_curve()`](fused_entropy_curve.md)), so no candidate
  matrices are rebuilt. `"saturation"` evaluates every integer `v` in
  `[vmin, vmax]` and selects the smallest one whose entropy reaches
  `tau` times the no-truncation entropy. `"entropy_elbow"` keeps the
  previous small-gain elbow heuristic on the candidate grid; its elbow
  is selected among interior grid points, so the smallest grid candidate
  and the no-truncation baseline cannot be selected directly (no
  truncation is recovered separately via the `v >= 0.8 * n` rule in the
  workflow).

- tau:

  Saturation fraction in (0, 1) used by `object = "saturation"`. Default
  `0.9`.

- early_stop:

  Logical; when `TRUE` and `object = "entropy_elbow"`, stop tuning once
  a stable small-gain elbow is reached. Default is `FALSE` so elbow
  selection is based on the full evaluated grid.

- elbow_rel_tol:

  Relative elbow threshold multiplier on max observed gain.

- elbow_abs_tol:

  Absolute lower bound for elbow threshold.

- elbow_min_points:

  Minimum evaluated points before early-stop is allowed.

- elbow_min_frac:

  Minimum evaluated fraction in `v` grid before early-stop.

- elbow_patience:

  Number of consecutive small-gain steps required to stop.

- elbow_smooth_window:

  Integer running-mean window used to smooth entropy gains before elbow
  detection.
