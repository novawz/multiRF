# mrf3 Stability Evaluation (No Refit)

Evaluate clustering stability by repeated subsampling/bootstrapping on
existing `mrf3_fit` outputs, without refitting forests.

## Usage

``` r
mrf3_stability(
  x,
  branches = c("specific_shared", "robust_clustering"),
  n_rep = 50L,
  sample_frac = 0.8,
  sample_mode = c("subsample", "bootstrap"),
  specific_k = NULL,
  robust_k = NULL,
  specific_method = c("auto", "PAM", "Spectral"),
  robust_method = c("auto", "PAM", "Spectral"),
  run_consensus = TRUE,
  consensus_method = c("PAM", "Spectral"),
  truth = NULL,
  seed = 529,
  verbose = TRUE
)
```

## Arguments

- x:

  An object from
  [`mrf3_fit()`](https://novawz.github.io/multiRF/reference/mrf3_fit.md).

- branches:

  Branches to evaluate. Supported: `"specific_shared"` and
  `"robust_clustering"`. Branches unavailable in `x` are skipped with a
  warning; an error is raised only if none of the requested branches are
  available.

- n_rep:

  Number of repeats.

- sample_frac:

  Sampling fraction per repeat.

- sample_mode:

  Sampling mode: `"subsample"` or `"bootstrap"`.

- specific_k:

  Optional `k` override for `"specific_shared"` branch.

- robust_k:

  Optional `k` override for `"robust_clustering"` branch.

- specific_method:

  Clustering backend for `"specific_shared"` repeats: `"auto"`,
  `"Spectral"`, or `"PAM"`.

- robust_method:

  Clustering backend for `"robust_clustering"` repeats: `"auto"`,
  `"Spectral"`, or `"PAM"`.

- run_consensus:

  Logical; whether to derive consensus clustering.

- consensus_method:

  Consensus clustering backend: `"Spectral"` or `"PAM"`.

- truth:

  Optional truth labels for external metrics. If names are provided,
  they are aligned to sample names.

- seed:

  Random seed.

- verbose:

  Logical; whether to print progress.

## Value

A list with per-branch repeat metrics and optional consensus results.
