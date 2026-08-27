# Simplified mrf3 Entry Point

A compact wrapper around [`mrf3_fit()`](mrf3_fit.md) for common use
cases. It keeps only the most frequently used parameters and forwards
advanced options through `...`. By default, this wrapper runs
shared/specific weighting and specific clustering branches.

## Usage

``` r
mrf3(
  dat.list,
  k = NULL,
  ntree = 500,
  top_v = NULL,
  samptype = c("swor", "swr"),
  main_clustering = c("similarity", "proximity", "enhanced_proximity"),
  filter_mode = c("auto", "none", "manual"),
  seed = 529,
  ...
)
```

## Arguments

- dat.list:

  A named list of omics matrices (samples in rows, features in columns).

- k:

  Optional cluster number used for both shared and specific clustering.
  If `NULL`, `k` is tuned in downstream clustering routines.

- ntree:

  Number of trees for RF fitting.

- top_v:

  Optional unified top-v cutoff. If set, applies to both `model_top_v`
  and `fused_top_v`.

- samptype:

  Sampling scheme passed through to [`mrf3_fit()`](mrf3_fit.md).

- main_clustering:

  Global clustering strategy applied consistently to shared, specific,
  and robust branches: `"similarity"` (default), `"proximity"`, or
  `"enhanced_proximity"`.

- filter_mode:

  Filtering mode passed to [`mrf3_fit()`](mrf3_fit.md).

- seed:

  Random seed.

- ...:

  Advanced options passed to [`mrf3_fit()`](mrf3_fit.md), e.g.
  `model_top_v_tune_args`, `fused_top_v_tune_args`, `clustering_args`,
  `run_imd`, `imd_args`, `run_variable_selection`,
  `variable_selection_args`, `cluster_imd_args`,
  `run_robust_clustering`, `compact_output`.

## Value

An object of class `"mrf3_fit"`.
