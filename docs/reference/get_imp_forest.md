# Get forest importance

Get forest importance

## Usage

``` r
get_imp_forest(
  mod,
  parallel = FALSE,
  robust = FALSE,
  calc = "Both",
  weighted = FALSE,
  use_depth = FALSE,
  normalized = FALSE,
  w = NULL,
  ytry = 1,
  cores = NULL,
  seed = -5
)
```

## Arguments

- mod:

  A fitted forest model from the native multiRF engine or the optional
  `randomForestSRC` fallback.

- parallel:

  Logical; whether to parallelize across trees.

- robust:

  Logical; whether to use robust matrix-based aggregation. Requires
  `calc = "Both"`.

- calc:

  Which importance side to compute: `"X"`, `"Y"`, or `"Both"`.

- weighted:

  Logical; whether to use weighted importance updates.

- use_depth:

  Logical; whether to aggregate non-zero depths instead of simple mean.

- normalized:

  Logical; whether to l2-normalize returned importance.

- w:

  Optional case weights. Currently unused by the post-hoc tree
  traversal; reserved for future use.

- ytry:

  Response sampling proportion used in node-level updates. Currently
  unused by the post-hoc tree traversal; reserved for future use.

- cores:

  Number of CPU cores used when `parallel = TRUE`; `NULL` (default) uses
  `parallel::detectCores() - 1`.

- seed:

  Random seed passed to stochastic components. Currently unused by the
  post-hoc tree traversal; reserved for future use.
