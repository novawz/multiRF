# Get multi-omics weights

Get multi-omics weights

## Usage

``` r
get_multi_weights(
  mod_list,
  dat.list,
  y = NULL,
  weighted = FALSE,
  use_depth = FALSE,
  robust = FALSE,
  parallel = FALSE,
  normalized = FALSE,
  calc = "Both",
  ytry = 1,
  w = NULL,
  cores = NULL,
  seed = -5,
  ...
)
```

## Arguments

- mod_list:

  A named list of fitted random forest models.

- dat.list:

  A named list of omics datasets used to align weights.

- y:

  Optional response vector (reserved for extensions).

- weighted:

  Logical; whether to use weighted importance updates.

- use_depth:

  Logical; whether to average depth-aware importances.

- robust:

  Logical; whether to use robust matrix-based aggregation. Requires
  `calc = "Both"`.

- parallel:

  Logical; whether to parallelize across models.

- normalized:

  Logical; whether to normalize the merged weights. The default is
  `FALSE`, preserving raw forest IMD on its 0-to-1 scale for variable
  selection.

- calc:

  Which importance side to compute: `"X"`, `"Y"`, or `"Both"`.

- ytry:

  Response sampling proportion used in node-level updates. Currently
  unused by the post-hoc tree traversal; reserved for future use.

- w:

  Optional case weights. Currently unused by the post-hoc tree
  traversal; reserved for future use.

- cores:

  Number of CPU cores used when `parallel = TRUE`; `NULL` (default) uses
  `parallel::detectCores() - 1`.

- seed:

  Random seed passed to stochastic components. Currently unused by the
  post-hoc tree traversal; reserved for future use.

- ...:

  Additional arguments for downstream helper functions.
