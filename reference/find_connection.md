# Find optimal directional connections from fitted RF models

Scores each directional RF model using modularity of the
diagonal-adjusted forest weight matrix (and optionally OOB normalized
MSE). Two-step selection: (1) drop the bottom `drop_bottom_q` fraction
of connections by rank-sum of quality metrics; (2) `select_one_per_pair`
keeps at most one direction per omics pair (lowest rank-sum).

## Usage

``` r
find_connection(
  mod.list,
  return_score = FALSE,
  drop_bottom_q = 0.2,
  select_one_per_pair = TRUE,
  compute_oob = FALSE,
  seed = 529L,
  ...
)
```

## Arguments

- mod.list:

  A named list of fitted RF models. Model names must follow
  `response_predictor` convention.

- return_score:

  Logical; whether to return the directional score matrix.

- drop_bottom_q:

  Proportion of models to drop based on quality rank-sum. Must be in
  `[0, 1)`. Default is `0.2`.

- select_one_per_pair:

  Logical; whether to keep at most one direction per omics pair among
  quality-filtered models.

- compute_oob:

  Logical; whether to compute OOB normalized MSE in addition to
  modularity. When `TRUE`, `quality_score = modularity - oob_nmse` and
  rank-sum uses both ranks; when `FALSE`, `quality_score = modularity`
  and `rank_sum = rank(-modularity)`. Default is `FALSE`.

- seed:

  Random seed used by Louvain modularity scoring.

- ...:

  Additional arguments (currently ignored).

## Value

A character vector of selected model names (`response_predictor`). If
`return_score = TRUE`, returns a list with `model_connection`,
`connect_list`, `score`, and model-level quality metrics.
