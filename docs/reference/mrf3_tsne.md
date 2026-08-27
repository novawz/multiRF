# compute tSNE embedding for MRF similarity matrix

compute tSNE embedding for MRF similarity matrix

## Usage

``` r
mrf3_tsne(
  mod,
  dims = 2,
  max_iter = 1000,
  learning_rate = 200,
  verbose = FALSE,
  tol = 1e-05,
  patience = 10,
  seed = NULL
)
```

## Arguments

- mod:

  mrf3 model

- dims:

  Number of dims

- max_iter:

  integer; Number of iterations (default: 1000)

- learning_rate:

  numeric; Learning rate (default: 200.0)

- verbose:

  Logical; whether to print optimization progress.

- tol:

  Convergence tolerance.

- patience:

  Early-stopping patience measured in iterations.

- seed:

  Optional integer seed used only to initialize the embedding.

## Value

Two dimensional tSNE embedding
