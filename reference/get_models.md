# Get fitted RF models

Extract the list of fitted random forest models.

## Usage

``` r
get_models(x, which = "init", ...)
```

## Arguments

- x:

  An `mrf3_fit` or `mrf3` object.

- which:

  For `mrf3_fit`: `"init"` (default) or `"refit"`
  (post-variable-selection refit, if available).

- ...:

  Additional arguments (unused).

## Value

A named list of rfsrc model objects, or `NULL`.
