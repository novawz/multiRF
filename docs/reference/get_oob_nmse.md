# OOB normalized MSE for a fitted forest model

Uses OOB forest weights to predict both xvar and yvar, then computes the
mean normalized MSE over all predictor and response coordinates. Lower
is better; when predictor and response dimensions differ, every
coordinate still receives equal weight.

## Usage

``` r
get_oob_nmse(mod)
```

## Arguments

- mod:

  A fitted model from `fit_forest`.

## Value

A single numeric value (normalized MSE).
