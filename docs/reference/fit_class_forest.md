# Fit a classification forest (native engine)

Uses the native multivariate regression forest on a one-hot encoded
response, then reconstructs training-set class probabilities and labels.
`err.rate` is the out-of-bag misclassification rate computed from OOB
forest weights (`NA` when OOB information is unavailable); samples with
zero forest weight receive `NA` probabilities and predictions.

## Usage

``` r
fit_class_forest(
  X,
  Y,
  ntree = 500L,
  mtry = NULL,
  nsplit = 10L,
  forest.wt = "all",
  proximity = c("all", "inbag", "oob", "none"),
  nodesize = 1L,
  max_depth = 0L,
  seed = 529L,
  samptype = c("swor", "swr"),
  nthread = getOption("multiRF.nthread", 0L),
  xvar.wt = NULL
)
```
