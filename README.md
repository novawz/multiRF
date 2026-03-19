# multiRF

`multiRF` provides methods for multi-omics data integration using random forests
and MRF-style weighting across modalities.

## Installation

```r
# after publishing to GitHub
remotes::install_github("noblegasss/multiRF")

# or install from a local checkout
pak::pak(".")
```

## Quick Start

```r
library(multiRF)
data("tcga_brca_data")  # loads tcga_brca and tcga_brca_clinical

# tcga_brca is a named list of matched omics matrices
names(tcga_brca)
#> [1] "gene"  "methy" "mirna"

# run the main workflow on the bundled TCGA BRCA example
fit <- mrf3(
  tcga_brca,
  k = 4,
  ntree = 100,
  filter_mode = "none",
  run_imd = TRUE,
  seed = 529
)

# inspect the fitted object
summary(fit)
table(fit$clusters)

# extract top IMD variables from each omics block
lapply(get_top_vars(fit, n = 10), head)
```

`mrf3()` is the recommended user-facing entry point. It is a thin wrapper
around `mrf3_fit()`, so advanced workflow arguments can still be passed through
`...` when needed.

```r
fit_full <- mrf3(
  tcga_brca,
  k = 4,
  ntree = 200,
  run_imd = TRUE,
  run_variable_selection = TRUE,
  run_robust_clustering = TRUE,
  variable_selection_args = list(method = "mixture"),
  model_top_v = 50,
  fused_top_v = 30,
  seed = 529
)
```

## Citation

If you use multiRF in your research, please cite:

> Zhang, W. et al. (2025). multiRF: An integrative framework for multi-omics data using random forests. *GigaScience*, 14, giaf148. [DOI:10.1093/gigascience/giaf148](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaf148/8374728)

## Main Functions

- `mrf3()`: recommended end-to-end user entry point. It wraps `mrf3_fit()` with simpler defaults for clustering, while still accepting advanced workflow arguments through `...`.
- `mrf3_fit()`: full stage-based workflow when you want the complete parameter surface exposed explicitly.
- `mrf3_vs()`: variable-selection stage that uses IMD weights to keep informative features, with optional model refitting.

## Survival Example

```r
# align clinical rows to fitted cluster labels via sample IDs
clinical_ids <- tcga_brca_clinical$sampleID
clusters <- get_clusters(fit)
common_ids <- intersect(names(clusters), clinical_ids)
clinical_sub <- tcga_brca_clinical[match(common_ids, clinical_ids), , drop = FALSE]

km <- plot_km(
  test_var = clusters[common_ids],
  time_var = "OS.time",
  event_var = "OS",
  pheno_mat = clinical_sub
)

km
```

## Bundled Data

- `tcga_brca`: TCGA BRCA multi-omics example data with `gene`, `methy`, and `mirna` blocks.
- `tcga_brca_clinical`: matched clinical metadata for the BRCA cohort.
