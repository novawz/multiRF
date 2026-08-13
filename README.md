# multiRF

**Fast multivariate random forests for multi-omics integration**

`multiRF` is an R package for integrating matched multi-omics datasets with
multivariate random forests (MRF). It fits directed forest models across omics
blocks, learns sample-by-sample similarity from shared terminal-node structure,
and decomposes the result into shared and omics-specific components for
clustering, variable selection, and visualization.

The package now uses a native C++ backend for multivariate regression,
unsupervised forests, forest weights, proximity matrices, and enhanced
proximity with sibling-leaf corrections. In practice, this gives a simpler
installation path and a faster MRF than the `randomForestSRC`-based MRF while keeping the same overall modeling logic.

Project website: [https://novawz.github.io/multiRF/](https://novawz.github.io/multiRF/)

## Installation

```r
remotes::install_github("novawz/multiRF")
```

The package compiles from source and requires a C++17 toolchain:

- macOS: Xcode Command Line Tools
- Windows: Rtools
- Linux: `g++` or `clang++`

OpenMP is recommended for parallel tree construction. `randomForestSRC` is not
required for the default workflow.

## Quick start

```r
library(multiRF)
data("tcga_brca_data")

names(tcga_brca)
#> [1] "gene"  "methy" "mirna"

fit <- mrf3(
  tcga_brca,
  k = 4,
  ntree = 100,
  filter_mode = "none",
  run_imd = TRUE,
  seed = 529
)

summary(fit)
table(get_clusters(fit))
get_top_vars(fit, n = 10)
```

`mrf3()` is the main user-facing entry point. It wraps the staged workflow in
`mrf3_fit()` and forwards advanced arguments through `...`. Both entry points
use automatic feature filtering by default; set `filter_mode = "none"` when
working with an already curated feature set, as in the example above.

The default similarity workflow retains all directed forests. It removes
their diagonal weights, row-normalizes after model top-v truncation, normalizes
modularity scores within each response block, and uniformly averages the
per-response matrices (Eqs. 6–8). A top-v value at least 80% of the sample size
is treated as no truncation. Shared and omics-specific similarities use
`S = W W^T` with a zero diagonal and spectral clustering by default.

## What it provides

- `mrf3()`: end-to-end workflow for fitting, reconstruction, and clustering
- `mrf3_fit()`: staged workflow with the full parameter surface exposed
- `mrf3_vs()`: variable selection from IMD weights
- `mrf3_stability()`: resampling-based cluster stability assessment
- `pairwise_imd()`: variable-level co-occurrence network analysis
- `plot_weights()`, `plot_tsne()`, `plot_umap()`, `plot_network()`: consistent, publication-oriented result visualizations
- `plot_km()`, `plot_circos()`: optional clinical and circular-network visualizations

## Clustering modes

```r
fit_sim <- mrf3(
  tcga_brca,
  k = 4,
  ntree = 100,
  main_clustering = "similarity", # default
  seed = 529
)

fit_prox <- mrf3(
  tcga_brca,
  k = 4,
  ntree = 100,
  main_clustering = "proximity",
  seed = 529
)

fit_enh <- mrf3(
  tcga_brca,
  k = 4,
  ntree = 100,
  main_clustering = "enhanced_proximity",
  seed = 529
)
```

## Full workflow example

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

Variable selection uses raw forest IMD on `[0, 1]`; `filter` tunes
`tau * sd(IMD)` using true OOB normalized
MSE across predictor and response coordinates; `mixture` includes a point mass
at zero; and `transformation` uses Eq. 16 with a Student-t reference
(`df = ntree - 1`). The newer connection-selection strategy is independent of
these feature-selection rules.

## Bundled data

- `tcga_brca`: TCGA BRCA example with `gene`, `methy`, and `mirna` blocks
- `tcga_brca_clinical`: matched clinical annotations including subtype and survival information

## Citation

If you use `multiRF` in your research, please cite:

> Zhang, W., Wang, L., Franzmann, E. J., and Chen, X. S. (2026). Multivariate Random Forests for Cross-Modal Multi-Omics Integration. *bioRxiv*. [doi:10.64898/2026.06.17.732933](https://doi.org/10.64898/2026.06.17.732933)

> Zhang, W. et al. (2025). An integrative multi-omics random forest framework for robust biomarker discovery. *GigaScience*, 14, giaf148. [doi:10.1093/gigascience/giaf148](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaf148/8374728)
