# multiRF

**Fast multivariate random forests for multi-omics integration**

`multiRF` integrates matched multi-omics datasets by fitting multivariate
random forests across pairs of data blocks, learning sample-by-sample
similarity structures directly from the data.  Two samples are considered
similar when trees consistently route them to the same terminal node,
yielding forest-weight and proximity matrices that capture cross-platform
relationships without predefined distance metrics.  These learned
similarities are decomposed into shared (cross-omics) and specific
(within-omics) components for downstream clustering, variable selection,
and visualization.

The computational backend is a native C++ engine (via Rcpp) that builds
multivariate regression and unsupervised forests with OpenMP thread-level
parallelism.  Forest weights, proximity matrices, and enhanced proximity
(with sibling-leaf corrections) are computed during tree construction in a
single pass.  On typical multi-omics datasets (500-1000 samples, 3
platforms, 300 trees), the native engine runs 20-30% faster than an
equivalent `randomForestSRC` pipeline while producing statistically
equivalent output.

## Multivariate splitting

Standard random forests split each node by minimizing impurity on a single
response.  `multiRF` uses a composite criterion over a vector-valued
response (Li and Xiao, 2011).  At each candidate split, the algorithm
evaluates a normalized between-group sum of squares across a random subset
of `ytry` response columns:

For a node with *n* samples and a candidate partition into left (*L*) and
right (*R*) children, define the centered left-child sum for response
column *j*:

    S_{L,j} = sum_{i in L} (Y_{ij} - Y_bar_j)

The node-level variance is sigma_j^2 = (1/n) sum_i (Y_{ij} - Y_bar_j)^2.
The standardized split statistic for column *j* is:

    Delta_j = S_{L,j}^2 / (n_L * sigma_j^2)  +  S_{R,j}^2 / (n_R * sigma_j^2)

The composite statistic averages over all informative columns (those with
nonzero variance):

    Delta = (1 / |J|) * sum_{j in J} Delta_j

The partition maximizing Delta is selected.  Normalizing each column by its
own variance ensures all responses contribute on a comparable scale.  This
matches the normalized composite splitting rule in `randomForestSRC`
(Ishwaran et al., 2008; Tang and Ishwaran, 2017).

## Installation

```r
remotes::install_github("noblegasss/multiRF")
```

The native C++ engine compiles from source and requires a C++17 compiler
(`Rtools` on Windows, `Xcode` CLI on macOS, `g++`/`clang++` on Linux).
OpenMP is optional but recommended.  `randomForestSRC` is not required.

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

`mrf3()` is the recommended entry point.  It wraps `mrf3_fit()` with
simpler defaults while accepting advanced arguments through `...`.

### Clustering backends

```r
# proximity-based clustering
fit_prox <- mrf3(tcga_brca, k = 4, main_clustering = "proximity",
                 ntree = 100, seed = 529)

# enhanced proximity (sibling-leaf corrections, computed in C++)
fit_enh <- mrf3(tcga_brca, k = 4, main_clustering = "enhanced_proximity",
                ntree = 100, seed = 529)
```

### Full workflow with variable selection

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

## Main functions

| Function | Purpose |
|----------|---------|
| `mrf3()` | End-to-end workflow (fit + cluster) |
| `mrf3_fit()` | Staged workflow with full parameter surface |
| `mrf3_vs()` | Variable selection using IMD weights |
| `mrf3_stability()` | Resampling-based stability assessment |
| `pairwise_imd()` | Variable-level co-occurrence network |
| `filter_omics()` | Pre-fitting feature filtering |

## Bundled data

- `tcga_brca`: TCGA BRCA multi-omics example with `gene`, `methy`, and
  `mirna` blocks.
- `tcga_brca_clinical`: matched clinical metadata including PAM50 subtypes
  and survival endpoints.

## Citation

If you use multiRF in your research, please cite:

> Zhang, W. et al. (2025). An integrative multi-omics random forest
> framework for robust biomarker discovery. *GigaScience*, 14, giaf148.
> [doi:10.1093/gigascience/giaf148](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaf148/8374728)

## References

- Li, M. and Xiao, Y. (2011). Multivariate random forests. *Wiley
  Interdisciplinary Reviews: Data Mining and Knowledge Discovery*, 1(1),
  80-87.
- Ishwaran, H. et al. (2008). Random survival forests. *Annals of Applied
  Statistics*, 2(3), 841-860.
- Tang, F. and Ishwaran, H. (2017). Random forest missing data algorithms.
  *Statistical Analysis and Data Mining*, 10(6), 363-377.
- Ishwaran, H. and Kogalur, U.B. randomForestSRC: multivariate splitting
  rule.
  [https://www.randomforestsrc.org/articles/mvsplit.html](https://www.randomforestsrc.org/articles/mvsplit.html)
