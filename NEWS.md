# multiRF 0.2.3

Maintenance release applying the fixes from a full-package audit (R CMD check
plus a verified line-by-line review of all R and C++ sources). No new features.

## Correctness fixes that change numerical results

- **IMD tree-network edge ratios**: `build_tree_network_cpp` assigned the wrong
  parent nodesize after backtracking, corrupting the `edge` ratio of
  second-child edges whose sibling subtree is internal (22–24% of all edges,
  median ~7x too small). An old-vs-new comparison on bit-identical fits showed
  no reachable IMD weight path consumed the corrupted values (the native
  engine's pre-computed IMD never uses the network; `weighted = FALSE` never
  reads `edge`; `weighted = TRUE, calc = "Y"` reads only the never-corrupted
  first-child edges; `calc = "X"/"Both"` crashed outright pre-fix), so IMD
  weights from earlier versions do NOT need re-running. The fix corrects the
  `edge` column of returned network data frames and, together with the
  aggregation fix below, makes `weighted = TRUE` usable for the first time.
- **Spectral clustering**: `spectral_cl()` now runs PAM on the row-normalized
  NJW embedding it computes (and returns), as the algorithm requires;
  previously it clustered the unnormalized eigenvectors.
- **OOB forest weights (sub-MRF path)**: `compute_oob_fw()` no longer lets a
  sample donate to its own OOB weight row; the diagonal is zero and rows are
  renormalized, matching the canonical C++ implementation.
- `cluster_metric_matrix()` computes the asymmetric purity metric in both
  directions instead of mirroring one triangle.
- Robust normalized importance used an L1 norm for Y but L2 for X due to a
  transposed expression; both now use L2.
- The unsupervised C++ forest now uses portable random draws (hand-rolled
  integer/real draws and Fisher–Yates shuffle over `std::mt19937`), making
  seeded results reproducible across platforms. Draw sequences differ from
  0.2.2, so seeded unsupervised results will change once.

## Crash fixes

- `mrf3_stability()` with default arguments no longer stops on fits without
  robust clustering (unavailable branches are skipped with a warning), and the
  robust branch no longer treats the `robust_clusters` list as a label vector
  (stability metrics were silently NA/garbage).
- `get_multi_weights(weighted = TRUE)` no longer errors (the `edge` column was
  dropped during aggregation).
- `mrf3_init(sub_mrf = TRUE, select_connection = TRUE, compute_oob = TRUE)`
  works: `compute_oob_forest_wt()` consumes a model's `forest.wt.oob` directly
  when present.
- Single-column response blocks no longer crash importance extraction
  (missing `drop = FALSE`).
- Legacy character-vector `connection` models no longer crash
  `get_reconstr_matrix()`.
- `pairwise_imd()`: duplicated feature names across omics blocks now raise an
  informative error instead of silently misattributing IMD mass; single-block
  and 1×1 edge cases are handled.
- Fixed a dangling-reference read after `nodes.push_back()` in the C++ tree
  builders (undefined behavior), and added missing `size_t` casts in the
  unsupervised n×n buffer indexing (overflow for n ≥ 46342).

## Performance (identical results, faster)

- Native engine speedups, all verified bit-identical across supervised,
  unsupervised, proximity, and weighted-sampling configurations: the split
  scan's response buffer is now sample-major (contiguous inner loop,
  1.4-1.9x), per-node sorted-index partitioning uses a per-tree arena
  instead of fresh allocations (up to 1.15x, more on wide blocks), and all
  per-node scratch buffers are hoisted out of the split search.
- On macOS, R's CRAN distribution bundles `libomp` but does not enable
  OpenMP for package compilation; installing with a user Makevars that sets
  `SHLIB_OPENMP_CXXFLAGS`/`SHLIB_OPENMP_LDFLAGS` (see `openmp.Makevars` in
  the project root) turns on the native engine's multithreaded tree
  building. Thread count does not affect results (verified identical for
  1/8/all threads); control it with `fit_forest(nthread = ...)` or
  `options(multiRF.nthread = ...)`.
- One parallel layer at a time: forests fitted inside forked workers
  (`parallel_connections`, `fit_sub_mrf(parallel = TRUE)`) now default to
  single-threaded OpenMP unless `multiRF.nthread` is set explicitly. This
  prevents core oversubscription and the fork-after-OpenMP hazard on
  Linux/libgomp; results are unchanged (verified identical).

## New: unified plot entry and pkgdown site

- New S3 method `plot(fit, type = ...)` for `mrf3_fit` objects: one entry
  point dispatching to all eight plot functions (`"tsne"` (default),
  `"umap"`, `"embed"`, `"network"`, `"circos"`, `"composition"`, `"km"`,
  `"weights"`), extracting the needed components from the fit and applying
  a consistent colour-blind-safe cluster palette across every plot type.
  All options of the underlying functions remain reachable via `...`.
- pkgdown site configured (`_pkgdown.yml`, grouped reference index, GitHub
  Actions workflow); `docs/` now holds the generated site.

## cluster_imd fast path (opt-in)

- `cluster_weighted_imd()` was repaired and vectorized, and `cluster_imd()`
  gained a `use_node_imd` argument (default `FALSE`): the node-score fast
  path is 119-883x faster than the per-tree depth traversal but measures a
  more global, split-score-weighted quantity that discriminates clusters
  less sharply (documented with measured agreement numbers), so the depth
  traversal remains the default estimator. `use_node_imd = NULL` picks the
  fast path automatically when models carry the needed statistics, with a
  message.

## Documentation accuracy

- The `method = "transformation"` selection rule is now described
  precisely (Eq. 16 t-score IMD): each feature's forest IMD is
  standardized against the block-wide mean forest IMD using the feature's
  across-tree standard error; the upper tail of a Student-t reference with
  `ntree - 1` degrees of freedom is kept at `level`, combined by strict
  majority vote across connected forests.

## Leaner internals

- Removed the unused exported alias `mrf3_full()` (use `mrf3_fit()`).
- `tcga_brca_clinical` pruned from 199 to 9 columns (IDs, PAM50 subtype,
  sample type, OS/PFI survival pairs, age, stage); `tcga_brca` unchanged.

- Removed ~600 lines of dead code identified by a call-graph audit
  (including the unused `em()`/`get_param()` mixture helpers and the
  `truncnorm` dependency, superseded connection-scoring helpers, orphaned
  tidy/fusion utilities, and a write-only C++ node field). No exported API
  was removed; `simulate_shared_specific_intersim()` now emits a formal
  deprecation warning.

## Other fixes and hardening

- Named `truth` vectors passed to `mrf3_stability()` are actually aligned by
  sample name (names were stripped before the check).
- `set.seed()` calls in `fit_sub_mrf()`, tuning, and `mrf3_stability()` no
  longer clobber the caller's global RNG state.
- Tuning with `parallel = TRUE` and `cores = NULL` now actually parallelizes.
- The same silent single-core bug is fixed in the IMD traversal
  (`get_imp_forest()`/`get_multi_weights()`), `cl_forest()`, and
  `fit_sub_mrf()`: with `parallel = TRUE`, `cores = NULL` now resolves to
  `parallel::detectCores() - 1` instead of 1. Results are identical; the
  documented parallelism just works now (measured 4-5x on the IMD stage).
- Exported documented API: `tune_model_top_v()`, `tune_fused_top_v()`,
  `get_reconstr_matrix()`, `get_shared_specific_weights()`,
  `simulate_intersim_shared_specific()`; registered `print.vs()` S3 method.
- Many small validation, error-message, and edge-case fixes across the
  clustering, tuning, and workflow layers (see the audit report).
- `DESCRIPTION`: added `parallel` to Imports; removed unused `tidyr`,
  `mclust`, `janitor`; `testthat (>= 3.1.7)`.
- Removed dead C++ fallback tree builders and stale in-tree build artifacts;
  the rendered vignette HTML moved to `docs/`.
- New regression-test file `test-review-fixes.R` locking in the fixes above.

# multiRF 0.2.2

- Added `plot_cluster_composition()` for auditable cluster-by-annotation
  composition heatmaps with counts, normalized percentages, and consistent
  package styling.

# multiRF 0.2.1

- **Fix**: Corrected the default `ytry` for supervised multivariate regression.

# multiRF 0.2.0

- Native C++ multivariate regression forest engine with OpenMP parallelism.
- Partition-based ranger-style tree building for faster splits.
- IMD importance computed during tree building (zero-cost).
- Unsupervised forest mode.

# multiRF 0.1.1

- Prepared package metadata for publishing (`DESCRIPTION`, license, and
  dependency declarations).
- Added a practical `README.md` with installation and usage instructions.
- Added `NEWS.md` to track user-facing changes.
- Hardened build settings to exclude development artifacts.
- Removed compiled artifacts from `src/` for clean source distribution.
