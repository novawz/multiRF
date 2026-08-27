# multiRF 0.2.3

- Added the unified `plot(fit, type = ...)` interface and pkgdown site.
- Fixed clustering, OOB weights, IMD, stability, and related edge cases.
- Improved C++ forest speed, parallel safety, and cross-platform
  reproducibility; seeded unsupervised results change once from 0.2.2.
- Added an optional fast `cluster_imd()` path and simplified internals,
  dependencies, tests, and documentation.
- Relicensed under GPL (>= 3) with RF-SRC provenance and citations.

# multiRF 0.2.2

- Added `plot_cluster_composition()` for auditable cluster-by-annotation
  composition heatmaps with counts, normalized percentages, and consistent
  package styling.

# multiRF 0.2.1

- **Fix**: Corrected the default `ytry` for supervised multivariate regression.

# multiRF 0.2.0

- C++ multivariate regression forest engine with OpenMP parallelism.
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
