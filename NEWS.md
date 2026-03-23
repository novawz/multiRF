# multiRF 0.2.1

- **Fix**: Changed default `ytry` for supervised multivariate regression from
  `min(qy, ceiling(px/3))` to `ceiling(qy/3)`. The previous default
  incorrectly used `px` (X columns) instead of `qy` (Y columns) to size the
  random Y subset, which could produce noisier splits when `px ≠ qy`.

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
