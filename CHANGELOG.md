# Changelog

## [Unreleased]

### Changed (Breaking)

- **`grid::Grid2D` storage layout** changed from x-major (`values[ix*ny+iy]`) to
  y-major (`values[iy*nx+ix]`) row-major order. Any code constructing a `Grid2D`
  with raw values or reading `Grid2D::values` directly must update the storage
  order.
- **`grid::Grid3D` storage layout** changed from `ix*(ny*nz)+iy*nz+iz` to
  `(iz*ny+iy)*nx+ix` (z-outermost, y-middle, x-innermost). Same caveat applies.
- **`grid::GridError` variants renamed**: `OutOfBounds` → `OutOfRange`;
  `axis` field type changed from `usize` to `&'static str` for clearer messages.
  `lo`/`hi` field names replace `min`/`max`.
- **`grid::Axis::validate_for_axis`** now takes `name: &'static str` instead of
  `axis: usize`.
- **`data::Provenance::cited`** now takes a single `source: &'static str` argument
  (previously two arguments).

### Added

- **`grid::AxisDirection`** is now `pub` and re-exported from the crate root, enabling
  callers to name it without importing the internal module.
- **`grid::algo` module** is now `pub`. Exposes direction-aware kernels:
  `validate_axis`, `linear_1d`, `bilinear`, `bilinear_unit`, `trilinear`,
  `trilinear_unit`.
- **`grid::Grid1D::interp_at_with`**, **`Grid2D::interp_at_with`**,
  **`Grid3D::interp_at_with`**: fallible versions that accept a pre-validated
  `ConstantRegion` hint to skip re-validation.
- **`grid::Grid2D::from_raw_row_major_y_descending`**: constructor for grids whose
  y-axis runs high-to-low (e.g. Leinert zodiacal data), storing internally with
  a reflection transform so queries use the normal ascending convention.
- **`grid::GridError::NonUniformStep`**: new variant signalled when a uniform axis
  has a non-uniform step after construction.
- **`grid::ConstantRegion`** re-exported from crate root.

- Added the initial `optica` crate with typed grids, spectra, participating-media
  building blocks, phase functions, and transport helpers.
