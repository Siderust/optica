# Changelog

All notable changes to `optica` are documented here.
Breaking changes are marked **Breaking** and each carries a migration note.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/) starting at 1.0;
pre-1.0 minor releases may contain breaking changes.

## [0.2.0] - 2026-06-21

### Added

- `#![cfg_attr(not(feature = "std"), no_std)]`: the crate now builds without `std`.
  Heap-using modules (`data`, `grid`, `phase`, `prelude`, `spectrum`) are gated
  behind the `alloc` feature; `medium`, `ray`, `scatter`, `transport` are always
  available.
- `serde` feature: derives `Serialize`/`Deserialize` on public data, error, and
  policy types (`Provenance`, `DataSource`, `Axis`, `AxisDirection`,
  `ConstantRegion`, `OutOfRange`, `GridError`, `Interpolation`, `SpectrumError`,
  `MieParams`, `ScatterError`, `OpticalCoefficientError`, `PhaseError`,
  `IntegrationOpts`, `IntegrationMethod`, `TransportError`). `TableSource`
  derives only `Serialize` (it borrows `&'static` data).
- `medium::InverseLength<U>` type alias (`Per<Ratio, U>`) for typed reciprocal
  lengths used by optical coefficients.
- `medium::OpticalCoefficientError` and `OpticalCoefficients::try_new`,
  `HomogeneousMedium::try_new`: fallible, validated constructors.
- `phase::PhaseError`, validated `HenyeyGreensteinPhaseFunction::try_new` and
  `DoubleHenyeyGreensteinPhaseFunction::try_new`.
- `scatter::ScatterError`, `MieParams::try_new`,
  `scatter::try_rayleigh_optical_depth_bodhaine99`,
  `scatter::try_mie_optical_depth`.
- `transport::IntegrationMethod` enum (`Midpoint`, `Trapezoidal`, `Simpson`,
  `GaussLegendre2`, `GaussLegendre4`) and `transport::TransportError`.
  `IntegrationOpts` now carries a `method` field. New
  `transport::try_integrate_optical_depth` validates inputs.
- `spectrum::SampledSpectrum`: `domain`, `contains`, `overlap_domain`,
  `resample_onto` (typed, takes `&[Quantity<X>]`), `resample_onto_raw` (raw
  escape hatch), `normalize_area` (returns `SampledSpectrum<X, Ratio>`),
  `normalize_area_raw` (raw escape hatch), `normalize_peak` (returns
  `SampledSpectrum<X, Ratio>`), `normalize_peak_raw` (raw escape hatch),
  `map_values_to` (typed closure), `map_values_raw` (raw escape hatch),
  `zip_with_to` (typed closure), `zip_with_raw` (raw escape hatch).
- `medium::TryMedium<C, F, U>` trait for fallible/tabulated media.
- `HomogeneousMedium::sigma_a()` and `sigma_s()` typed getters returning
  `Quantity<InverseLength<U>>`.
- `grid::Axis::bounds`, `Grid2D::x_bounds`/`y_bounds`/`domain`,
  `Grid3D::x_bounds`/`y_bounds`/`z_bounds`/`domain`.
- `spectrum::SpectrumError::InvalidValue` variant. Enum marked `#[non_exhaustive]`.
- `tests/serde_roundtrip.rs`: JSON round-trip smoke tests for serialized types.

### Changed (Breaking)

- `HomogeneousMedium::new`, `OpticalCoefficients::new`,
  `HenyeyGreensteinPhaseFunction { g: … }` literal construction,
  `DoubleHenyeyGreensteinPhaseFunction { … }` literal construction,
  `MieParams::new`, `rayleigh_optical_depth_bodhaine99`, and `mie_optical_depth`
  are removed. Use the corresponding `try_new`/`try_*` variants.
- `PhaseModel::HenyeyGreenstein { g }` enum variant is now
  `PhaseModel::HenyeyGreenstein(HenyeyGreensteinPhaseFunction)`; the wrapped
  function carries its own validated asymmetry parameter.
- The fake unit markers `AbsorptionCoeff`, `ScatteringCoeff`, and
  `ExtinctionCoeff` are removed. Coefficients now use the real reciprocal-length
  unit `InverseLength<U> = Per<Ratio, U>`.
- `IntegrationOpts { n_steps }` literal construction must now supply `method`;
  use `IntegrationOpts::default()` or `IntegrationOpts::new(n_steps, method)`.
- `SpectrumError` is now `#[non_exhaustive]`.
- `HenyeyGreensteinPhaseFunction::try_new` and
  `DoubleHenyeyGreensteinPhaseFunction::try_new` now reject `g = ±1` (open
  interval `(-1, 1)` strictly). The HG formula degenerates to a Dirac delta at
  those limits.
- `SampledSpectrum::normalize_area` returns `SampledSpectrum<X, Ratio>` (not
  `Self`); `normalize_peak` likewise.
- `SampledSpectrum::map_values` renamed to `map_values_raw`; `zip_with` renamed
  to `zip_with_raw`; `resample_onto(&[f64])` renamed to `resample_onto_raw`.
- `HomogeneousMedium` internal fields now stored as `Quantity<InverseLength<U>>`.

## [0.1.0] — 2026-05-26

Initial release preparation.

### Added

- Added the initial `optica` crate with typed grids, spectra, participating-media
  building blocks, phase functions, and transport helpers.
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
- Three runnable examples under `examples/`: `spectrum_ascii`, `grid2d_interp`,
  `optical_depth`.
- GitHub CI workflow (`.github/workflows/ci.yml`) and publish workflow
  (`.github/workflows/publish.yml`).

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

[Unreleased]: https://github.com/Siderust/optica/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Siderust/optica/releases/tag/v0.1.0
