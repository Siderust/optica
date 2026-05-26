# optica

[![Crates.io](https://img.shields.io/crates/v/optica)](https://crates.io/crates/optica)
[![docs.rs](https://img.shields.io/docsrs/optica)](https://docs.rs/optica)
[![CI](https://github.com/Siderust/optica/actions/workflows/ci.yml/badge.svg)](https://github.com/Siderust/optica/actions/workflows/ci.yml)
[![License: AGPL-3.0-only](https://img.shields.io/badge/license-AGPL--3.0--only-blue)](LICENSE)

**Fast, typed participating-media and optics foundations for Rust.**

`optica` provides the core building blocks for radiative-transfer and
optical-depth workloads: typed rays, optical coefficients, phase functions,
sampled spectra, interpolation tables, and Beer–Lambert integration.
It is intentionally domain-agnostic — no astronomy policies, planetary
constants, ephemerides, or observatory presets.

> **Pre-1.0 API:** Breaking changes may occur in minor releases.
> Each breaking change is documented in [CHANGELOG.md](CHANGELOG.md).

## Scope

- Typed rays and ray segments (affn geometry, phantom length unit)
- Optical coefficients: absorption `σ_a`, scattering `σ_s`, extinction `σ_t`, SSA `ω₀`
- Scattering phase functions: Rayleigh, Henyey–Greenstein, double-HG, tabulated
- Rayleigh and Ångström (Mie) optical-depth formulas (caller-supplied parameters)
- Sampled spectra: 1-D tables with linear, nearest, step, and cubic-spline interpolation
- Two-column ASCII spectrum loader (whitespace or CSV, unit-scale factors, provenance)
- 1-D, 2-D, and 3-D typed interpolation grids (ascending and descending axes)
- Optical-depth integration over a ray segment: midpoint, trapezoidal, Simpson, and Gauss-Legendre rules
- Beer–Lambert transmittance and van Rhijn path-length factor
- Provenance metadata for tabulated inputs and generated products

## Non-goals

- Astronomical bodies, ephemerides, or site-specific observatory constants
- Atmospheric profiles or specific planetary atmospheric models
- Radiometric source models (solar spectra, thermal emission)
- GPU or SIMD acceleration (the crate prioritises correctness and `no_std` support)

## Installation

```toml
[dependencies]
optica = "0.1"
```

Default features include `std`. Opt out for `no_std` + alloc:

```toml
[dependencies]
optica = { version = "0.1", default-features = false, features = ["alloc"] }
```

Enable Serde serialization:

```toml
[dependencies]
optica = { version = "0.1", features = ["serde"] }
```

## Feature flags

| Feature | Default | Description |
|---------|---------|-------------|
| `std`   | ✓       | Enables `std`-dependent helpers (ASCII/string parsing) and implies `alloc`. |
| `alloc` |         | Enables heap-backed types (`Vec`/`Box`/`String`) for `no_std` targets. With this off, only `medium`, `ray`, `scatter`, and `transport` are compiled. |
| `serde` |         | Derives `Serialize`/`Deserialize` on the public data, error, and policy types. `TableSource` derives only `Serialize` because it borrows static slices. |
| `astro` |         | Reserved for future astronomy-specific adapters; currently a no-op. |

## Serde support

When the `serde` feature is enabled, the following types derive `Serialize` and `Deserialize`:
`AxisDirection`, `OutOfRange`, `Provenance`, `TableSource` (`Serialize` only),
`Interpolation`, `SpectrumError`, `OpticalCoefficientError`, `PhaseError`,
`ScatterError`, `TransportError`, `IntegrationMethod`, `IntegrationOpts`,
`MieParams`.

The following types **do not** derive serde:
`SampledSpectrum`, `Grid1D`, `Grid2D`, `Grid3D`, `HomogeneousMedium`,
`HenyeyGreensteinPhaseFunction`, `DoubleHenyeyGreensteinPhaseFunction`,
`RayleighPhaseFunction`, `PhaseModel`, `PhaseTable`.
These containers hold variable-length data or opaque function pointers; their
serialization format is left to the caller.

## Examples

### Load a two-column ASCII spectrum

```rust
use optica::spectrum::{Interpolation, loaders::ascii};
use optica::grid::OutOfRange;
use qtty::unit::{Nanometer, Ratio};

let data = "# wavelength(nm) transmittance\n400.0 0.0\n550.0 0.85\n700.0 0.92\n";
let spectrum = ascii::two_column::<Nanometer, Ratio>(
    data, 1.0, 1.0,
    Interpolation::Linear,
    OutOfRange::ClampToEndpoints,
    None,
).unwrap();

use qtty::Quantity;
let t = spectrum.interp_at(Quantity::<Nanometer>::new(550.0));
assert!((t.value() - 0.85).abs() < 1e-12);
```

### Interpolate a 2-D table

```rust
use optica::grid::{Grid2D, OutOfRange};
use qtty::{Quantity, unit::{Nanometer, Radian, Ratio}};

// Row-major: row 0 (y=0°): V(400nm)=0.1, V(700nm)=0.4
//            row 1 (y=1°): V(400nm)=0.5, V(700nm)=0.9
let grid = Grid2D::<Nanometer, Radian, Ratio>::from_raw_row_major(
    &[400.0, 700.0],
    &[0.0, 1.0],
    &[0.1, 0.4, 0.5, 0.9],
    OutOfRange::ClampToEndpoints,
).unwrap();

let v = grid.interp_at(Quantity::new(550.0), Quantity::new(0.5));
assert!((v.value() - 0.475).abs() < 1e-12);
```

### Integrate optical depth along a ray

```rust
use affn::{CartesianDirection, Position, ReferenceCenter, ReferenceFrame};
use optica::medium::HomogeneousMedium;
use optica::ray::{Ray, RaySegment};
use optica::transport::{integrate_optical_depth, IntegrationMethod, IntegrationOpts};
use qtty::length::{Kilometers, Nanometers};
use qtty::unit::Kilometer;

#[derive(Debug, Copy, Clone)]
struct Origin;
impl ReferenceCenter for Origin {
    type Params = ();
    fn center_name() -> &'static str { "Origin" }
}

#[derive(Debug, Copy, Clone)]
struct LocalFrame;
impl ReferenceFrame for LocalFrame {
    fn frame_name() -> &'static str { "Local" }
}

// σ_a = 0.1 km⁻¹, σ_s = 0.2 km⁻¹ → σ_t = 0.3 km⁻¹
let medium = HomogeneousMedium::<Kilometer>::try_new(0.1, 0.2).unwrap();
let ray = Ray::new(
    Position::<Origin, LocalFrame, Kilometer>::new(0.0, 0.0, 0.0),
    CartesianDirection::<LocalFrame>::new(0.0, 0.0, 1.0),
);
let tau = integrate_optical_depth(
    &medium,
    &ray,
    RaySegment::new(Kilometers::new(0.0), Kilometers::new(10.0)),
    Nanometers::new(550.0),
    IntegrationOpts::new(32, IntegrationMethod::Midpoint),
);
// τ = σ_t × distance = 0.3 × 10 = 3.0
assert!((tau.value() - 3.0).abs() < 1e-12);
```

## Relationship with `siderust`

[`siderust`](https://crates.io/crates/siderust) builds domain-specific
astronomy and atmospheric science on top of `optica`. It adds:

- Concrete atmospheric profiles (ozone, Rayleigh, aerosol models)
- Observatory and site-specific presets
- Ephemerides and body constants
- Astronomical epoch and coordinate types

`optica` is the generic foundation; `siderust` is the astronomy adapter.
Users who only need participating-media primitives can depend on `optica`
directly without pulling in `siderust`.

## License

AGPL-3.0-only — see [LICENSE](LICENSE).
