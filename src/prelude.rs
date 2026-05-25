// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Convenient re-exports for common `optica` workflows.
//!
//! # Examples
//!
//! ```rust
//! use optica::prelude::*;
//! use qtty::length::Nanometers;
//!
//! let spectrum = SampledSpectrum::<qtty::unit::Nanometer, qtty::unit::Ratio>::from_sorted(
//!     &[450.0, 550.0],
//!     &[0.2, 0.4],
//!     Interpolation::Linear,
//!     OutOfRange::ClampToEndpoints,
//! )
//! .unwrap();
//!
//! let sample = spectrum.interp_at(Nanometers::new(500.0));
//! assert!((sample.value() - 0.3).abs() < 1e-12);
//! ```

pub use crate::data::{DataSource, Provenance, TableSource};
pub use crate::grid::{Axis, Grid1D, Grid2D, Grid3D, GridError, OutOfRange};
pub use crate::medium::{
    AbsorptionCoeff, ExtinctionCoeff, HomogeneousMedium, Medium, OpticalCoefficients,
    ScatteringCoeff,
};
pub use crate::phase::{
    rayleigh_phase, DoubleHenyeyGreensteinPhaseFunction, HenyeyGreensteinPhaseFunction,
    PhaseFunction, PhaseModel, PhaseTable, RayleighPhaseFunction, ScatteringFactor,
};
pub use crate::ray::{Ray, RaySegment};
pub use crate::scatter::{mie_optical_depth, rayleigh_optical_depth_bodhaine99, MieParams};
pub use crate::spectrum::{Interpolation, SampledSpectrum};
pub use crate::transport::{
    integrate_optical_depth, transmittance, van_rhijn_factor, IntegrationOpts,
};
