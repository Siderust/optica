// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Sampled-spectrum container, interpolation, integration, and loaders.
//!
//! ## Scientific scope
//!
//! Sampled spectra represent continuous optical functions at discrete
//! wavelength (or frequency) samples: passbands, source spectra, atmospheric
//! transmittance curves, and phase function tables.
//!
//! ## Technical scope
//!
//! [`SampledSpectrum<X, Y>`] stores axis and value data as `Box<[f64]>` slices
//! for allocation-free evaluation. Five interpolation modes are supported via
//! [`Interpolation`]. Typed trapezoidal integration is available via the
//! methods on `SampledSpectrum`; raw kernels are re-exported from [`integrate`].
//!
//! # References
//!
//! - Atkinson, K. E. (1989). *An Introduction to Numerical Analysis*,
//!   2nd ed., §5.2 (composite trapezoidal rule). John Wiley & Sons.

pub mod algo;
pub mod error;
pub mod integrate;
pub mod loaders;
mod sampled;

pub use crate::grid::OutOfRange;
pub use error::SpectrumError;
pub use sampled::SampledSpectrum;

/// Interpolation mode for sampled spectra.
///
/// # Examples
///
/// ```rust
/// use optica::spectrum::Interpolation;
///
/// let mode = Interpolation::Linear;
/// assert!(matches!(mode, Interpolation::Linear));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interpolation {
    /// Piecewise-linear interpolation.
    ///
    /// Matches `numpy.interp` exactly when out-of-range is
    /// [`OutOfRange::ClampToEndpoints`].
    Linear,
    /// Nearest-neighbour interpolation.
    ///
    /// Ties (exact midpoint) resolve to the lower index, matching
    /// `scipy.interpolate.interp1d(kind="nearest")`.
    Nearest,
    /// Left-continuous piecewise-constant ("previous") step interpolation.
    PiecewiseConstantLeft,
    /// Right-continuous piecewise-constant ("next") step interpolation.
    PiecewiseConstantRight,
    /// Natural cubic spline.
    ///
    /// Coefficients are precomputed once at construction time.
    CubicSpline,
}
