// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Participating-medium abstractions and optical coefficients.
//!
//! This module keeps the public API typed over geometry (`affn`) and wavelength
//! (`qtty`) while representing coefficient magnitudes as raw `f64` values in the
//! implied reciprocal path-length unit.

use core::marker::PhantomData;

use affn::{Position, ReferenceCenter, ReferenceFrame};
use qtty::dimensionless::Albedos;
use qtty::length::LengthUnit;
use qtty::{Dimensionless, Quantity, Unit};

/// Unit marker for absorption coefficient `[1 / length]`.
///
/// # Examples
///
/// ```rust
/// use optica::medium::AbsorptionCoeff;
/// use qtty::Quantity;
///
/// let sigma = Quantity::<AbsorptionCoeff>::new(0.1);
/// assert_eq!(sigma.value(), 0.1);
/// ```
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct AbsorptionCoeff;

impl Unit for AbsorptionCoeff {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "m⁻¹";
}

/// Unit marker for scattering coefficient `[1 / length]`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ScatteringCoeff;

impl Unit for ScatteringCoeff {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "m⁻¹";
}

/// Unit marker for extinction coefficient `[1 / length]`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ExtinctionCoeff;

impl Unit for ExtinctionCoeff {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "m⁻¹";
}

/// Per-wavelength optical coefficients at a point in a medium.
///
/// The type parameter `U` tracks the implied length unit so that `sigma_t.value()`
/// is interpreted in units of `1 / U`, consistent with path lengths expressed in `U`.
///
/// # Examples
///
/// ```rust
/// use optica::medium::OpticalCoefficients;
/// use qtty::unit::Kilometer;
///
/// let coeffs = OpticalCoefficients::<Kilometer>::new(0.1, 0.2);
/// assert!((coeffs.sigma_t.value() - 0.3).abs() < 1e-12);
/// assert!((coeffs.ssa.value() - (2.0 / 3.0)).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpticalCoefficients<U: LengthUnit> {
    /// Absorption coefficient `σ_a`.
    pub sigma_a: Quantity<AbsorptionCoeff>,
    /// Scattering coefficient `σ_s`.
    pub sigma_s: Quantity<ScatteringCoeff>,
    /// Extinction coefficient `σ_t = σ_a + σ_s`.
    pub sigma_t: Quantity<ExtinctionCoeff>,
    /// Single-scattering albedo `ω₀ = σ_s / σ_t`.
    pub ssa: Albedos,
    _unit: PhantomData<U>,
}

impl<U: LengthUnit> OpticalCoefficients<U> {
    /// Builds optical coefficients from absorption and scattering magnitudes.
    #[must_use]
    pub fn new(sigma_a: f64, sigma_s: f64) -> Self {
        let sigma_t = sigma_a + sigma_s;
        let ssa = if sigma_t == 0.0 { 0.0 } else { sigma_s / sigma_t };
        Self {
            sigma_a: Quantity::new(sigma_a),
            sigma_s: Quantity::new(sigma_s),
            sigma_t: Quantity::new(sigma_t),
            ssa: Albedos::new(ssa),
            _unit: PhantomData,
        }
    }

    /// Returns a fully transparent medium state.
    #[must_use]
    pub fn transparent() -> Self {
        Self::new(0.0, 0.0)
    }
}

/// A participating medium mapping position and wavelength to optical coefficients.
///
/// # Examples
///
/// ```rust
/// use affn::Position;
/// use affn::{ReferenceCenter, ReferenceFrame};
/// use optica::medium::{HomogeneousMedium, Medium};
/// use qtty::length::Nanometers;
/// use qtty::unit::Kilometer;
///
/// #[derive(Debug, Copy, Clone)]
/// struct Center;
/// impl ReferenceCenter for Center {
///     type Params = ();
///     fn center_name() -> &'static str { "Center" }
/// }
///
/// #[derive(Debug, Copy, Clone)]
/// struct Frame;
/// impl ReferenceFrame for Frame {
///     fn frame_name() -> &'static str { "Frame" }
/// }
///
/// let medium = HomogeneousMedium::<Kilometer>::new(0.1, 0.2);
/// let coeffs = medium.coefficients(
///     Position::<Center, Frame, Kilometer>::new(0.0, 0.0, 0.0),
///     Nanometers::new(550.0),
/// );
/// assert!((coeffs.sigma_t.value() - 0.3).abs() < 1e-12);
/// ```
pub trait Medium<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> {
    /// Returns optical coefficients at a position and wavelength.
    fn coefficients(&self, p: Position<C, F, U>, wavelength: qtty::length::Nanometers) -> OpticalCoefficients<U>;
}

/// Spatially uniform medium with wavelength-independent coefficients.
///
/// # Examples
///
/// ```rust
/// use optica::medium::HomogeneousMedium;
/// use qtty::unit::Kilometer;
///
/// let medium = HomogeneousMedium::<Kilometer>::new(0.1, 0.2);
/// assert_eq!(medium, HomogeneousMedium::<Kilometer>::new(0.1, 0.2));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HomogeneousMedium<U: LengthUnit> {
    sigma_a: f64,
    sigma_s: f64,
    _phantom: PhantomData<U>,
}

impl<U: LengthUnit> HomogeneousMedium<U> {
    /// Creates a homogeneous medium from constant absorption and scattering coefficients.
    #[must_use]
    pub fn new(sigma_a: f64, sigma_s: f64) -> Self {
        Self {
            sigma_a,
            sigma_s,
            _phantom: PhantomData,
        }
    }

    /// Creates a transparent homogeneous medium.
    #[must_use]
    pub fn transparent() -> Self {
        Self::new(0.0, 0.0)
    }
}

impl<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> Medium<C, F, U> for HomogeneousMedium<U> {
    fn coefficients(&self, _p: Position<C, F, U>, _wavelength: qtty::length::Nanometers) -> OpticalCoefficients<U> {
        OpticalCoefficients::new(self.sigma_a, self.sigma_s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transparent_coefficients_are_zero() {
        let coeffs = OpticalCoefficients::<qtty::unit::Kilometer>::transparent();
        assert_eq!(coeffs.sigma_t.value(), 0.0);
        assert_eq!(coeffs.ssa.value(), 0.0);
    }
}
