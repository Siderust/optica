// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Scattering phase functions for participating-media models.
//!
//! The built-in models here are wavelength-independent analytic kernels plus a
//! typed tabulated adapter for wavelength-dependent phase tables.
//!
//! ## References
//!
//! - Chandrasekhar, *Radiative Transfer*.
//! - Henyey & Greenstein (1941).

use qtty::angular::Radians;
use qtty::length::Nanometers;
use qtty::unit::{Nanometer, Radian};
use qtty::{Dimensionless, Quantity, Unit};

/// Dimensionless value marker for phase-function values.
///
/// # Examples
///
/// ```rust
/// use optica::phase::ScatteringFactor;
/// use qtty::Quantity;
///
/// let value = Quantity::<ScatteringFactor>::new(0.25);
/// assert_eq!(value.value(), 0.25);
/// ```
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ScatteringFactor;

impl Unit for ScatteringFactor {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}

/// A wavelength-dependent scattering phase function.
///
/// # Examples
///
/// ```rust
/// use optica::phase::{PhaseFunction, RayleighPhaseFunction};
/// use qtty::angular::Radians;
/// use qtty::length::Nanometers;
///
/// let model = RayleighPhaseFunction;
/// let value = model.phase(Nanometers::new(550.0), Radians::new(0.0));
/// assert!(value.value() > 0.0);
/// ```
pub trait PhaseFunction {
    /// Evaluates the phase function at a wavelength and scattering angle.
    fn phase(&self, wavelength: Nanometers, theta: Radians) -> Quantity<ScatteringFactor>;
}

/// Built-in phase models for allocation-free hot-loop use.
///
/// # Examples
///
/// ```rust
/// use optica::phase::{PhaseFunction, PhaseModel};
/// use qtty::angular::Radians;
/// use qtty::length::Nanometers;
///
/// let model = PhaseModel::HenyeyGreenstein { g: 0.3 };
/// let value = model.phase(Nanometers::new(550.0), Radians::new(0.2));
/// assert!(value.value().is_finite());
/// ```
#[derive(Debug, Clone, Copy)]
pub enum PhaseModel<'a> {
    /// The analytic Rayleigh phase function.
    Rayleigh,
    /// The analytic Henyey-Greenstein model.
    HenyeyGreenstein {
        /// Asymmetry parameter.
        g: f64,
    },
    /// Weighted sum of two Henyey-Greenstein lobes.
    DoubleHenyeyGreenstein {
        /// First asymmetry parameter.
        g1: f64,
        /// Second asymmetry parameter.
        g2: f64,
        /// Weight applied to the first lobe.
        weight: f64,
    },
    /// A tabulated wavelength-angle phase table.
    Tabulated(&'a PhaseTable),
}

impl<'a> PhaseFunction for PhaseModel<'a> {
    fn phase(&self, wavelength: Nanometers, theta: Radians) -> Quantity<ScatteringFactor> {
        match *self {
            Self::Rayleigh => RayleighPhaseFunction.phase(wavelength, theta),
            Self::HenyeyGreenstein { g } => HenyeyGreensteinPhaseFunction { g }.phase(wavelength, theta),
            Self::DoubleHenyeyGreenstein { g1, g2, weight } => {
                DoubleHenyeyGreensteinPhaseFunction { g1, g2, weight }.phase(wavelength, theta)
            }
            Self::Tabulated(table) => table.interp_at(wavelength, theta),
        }
    }
}

/// Phase table type: 2-D grid over wavelength × scattering angle.
pub type PhaseTable = crate::grid::Grid2D<Nanometer, Radian, ScatteringFactor>;

/// Rayleigh phase function: `P(θ) = (3 / 16π) (1 + cos²θ)`.
///
/// # Examples
///
/// ```rust
/// use optica::phase::{PhaseFunction, RayleighPhaseFunction};
/// use qtty::angular::Radians;
/// use qtty::length::Nanometers;
///
/// let value = RayleighPhaseFunction.phase(Nanometers::new(550.0), Radians::new(0.0));
/// assert!(value.value() > 0.0);
/// ```
#[derive(Debug, Clone, Copy, Default)]
pub struct RayleighPhaseFunction;

impl PhaseFunction for RayleighPhaseFunction {
    fn phase(&self, _wavelength: Nanometers, theta: Radians) -> Quantity<ScatteringFactor> {
        rayleigh_phase(theta)
    }
}

/// Henyey-Greenstein phase function.
///
/// # Examples
///
/// ```rust
/// use optica::phase::{HenyeyGreensteinPhaseFunction, PhaseFunction};
/// use qtty::angular::Radians;
/// use qtty::length::Nanometers;
///
/// let hg = HenyeyGreensteinPhaseFunction { g: 0.2 };
/// let value = hg.phase(Nanometers::new(550.0), Radians::new(0.5));
/// assert!(value.value() > 0.0);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct HenyeyGreensteinPhaseFunction {
    /// Asymmetry parameter.
    pub g: f64,
}

impl PhaseFunction for HenyeyGreensteinPhaseFunction {
    fn phase(&self, _wavelength: Nanometers, theta: Radians) -> Quantity<ScatteringFactor> {
        let cos_theta = theta.value().cos();
        let g2 = self.g * self.g;
        let denom = 1.0 + g2 - 2.0 * self.g * cos_theta;
        Quantity::new((1.0 - g2) / (4.0 * core::f64::consts::PI * denom.powf(1.5)))
    }
}

/// Double Henyey-Greenstein phase function.
///
/// # Examples
///
/// ```rust
/// use optica::phase::{DoubleHenyeyGreensteinPhaseFunction, PhaseFunction};
/// use qtty::angular::Radians;
/// use qtty::length::Nanometers;
///
/// let model = DoubleHenyeyGreensteinPhaseFunction {
///     g1: 0.7,
///     g2: -0.3,
///     weight: 0.8,
/// };
/// let value = model.phase(Nanometers::new(550.0), Radians::new(0.5));
/// assert!(value.value() > 0.0);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct DoubleHenyeyGreensteinPhaseFunction {
    /// First asymmetry parameter.
    pub g1: f64,
    /// Second asymmetry parameter.
    pub g2: f64,
    /// Weight applied to the first lobe.
    pub weight: f64,
}

impl PhaseFunction for DoubleHenyeyGreensteinPhaseFunction {
    fn phase(&self, wavelength: Nanometers, theta: Radians) -> Quantity<ScatteringFactor> {
        let hg1 = HenyeyGreensteinPhaseFunction { g: self.g1 }.phase(wavelength, theta);
        let hg2 = HenyeyGreensteinPhaseFunction { g: self.g2 }.phase(wavelength, theta);
        Quantity::new(self.weight * hg1.value() + (1.0 - self.weight) * hg2.value())
    }
}

/// Computes the Rayleigh phase function directly.
///
/// # Examples
///
/// ```rust
/// use optica::phase::rayleigh_phase;
/// use qtty::angular::Radians;
///
/// let value = rayleigh_phase(Radians::new(0.0));
/// assert!(value.value() > 0.0);
/// ```
#[must_use]
pub fn rayleigh_phase(theta: Radians) -> Quantity<ScatteringFactor> {
    let cosine = theta.value().cos();
    Quantity::<ScatteringFactor>::new(3.0 / (16.0 * core::f64::consts::PI) * (1.0 + cosine * cosine))
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn rayleigh_matches_closed_form() {
        let theta = Radians::new(0.0);
        let expected = 3.0 / (8.0 * core::f64::consts::PI);
        assert_relative_eq!(rayleigh_phase(theta).value(), expected, epsilon = 1.0e-12);
    }

    #[test]
    fn isotropic_hg_matches_one_over_four_pi() {
        let hg = HenyeyGreensteinPhaseFunction { g: 0.0 };
        let value = hg.phase(Nanometers::new(550.0), Radians::new(1.0));
        assert_relative_eq!(value.value(), 1.0 / (4.0 * core::f64::consts::PI), epsilon = 1.0e-12);
    }
}
