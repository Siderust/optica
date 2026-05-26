// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Scattering phase functions for participating-media models.
//!
//! Analytic kernels validate their parameters via `try_new` constructors that
//! return a [`PhaseError`] on invalid input.
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

/// Errors produced when validating phase-function parameters.
///
/// # Examples
///
/// ```rust
/// use optica::phase::PhaseError;
///
/// let err = PhaseError::AsymmetryOutOfRange { field: "g", value: 1.5 };
/// assert!(err.to_string().contains("must lie in (-1, 1)"));
/// ```
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PhaseError {
    /// Parameter was non-finite (NaN or ±∞).
    #[error("{field} must be finite (got {value})")]
    NonFinite {
        /// Offending field name.
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
    /// Asymmetry parameter outside `(-1, 1)`.
    ///
    /// The values `g = ±1` are not accepted because the Henyey-Greenstein
    /// formula degenerates to a Dirac delta at those limits.
    #[error("{field} must lie in (-1, 1) (got {value})")]
    AsymmetryOutOfRange {
        /// Offending field name.
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
    /// Mixing weight outside `[0, 1]`.
    #[error("{field} must lie in [0, 1] (got {value})")]
    WeightOutOfRange {
        /// Offending field name.
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
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
/// The variants carry already-validated kernels constructed via the matching
/// `try_new` helpers, e.g. [`HenyeyGreensteinPhaseFunction::try_new`].
///
/// # Examples
///
/// ```rust
/// use optica::phase::{HenyeyGreensteinPhaseFunction, PhaseFunction, PhaseModel};
/// use qtty::angular::Radians;
/// use qtty::length::Nanometers;
///
/// let hg = HenyeyGreensteinPhaseFunction::try_new(0.3).unwrap();
/// let model = PhaseModel::HenyeyGreenstein(hg);
/// let value = model.phase(Nanometers::new(550.0), Radians::new(0.2));
/// assert!(value.value().is_finite());
/// ```
#[derive(Debug, Clone, Copy)]
pub enum PhaseModel<'a> {
    /// The analytic Rayleigh phase function.
    Rayleigh,
    /// The analytic Henyey-Greenstein model.
    HenyeyGreenstein(HenyeyGreensteinPhaseFunction),
    /// Weighted sum of two Henyey-Greenstein lobes.
    DoubleHenyeyGreenstein(DoubleHenyeyGreensteinPhaseFunction),
    /// A tabulated wavelength-angle phase table.
    Tabulated(&'a PhaseTable),
}

impl<'a> PhaseFunction for PhaseModel<'a> {
    fn phase(&self, wavelength: Nanometers, theta: Radians) -> Quantity<ScatteringFactor> {
        match self {
            Self::Rayleigh => RayleighPhaseFunction.phase(wavelength, theta),
            Self::HenyeyGreenstein(hg) => hg.phase(wavelength, theta),
            Self::DoubleHenyeyGreenstein(dhg) => dhg.phase(wavelength, theta),
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

/// Henyey-Greenstein phase function with a validated asymmetry parameter.
///
/// # Examples
///
/// ```rust
/// use optica::phase::{HenyeyGreensteinPhaseFunction, PhaseFunction};
/// use qtty::angular::Radians;
/// use qtty::length::Nanometers;
///
/// let hg = HenyeyGreensteinPhaseFunction::try_new(0.2).unwrap();
/// let value = hg.phase(Nanometers::new(550.0), Radians::new(0.5));
/// assert!(value.value() > 0.0);
/// assert!(HenyeyGreensteinPhaseFunction::try_new(1.5).is_err());
/// assert!(HenyeyGreensteinPhaseFunction::try_new(1.0).is_err());  // g=1 degenerates
/// ```
#[derive(Debug, Clone, Copy)]
pub struct HenyeyGreensteinPhaseFunction {
    g: f64,
}

impl HenyeyGreensteinPhaseFunction {
    /// Constructs a Henyey-Greenstein kernel after validating `g ∈ (-1, 1)`.
    ///
    /// The values `g = ±1` are rejected because at those limits the HG formula
    /// `(1 − g²) / (1 + g² − 2g cos θ)^(3/2)` collapses to a Dirac delta.
    ///
    /// # Errors
    ///
    /// Returns [`PhaseError`] when `g` is non-finite or outside `(-1, 1)`.
    pub fn try_new(g: f64) -> Result<Self, PhaseError> {
        check_asymmetry("g", g)?;
        Ok(Self { g })
    }

    /// Returns the asymmetry parameter `g`.
    #[must_use]
    pub fn g(&self) -> f64 {
        self.g
    }
}

impl PhaseFunction for HenyeyGreensteinPhaseFunction {
    fn phase(&self, _wavelength: Nanometers, theta: Radians) -> Quantity<ScatteringFactor> {
        let cos_theta = theta.value().cos();
        let g2 = self.g * self.g;
        let denom = 1.0 + g2 - 2.0 * self.g * cos_theta;
        Quantity::new((1.0 - g2) / (4.0 * core::f64::consts::PI * denom.powf(1.5)))
    }
}

/// Double Henyey-Greenstein phase function with validated parameters.
///
/// # Examples
///
/// ```rust
/// use optica::phase::{DoubleHenyeyGreensteinPhaseFunction, PhaseFunction};
/// use qtty::angular::Radians;
/// use qtty::length::Nanometers;
///
/// let model = DoubleHenyeyGreensteinPhaseFunction::try_new(0.7, -0.3, 0.8).unwrap();
/// let value = model.phase(Nanometers::new(550.0), Radians::new(0.5));
/// assert!(value.value() > 0.0);
/// assert!(DoubleHenyeyGreensteinPhaseFunction::try_new(0.7, -0.3, 1.2).is_err());
/// ```
#[derive(Debug, Clone, Copy)]
pub struct DoubleHenyeyGreensteinPhaseFunction {
    g1: f64,
    g2: f64,
    weight: f64,
}

impl DoubleHenyeyGreensteinPhaseFunction {
    /// Constructs a validated double-HG kernel.
    ///
    /// # Errors
    ///
    /// Returns [`PhaseError`] when either asymmetry lies outside `(-1, 1)` or
    /// the weight lies outside `[0, 1]`. Non-finite inputs are also rejected.
    pub fn try_new(g1: f64, g2: f64, weight: f64) -> Result<Self, PhaseError> {
        check_asymmetry("g1", g1)?;
        check_asymmetry("g2", g2)?;
        check_weight("weight", weight)?;
        Ok(Self { g1, g2, weight })
    }

    /// Returns the first lobe asymmetry parameter.
    #[must_use]
    pub fn g1(&self) -> f64 {
        self.g1
    }

    /// Returns the second lobe asymmetry parameter.
    #[must_use]
    pub fn g2(&self) -> f64 {
        self.g2
    }

    /// Returns the mixing weight applied to the first lobe.
    #[must_use]
    pub fn weight(&self) -> f64 {
        self.weight
    }
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
    Quantity::<ScatteringFactor>::new(
        3.0 / (16.0 * core::f64::consts::PI) * (1.0 + cosine * cosine),
    )
}

fn check_asymmetry(field: &'static str, value: f64) -> Result<(), PhaseError> {
    if !value.is_finite() {
        return Err(PhaseError::NonFinite { field, value });
    }
    if value <= -1.0 || value >= 1.0 {
        return Err(PhaseError::AsymmetryOutOfRange { field, value });
    }
    Ok(())
}

fn check_weight(field: &'static str, value: f64) -> Result<(), PhaseError> {
    if !value.is_finite() {
        return Err(PhaseError::NonFinite { field, value });
    }
    if !(0.0..=1.0).contains(&value) {
        return Err(PhaseError::WeightOutOfRange { field, value });
    }
    Ok(())
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
        let hg = HenyeyGreensteinPhaseFunction::try_new(0.0).unwrap();
        let value = hg.phase(Nanometers::new(550.0), Radians::new(1.0));
        assert_relative_eq!(
            value.value(),
            1.0 / (4.0 * core::f64::consts::PI),
            epsilon = 1.0e-12
        );
    }

    #[test]
    fn hg_rejects_out_of_range_g() {
        assert!(matches!(
            HenyeyGreensteinPhaseFunction::try_new(1.5),
            Err(PhaseError::AsymmetryOutOfRange { .. })
        ));
        assert!(matches!(
            HenyeyGreensteinPhaseFunction::try_new(1.0),
            Err(PhaseError::AsymmetryOutOfRange { .. })
        ));
        assert!(matches!(
            HenyeyGreensteinPhaseFunction::try_new(-1.0),
            Err(PhaseError::AsymmetryOutOfRange { .. })
        ));
        assert!(matches!(
            HenyeyGreensteinPhaseFunction::try_new(f64::NAN),
            Err(PhaseError::NonFinite { .. })
        ));
    }

    #[test]
    fn double_hg_rejects_out_of_range_weight() {
        assert!(matches!(
            DoubleHenyeyGreensteinPhaseFunction::try_new(0.5, -0.5, 2.0),
            Err(PhaseError::WeightOutOfRange { .. })
        ));
    }
}
