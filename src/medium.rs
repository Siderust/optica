// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Participating-medium abstractions and optical coefficients.
//!
//! Coefficients are typed as proper reciprocal-length quantities via
//! [`InverseLength<U>`], so a coefficient times a length integrates dimensionally
//! into a dimensionless optical depth automatically.

use core::marker::PhantomData;

use affn::{Position, ReferenceCenter, ReferenceFrame};
use qtty::dimensionless::{Albedos, Ratio};
use qtty::length::LengthUnit;
use qtty::unit::Per;
use qtty::Quantity;

/// Reciprocal-length unit alias parameterised by the chosen length unit `U`.
///
/// `Quantity<InverseLength<U>> · Quantity<U>` reduces to a dimensionless
/// quantity, matching the physical convention that `σ · ℓ = τ`.
///
/// # Examples
///
/// ```rust
/// use optica::medium::InverseLength;
/// use qtty::unit::Kilometer;
/// use qtty::Quantity;
///
/// let sigma = Quantity::<InverseLength<Kilometer>>::new(0.05);
/// assert!((sigma.value() - 0.05).abs() < 1e-12);
/// ```
pub type InverseLength<U> = Per<Ratio, U>;

/// Errors produced when validating optical-coefficient inputs.
///
/// # Examples
///
/// ```rust
/// use optica::medium::OpticalCoefficientError;
///
/// let error = OpticalCoefficientError::Negative {
///     field: "sigma_a",
///     value: -1.0,
/// };
/// assert!(error.to_string().contains("must be non-negative"));
/// ```
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum OpticalCoefficientError {
    /// A coefficient was non-finite (NaN or ±∞).
    #[error("{field} must be finite (got {value})")]
    NonFinite {
        /// Name of the offending field.
        field: &'static str,
        /// The supplied value.
        value: f64,
    },
    /// A coefficient was negative.
    #[error("{field} must be non-negative (got {value})")]
    Negative {
        /// Name of the offending field.
        field: &'static str,
        /// The supplied value.
        value: f64,
    },
}

/// Per-wavelength optical coefficients at a point in a medium.
///
/// The type parameter `U` is the path-length unit. All fields are typed as
/// [`Quantity<InverseLength<U>>`] so that `sigma_t * ds` is dimensionally a
/// dimensionless optical depth.
///
/// # Examples
///
/// ```rust
/// use optica::medium::OpticalCoefficients;
/// use qtty::unit::Kilometer;
///
/// let coeffs = OpticalCoefficients::<Kilometer>::try_new(0.1, 0.2).unwrap();
/// assert!((coeffs.sigma_t.value() - 0.3).abs() < 1e-12);
/// assert!((coeffs.ssa.value() - (2.0 / 3.0)).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpticalCoefficients<U: LengthUnit> {
    /// Absorption coefficient `σ_a` in `1 / U`.
    pub sigma_a: Quantity<InverseLength<U>>,
    /// Scattering coefficient `σ_s` in `1 / U`.
    pub sigma_s: Quantity<InverseLength<U>>,
    /// Extinction coefficient `σ_t = σ_a + σ_s` in `1 / U`.
    pub sigma_t: Quantity<InverseLength<U>>,
    /// Single-scattering albedo `ω₀ = σ_s / σ_t`.
    pub ssa: Albedos,
    _unit: PhantomData<U>,
}

impl<U: LengthUnit> OpticalCoefficients<U> {
    /// Constructs validated optical coefficients from absorption and scattering magnitudes.
    ///
    /// # Errors
    ///
    /// Returns [`OpticalCoefficientError`] when either input is non-finite or negative.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::medium::OpticalCoefficients;
    /// use qtty::unit::Kilometer;
    ///
    /// let ok = OpticalCoefficients::<Kilometer>::try_new(0.0, 0.0).unwrap();
    /// assert_eq!(ok.sigma_t.value(), 0.0);
    /// assert!(OpticalCoefficients::<Kilometer>::try_new(-0.1, 0.0).is_err());
    /// ```
    pub fn try_new(sigma_a: f64, sigma_s: f64) -> Result<Self, OpticalCoefficientError> {
        check_coeff("sigma_a", sigma_a)?;
        check_coeff("sigma_s", sigma_s)?;
        let sigma_t = sigma_a + sigma_s;
        let ssa = if sigma_t == 0.0 {
            0.0
        } else {
            sigma_s / sigma_t
        };
        Ok(Self {
            sigma_a: Quantity::new(sigma_a),
            sigma_s: Quantity::new(sigma_s),
            sigma_t: Quantity::new(sigma_t),
            ssa: Albedos::new(ssa),
            _unit: PhantomData,
        })
    }

    /// Returns the fully transparent state (`σ_a = σ_s = 0`).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::medium::OpticalCoefficients;
    /// use qtty::unit::Kilometer;
    ///
    /// let coeffs = OpticalCoefficients::<Kilometer>::transparent();
    /// assert_eq!(coeffs.sigma_t.value(), 0.0);
    /// ```
    #[must_use]
    pub fn transparent() -> Self {
        // Safe by construction.
        Self::try_new(0.0, 0.0).expect("zero coefficients are always valid")
    }
}

fn check_coeff(field: &'static str, value: f64) -> Result<(), OpticalCoefficientError> {
    if !value.is_finite() {
        return Err(OpticalCoefficientError::NonFinite { field, value });
    }
    if value < 0.0 {
        return Err(OpticalCoefficientError::Negative { field, value });
    }
    Ok(())
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
/// let medium = HomogeneousMedium::<Kilometer>::try_new(0.1, 0.2).unwrap();
/// let coeffs = medium.coefficients(
///     Position::<Center, Frame, Kilometer>::new(0.0, 0.0, 0.0),
///     Nanometers::new(550.0),
/// );
/// assert!((coeffs.sigma_t.value() - 0.3).abs() < 1e-12);
/// ```
pub trait Medium<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> {
    /// Returns optical coefficients at a position and wavelength.
    fn coefficients(
        &self,
        p: Position<C, F, U>,
        wavelength: qtty::length::Nanometers,
    ) -> OpticalCoefficients<U>;
}

/// Spatially uniform medium with wavelength-independent coefficients.
///
/// # Examples
///
/// ```rust
/// use optica::medium::HomogeneousMedium;
/// use qtty::unit::Kilometer;
///
/// let medium = HomogeneousMedium::<Kilometer>::try_new(0.1, 0.2).unwrap();
/// assert_eq!(medium, HomogeneousMedium::<Kilometer>::try_new(0.1, 0.2).unwrap());
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HomogeneousMedium<U: LengthUnit> {
    sigma_a: f64,
    sigma_s: f64,
    _phantom: PhantomData<U>,
}

impl<U: LengthUnit> HomogeneousMedium<U> {
    /// Constructs a homogeneous medium from constant validated coefficients.
    ///
    /// # Errors
    ///
    /// Returns [`OpticalCoefficientError`] when either input is non-finite or negative.
    pub fn try_new(sigma_a: f64, sigma_s: f64) -> Result<Self, OpticalCoefficientError> {
        check_coeff("sigma_a", sigma_a)?;
        check_coeff("sigma_s", sigma_s)?;
        Ok(Self {
            sigma_a,
            sigma_s,
            _phantom: PhantomData,
        })
    }

    /// Creates a transparent homogeneous medium.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::medium::HomogeneousMedium;
    /// use qtty::unit::Kilometer;
    ///
    /// let medium = HomogeneousMedium::<Kilometer>::transparent();
    /// assert_eq!(medium, HomogeneousMedium::<Kilometer>::try_new(0.0, 0.0).unwrap());
    /// ```
    #[must_use]
    pub fn transparent() -> Self {
        Self::try_new(0.0, 0.0).expect("zero coefficients are always valid")
    }
}

impl<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> Medium<C, F, U>
    for HomogeneousMedium<U>
{
    fn coefficients(
        &self,
        _p: Position<C, F, U>,
        _wavelength: qtty::length::Nanometers,
    ) -> OpticalCoefficients<U> {
        // Invariants enforced at construction time.
        OpticalCoefficients::try_new(self.sigma_a, self.sigma_s)
            .expect("HomogeneousMedium fields are validated at construction")
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

    #[test]
    fn try_new_rejects_negative_inputs() {
        let r = OpticalCoefficients::<qtty::unit::Kilometer>::try_new(-1.0, 0.0);
        assert!(matches!(r, Err(OpticalCoefficientError::Negative { .. })));
    }

    #[test]
    fn try_new_rejects_nan_inputs() {
        let r = OpticalCoefficients::<qtty::unit::Kilometer>::try_new(f64::NAN, 0.0);
        assert!(matches!(r, Err(OpticalCoefficientError::NonFinite { .. })));
    }

    #[test]
    fn homogeneous_medium_rejects_negative() {
        assert!(HomogeneousMedium::<qtty::unit::Kilometer>::try_new(0.0, -0.1).is_err());
    }
}
