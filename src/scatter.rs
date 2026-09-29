// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic single-scattering depth parameterizations.
//!
//! This module contains atmosphere-agnostic formulas. Callers must supply all
//! physical parameters explicitly; there are no Earth defaults or site presets.
//!
//! ## References
//!
//! - Bodhaine et al. (1999) for the Rayleigh approximation implemented here.
//! - Ångström aerosol optical-depth power law.

use qtty::dimensionless::OpticalDepths;
use qtty::length::{Kilometers, Nanometers};
use qtty::pressure::Hectopascals;
use qtty::unit::Micrometer;
#[cfg(not(feature = "std"))]
use qtty::Real;

/// Errors produced when validating scattering inputs.
///
/// # Examples
///
/// ```rust
/// use optica::scatter::ScatterError;
///
/// let err = ScatterError::NonPositive { field: "wavelength", value: 0.0 };
/// assert!(err.to_string().contains("must be strictly positive"));
/// ```
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ScatterError {
    /// Input is non-finite.
    #[error("{field} must be finite (got {value})")]
    NonFinite {
        /// Offending field.
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
    /// Input is not strictly positive.
    #[error("{field} must be strictly positive (got {value})")]
    NonPositive {
        /// Offending field.
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
    /// Input is negative.
    #[error("{field} must be non-negative (got {value})")]
    Negative {
        /// Offending field.
        field: &'static str,
        /// Supplied value.
        value: f64,
    },
}

/// Rayleigh optical depth using the Bodhaine et al. (1999) approximation.
///
/// Performs input validation and returns a fallible result. All environmental
/// inputs are caller supplied; the helper bakes in no planetary constants.
///
/// # Errors
///
/// Returns [`ScatterError`] when any input is non-finite, the wavelength or
/// scale height is non-positive, or the surface pressure is negative.
///
/// # Examples
///
/// ```rust
/// use optica::scatter::try_rayleigh_optical_depth_bodhaine99;
/// use qtty::length::{Kilometers, Nanometers};
/// use qtty::pressure::Hectopascals;
///
/// let tau = try_rayleigh_optical_depth_bodhaine99(
///     Nanometers::new(550.0),
///     Hectopascals::new(1013.25),
///     Kilometers::new(0.0),
///     Kilometers::new(8.0),
/// )
/// .unwrap();
/// assert!(tau.value() > 0.0);
/// ```
pub fn try_rayleigh_optical_depth_bodhaine99(
    wavelength: Nanometers,
    surface_pressure: Hectopascals,
    observer_altitude: Kilometers,
    scale_height: Kilometers,
) -> Result<OpticalDepths, ScatterError> {
    let lam_um = wavelength.to::<Micrometer>().value();
    let p_atm = surface_pressure.value() / 1013.25;
    let h_km = observer_altitude.value();
    let scale_km = scale_height.value();

    check_finite("wavelength", lam_um)?;
    check_positive("wavelength", lam_um)?;
    check_finite("surface_pressure", surface_pressure.value())?;
    check_non_negative("surface_pressure", surface_pressure.value())?;
    check_finite("observer_altitude", h_km)?;
    check_finite("scale_height", scale_km)?;
    check_positive("scale_height", scale_km)?;

    let lam_sq = lam_um * lam_um;
    let tau_sea = 0.002_152_0 * (1.045_599_6 - 341.290_61 / lam_sq - 0.902_308_50 * lam_sq)
        / (1.0 + 0.002_705_988_9 / lam_sq - 85.968_563 * lam_sq);
    Ok(OpticalDepths::new(
        p_atm * tau_sea * (-(h_km / scale_km)).exp(),
    ))
}

/// Parameters for the Ångström power-law aerosol optical depth.
///
/// Construct via [`MieParams::try_new`] to enforce non-negative `tau0`, a
/// finite Ångström exponent, and a strictly positive reference wavelength.
///
/// # Examples
///
/// ```rust
/// use optica::scatter::MieParams;
///
/// let params = MieParams::try_new(0.2, -1.3, 550.0).unwrap();
/// assert_eq!(params.lambda_ref.value(), 550.0);
/// assert!(MieParams::try_new(-0.1, -1.0, 500.0).is_err());
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MieParams {
    /// Optical depth at the reference wavelength.
    pub tau0: OpticalDepths,
    /// Ångström exponent.
    pub alpha: f64,
    /// Reference wavelength.
    pub lambda_ref: Nanometers,
}

impl MieParams {
    /// Constructs validated parameters from scalar values.
    ///
    /// # Errors
    ///
    /// Returns [`ScatterError`] when `tau0` is non-finite or negative, when
    /// `alpha` is non-finite, or when `lambda_ref_nm` is non-finite or
    /// non-positive.
    pub fn try_new(tau0: f64, alpha: f64, lambda_ref_nm: f64) -> Result<Self, ScatterError> {
        check_finite("tau0", tau0)?;
        check_non_negative("tau0", tau0)?;
        check_finite("alpha", alpha)?;
        check_finite("lambda_ref", lambda_ref_nm)?;
        check_positive("lambda_ref", lambda_ref_nm)?;
        Ok(Self {
            tau0: OpticalDepths::new(tau0),
            alpha,
            lambda_ref: Nanometers::new(lambda_ref_nm),
        })
    }
}

/// Aerosol optical depth: `τ(λ) = τ₀ (λ / λ_ref)^α`.
///
/// # Errors
///
/// Returns [`ScatterError`] when `wavelength` is non-finite or non-positive.
///
/// # Examples
///
/// ```rust
/// use optica::scatter::{try_mie_optical_depth, MieParams};
/// use qtty::length::Nanometers;
///
/// let params = MieParams::try_new(0.2, -1.0, 500.0).unwrap();
/// let tau = try_mie_optical_depth(&params, Nanometers::new(1000.0)).unwrap();
/// assert!((tau.value() - 0.1).abs() < 1e-12);
/// ```
pub fn try_mie_optical_depth(
    params: &MieParams,
    wavelength: Nanometers,
) -> Result<OpticalDepths, ScatterError> {
    let lambda = wavelength.value();
    check_finite("wavelength", lambda)?;
    check_positive("wavelength", lambda)?;
    let lambda_ref = params.lambda_ref.value();
    Ok(OpticalDepths::new(
        params.tau0.value() * (lambda / lambda_ref).powf(params.alpha),
    ))
}

fn check_finite(field: &'static str, value: f64) -> Result<(), ScatterError> {
    if !value.is_finite() {
        return Err(ScatterError::NonFinite { field, value });
    }
    Ok(())
}

fn check_positive(field: &'static str, value: f64) -> Result<(), ScatterError> {
    if value <= 0.0 {
        return Err(ScatterError::NonPositive { field, value });
    }
    Ok(())
}

fn check_non_negative(field: &'static str, value: f64) -> Result<(), ScatterError> {
    if value < 0.0 {
        return Err(ScatterError::Negative { field, value });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn mie_power_law_matches_reference_ratio() {
        let params = MieParams::try_new(0.2, -1.0, 500.0).unwrap();
        let tau = try_mie_optical_depth(&params, Nanometers::new(1000.0)).unwrap();
        assert_relative_eq!(tau.value(), 0.1, epsilon = 1.0e-12);
    }

    #[test]
    fn mie_rejects_negative_tau0() {
        assert!(matches!(
            MieParams::try_new(-0.1, -1.0, 500.0),
            Err(ScatterError::Negative { .. })
        ));
    }

    #[test]
    fn rayleigh_validates_inputs() {
        let r = try_rayleigh_optical_depth_bodhaine99(
            Nanometers::new(550.0),
            Hectopascals::new(1013.25),
            Kilometers::new(0.0),
            Kilometers::new(0.0),
        );
        assert!(matches!(r, Err(ScatterError::NonPositive { .. })));
    }

    #[test]
    fn rayleigh_succeeds_for_reasonable_inputs() {
        let tau = try_rayleigh_optical_depth_bodhaine99(
            Nanometers::new(550.0),
            Hectopascals::new(1013.25),
            Kilometers::new(0.0),
            Kilometers::new(8.0),
        )
        .unwrap();
        assert!(tau.value() > 0.0 && tau.value() < 1.0);
    }
}
