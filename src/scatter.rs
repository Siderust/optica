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

/// Rayleigh optical depth using the Bodhaine et al. (1999) approximation.
///
/// This helper is generic in the sense that all environmental inputs are caller
/// supplied. It does not bake in any planetary constants or observatory presets.
///
/// # Examples
///
/// ```rust
/// use optica::scatter::rayleigh_optical_depth_bodhaine99;
/// use qtty::length::{Kilometers, Nanometers};
/// use qtty::pressure::Hectopascals;
///
/// let tau = rayleigh_optical_depth_bodhaine99(
///     Nanometers::new(550.0),
///     Hectopascals::new(1013.25),
///     Kilometers::new(0.0),
///     Kilometers::new(8.0),
/// );
/// assert!(tau.value() > 0.0);
/// ```
#[must_use]
pub fn rayleigh_optical_depth_bodhaine99(
    wavelength: Nanometers,
    surface_pressure: Hectopascals,
    observer_altitude: Kilometers,
    scale_height: Kilometers,
) -> OpticalDepths {
    let lam_um = wavelength.to::<Micrometer>().value();
    let p_atm = surface_pressure.value() / 1013.25;
    let h_km = observer_altitude.value();
    let scale_km = scale_height.value();
    if !lam_um.is_finite()
        || lam_um <= 0.0
        || !p_atm.is_finite()
        || !h_km.is_finite()
        || !scale_km.is_finite()
        || scale_km <= 0.0
    {
        return OpticalDepths::new(f64::NAN);
    }

    let lam_sq = lam_um * lam_um;
    let tau_sea = 0.002_152_0 * (1.045_599_6 - 341.290_61 / lam_sq - 0.902_308_50 * lam_sq)
        / (1.0 + 0.002_705_988_9 / lam_sq - 85.968_563 * lam_sq);
    OpticalDepths::new(p_atm * tau_sea * (-(h_km / scale_km)).exp())
}

/// Parameters for the Ångström power-law aerosol optical depth.
///
/// # Examples
///
/// ```rust
/// use optica::scatter::MieParams;
///
/// let params = MieParams::new(0.2, -1.3, 550.0);
/// assert_eq!(params.lambda_ref.value(), 550.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MieParams {
    /// Optical depth at the reference wavelength.
    pub tau0: OpticalDepths,
    /// Ångström exponent.
    pub alpha: f64,
    /// Reference wavelength.
    pub lambda_ref: Nanometers,
}

impl MieParams {
    /// Constructs parameters from explicit scalar values.
    #[must_use]
    pub fn new(tau0: f64, alpha: f64, lambda_ref_nm: f64) -> Self {
        Self {
            tau0: OpticalDepths::new(tau0),
            alpha,
            lambda_ref: Nanometers::new(lambda_ref_nm),
        }
    }
}

/// Aerosol optical depth: `τ(λ) = τ₀ (λ / λ_ref)^α`.
///
/// # Examples
///
/// ```rust
/// use optica::scatter::{mie_optical_depth, MieParams};
/// use qtty::length::Nanometers;
///
/// let params = MieParams::new(0.2, -1.0, 500.0);
/// let tau = mie_optical_depth(&params, Nanometers::new(1000.0));
/// assert!((tau.value() - 0.1).abs() < 1e-12);
/// ```
#[must_use]
pub fn mie_optical_depth(params: &MieParams, wavelength: Nanometers) -> OpticalDepths {
    let lambda = wavelength.value();
    let lambda_ref = params.lambda_ref.value();
    if !lambda.is_finite() || !lambda_ref.is_finite() || lambda <= 0.0 || lambda_ref <= 0.0 {
        return OpticalDepths::new(f64::NAN);
    }
    OpticalDepths::new(params.tau0.value() * (lambda / lambda_ref).powf(params.alpha))
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn mie_power_law_matches_reference_ratio() {
        let params = MieParams::new(0.2, -1.0, 500.0);
        let tau = mie_optical_depth(&params, Nanometers::new(1000.0));
        assert_relative_eq!(tau.value(), 0.1, epsilon = 1.0e-12);
    }
}
