// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Optical-depth integration and transport helpers.
//!
//! The numerical integration helper uses a midpoint rule over typed rays and
//! medium queries. This keeps public APIs fully typed while keeping the inner loop
//! on raw `f64` scalars.

use affn::{ReferenceCenter, ReferenceFrame};
use qtty::angular::Radians;
use qtty::dimensionless::{OpticalDepths, Transmittances};
use qtty::length::{Kilometers, LengthUnit, Nanometers};
use qtty::Quantity;

use crate::medium::Medium;
use crate::ray::{Ray, RaySegment};

/// Beer-Lambert transmittance: `T = exp(-τ)`.
///
/// # Examples
///
/// ```rust
/// use optica::transport::transmittance;
/// use qtty::dimensionless::OpticalDepths;
///
/// let transmission = transmittance(OpticalDepths::new(0.0));
/// assert_eq!(transmission.value(), 1.0);
/// ```
#[must_use]
pub fn transmittance(tau: OpticalDepths) -> Transmittances {
    Transmittances::new((-tau.value()).exp())
}

/// Options for numerical optical-depth integration.
///
/// # Examples
///
/// ```rust
/// use optica::transport::IntegrationOpts;
///
/// assert_eq!(IntegrationOpts::default().n_steps, 64);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IntegrationOpts {
    /// Number of midpoint samples used along the segment.
    pub n_steps: usize,
}

impl Default for IntegrationOpts {
    fn default() -> Self {
        Self { n_steps: 64 }
    }
}

/// Integrates optical depth along a ray segment with the midpoint rule.
///
/// `τ = ∫ σ_t(p(t)) dt` from `t_min` to `t_max`.
///
/// # Examples
///
/// ```rust
/// use affn::{CartesianDirection, Position, ReferenceCenter, ReferenceFrame};
/// use optica::medium::HomogeneousMedium;
/// use optica::ray::{Ray, RaySegment};
/// use optica::transport::{integrate_optical_depth, IntegrationOpts};
/// use qtty::length::{Kilometers, Nanometers};
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
/// let ray = Ray::new(
///     Position::<Center, Frame, Kilometer>::new(0.0, 0.0, 0.0),
///     CartesianDirection::<Frame>::new(0.0, 0.0, 1.0),
/// );
/// let tau = integrate_optical_depth(
///     &medium,
///     &ray,
///     RaySegment::new(Kilometers::new(0.0), Kilometers::new(10.0)),
///     Nanometers::new(550.0),
///     IntegrationOpts { n_steps: 16 },
/// );
/// assert!((tau.value() - 3.0).abs() < 1e-12);
/// ```
#[must_use]
pub fn integrate_optical_depth<M, C, F, U>(
    medium: &M,
    ray: &Ray<C, F, U>,
    segment: RaySegment<U>,
    wavelength: Nanometers,
    opts: IntegrationOpts,
) -> OpticalDepths
where
    M: Medium<C, F, U>,
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit + Copy,
{
    let n = opts.n_steps.max(1);
    let dt = (segment.t_max.value() - segment.t_min.value()) / n as f64;
    let mut tau = 0.0;
    for i in 0..n {
        let t = segment.t_min.value() + (i as f64 + 0.5) * dt;
        let disp = ray.direction * Quantity::<U>::new(t);
        let p = ray.origin.clone() + disp;
        let coeffs = medium.coefficients(p, wavelength);
        tau += coeffs.sigma_t.value() * dt;
    }
    OpticalDepths::new(tau)
}

/// Generic van Rhijn shell path-length factor for a spherical emitting layer.
///
/// `V(z, h) = 1 / sqrt(1 - (R / (R + h) * sin(z))^2)`
///
/// # Examples
///
/// ```rust
/// use optica::transport::van_rhijn_factor;
/// use qtty::angular::Radians;
/// use qtty::length::Kilometers;
///
/// let factor = van_rhijn_factor(Radians::new(0.0), Kilometers::new(90.0), Kilometers::new(6371.0));
/// assert_eq!(factor, 1.0);
/// ```
#[must_use]
pub fn van_rhijn_factor(zenith: Radians, emission_height: Kilometers, body_radius: Kilometers) -> f64 {
    let r = body_radius.value();
    let h = emission_height.value();
    if !zenith.value().is_finite() || !h.is_finite() || !r.is_finite() || h <= 0.0 || r <= 0.0 {
        return f64::NAN;
    }
    let ratio = r / (r + h);
    let sine = zenith.value().sin();
    let inner = 1.0 - (ratio * sine) * (ratio * sine);
    if inner <= 0.0 {
        f64::INFINITY
    } else {
        inner.sqrt().recip()
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;
    use affn::{CartesianDirection, Position};
    use crate::medium::HomogeneousMedium;

    #[derive(Debug, Copy, Clone)]
    struct Center;

    impl ReferenceCenter for Center {
        type Params = ();
        fn center_name() -> &'static str {
            "Center"
        }
    }

    #[derive(Debug, Copy, Clone)]
    struct Frame;

    impl ReferenceFrame for Frame {
        fn frame_name() -> &'static str {
            "Frame"
        }
    }

    #[test]
    fn integrates_homogeneous_medium_exactly() {
        let medium = HomogeneousMedium::<qtty::unit::Kilometer>::new(0.1, 0.2);
        let ray = Ray::new(
            Position::<Center, Frame, qtty::unit::Kilometer>::new(0.0, 0.0, 0.0),
            CartesianDirection::<Frame>::new(0.0, 0.0, 1.0),
        );
        let tau = integrate_optical_depth(
            &medium,
            &ray,
            RaySegment::new(Kilometers::new(0.0), Kilometers::new(10.0)),
            Nanometers::new(550.0),
            IntegrationOpts { n_steps: 8 },
        );
        assert_relative_eq!(tau.value(), 3.0, epsilon = 1.0e-12);
    }
}
