// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Optical-depth integration and transport helpers.
//!
//! The integrator dispatches on [`IntegrationMethod`] and validates inputs via
//! [`try_integrate_optical_depth`]. The inner loop runs on raw `f64` scalars
//! while all public surfaces remain typed over [`affn`] and [`qtty`].

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

/// Numerical integration kernels supported by [`try_integrate_optical_depth`].
///
/// Choose by the smoothness of `σ_t(t)` along the segment:
///
/// - [`Midpoint`](Self::Midpoint): first-order, robust for piecewise data.
/// - [`Trapezoidal`](Self::Trapezoidal): second-order; same node count as midpoint.
/// - [`Simpson`](Self::Simpson): fourth-order on smooth integrands. Requires
///   `n_steps` to be even; odd values are rounded up.
/// - [`GaussLegendre2`](Self::GaussLegendre2): 2-point Gauss-Legendre quadrature
///   applied on each sub-interval; exact for polynomials of degree ≤ 3 per cell.
/// - [`GaussLegendre4`](Self::GaussLegendre4): 4-point Gauss-Legendre quadrature
///   per sub-interval; exact for polynomials of degree ≤ 7 per cell.
///
/// # Examples
///
/// ```rust
/// use optica::transport::IntegrationMethod;
///
/// assert_eq!(IntegrationMethod::default(), IntegrationMethod::Midpoint);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum IntegrationMethod {
    /// Composite midpoint rule (default).
    #[default]
    Midpoint,
    /// Composite trapezoidal rule.
    Trapezoidal,
    /// Composite Simpson's rule (requires even `n_steps`).
    Simpson,
    /// 2-point Gauss-Legendre quadrature per sub-interval.
    GaussLegendre2,
    /// 4-point Gauss-Legendre quadrature per sub-interval.
    GaussLegendre4,
}

/// Options for numerical optical-depth integration.
///
/// # Examples
///
/// ```rust
/// use optica::transport::{IntegrationMethod, IntegrationOpts};
///
/// let opts = IntegrationOpts::default();
/// assert_eq!(opts.n_steps, 64);
/// assert_eq!(opts.method, IntegrationMethod::Midpoint);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct IntegrationOpts {
    /// Number of sub-intervals along the segment.
    pub n_steps: usize,
    /// Quadrature kernel.
    pub method: IntegrationMethod,
}

impl Default for IntegrationOpts {
    fn default() -> Self {
        Self {
            n_steps: 64,
            method: IntegrationMethod::Midpoint,
        }
    }
}

impl IntegrationOpts {
    /// Constructs options with an explicit method.
    #[must_use]
    pub const fn new(n_steps: usize, method: IntegrationMethod) -> Self {
        Self { n_steps, method }
    }
}

/// Errors produced by [`try_integrate_optical_depth`].
///
/// # Examples
///
/// ```rust
/// use optica::transport::TransportError;
///
/// let err = TransportError::NonPositiveSteps;
/// assert!(err.to_string().contains("n_steps"));
/// ```
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TransportError {
    /// `n_steps` was zero.
    #[error("n_steps must be ≥ 1")]
    NonPositiveSteps,
    /// The integration segment was empty or descending (`t_max < t_min`).
    #[error("segment must satisfy t_max ≥ t_min (got [{lo}, {hi}])")]
    InvalidSegment {
        /// Lower bound supplied.
        lo: f64,
        /// Upper bound supplied.
        hi: f64,
    },
    /// One of the segment endpoints was non-finite.
    #[error("segment endpoints must be finite (got [{lo}, {hi}])")]
    NonFiniteSegment {
        /// Lower bound supplied.
        lo: f64,
        /// Upper bound supplied.
        hi: f64,
    },
}

/// Integrates optical depth along a ray segment using [`IntegrationOpts::default`].
///
/// Convenience wrapper around [`try_integrate_optical_depth`] that panics on
/// invalid input; intended for hot paths where `n_steps ≥ 1` and the segment
/// is finite by construction. For caller-supplied numerical configuration prefer
/// the fallible variant.
///
/// # Panics
///
/// Panics when validation by [`try_integrate_optical_depth`] fails.
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
/// let medium = HomogeneousMedium::<Kilometer>::try_new(0.1, 0.2).unwrap();
/// let ray = Ray::new(
///     Position::<Center, Frame, Kilometer>::new(0.0, 0.0, 0.0),
///     CartesianDirection::<Frame>::new(0.0, 0.0, 1.0),
/// );
/// let tau = integrate_optical_depth(
///     &medium,
///     &ray,
///     RaySegment::new(Kilometers::new(0.0), Kilometers::new(10.0)),
///     Nanometers::new(550.0),
///     IntegrationOpts { n_steps: 16, ..Default::default() },
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
    try_integrate_optical_depth(medium, ray, segment, wavelength, opts)
        .expect("validated inputs satisfy try_integrate_optical_depth contract")
}

/// Fallible variant of [`integrate_optical_depth`] that validates inputs.
///
/// # Errors
///
/// Returns [`TransportError`] when `opts.n_steps == 0`, the segment endpoints
/// are non-finite, or `t_max < t_min`.
pub fn try_integrate_optical_depth<M, C, F, U>(
    medium: &M,
    ray: &Ray<C, F, U>,
    segment: RaySegment<U>,
    wavelength: Nanometers,
    opts: IntegrationOpts,
) -> Result<OpticalDepths, TransportError>
where
    M: Medium<C, F, U>,
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit + Copy,
{
    let lo = segment.t_min.value();
    let hi = segment.t_max.value();
    if !lo.is_finite() || !hi.is_finite() {
        return Err(TransportError::NonFiniteSegment { lo, hi });
    }
    if hi < lo {
        return Err(TransportError::InvalidSegment { lo, hi });
    }
    if opts.n_steps == 0 {
        return Err(TransportError::NonPositiveSteps);
    }
    let sample = |t: f64| -> f64 {
        let disp = ray.direction * Quantity::<U>::new(t);
        let p = ray.origin.clone() + disp;
        medium.coefficients(p, wavelength).sigma_t.value()
    };

    let tau = match opts.method {
        IntegrationMethod::Midpoint => composite_midpoint(lo, hi, opts.n_steps, sample),
        IntegrationMethod::Trapezoidal => composite_trapezoidal(lo, hi, opts.n_steps, sample),
        IntegrationMethod::Simpson => {
            let n = if opts.n_steps % 2 == 0 {
                opts.n_steps
            } else {
                opts.n_steps + 1
            };
            composite_simpson(lo, hi, n, sample)
        }
        IntegrationMethod::GaussLegendre2 => {
            composite_gauss(lo, hi, opts.n_steps, &GAUSS2_NODES, &GAUSS2_WEIGHTS, sample)
        }
        IntegrationMethod::GaussLegendre4 => {
            composite_gauss(lo, hi, opts.n_steps, &GAUSS4_NODES, &GAUSS4_WEIGHTS, sample)
        }
    };
    Ok(OpticalDepths::new(tau))
}

fn composite_midpoint(lo: f64, hi: f64, n: usize, f: impl Fn(f64) -> f64) -> f64 {
    let dt = (hi - lo) / n as f64;
    let mut acc = 0.0;
    for i in 0..n {
        let t = lo + (i as f64 + 0.5) * dt;
        acc += f(t) * dt;
    }
    acc
}

fn composite_trapezoidal(lo: f64, hi: f64, n: usize, f: impl Fn(f64) -> f64) -> f64 {
    let dt = (hi - lo) / n as f64;
    let mut acc = 0.5 * (f(lo) + f(hi)) * dt;
    for i in 1..n {
        let t = lo + i as f64 * dt;
        acc += f(t) * dt;
    }
    acc
}

fn composite_simpson(lo: f64, hi: f64, n: usize, f: impl Fn(f64) -> f64) -> f64 {
    debug_assert!(n % 2 == 0 && n >= 2);
    let dt = (hi - lo) / n as f64;
    let mut acc = f(lo) + f(hi);
    for i in 1..n {
        let t = lo + i as f64 * dt;
        let w = if i % 2 == 0 { 2.0 } else { 4.0 };
        acc += w * f(t);
    }
    acc * dt / 3.0
}

fn composite_gauss(
    lo: f64,
    hi: f64,
    n: usize,
    nodes: &[f64],
    weights: &[f64],
    f: impl Fn(f64) -> f64,
) -> f64 {
    let dt = (hi - lo) / n as f64;
    let half = 0.5 * dt;
    let mut acc = 0.0;
    for i in 0..n {
        let center = lo + (i as f64 + 0.5) * dt;
        for (&x, &w) in nodes.iter().zip(weights.iter()) {
            acc += w * f(center + half * x);
        }
    }
    acc * half
}

// 2-point Gauss-Legendre on [-1, 1].
const GAUSS2_NODES: [f64; 2] = [-0.577_350_269_189_625_8, 0.577_350_269_189_625_8];
const GAUSS2_WEIGHTS: [f64; 2] = [1.0, 1.0];

// 4-point Gauss-Legendre on [-1, 1].
const GAUSS4_NODES: [f64; 4] = [
    -0.861_136_311_594_052_6,
    -0.339_981_043_584_856_3,
    0.339_981_043_584_856_3,
    0.861_136_311_594_052_6,
];
const GAUSS4_WEIGHTS: [f64; 4] = [
    0.347_854_845_137_453_8,
    0.652_145_154_862_546_3,
    0.652_145_154_862_546_3,
    0.347_854_845_137_453_8,
];

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
pub fn van_rhijn_factor(
    zenith: Radians,
    emission_height: Kilometers,
    body_radius: Kilometers,
) -> f64 {
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
    use crate::medium::HomogeneousMedium;
    use affn::{CartesianDirection, Position};

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

    fn setup() -> (
        HomogeneousMedium<qtty::unit::Kilometer>,
        Ray<Center, Frame, qtty::unit::Kilometer>,
    ) {
        let medium = HomogeneousMedium::<qtty::unit::Kilometer>::try_new(0.1, 0.2).unwrap();
        let ray = Ray::new(
            Position::<Center, Frame, qtty::unit::Kilometer>::new(0.0, 0.0, 0.0),
            CartesianDirection::<Frame>::new(0.0, 0.0, 1.0),
        );
        (medium, ray)
    }

    #[test]
    fn integrates_homogeneous_medium_exactly_midpoint() {
        let (medium, ray) = setup();
        let tau = integrate_optical_depth(
            &medium,
            &ray,
            RaySegment::new(Kilometers::new(0.0), Kilometers::new(10.0)),
            Nanometers::new(550.0),
            IntegrationOpts::new(8, IntegrationMethod::Midpoint),
        );
        assert_relative_eq!(tau.value(), 3.0, epsilon = 1.0e-12);
    }

    #[test]
    fn integrates_with_trapezoidal_and_simpson_and_gauss() {
        let (medium, ray) = setup();
        for method in [
            IntegrationMethod::Trapezoidal,
            IntegrationMethod::Simpson,
            IntegrationMethod::GaussLegendre2,
            IntegrationMethod::GaussLegendre4,
        ] {
            let tau = integrate_optical_depth(
                &medium,
                &ray,
                RaySegment::new(Kilometers::new(0.0), Kilometers::new(10.0)),
                Nanometers::new(550.0),
                IntegrationOpts::new(8, method),
            );
            assert_relative_eq!(tau.value(), 3.0, epsilon = 1.0e-10);
        }
    }

    #[test]
    fn rejects_zero_steps() {
        let (medium, ray) = setup();
        let r = try_integrate_optical_depth(
            &medium,
            &ray,
            RaySegment::new(Kilometers::new(0.0), Kilometers::new(1.0)),
            Nanometers::new(550.0),
            IntegrationOpts::new(0, IntegrationMethod::Midpoint),
        );
        assert!(matches!(r, Err(TransportError::NonPositiveSteps)));
    }

    #[test]
    fn rejects_descending_segment_when_construct_manually() {
        // RaySegment::new auto-orders, so we construct directly to bypass.
        let (medium, ray) = setup();
        let bad = RaySegment {
            t_min: Kilometers::new(1.0),
            t_max: Kilometers::new(0.0),
        };
        let r = try_integrate_optical_depth(
            &medium,
            &ray,
            bad,
            Nanometers::new(550.0),
            IntegrationOpts::default(),
        );
        assert!(matches!(r, Err(TransportError::InvalidSegment { .. })));
    }
}
