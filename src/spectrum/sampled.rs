// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Allocation-free sampled spectrum container.
//!
//! ## Scientific scope
//!
//! A sampled spectrum is a discrete representation of a physical quantity
//! as a function of wavelength: measured passbands, solar reference spectra,
//! airglow continua, and ozone transmittance tables are all represented
//! as `SampledSpectrum`.
//!
//! ## Technical scope
//!
//! `SampledSpectrum<X, Y>` stores axis and value arrays as `Box<[f64]>`.
//! Interpolation and integration operate on raw slices without allocation;
//! the only heap allocation after construction is the optional cubic spline
//! coefficient buffer, precomputed once when `Interpolation::CubicSpline`
//! is chosen.

use core::marker::PhantomData;

use qtty::{DimMul, Dimension, Prod, Quantity, Unit};

use crate::data::Provenance;
use crate::grid::OutOfRange;
use crate::spectrum::algo::{self, CubicSplineCoeffs};
use crate::spectrum::{Interpolation, SpectrumError};

/// A sampled spectral function over axis unit `X` with value unit `Y`.
///
/// # Examples
///
/// ```rust
/// use optica::grid::OutOfRange;
/// use optica::spectrum::{Interpolation, SampledSpectrum};
/// use qtty::unit::{Nanometer, Ratio};
///
/// let xs = vec![400.0_f64, 500.0, 600.0, 700.0];
/// let ys = vec![0.0_f64, 0.5, 1.0, 0.5];
/// let s = SampledSpectrum::<Nanometer, Ratio>::from_raw(
///     xs,
///     ys,
///     Interpolation::Linear,
///     OutOfRange::ClampToEndpoints,
///     None,
/// )
/// .unwrap();
/// assert_eq!(s.len(), 4);
/// ```
#[derive(Debug, Clone)]
pub struct SampledSpectrum<X: Unit, Y: Unit> {
    xs: Box<[f64]>,
    ys: Box<[f64]>,
    interp: Interpolation,
    oor: OutOfRange,
    spline: Option<CubicSplineCoeffs>,
    provenance: Option<Provenance>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<X: Unit, Y: Unit> SampledSpectrum<X, Y> {
    /// Construct from owned `Vec`s with full validation.
    ///
    /// Validates length match, at least 2 samples, and strict monotonic increase
    /// of the x-axis. For [`Interpolation::CubicSpline`], precomputes coefficients.
    pub fn from_raw(
        xs: Vec<f64>,
        ys: Vec<f64>,
        interp: Interpolation,
        oor: OutOfRange,
        provenance: Option<Provenance>,
    ) -> Result<Self, SpectrumError> {
        algo::validate(&xs, &ys)?;
        let spline = if matches!(interp, Interpolation::CubicSpline) {
            Some(CubicSplineCoeffs::natural(&xs, &ys)?)
        } else {
            None
        };
        Ok(Self {
            xs: xs.into_boxed_slice(),
            ys: ys.into_boxed_slice(),
            interp,
            oor,
            spline,
            provenance,
            _x: PhantomData,
            _y: PhantomData,
        })
    }

    /// Construct from borrowed slices (copies data into `Box<[f64]>`).
    pub fn from_sorted(
        xs: &[f64],
        ys: &[f64],
        interp: Interpolation,
        oor: OutOfRange,
    ) -> Result<Self, SpectrumError> {
        Self::from_raw(xs.to_vec(), ys.to_vec(), interp, oor, None)
    }

    /// Raw x-axis values (zero-copy).
    #[inline]
    pub fn xs_raw(&self) -> &[f64] {
        &self.xs
    }

    /// Raw y-axis values (zero-copy).
    #[inline]
    pub fn ys_raw(&self) -> &[f64] {
        &self.ys
    }

    /// Number of samples.
    #[inline]
    pub fn len(&self) -> usize {
        self.xs.len()
    }

    /// Whether the spectrum stores no samples.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.xs.is_empty()
    }

    /// Interpolation mode.
    #[inline]
    pub fn interpolation(&self) -> Interpolation {
        self.interp
    }

    /// Out-of-range policy.
    #[inline]
    pub fn out_of_range(&self) -> OutOfRange {
        self.oor
    }

    /// Provenance metadata, if any.
    pub fn provenance(&self) -> Option<&Provenance> {
        self.provenance.as_ref()
    }

    /// Replace the provenance record.
    pub fn set_provenance(&mut self, provenance: Option<Provenance>) {
        self.provenance = provenance;
    }

    /// Builder-style provenance attachment.
    #[must_use]
    pub fn with_provenance(mut self, provenance: Provenance) -> Self {
        self.provenance = Some(provenance);
        self
    }

    /// Infallible interpolation at `x`.
    ///
    /// Out-of-range queries are always clamped to the nearest endpoint,
    /// regardless of the stored [`OutOfRange`] policy. Use [`Self::try_interp_at`]
    /// to enforce the policy.
    #[inline]
    pub fn interp_at(&self, x: Quantity<X>) -> Quantity<Y> {
        let raw = self
            .eval(x.value(), OutOfRange::ClampToEndpoints)
            .expect("ClampToEndpoints never errors");
        Quantity::<Y>::new(raw)
    }

    /// Fallible interpolation at `x`, respecting the stored [`OutOfRange`] policy.
    #[inline]
    pub fn try_interp_at(&self, x: Quantity<X>) -> Result<Quantity<Y>, SpectrumError> {
        let raw = self.eval(x.value(), self.oor)?;
        Ok(Quantity::<Y>::new(raw))
    }

    /// Trapezoidal integral over the full sampled domain.
    pub fn integrate(&self) -> Quantity<Prod<Y, X>>
    where
        Y::Dim: DimMul<X::Dim>,
        <Y::Dim as DimMul<X::Dim>>::Output: Dimension,
        Prod<Y, X>: Unit,
    {
        Quantity::<Prod<Y, X>>::new(algo::trapz(&self.xs, &self.ys))
    }

    /// Trapezoidal integral restricted to `[lo, hi]`.
    pub fn integrate_range(&self, lo: Quantity<X>, hi: Quantity<X>) -> Quantity<Prod<Y, X>>
    where
        Y::Dim: DimMul<X::Dim>,
        <Y::Dim as DimMul<X::Dim>>::Output: Dimension,
        Prod<Y, X>: Unit,
    {
        Quantity::<Prod<Y, X>>::new(algo::trapz_range(
            &self.xs,
            &self.ys,
            lo.value(),
            hi.value(),
        ))
    }

    /// Trapezoidal integral of `self(x) · weight(x)` over `weight`'s grid.
    pub fn integrate_weighted<Yw: Unit>(
        &self,
        weight: &SampledSpectrum<X, Yw>,
    ) -> Quantity<Prod<Prod<Y, Yw>, X>>
    where
        Y::Dim: DimMul<Yw::Dim>,
        <Y::Dim as DimMul<Yw::Dim>>::Output: DimMul<X::Dim>,
        <<Y::Dim as DimMul<Yw::Dim>>::Output as DimMul<X::Dim>>::Output: Dimension,
        Prod<Y, Yw>: Unit,
        Prod<Prod<Y, Yw>, X>: Unit,
    {
        Quantity::<Prod<Prod<Y, Yw>, X>>::new(algo::trapz_weighted(
            &self.xs,
            &self.ys,
            weight.xs_raw(),
            weight.ys_raw(),
        ))
    }

    fn eval(&self, x: f64, oor: OutOfRange) -> Result<f64, SpectrumError> {
        match self.interp {
            Interpolation::CubicSpline => {
                let coeffs = self
                    .spline
                    .as_ref()
                    .expect("spline precomputed at construction");
                algo::interp_cubic_spline(&self.xs, &self.ys, coeffs, x, oor)
            }
            _ => algo::interp(&self.xs, &self.ys, x, self.interp, oor),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use qtty::unit::{Nanometer, Ratio};

    fn make_linear() -> SampledSpectrum<Nanometer, Ratio> {
        let xs: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 * x + 1.0).collect();
        SampledSpectrum::from_raw(
            xs,
            ys,
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
        .unwrap()
    }

    #[test]
    fn interp_at_midpoint() {
        let s = make_linear();
        let val: f64 = s.interp_at(Quantity::<Nanometer>::new(2.5)).value();
        assert_abs_diff_eq!(val, 6.0, epsilon = 1e-12);
    }

    #[test]
    fn interp_at_clamps_oob_when_oor_is_error() {
        let xs = vec![0.0_f64, 1.0, 2.0];
        let ys = vec![1.0_f64, 2.0, 3.0];
        let s = SampledSpectrum::<Nanometer, Ratio>::from_raw(
            xs,
            ys,
            Interpolation::Linear,
            OutOfRange::Error,
            None,
        )
        .unwrap();
        let val = s.interp_at(Quantity::<Nanometer>::new(-5.0)).value();
        assert_abs_diff_eq!(val, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn try_interp_at_errors_when_oor_is_error() {
        let xs = vec![0.0_f64, 1.0, 2.0];
        let ys = vec![1.0_f64, 2.0, 3.0];
        let s = SampledSpectrum::<Nanometer, Ratio>::from_raw(
            xs,
            ys,
            Interpolation::Linear,
            OutOfRange::Error,
            None,
        )
        .unwrap();
        assert!(s.try_interp_at(Quantity::<Nanometer>::new(-5.0)).is_err());
    }

    #[test]
    fn xs_raw_ys_raw_zero_copy() {
        let s = make_linear();
        assert_eq!(s.xs_raw().len(), 5);
        assert_eq!(s.ys_raw().len(), 5);
    }

    #[test]
    fn from_sorted_works() {
        let xs = [0.0_f64, 1.0, 2.0];
        let ys = [1.0_f64, 2.0, 3.0];
        let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
            &xs,
            &ys,
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();
        assert_eq!(s.xs_raw(), &[0.0, 1.0, 2.0]);
    }
}
