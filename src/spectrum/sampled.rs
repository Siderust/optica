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

use alloc::{boxed::Box, vec::Vec};

use qtty::unit::Ratio;
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

    /// Inclusive sampled domain as `(x_min, x_max)`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[400.0, 500.0, 600.0],
    ///     &[0.1, 0.2, 0.3],
    ///     Interpolation::Linear,
    ///     OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let (lo, hi) = s.domain();
    /// assert_eq!(lo.value(), 400.0);
    /// assert_eq!(hi.value(), 600.0);
    /// ```
    #[inline]
    pub fn domain(&self) -> (Quantity<X>, Quantity<X>) {
        (
            Quantity::<X>::new(self.xs[0]),
            Quantity::<X>::new(self.xs[self.xs.len() - 1]),
        )
    }

    /// Returns whether `x` is inside the inclusive sampled domain.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::Quantity;
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[400.0, 500.0],
    ///     &[0.1, 0.2],
    ///     Interpolation::Linear,
    ///     OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// assert!(s.contains(Quantity::<Nanometer>::new(450.0)));
    /// assert!(!s.contains(Quantity::<Nanometer>::new(700.0)));
    /// ```
    #[inline]
    pub fn contains(&self, x: Quantity<X>) -> bool {
        let v = x.value();
        v >= self.xs[0] && v <= self.xs[self.xs.len() - 1]
    }

    /// Returns the intersection of two sampled domains as `(lo, hi)`, or `None`
    /// when they do not overlap.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let a = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[400.0, 600.0], &[0.1, 0.2],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let b = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[500.0, 700.0], &[0.3, 0.4],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let (lo, hi) = a.overlap_domain(&b).unwrap();
    /// assert_eq!(lo.value(), 500.0);
    /// assert_eq!(hi.value(), 600.0);
    /// ```
    #[must_use]
    pub fn overlap_domain<Yo: Unit>(
        &self,
        other: &SampledSpectrum<X, Yo>,
    ) -> Option<(Quantity<X>, Quantity<X>)> {
        let (alo, ahi) = self.domain();
        let (blo, bhi) = other.domain();
        let lo = if alo.value() >= blo.value() { alo } else { blo };
        let hi = if ahi.value() <= bhi.value() { ahi } else { bhi };
        if lo.value() <= hi.value() {
            Some((lo, hi))
        } else {
            None
        }
    }

    /// Resamples `self` onto an explicit set of typed x-locations.
    ///
    /// Each input point is interpolated using the current [`Interpolation`]
    /// kernel; the produced spectrum uses [`Interpolation::Linear`] with
    /// [`OutOfRange::ClampToEndpoints`] for downstream querying. Provenance is
    /// not propagated.
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError`] when the requested grid violates monotonicity
    /// or has fewer than two points, or when the source out-of-range policy
    /// rejects an evaluation.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::{Quantity, unit::{Nanometer, Ratio}};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[400.0, 600.0], &[0.0, 1.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let grid = [Quantity::<Nanometer>::new(450.0), Quantity::<Nanometer>::new(550.0)];
    /// let r = s.resample_onto(&grid).unwrap();
    /// assert!((r.ys_raw()[0] - 0.25).abs() < 1e-12);
    /// assert!((r.ys_raw()[1] - 0.75).abs() < 1e-12);
    /// ```
    pub fn resample_onto(&self, xs: &[Quantity<X>]) -> Result<Self, SpectrumError> {
        let raw: Vec<f64> = xs.iter().map(|q| q.value()).collect();
        self.resample_onto_raw(&raw)
    }

    /// Resamples `self` onto an explicit set of raw x-values (in axis units).
    ///
    /// This is an escape hatch for callers that already hold raw `f64` grids.
    /// Prefer [`resample_onto`](Self::resample_onto) for unit-safe code.
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError`] when the requested grid violates monotonicity
    /// or has fewer than two points, or when the source out-of-range policy
    /// rejects an evaluation.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[400.0, 600.0], &[0.0, 1.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let r = s.resample_onto_raw(&[450.0, 550.0]).unwrap();
    /// assert!((r.ys_raw()[0] - 0.25).abs() < 1e-12);
    /// assert!((r.ys_raw()[1] - 0.75).abs() < 1e-12);
    /// ```
    pub fn resample_onto_raw(&self, xs: &[f64]) -> Result<Self, SpectrumError> {
        let mut ys = Vec::with_capacity(xs.len());
        for &x in xs {
            ys.push(self.eval(x, self.oor)?);
        }
        Self::from_raw(
            xs.to_vec(),
            ys,
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
            None,
        )
    }

    /// Returns a new spectrum scaled so that its trapezoidal integral over the
    /// sampled domain equals 1 (in raw axis units).
    ///
    /// The returned spectrum has value unit [`Ratio`](qtty::unit::Ratio)
    /// regardless of `Y`, because dividing every sample by the area strips the
    /// dimensional meaning of the original values.  To preserve `Y` (at the
    /// cost of dimensional unsoundness), see [`normalize_area_raw`](Self::normalize_area_raw).
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError::InvalidValue`] when the integral is zero or
    /// non-finite.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[0.0, 1.0, 2.0], &[1.0, 1.0, 1.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let n: SampledSpectrum<Nanometer, Ratio> = s.normalize_area().unwrap();
    /// // ∫ 1 dx over [0, 2] = 2, so each sample becomes 0.5.
    /// assert!((n.ys_raw()[0] - 0.5).abs() < 1e-12);
    /// ```
    pub fn normalize_area(&self) -> Result<SampledSpectrum<X, Ratio>, SpectrumError> {
        let area = algo::trapz(&self.xs, &self.ys);
        if !area.is_finite() || area == 0.0 {
            return Err(SpectrumError::InvalidValue {
                what: alloc::string::String::from("integrated area must be finite and non-zero"),
            });
        }
        let inv = area.recip();
        let ys: Vec<f64> = self.ys.iter().map(|&y| y * inv).collect();
        let xs: Vec<f64> = self.xs.to_vec();
        SampledSpectrum::<X, Ratio>::from_raw(
            xs,
            ys,
            self.interp,
            self.oor,
            self.provenance.clone(),
        )
    }

    /// Returns a new spectrum scaled so that its trapezoidal integral equals 1,
    /// preserving the original value unit `Y`.
    ///
    /// This is an escape hatch for callers that need to keep the nominal `Y`
    /// type despite the unit-semantic change that normalisation produces.
    /// Prefer [`normalize_area`](Self::normalize_area) for dimensionally sound
    /// code.
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError::InvalidValue`] when the integral is zero or
    /// non-finite.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[0.0, 1.0, 2.0], &[1.0, 1.0, 1.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let n = s.normalize_area_raw().unwrap();
    /// assert!((n.ys_raw()[0] - 0.5).abs() < 1e-12);
    /// ```
    pub fn normalize_area_raw(&self) -> Result<Self, SpectrumError> {
        let area = algo::trapz(&self.xs, &self.ys);
        if !area.is_finite() || area == 0.0 {
            return Err(SpectrumError::InvalidValue {
                what: alloc::string::String::from("integrated area must be finite and non-zero"),
            });
        }
        let inv = area.recip();
        let ys: Vec<f64> = self.ys.iter().map(|&y| y * inv).collect();
        let xs: Vec<f64> = self.xs.to_vec();
        Self::from_raw(xs, ys, self.interp, self.oor, self.provenance.clone())
    }

    /// Returns a new spectrum scaled so that its maximum sample equals 1.
    ///
    /// The returned spectrum has value unit [`Ratio`](qtty::unit::Ratio)
    /// because peak normalisation divides every sample by the peak value,
    /// which carries the same unit `Y`, yielding a dimensionless relative shape.
    /// To preserve `Y` (at the cost of dimensional unsoundness), see
    /// [`normalize_peak_raw`](Self::normalize_peak_raw).
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError::InvalidValue`] when the maximum sample is zero
    /// or non-finite.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[0.0, 1.0, 2.0], &[1.0, 4.0, 2.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let n: SampledSpectrum<Nanometer, Ratio> = s.normalize_peak().unwrap();
    /// assert!((n.ys_raw()[1] - 1.0).abs() < 1e-12);
    /// ```
    pub fn normalize_peak(&self) -> Result<SampledSpectrum<X, Ratio>, SpectrumError> {
        let peak = self.ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !peak.is_finite() || peak == 0.0 {
            return Err(SpectrumError::InvalidValue {
                what: alloc::string::String::from("peak value must be finite and non-zero"),
            });
        }
        let inv = peak.recip();
        let ys: Vec<f64> = self.ys.iter().map(|&y| y * inv).collect();
        let xs: Vec<f64> = self.xs.to_vec();
        SampledSpectrum::<X, Ratio>::from_raw(
            xs,
            ys,
            self.interp,
            self.oor,
            self.provenance.clone(),
        )
    }

    /// Returns a new spectrum scaled so that its maximum sample equals 1,
    /// preserving the original value unit `Y`.
    ///
    /// This is an escape hatch for callers that need to keep the nominal `Y`
    /// type despite the unit-semantic change. Prefer
    /// [`normalize_peak`](Self::normalize_peak) for dimensionally sound code.
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError::InvalidValue`] when the maximum sample is zero
    /// or non-finite.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[0.0, 1.0, 2.0], &[1.0, 4.0, 2.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let n = s.normalize_peak_raw().unwrap();
    /// assert!((n.ys_raw()[1] - 1.0).abs() < 1e-12);
    /// ```
    pub fn normalize_peak_raw(&self) -> Result<Self, SpectrumError> {
        let peak = self.ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !peak.is_finite() || peak == 0.0 {
            return Err(SpectrumError::InvalidValue {
                what: alloc::string::String::from("peak value must be finite and non-zero"),
            });
        }
        let inv = peak.recip();
        let ys: Vec<f64> = self.ys.iter().map(|&y| y * inv).collect();
        let xs: Vec<f64> = self.xs.to_vec();
        Self::from_raw(xs, ys, self.interp, self.oor, self.provenance.clone())
    }

    /// Maps every value through a typed closure, producing a spectrum with a
    /// new value unit `Y2`.
    ///
    /// This is the dimensionally safe variant: the closure receives a
    /// [`Quantity<Y>`] and returns a [`Quantity<Y2>`], so the compiler enforces
    /// that the output unit is intentional.
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError`] when validation of the transformed values fails.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::{Quantity, unit::{Nanometer, Ratio}};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[400.0, 500.0], &[0.1, 0.2],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let s2: SampledSpectrum<Nanometer, Ratio> =
    ///     s.map_values_to(|y: Quantity<Ratio>| Quantity::<Ratio>::new(y.value() * 10.0))
    ///      .unwrap();
    /// assert_eq!(s2.ys_raw(), &[1.0, 2.0]);
    /// ```
    pub fn map_values_to<Y2: Unit, F>(&self, f: F) -> Result<SampledSpectrum<X, Y2>, SpectrumError>
    where
        F: Fn(Quantity<Y>) -> Quantity<Y2>,
    {
        let ys: Vec<f64> = self
            .ys
            .iter()
            .map(|&y| f(Quantity::<Y>::new(y)).value())
            .collect();
        let xs: Vec<f64> = self.xs.to_vec();
        SampledSpectrum::<X, Y2>::from_raw(xs, ys, self.interp, self.oor, self.provenance.clone())
    }

    /// Maps every value through a raw `f64` closure, preserving x-axis and metadata.
    ///
    /// This is an escape hatch for callers that need low-level control or that
    /// are already working in raw (unit-normalised) values.  Prefer
    /// [`map_values_to`](Self::map_values_to) for dimensionally safe code.
    ///
    /// `f` must be deterministic and finite-preserving; non-finite outputs will
    /// be rejected by validation.
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError`] when validation of the transformed values fails.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let s = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[400.0, 500.0], &[0.1, 0.2],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let s2 = s.map_values_raw(|y| y * 10.0).unwrap();
    /// assert_eq!(s2.ys_raw(), &[1.0, 2.0]);
    /// ```
    pub fn map_values_raw<F>(&self, f: F) -> Result<Self, SpectrumError>
    where
        F: Fn(f64) -> f64,
    {
        let ys: Vec<f64> = self.ys.iter().map(|&y| f(y)).collect();
        let xs: Vec<f64> = self.xs.to_vec();
        Self::from_raw(xs, ys, self.interp, self.oor, self.provenance.clone())
    }

    /// Combines two spectra sample-by-sample using a typed closure, producing a
    /// spectrum with a new value unit `Yout`.
    ///
    /// `other` is resampled onto `self`'s x-axis before the closure is applied.
    /// The output uses `self`'s interpolation and out-of-range policy.
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError`] when `other` cannot be evaluated at a point in
    /// `self`'s grid under its own out-of-range policy.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::{Quantity, unit::{Nanometer, Ratio}};
    ///
    /// let a = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[0.0, 1.0, 2.0], &[1.0, 2.0, 3.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let b = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[0.0, 2.0], &[10.0, 20.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let c: SampledSpectrum<Nanometer, Ratio> = a.zip_with_to(&b, |ya, yb| {
    ///     Quantity::<Ratio>::new(ya.value() + yb.value())
    /// }).unwrap();
    /// assert_eq!(c.ys_raw(), &[11.0, 17.0, 23.0]);
    /// ```
    pub fn zip_with_to<Yo: Unit, Yout: Unit, F>(
        &self,
        other: &SampledSpectrum<X, Yo>,
        f: F,
    ) -> Result<SampledSpectrum<X, Yout>, SpectrumError>
    where
        F: Fn(Quantity<Y>, Quantity<Yo>) -> Quantity<Yout>,
    {
        let mut ys = Vec::with_capacity(self.xs.len());
        for (&x, &ya) in self.xs.iter().zip(self.ys.iter()) {
            let yb = other.eval(x, other.oor)?;
            ys.push(f(Quantity::<Y>::new(ya), Quantity::<Yo>::new(yb)).value());
        }
        let xs: Vec<f64> = self.xs.to_vec();
        SampledSpectrum::<X, Yout>::from_raw(xs, ys, self.interp, self.oor, self.provenance.clone())
    }

    /// Combines two spectra sample-by-sample using a raw `f64` closure.
    ///
    /// This is an escape hatch for callers already working in raw values.
    /// Prefer [`zip_with_to`](Self::zip_with_to) for dimensionally safe code.
    ///
    /// `other` is resampled onto `self`'s x-axis; the output uses `self`'s
    /// interpolation and out-of-range policy.
    ///
    /// # Errors
    ///
    /// Returns [`SpectrumError`] when `other` cannot be evaluated at a point in
    /// `self`'s grid under its own out-of-range policy.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::OutOfRange;
    /// use optica::spectrum::{Interpolation, SampledSpectrum};
    /// use qtty::unit::{Nanometer, Ratio};
    ///
    /// let a = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[0.0, 1.0, 2.0], &[1.0, 2.0, 3.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let b = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
    ///     &[0.0, 2.0], &[10.0, 20.0],
    ///     Interpolation::Linear, OutOfRange::ClampToEndpoints,
    /// ).unwrap();
    /// let c = a.zip_with_raw(&b, |ya, yb| ya + yb).unwrap();
    /// assert_eq!(c.ys_raw(), &[11.0, 17.0, 23.0]);
    /// ```
    pub fn zip_with_raw<Yo: Unit, F>(
        &self,
        other: &SampledSpectrum<X, Yo>,
        f: F,
    ) -> Result<Self, SpectrumError>
    where
        F: Fn(f64, f64) -> f64,
    {
        let mut ys = Vec::with_capacity(self.xs.len());
        for (&x, &ya) in self.xs.iter().zip(self.ys.iter()) {
            let yb = other.eval(x, other.oor)?;
            ys.push(f(ya, yb));
        }
        let xs: Vec<f64> = self.xs.to_vec();
        Self::from_raw(xs, ys, self.interp, self.oor, self.provenance.clone())
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
