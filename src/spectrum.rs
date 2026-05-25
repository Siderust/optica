// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Sampled spectra with allocation-free interpolation.

use core::marker::PhantomData;

use qtty::{Quantity, Unit};

use crate::data::Provenance;
use crate::grid::algo::lerp;
use crate::grid::{Axis, GridError, OutOfRange};

/// Interpolation strategy for [`SampledSpectrum`].
///
/// # Examples
///
/// ```rust
/// use optica::spectrum::Interpolation;
///
/// assert_ne!(Interpolation::Linear, Interpolation::NearestNeighbor);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interpolation {
    /// Piecewise linear interpolation.
    Linear,
    /// Nearest-neighbor sampling.
    NearestNeighbor,
}

/// A performance-focused sampled spectrum.
///
/// The spectrum stores both the sample coordinates and values in contiguous
/// `Box<[f64]>` buffers. Query methods perform no heap allocation.
///
/// # Examples
///
/// ```rust
/// use optica::grid::OutOfRange;
/// use optica::spectrum::{Interpolation, SampledSpectrum};
/// use qtty::length::Nanometers;
///
/// let spectrum = SampledSpectrum::<qtty::unit::Nanometer, qtty::unit::Ratio>::from_sorted(
///     &[400.0, 500.0],
///     &[0.2, 0.4],
///     Interpolation::Linear,
///     OutOfRange::ClampToEndpoints,
/// )
/// .unwrap();
///
/// let sample = spectrum.interp_at(Nanometers::new(450.0));
/// assert!((sample.value() - 0.3).abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct SampledSpectrum<X: Unit, Y: Unit> {
    axis: Axis,
    values: Box<[f64]>,
    interp: Interpolation,
    out_of_range: OutOfRange,
    provenance: Option<Provenance>,
    _phantom: PhantomData<(X, Y)>,
}

impl<X: Unit, Y: Unit> SampledSpectrum<X, Y> {
    /// Builds a spectrum from sorted sample coordinates and values.
    pub fn from_sorted(
        xs: &[f64],
        ys: &[f64],
        interp: Interpolation,
        oor: OutOfRange,
    ) -> Result<Self, GridError> {
        let axis = Axis::non_uniform(xs.to_vec().into_boxed_slice())?;
        if axis.len() != ys.len() {
            return Err(GridError::ShapeMismatch {
                expected: axis.len(),
                got: ys.len(),
            });
        }
        Ok(Self {
            axis,
            values: ys.to_vec().into_boxed_slice(),
            interp,
            out_of_range: oor,
            provenance: None,
            _phantom: PhantomData,
        })
    }

    /// Samples the spectrum at `x` using the configured interpolation mode.
    #[must_use]
    pub fn interp_at(&self, x: Quantity<X>) -> Quantity<Y> {
        match self.locate_query(x.value(), false) {
            Ok(Some((low, t))) => Quantity::new(self.interpolate(low, t)),
            Ok(None) => Quantity::zero(),
            Err(_) => {
                let (low, t) = self.axis.locate(x.value());
                Quantity::new(self.interpolate(low, t))
            }
        }
    }

    /// Samples the spectrum at `x`, returning an error when requested.
    pub fn try_interp_at(&self, x: Quantity<X>) -> Result<Quantity<Y>, GridError> {
        match self.locate_query(x.value(), true)? {
            Some((low, t)) => Ok(Quantity::new(self.interpolate(low, t))),
            None => Ok(Quantity::zero()),
        }
    }

    /// Returns the raw sample coordinates without allocation.
    #[must_use]
    pub fn xs_raw(&self) -> &[f64] {
        self.axis
            .raw_points()
            .expect("SampledSpectrum stores explicit sample coordinates")
    }

    /// Returns the raw sample values without allocation.
    #[must_use]
    pub fn ys_raw(&self) -> &[f64] {
        &self.values
    }

    /// Returns the number of samples in the spectrum.
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Returns `true` when the spectrum stores no samples.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Integrates the spectrum with the trapezoidal rule over the raw sample axis.
    #[must_use]
    pub fn integrate(&self) -> f64 {
        let xs = self.xs_raw();
        xs.windows(2)
            .zip(self.values.windows(2))
            .map(|(xw, yw)| (xw[1] - xw[0]) * 0.5 * (yw[0] + yw[1]))
            .sum()
    }

    /// Attaches provenance metadata.
    #[must_use]
    pub fn with_provenance(mut self, provenance: Provenance) -> Self {
        self.provenance = Some(provenance);
        self
    }

    /// Returns the attached provenance metadata, if any.
    #[must_use]
    pub fn provenance(&self) -> Option<&Provenance> {
        self.provenance.as_ref()
    }

    fn interpolate(&self, low: usize, t: f64) -> f64 {
        match self.interp {
            Interpolation::Linear => lerp(self.values[low], self.values[low + 1], t),
            Interpolation::NearestNeighbor => {
                let idx = if t < 0.5 { low } else { low + 1 };
                self.values[idx]
            }
        }
    }

    fn locate_query(&self, x: f64, strict_error: bool) -> Result<Option<(usize, f64)>, GridError> {
        if self.axis.contains(x) {
            return Ok(Some(self.axis.locate(x)));
        }
        match self.out_of_range {
            OutOfRange::ClampToEndpoints => Ok(Some(self.axis.locate(x))),
            OutOfRange::Zero => Ok(None),
            OutOfRange::Error if strict_error => {
                let (lo, hi) = self.axis.bounds();
                Err(GridError::OutOfRange {
                    axis: "wavelength",
                    value: x,
                    lo,
                    hi,
                })
            }
            OutOfRange::Error => Ok(Some(self.axis.locate(x))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::unit::{Nanometer, Ratio};

    #[test]
    fn trapz_matches_linear_ramp() {
        let spectrum = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
            &[0.0, 1.0, 2.0],
            &[0.0, 1.0, 2.0],
            Interpolation::Linear,
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();

        assert_eq!(spectrum.integrate(), 2.0);
    }

    #[test]
    fn nearest_neighbor_uses_closest_sample() {
        let spectrum = SampledSpectrum::<Nanometer, Ratio>::from_sorted(
            &[400.0, 500.0],
            &[2.0, 4.0],
            Interpolation::NearestNeighbor,
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();

        let value = spectrum.interp_at(Quantity::<Nanometer>::new(490.0));
        assert_eq!(value.value(), 4.0);
    }
}
