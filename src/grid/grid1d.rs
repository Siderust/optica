// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! One-dimensional typed interpolation tables.

use core::marker::PhantomData;

use qtty::{Quantity, Unit};

use crate::data::Provenance;
use crate::grid::algo::lerp;
use crate::grid::{Axis, GridError, OutOfRange};

/// One-dimensional lookup table with typed axes and values.
///
/// Axes and values are stored as contiguous `Box<[f64]>`, so `interp_at` performs
/// no heap allocation.
///
/// # Examples
///
/// ```rust
/// use optica::grid::{Grid1D, OutOfRange};
/// use qtty::{Quantity, unit::{Nanometer, Ratio}};
///
/// let grid = Grid1D::<Nanometer, Ratio>::from_sorted(
///     &[400.0, 500.0],
///     &[1.0, 3.0],
///     OutOfRange::ClampToEndpoints,
/// )
/// .unwrap();
///
/// let value = grid.interp_at(Quantity::<Nanometer>::new(450.0));
/// assert!((value.value() - 2.0).abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct Grid1D<X: Unit, V: Unit> {
    axis: Axis,
    values: Box<[f64]>,
    out_of_range: OutOfRange,
    provenance: Option<Provenance>,
    _phantom: PhantomData<(X, V)>,
}

impl<X: Unit, V: Unit> Grid1D<X, V> {
    /// Builds a validated 1-D grid from sorted sample coordinates and values.
    pub fn from_sorted(xs: &[f64], ys: &[f64], oor: OutOfRange) -> Result<Self, GridError> {
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
            out_of_range: oor,
            provenance: None,
            _phantom: PhantomData,
        })
    }

    /// Builds a validated uniformly spaced 1-D grid.
    pub fn uniform(start: f64, step: f64, ys: &[f64], oor: OutOfRange) -> Result<Self, GridError> {
        let axis = Axis::uniform(start, step, ys.len())?;
        Ok(Self {
            axis,
            values: ys.to_vec().into_boxed_slice(),
            out_of_range: oor,
            provenance: None,
            _phantom: PhantomData,
        })
    }

    /// Interpolates a value at `x` using the grid's out-of-range policy.
    ///
    /// When `OutOfRange::Error` is configured, this infallible method clamps to the
    /// nearest endpoint. Use [`try_interp_at`](Self::try_interp_at) to detect that case.
    #[must_use]
    pub fn interp_at(&self, x: Quantity<X>) -> Quantity<V> {
        match self.locate_query(x.value(), false) {
            Ok(Some((low, t))) => Quantity::new(lerp(self.values[low], self.values[low + 1], t)),
            Ok(None) => Quantity::zero(),
            Err(_) => {
                let (low, t) = self.axis.locate(x.value());
                Quantity::new(lerp(self.values[low], self.values[low + 1], t))
            }
        }
    }

    /// Interpolates a value at `x`, returning an error when `OutOfRange::Error` is active.
    pub fn try_interp_at(&self, x: Quantity<X>) -> Result<Quantity<V>, GridError> {
        match self.locate_query(x.value(), true)? {
            Some((low, t)) => Ok(Quantity::new(lerp(
                self.values[low],
                self.values[low + 1],
                t,
            ))),
            None => Ok(Quantity::zero()),
        }
    }

    /// Interpolates a value at `x`, overriding the stored out-of-range policy for this call.
    pub fn interp_at_with(
        &self,
        x: Quantity<X>,
        oor: OutOfRange,
    ) -> Result<Quantity<V>, GridError> {
        let xv = x.value();
        if self.axis.contains(xv) {
            let (low, t) = self.axis.locate(xv);
            return Ok(Quantity::new(lerp(
                self.values[low],
                self.values[low + 1],
                t,
            )));
        }
        match oor {
            OutOfRange::ClampToEndpoints => {
                let (low, t) = self.axis.locate(xv);
                Ok(Quantity::new(lerp(
                    self.values[low],
                    self.values[low + 1],
                    t,
                )))
            }
            OutOfRange::Zero => Ok(Quantity::zero()),
            OutOfRange::Error => {
                let (lo, hi) = self.axis.bounds();
                Err(GridError::OutOfRange {
                    axis: "x",
                    value: xv,
                    lo,
                    hi,
                })
            }
        }
    }

    /// Returns the number of samples in the grid.
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Returns `true` when the grid stores no values.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
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
                    axis: "x",
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
    fn linear_interpolation_works() {
        let grid = Grid1D::<Nanometer, Ratio>::from_sorted(
            &[400.0, 500.0],
            &[2.0, 4.0],
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();

        let value = grid.interp_at(Quantity::<Nanometer>::new(450.0));
        assert_eq!(value.value(), 3.0);
    }

    #[test]
    fn zero_policy_returns_zero() {
        let grid =
            Grid1D::<Nanometer, Ratio>::from_sorted(&[400.0, 500.0], &[2.0, 4.0], OutOfRange::Zero)
                .unwrap();

        let value = grid.interp_at(Quantity::<Nanometer>::new(300.0));
        assert_eq!(value.value(), 0.0);
    }
}
