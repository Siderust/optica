// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Two-dimensional typed interpolation tables.

use core::marker::PhantomData;

use qtty::{Quantity, Unit};

use crate::data::Provenance;
use crate::grid::algo::bilerp;
use crate::grid::{Axis, GridError, OutOfRange};

/// Two-dimensional lookup table over typed `x` and `y` axes.
///
/// Values are stored in row-major order as `values[ix * ny + iy]`.
///
/// # Examples
///
/// ```rust
/// use optica::grid::{Grid2D, OutOfRange};
/// use qtty::{Quantity, unit::{Nanometer, Radian, Ratio}};
///
/// let grid = Grid2D::<Nanometer, Radian, Ratio>::from_raw_row_major(
///     &[400.0, 500.0],
///     &[0.0, 1.0],
///     &[1.0, 2.0, 3.0, 4.0],
///     OutOfRange::ClampToEndpoints,
/// )
/// .unwrap();
///
/// let value = grid.interp_at(Quantity::<Nanometer>::new(450.0), Quantity::<Radian>::new(0.5));
/// assert!((value.value() - 2.5).abs() < 1e-12);
/// ```
#[derive(Debug, Clone)]
pub struct Grid2D<X: Unit, Y: Unit, V: Unit> {
    x_axis: Axis,
    y_axis: Axis,
    values: Box<[f64]>,
    out_of_range: OutOfRange,
    provenance: Option<Provenance>,
    _phantom: PhantomData<(X, Y, V)>,
}

impl<X: Unit, Y: Unit, V: Unit> Grid2D<X, Y, V> {
    /// Builds a validated 2-D grid from sorted axes and row-major values.
    pub fn from_raw_row_major(
        xs: &[f64],
        ys: &[f64],
        values: &[f64],
        oor: OutOfRange,
    ) -> Result<Self, GridError> {
        let x_axis = Axis::NonUniform(xs.to_vec().into_boxed_slice());
        x_axis.validate_for_axis(0)?;
        let y_axis = Axis::NonUniform(ys.to_vec().into_boxed_slice());
        y_axis.validate_for_axis(1)?;
        let expected = x_axis.len() * y_axis.len();
        if expected != values.len() {
            return Err(GridError::ShapeMismatch {
                expected,
                got: values.len(),
            });
        }
        Ok(Self {
            x_axis,
            y_axis,
            values: values.to_vec().into_boxed_slice(),
            out_of_range: oor,
            provenance: None,
            _phantom: PhantomData,
        })
    }

    /// Interpolates a value at `(x, y)` using the configured out-of-range policy.
    #[must_use]
    pub fn interp_at(&self, x: Quantity<X>, y: Quantity<Y>) -> Quantity<V> {
        match self.locate_query(x.value(), y.value(), false) {
            Ok(Some((ix, tx, iy, ty))) => Quantity::new(self.interpolate(ix, tx, iy, ty)),
            Ok(None) => Quantity::zero(),
            Err(_) => {
                let (ix, tx) = self.x_axis.locate(x.value());
                let (iy, ty) = self.y_axis.locate(y.value());
                Quantity::new(self.interpolate(ix, tx, iy, ty))
            }
        }
    }

    /// Interpolates a value at `(x, y)`, returning an error when requested.
    pub fn try_interp_at(&self, x: Quantity<X>, y: Quantity<Y>) -> Result<Quantity<V>, GridError> {
        match self.locate_query(x.value(), y.value(), true)? {
            Some((ix, tx, iy, ty)) => Ok(Quantity::new(self.interpolate(ix, tx, iy, ty))),
            None => Ok(Quantity::zero()),
        }
    }

    /// Returns the number of samples on the `x` axis.
    #[must_use]
    pub fn nx(&self) -> usize {
        self.x_axis.len()
    }

    /// Returns the number of samples on the `y` axis.
    #[must_use]
    pub fn ny(&self) -> usize {
        self.y_axis.len()
    }

    /// Returns the total number of stored values.
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

    fn interpolate(&self, ix: usize, tx: f64, iy: usize, ty: f64) -> f64 {
        let ny = self.ny();
        let v00 = self.values[ix * ny + iy];
        let v01 = self.values[ix * ny + (iy + 1)];
        let v10 = self.values[(ix + 1) * ny + iy];
        let v11 = self.values[(ix + 1) * ny + (iy + 1)];
        bilerp(v00, v01, v10, v11, tx, ty)
    }

    fn locate_query(
        &self,
        x: f64,
        y: f64,
        strict_error: bool,
    ) -> Result<Option<(usize, f64, usize, f64)>, GridError> {
        let x_in = self.x_axis.contains(x);
        let y_in = self.y_axis.contains(y);
        if x_in && y_in {
            let (ix, tx) = self.x_axis.locate(x);
            let (iy, ty) = self.y_axis.locate(y);
            return Ok(Some((ix, tx, iy, ty)));
        }

        match self.out_of_range {
            OutOfRange::ClampToEndpoints => {
                let (ix, tx) = self.x_axis.locate(x);
                let (iy, ty) = self.y_axis.locate(y);
                Ok(Some((ix, tx, iy, ty)))
            }
            OutOfRange::Zero => Ok(None),
            OutOfRange::Error if strict_error => {
                if !x_in {
                    let (min, max) = self.x_axis.bounds();
                    return Err(GridError::OutOfRange {
                        axis: 0,
                        value: x,
                        min,
                        max,
                    });
                }
                let (min, max) = self.y_axis.bounds();
                Err(GridError::OutOfRange {
                    axis: 1,
                    value: y,
                    min,
                    max,
                })
            }
            OutOfRange::Error => {
                let (ix, tx) = self.x_axis.locate(x);
                let (iy, ty) = self.y_axis.locate(y);
                Ok(Some((ix, tx, iy, ty)))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::unit::{Nanometer, Radian, Ratio};

    #[test]
    fn bilinear_interpolation_works() {
        let grid = Grid2D::<Nanometer, Radian, Ratio>::from_raw_row_major(
            &[400.0, 500.0],
            &[0.0, 1.0],
            &[1.0, 2.0, 3.0, 4.0],
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();

        let value = grid.interp_at(Quantity::<Nanometer>::new(450.0), Quantity::<Radian>::new(0.5));
        assert_eq!(value.value(), 2.5);
    }
}
