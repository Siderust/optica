// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Three-dimensional typed interpolation tables.

use core::marker::PhantomData;

use qtty::{Quantity, Unit};

use crate::data::Provenance;
use crate::grid::algo::trilerp;
use crate::grid::{Axis, GridError, OutOfRange};

/// Three-dimensional lookup table over typed `x`, `y`, and `z` axes.
///
/// Values are stored in row-major order as `values[ix * (ny * nz) + iy * nz + iz]`.
///
/// # Examples
///
/// ```rust
/// use optica::grid::{Grid3D, OutOfRange};
/// use qtty::{Quantity, unit::{Kilometer, Nanometer, Radian, Ratio}};
///
/// let grid = Grid3D::<Nanometer, Radian, Kilometer, Ratio>::from_raw_row_major(
///     &[400.0, 500.0],
///     &[0.0, 1.0],
///     &[0.0, 2.0],
///     &[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0],
///     OutOfRange::ClampToEndpoints,
/// )
/// .unwrap();
///
/// let value = grid.interp_at(
///     Quantity::<Nanometer>::new(450.0),
///     Quantity::<Radian>::new(0.5),
///     Quantity::<Kilometer>::new(1.0),
/// );
/// assert_eq!(value.value(), 7.0);
/// ```
#[derive(Debug, Clone)]
pub struct Grid3D<X: Unit, Y: Unit, Z: Unit, V: Unit> {
    x_axis: Axis,
    y_axis: Axis,
    z_axis: Axis,
    values: Box<[f64]>,
    out_of_range: OutOfRange,
    provenance: Option<Provenance>,
    _phantom: PhantomData<(X, Y, Z, V)>,
}

type Locate3 = (usize, f64, usize, f64, usize, f64);

impl<X: Unit, Y: Unit, Z: Unit, V: Unit> Grid3D<X, Y, Z, V> {
    /// Builds a validated 3-D grid from sorted axes and row-major values.
    pub fn from_raw_row_major(
        xs: &[f64],
        ys: &[f64],
        zs: &[f64],
        values: &[f64],
        oor: OutOfRange,
    ) -> Result<Self, GridError> {
        let x_axis = Axis::NonUniform(xs.to_vec().into_boxed_slice());
        x_axis.validate_for_axis(0)?;
        let y_axis = Axis::NonUniform(ys.to_vec().into_boxed_slice());
        y_axis.validate_for_axis(1)?;
        let z_axis = Axis::NonUniform(zs.to_vec().into_boxed_slice());
        z_axis.validate_for_axis(2)?;
        let expected = x_axis.len() * y_axis.len() * z_axis.len();
        if expected != values.len() {
            return Err(GridError::ShapeMismatch {
                expected,
                got: values.len(),
            });
        }
        Ok(Self {
            x_axis,
            y_axis,
            z_axis,
            values: values.to_vec().into_boxed_slice(),
            out_of_range: oor,
            provenance: None,
            _phantom: PhantomData,
        })
    }

    /// Interpolates a value at `(x, y, z)` using the configured out-of-range policy.
    #[must_use]
    pub fn interp_at(&self, x: Quantity<X>, y: Quantity<Y>, z: Quantity<Z>) -> Quantity<V> {
        match self.locate_query(x.value(), y.value(), z.value(), false) {
            Ok(Some((ix, tx, iy, ty, iz, tz))) => {
                Quantity::new(self.interpolate(ix, tx, iy, ty, iz, tz))
            }
            Ok(None) => Quantity::zero(),
            Err(_) => {
                let (ix, tx) = self.x_axis.locate(x.value());
                let (iy, ty) = self.y_axis.locate(y.value());
                let (iz, tz) = self.z_axis.locate(z.value());
                Quantity::new(self.interpolate(ix, tx, iy, ty, iz, tz))
            }
        }
    }

    /// Interpolates a value at `(x, y, z)`, returning an error when requested.
    pub fn try_interp_at(
        &self,
        x: Quantity<X>,
        y: Quantity<Y>,
        z: Quantity<Z>,
    ) -> Result<Quantity<V>, GridError> {
        match self.locate_query(x.value(), y.value(), z.value(), true)? {
            Some((ix, tx, iy, ty, iz, tz)) => Ok(Quantity::new(self.interpolate(ix, tx, iy, ty, iz, tz))),
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

    /// Returns the number of samples on the `z` axis.
    #[must_use]
    pub fn nz(&self) -> usize {
        self.z_axis.len()
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

    fn interpolate(&self, ix: usize, tx: f64, iy: usize, ty: f64, iz: usize, tz: f64) -> f64 {
        let ny = self.ny();
        let nz = self.nz();
        let base0 = ix * (ny * nz);
        let base1 = (ix + 1) * (ny * nz);
        let row0 = iy * nz;
        let row1 = (iy + 1) * nz;
        let v000 = self.values[base0 + row0 + iz];
        let v001 = self.values[base0 + row0 + (iz + 1)];
        let v010 = self.values[base0 + row1 + iz];
        let v011 = self.values[base0 + row1 + (iz + 1)];
        let v100 = self.values[base1 + row0 + iz];
        let v101 = self.values[base1 + row0 + (iz + 1)];
        let v110 = self.values[base1 + row1 + iz];
        let v111 = self.values[base1 + row1 + (iz + 1)];
        trilerp(v000, v001, v010, v011, v100, v101, v110, v111, tx, ty, tz)
    }

    fn locate_query(
        &self,
        x: f64,
        y: f64,
        z: f64,
        strict_error: bool,
    ) -> Result<Option<Locate3>, GridError> {
        let x_in = self.x_axis.contains(x);
        let y_in = self.y_axis.contains(y);
        let z_in = self.z_axis.contains(z);
        if x_in && y_in && z_in {
            let (ix, tx) = self.x_axis.locate(x);
            let (iy, ty) = self.y_axis.locate(y);
            let (iz, tz) = self.z_axis.locate(z);
            return Ok(Some((ix, tx, iy, ty, iz, tz)));
        }

        match self.out_of_range {
            OutOfRange::ClampToEndpoints => {
                let (ix, tx) = self.x_axis.locate(x);
                let (iy, ty) = self.y_axis.locate(y);
                let (iz, tz) = self.z_axis.locate(z);
                Ok(Some((ix, tx, iy, ty, iz, tz)))
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
                if !y_in {
                    let (min, max) = self.y_axis.bounds();
                    return Err(GridError::OutOfRange {
                        axis: 1,
                        value: y,
                        min,
                        max,
                    });
                }
                let (min, max) = self.z_axis.bounds();
                Err(GridError::OutOfRange {
                    axis: 2,
                    value: z,
                    min,
                    max,
                })
            }
            OutOfRange::Error => {
                let (ix, tx) = self.x_axis.locate(x);
                let (iy, ty) = self.y_axis.locate(y);
                let (iz, tz) = self.z_axis.locate(z);
                Ok(Some((ix, tx, iy, ty, iz, tz)))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::unit::{Kilometer, Nanometer, Radian, Ratio};

    #[test]
    fn trilinear_interpolation_works() {
        let grid = Grid3D::<Nanometer, Radian, Kilometer, Ratio>::from_raw_row_major(
            &[400.0, 500.0],
            &[0.0, 1.0],
            &[0.0, 2.0],
            &[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0],
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();

        let value = grid.interp_at(
            Quantity::<Nanometer>::new(450.0),
            Quantity::<Radian>::new(0.5),
            Quantity::<Kilometer>::new(1.0),
        );
        assert_eq!(value.value(), 7.0);
    }
}
