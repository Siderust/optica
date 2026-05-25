// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Three-dimensional typed interpolation tables.

use core::marker::PhantomData;

use qtty::{Quantity, Unit};

use crate::data::Provenance;
use crate::grid::algo;
use crate::grid::{GridError, OutOfRange};

/// Three-dimensional lookup table over typed `x`, `y`, and `z` axes.
///
/// Values are stored in z-outermost, y-middle, x-innermost order:
/// `values[(iz * ny + iy) * nx + ix]`.
///
/// # Examples
///
/// ```rust
/// use optica::grid::{Grid3D, OutOfRange};
/// use qtty::{Quantity, unit::{Kilometer, Nanometer, Radian, Ratio}};
///
/// // 2×2×2 grid; storage: (iz*2+iy)*2+ix
/// // iz=0,iy=0: [0, 2]; iz=0,iy=1: [4, 6]; iz=1,iy=0: [8, 10]; iz=1,iy=1: [12, 14]
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
    x_axis: crate::grid::Axis,
    y_axis: crate::grid::Axis,
    z_axis: crate::grid::Axis,
    values: Box<[f64]>,
    out_of_range: OutOfRange,
    provenance: Option<Provenance>,
    _phantom: PhantomData<(X, Y, Z, V)>,
}

type Locate3 = (usize, f64, usize, f64, usize, f64);

impl<X: Unit, Y: Unit, Z: Unit, V: Unit> Grid3D<X, Y, Z, V> {
    /// Builds a validated 3-D grid from ascending axes and z-outermost values.
    ///
    /// Storage convention: `values[(iz * ny + iy) * nx + ix]`.
    ///
    /// # Errors
    ///
    /// Returns [`GridError`] when any axis is invalid or the values length does not match
    /// `xs.len() * ys.len() * zs.len()`.
    pub fn from_raw_row_major(
        xs: &[f64],
        ys: &[f64],
        zs: &[f64],
        values: &[f64],
        oor: OutOfRange,
    ) -> Result<Self, GridError> {
        let x_axis = crate::grid::Axis::NonUniform(xs.to_vec().into_boxed_slice());
        x_axis.validate_for_axis("x")?;
        let y_axis = crate::grid::Axis::NonUniform(ys.to_vec().into_boxed_slice());
        y_axis.validate_for_axis("y")?;
        let z_axis = crate::grid::Axis::NonUniform(zs.to_vec().into_boxed_slice());
        z_axis.validate_for_axis("z")?;
        let expected = x_axis.len() * y_axis.len() * z_axis.len();
        if expected != values.len() {
            return Err(GridError::ShapeMismatch { expected, got: values.len() });
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
            Some((ix, tx, iy, ty, iz, tz)) => {
                Ok(Quantity::new(self.interpolate(ix, tx, iy, ty, iz, tz)))
            }
            None => Ok(Quantity::zero()),
        }
    }

    /// Interpolates a value at `(x, y, z)`, overriding the stored out-of-range policy per axis.
    pub fn interp_at_with(
        &self,
        x: Quantity<X>,
        y: Quantity<Y>,
        z: Quantity<Z>,
        oor_x: OutOfRange,
        oor_y: OutOfRange,
        oor_z: OutOfRange,
    ) -> Result<Quantity<V>, GridError> {
        let xv = x.value();
        let yv = y.value();
        let zv = z.value();
        let (x_lo, x_hi) = self.x_axis.bounds();
        let (y_lo, y_hi) = self.y_axis.bounds();
        let (z_lo, z_hi) = self.z_axis.bounds();
        if !algo::check_oor(xv, x_lo, x_hi, oor_x, "x")? {
            return Ok(Quantity::zero());
        }
        if !algo::check_oor(yv, y_lo, y_hi, oor_y, "y")? {
            return Ok(Quantity::zero());
        }
        if !algo::check_oor(zv, z_lo, z_hi, oor_z, "z")? {
            return Ok(Quantity::zero());
        }
        let (ix, tx) = self.x_axis.locate(xv);
        let (iy, ty) = self.y_axis.locate(yv);
        let (iz, tz) = self.z_axis.locate(zv);
        Ok(Quantity::new(self.interpolate(ix, tx, iy, ty, iz, tz)))
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
        let nx = self.nx();
        let ny = self.ny();
        let idx = |iz: usize, iy: usize, ix: usize| (iz * ny + iy) * nx + ix;
        let v000 = self.values[idx(iz, iy, ix)];
        let v100 = self.values[idx(iz, iy, ix + 1)];
        let v010 = self.values[idx(iz, iy + 1, ix)];
        let v110 = self.values[idx(iz, iy + 1, ix + 1)];
        let v001 = self.values[idx(iz + 1, iy, ix)];
        let v101 = self.values[idx(iz + 1, iy, ix + 1)];
        let v011 = self.values[idx(iz + 1, iy + 1, ix)];
        let v111 = self.values[idx(iz + 1, iy + 1, ix + 1)];
        algo::trilinear_unit(v000, v100, v010, v110, v001, v101, v011, v111, tx, ty, tz)
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
                    let (lo, hi) = self.x_axis.bounds();
                    return Err(GridError::OutOfRange { axis: "x", value: x, lo, hi });
                }
                if !y_in {
                    let (lo, hi) = self.y_axis.bounds();
                    return Err(GridError::OutOfRange { axis: "y", value: y, lo, hi });
                }
                let (lo, hi) = self.z_axis.bounds();
                Err(GridError::OutOfRange { axis: "z", value: z, lo, hi })
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

    /// Trilinear midpoint — passes regardless of storage layout.
    #[test]
    fn trilinear_midpoint_works() {
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

    /// Discriminative test: verifies `(iz*ny+iy)*nx+ix` storage.
    ///
    /// With `values[0..7] = [0,2,4,6,8,10,12,14]` and storage `(iz*2+iy)*2+ix`:
    ///   iz=0,iy=0: V(400)=0, V(500)=2
    ///   iz=0,iy=1: V(400)=4, V(500)=6
    ///   iz=1,iy=0: V(400)=8, V(500)=10
    ///   iz=1,iy=1: V(400)=12, V(500)=14
    /// At (400, 0, 0): must return 0; at (500, 0, 0): must return 2.
    #[test]
    fn storage_is_z_outer_y_mid_x_inner() {
        let grid = Grid3D::<Nanometer, Radian, Kilometer, Ratio>::from_raw_row_major(
            &[400.0, 500.0],
            &[0.0, 1.0],
            &[0.0, 2.0],
            &[0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0],
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();

        // iz=0, iy=0: row [0, 2] → at x=400 → 0, at x=500 → 2
        assert_eq!(
            grid.interp_at(Quantity::new(400.0), Quantity::new(0.0), Quantity::new(0.0))
                .value(),
            0.0
        );
        assert_eq!(
            grid.interp_at(Quantity::new(500.0), Quantity::new(0.0), Quantity::new(0.0))
                .value(),
            2.0
        );
        // iz=0, iy=1: row [4, 6] → at x=400 → 4
        assert_eq!(
            grid.interp_at(Quantity::new(400.0), Quantity::new(1.0), Quantity::new(0.0))
                .value(),
            4.0
        );
        // iz=1, iy=0: row [8, 10] → at x=400 → 8
        assert_eq!(
            grid.interp_at(Quantity::new(400.0), Quantity::new(0.0), Quantity::new(2.0))
                .value(),
            8.0
        );
    }
}
