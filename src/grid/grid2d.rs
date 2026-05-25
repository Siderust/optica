// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Two-dimensional typed interpolation tables.

use core::marker::PhantomData;

use qtty::{Quantity, Unit};

use crate::data::Provenance;
use crate::grid::algo;
use crate::grid::{AxisDirection, GridError, OutOfRange};

/// A constant-value region that short-circuits interpolation.
///
/// When a query point `(x, y)` falls within the bounds defined by this region,
/// `Grid2D` returns `value` directly without interpolating the underlying table.
///
/// # Examples
///
/// ```rust
/// use optica::grid::ConstantRegion;
///
/// let region = ConstantRegion::lower_corner(10.0, 5.0, 0.0);
/// assert!(region.contains(5.0, 3.0));
/// assert!(!region.contains(5.0, 6.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantRegion {
    /// When set, the region applies only when `x <= x_upper_bound`.
    pub x_upper_bound: Option<f64>,
    /// When set, the region applies only when `y <= y_upper_bound`.
    pub y_upper_bound: Option<f64>,
    /// Value returned when the region matches.
    pub value: f64,
}

impl ConstantRegion {
    /// Creates a lower-corner constant region active when `x ≤ x_max` **and** `y ≤ y_max`.
    #[must_use]
    pub fn lower_corner(x_max: f64, y_max: f64, value: f64) -> Self {
        Self {
            x_upper_bound: Some(x_max),
            y_upper_bound: Some(y_max),
            value,
        }
    }

    /// Returns `true` when `(x, y)` falls within this region.
    #[must_use]
    pub fn contains(&self, x: f64, y: f64) -> bool {
        let x_ok = self.x_upper_bound.is_none_or(|xb| x <= xb);
        let y_ok = self.y_upper_bound.is_none_or(|yb| y <= yb);
        x_ok && y_ok
    }
}

/// Two-dimensional lookup table over typed `x` and `y` axes.
///
/// Values are stored in y-major (row-major) order: `values[iy * nx + ix]`, so
/// each contiguous row corresponds to a fixed `y` value varying over `x`.
///
/// Both ascending and descending axes are supported. The y-descending
/// constructor ([`from_raw_row_major_y_descending`]) accepts Leinert-style
/// tables whose first axis is the highest physical value without requiring
/// callers to reorder the data.
///
/// # Examples
///
/// ```rust
/// use optica::grid::{Grid2D, OutOfRange};
/// use qtty::{Quantity, unit::{Nanometer, Radian, Ratio}};
///
/// // Row-major: row 0 (y=0.0): V(400)=1, V(500)=2; row 1 (y=1.0): V(400)=3, V(500)=4
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
///
/// [`from_raw_row_major_y_descending`]: Self::from_raw_row_major_y_descending
#[derive(Debug, Clone)]
pub struct Grid2D<X: Unit, Y: Unit, V: Unit> {
    xs: Box<[f64]>,
    /// Internal y-axis, always ascending (either as-supplied or reversed).
    ys: Box<[f64]>,
    /// y-major flat storage: `values[iy * nx + ix]`.
    values: Box<[f64]>,
    nx: usize,
    ny: usize,
    dir_x: AxisDirection,
    /// Direction of the internal ys (always `Ascending`; descending inputs are normalized).
    dir_y: AxisDirection,
    /// For y-descending inputs: `ys_desc[0] + ys_desc[ny-1]`.
    y_reflect_offset: Option<f64>,
    regions: Vec<ConstantRegion>,
    out_of_range: OutOfRange,
    provenance: Option<Provenance>,
    _phantom: PhantomData<(X, Y, V)>,
}

impl<X: Unit, Y: Unit, V: Unit> Grid2D<X, Y, V> {
    /// Builds a validated 2-D grid from sorted axes and y-major row-major values.
    ///
    /// Both axes may be ascending or descending. Values are stored as
    /// `values[iy * nx + ix]`: the first `nx` entries are the first row (`y = ys[0]`).
    ///
    /// # Errors
    ///
    /// Returns [`GridError`] when either axis is invalid or the values length does not
    /// match `xs.len() * ys.len()`.
    pub fn from_raw_row_major(
        xs: &[f64],
        ys: &[f64],
        values: &[f64],
        oor: OutOfRange,
    ) -> Result<Self, GridError> {
        let dir_x = algo::validate_axis("x", xs)?;
        let dir_y = algo::validate_axis("y", ys)?;
        let nx = xs.len();
        let ny = ys.len();
        let expected = nx * ny;
        if expected != values.len() {
            return Err(GridError::ShapeMismatch { expected, got: values.len() });
        }
        Ok(Self {
            xs: xs.to_vec().into_boxed_slice(),
            ys: ys.to_vec().into_boxed_slice(),
            values: values.to_vec().into_boxed_slice(),
            nx,
            ny,
            dir_x,
            dir_y,
            y_reflect_offset: None,
            regions: Vec::new(),
            out_of_range: oor,
            provenance: None,
            _phantom: PhantomData,
        })
    }

    /// Builds a validated 2-D grid where the y-axis is strictly descending and uniform.
    ///
    /// This constructor is intended for Leinert-style tables that store the highest
    /// physical y-value first. It accepts a y-descending axis and keeps the value rows in
    /// their original order by applying a reflection transform at query time:
    /// `y_internal = (ys_desc[0] + ys_desc[ny-1]) − y_query`.
    ///
    /// The y-axis must be strictly descending **and** uniformly spaced (equal absolute
    /// step between all consecutive pairs, compared bit-for-bit in IEEE 754).
    ///
    /// # Errors
    ///
    /// Returns [`GridError`] when the y-axis is invalid, not uniform, or the value count
    /// does not match `xs.len() * ys_desc.len()`.
    pub fn from_raw_row_major_y_descending(
        xs: &[f64],
        ys_desc: &[f64],
        values: &[f64],
    ) -> Result<Self, GridError> {
        let dir_x = algo::validate_axis("x", xs)?;
        let nx = xs.len();
        let ny = ys_desc.len();
        if ny < 2 {
            return Err(GridError::TooFewSamples { axis: "y", len: ny });
        }
        for (i, &v) in ys_desc.iter().enumerate() {
            if !v.is_finite() {
                return Err(GridError::NonFinite { axis: "y", index: i });
            }
        }
        // Validate strictly descending
        for i in 1..ny {
            if ys_desc[i] >= ys_desc[i - 1] {
                return Err(GridError::NotMonotonic { axis: "y", at_index: i });
            }
        }
        // Validate uniform step (bit-exact IEEE 754 comparison matching siderust parity)
        let step = ys_desc[0] - ys_desc[1];
        for i in 1..ny {
            let got = ys_desc[i - 1] - ys_desc[i];
            if got != step {
                return Err(GridError::NonUniformStep { axis: "y", expected: step, got });
            }
        }
        let expected = nx * ny;
        if expected != values.len() {
            return Err(GridError::ShapeMismatch { expected, got: values.len() });
        }
        // Reverse y-axis to ascending but keep rows in their original descending order.
        // The reflection y_internal = (ys_desc[0] + ys_desc[ny-1]) - y_user maps user
        // coordinates back to the original row indices without flipping the value array.
        let y_reflect_offset = ys_desc[0] + ys_desc[ny - 1];
        let ys_asc: Box<[f64]> = ys_desc.iter().copied().rev().collect();
        Ok(Self {
            xs: xs.to_vec().into_boxed_slice(),
            ys: ys_asc,
            values: values.to_vec().into_boxed_slice(),
            nx,
            ny,
            dir_x,
            dir_y: AxisDirection::Ascending,
            y_reflect_offset: Some(y_reflect_offset),
            regions: Vec::new(),
            out_of_range: OutOfRange::ClampToEndpoints,
            provenance: None,
            _phantom: PhantomData,
        })
    }

    /// Interpolates a value at `(x, y)` using the stored out-of-range policy.
    ///
    /// When `OutOfRange::Error` is configured, this infallible method clamps to the
    /// nearest endpoint. Use [`try_interp_at`](Self::try_interp_at) to detect that case.
    #[must_use]
    pub fn interp_at(&self, x: Quantity<X>, y: Quantity<Y>) -> Quantity<V> {
        self.eval(x.value(), y.value(), self.out_of_range, self.out_of_range)
            .map(Quantity::new)
            .unwrap_or_else(|_| {
                // OOR::Error configured but we're infallible: clamp
                Quantity::new(self.eval_clamped(x.value(), y.value()))
            })
    }

    /// Interpolates a value at `(x, y)`, returning an error when `OutOfRange::Error` is active.
    pub fn try_interp_at(&self, x: Quantity<X>, y: Quantity<Y>) -> Result<Quantity<V>, GridError> {
        Ok(Quantity::new(self.eval(x.value(), y.value(), self.out_of_range, self.out_of_range)?))
    }

    /// Interpolates a value at `(x, y)`, overriding the stored out-of-range policy per axis.
    pub fn interp_at_with(
        &self,
        x: Quantity<X>,
        y: Quantity<Y>,
        oor_x: OutOfRange,
        oor_y: OutOfRange,
    ) -> Result<Quantity<V>, GridError> {
        Ok(Quantity::new(self.eval(x.value(), y.value(), oor_x, oor_y)?))
    }

    /// Appends a constant-value region.
    #[must_use]
    pub fn with_constant_region(mut self, region: ConstantRegion) -> Self {
        self.regions.push(region);
        self
    }

    /// Attaches provenance metadata.
    #[must_use]
    pub fn with_provenance(mut self, provenance: Provenance) -> Self {
        self.provenance = Some(provenance);
        self
    }

    /// Returns the number of samples on the `x` axis.
    #[must_use]
    pub fn nx(&self) -> usize {
        self.nx
    }

    /// Returns the number of samples on the `y` axis.
    #[must_use]
    pub fn ny(&self) -> usize {
        self.ny
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

    /// Returns the attached provenance metadata, if any.
    #[must_use]
    pub fn provenance(&self) -> Option<&Provenance> {
        self.provenance.as_ref()
    }

    fn eval(&self, xv: f64, yv: f64, oor_x: OutOfRange, oor_y: OutOfRange) -> Result<f64, GridError> {
        // Constant regions are checked in the user's original (possibly reflected) y coordinates
        for region in &self.regions {
            if region.contains(xv, yv) {
                return Ok(region.value);
            }
        }

        let yv_internal = match self.y_reflect_offset {
            Some(offset) => offset - yv,
            None => yv,
        };

        let (x_lo, x_hi) = algo::axis_range(&self.xs, self.dir_x);
        let (y_lo, y_hi) = algo::axis_range(&self.ys, self.dir_y);

        if !algo::check_oor(xv, x_lo, x_hi, oor_x, "x")? {
            return Ok(0.0);
        }
        if !algo::check_oor(yv_internal, y_lo, y_hi, oor_y, "y")? {
            return Ok(0.0);
        }

        let (ix0, tx) = algo::locate_dir(&self.xs, xv, self.dir_x);
        let (iy0, ty) = algo::locate_dir(&self.ys, yv_internal, self.dir_y);
        let nx = self.nx;
        let f00 = self.values[iy0 * nx + ix0];
        let f10 = self.values[iy0 * nx + (ix0 + 1)];
        let f01 = self.values[(iy0 + 1) * nx + ix0];
        let f11 = self.values[(iy0 + 1) * nx + (ix0 + 1)];
        Ok(algo::bilinear_unit(f00, f10, f01, f11, tx, ty))
    }

    fn eval_clamped(&self, xv: f64, yv: f64) -> f64 {
        for region in &self.regions {
            if region.contains(xv, yv) {
                return region.value;
            }
        }
        let yv_internal = match self.y_reflect_offset {
            Some(offset) => offset - yv,
            None => yv,
        };
        let (ix0, tx) = algo::locate_dir(&self.xs, xv, self.dir_x);
        let (iy0, ty) = algo::locate_dir(&self.ys, yv_internal, self.dir_y);
        let nx = self.nx;
        let ix1 = (ix0 + 1).min(self.nx - 1);
        let iy1 = (iy0 + 1).min(self.ny - 1);
        let f00 = self.values[iy0 * nx + ix0];
        let f10 = self.values[iy0 * nx + ix1];
        let f01 = self.values[iy1 * nx + ix0];
        let f11 = self.values[iy1 * nx + ix1];
        algo::bilinear_unit(f00, f10, f01, f11, tx, ty)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::unit::{Nanometer, Radian, Ratio};

    /// Bilinear midpoint — passes regardless of storage layout.
    #[test]
    fn bilinear_midpoint_works() {
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

    /// Discriminative test: verifies y-major storage `values[iy*nx+ix]`.
    ///
    /// With `values = [1, 2, 3, 4]` and y-major storage:
    ///   row 0 (y=0): V(400)=1, V(500)=2 → at (450, 0) = 1.5
    ///   row 1 (y=1): V(400)=3, V(500)=4
    /// x-major storage would give V(400,0)=1, V(400,1)=2, V(500,0)=3 → at (450,0) = 2.0.
    #[test]
    fn storage_is_y_major() {
        let grid = Grid2D::<Nanometer, Radian, Ratio>::from_raw_row_major(
            &[400.0, 500.0],
            &[0.0, 1.0],
            &[1.0, 2.0, 3.0, 4.0],
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();

        // At y=0 (first row): lerp(V(400,0)=1, V(500,0)=2, 0.5) = 1.5
        let v = grid.interp_at(Quantity::<Nanometer>::new(450.0), Quantity::<Radian>::new(0.0));
        assert!((v.value() - 1.5).abs() < 1e-12, "expected 1.5, got {}", v.value());

        // At y=1 (second row): lerp(V(400,1)=3, V(500,1)=4, 0.5) = 3.5
        let v = grid.interp_at(Quantity::<Nanometer>::new(450.0), Quantity::<Radian>::new(1.0));
        assert!((v.value() - 3.5).abs() < 1e-12, "expected 3.5, got {}", v.value());
    }

    #[test]
    fn corner_values_exact() {
        let grid = Grid2D::<Nanometer, Radian, Ratio>::from_raw_row_major(
            &[400.0, 500.0],
            &[0.0, 1.0],
            &[10.0, 20.0, 30.0, 40.0],
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();

        // V(400,0)=10, V(500,0)=20, V(400,1)=30, V(500,1)=40
        assert_eq!(grid.interp_at(Quantity::new(400.0), Quantity::new(0.0)).value(), 10.0);
        assert_eq!(grid.interp_at(Quantity::new(500.0), Quantity::new(0.0)).value(), 20.0);
        assert_eq!(grid.interp_at(Quantity::new(400.0), Quantity::new(1.0)).value(), 30.0);
        assert_eq!(grid.interp_at(Quantity::new(500.0), Quantity::new(1.0)).value(), 40.0);
    }

    #[test]
    fn y_descending_matches_ascending_equivalent() {
        // Build a descending grid with the same physical values as an ascending grid.
        // ys_desc = [1.0, 0.0]; row 0 (y=1.0): [10, 20]; row 1 (y=0.0): [30, 40]
        let desc = Grid2D::<Nanometer, Radian, Ratio>::from_raw_row_major_y_descending(
            &[400.0, 500.0],
            &[1.0, 0.0],
            &[10.0, 20.0, 30.0, 40.0],
        )
        .unwrap();

        // The ascending equivalent:
        // ys_asc = [0.0, 1.0]; row 0 (y=0.0): [30, 40]; row 1 (y=1.0): [10, 20]
        let asc = Grid2D::<Nanometer, Radian, Ratio>::from_raw_row_major(
            &[400.0, 500.0],
            &[0.0, 1.0],
            &[30.0, 40.0, 10.0, 20.0],
            OutOfRange::ClampToEndpoints,
        )
        .unwrap();

        for &xv in &[400.0_f64, 430.0, 450.0, 480.0, 500.0] {
            for &yv in &[0.0_f64, 0.25, 0.5, 0.75, 1.0] {
                let vd = desc.interp_at(Quantity::new(xv), Quantity::new(yv)).value();
                let va = asc.interp_at(Quantity::new(xv), Quantity::new(yv)).value();
                assert!(
                    (vd - va).abs() < 1e-10,
                    "mismatch at ({xv},{yv}): desc={vd}, asc={va}"
                );
            }
        }
    }

    #[test]
    fn constant_region_short_circuits() {
        let grid = Grid2D::<Nanometer, Radian, Ratio>::from_raw_row_major(
            &[400.0, 500.0],
            &[0.0, 1.0],
            &[1.0, 2.0, 3.0, 4.0],
            OutOfRange::ClampToEndpoints,
        )
        .unwrap()
        .with_constant_region(ConstantRegion::lower_corner(420.0, 0.3, 0.0));

        let v = grid.interp_at(Quantity::new(410.0), Quantity::new(0.1)).value();
        assert_eq!(v, 0.0, "constant region should return 0.0");

        let v = grid.interp_at(Quantity::new(450.0), Quantity::new(0.5)).value();
        assert_ne!(v, 0.0, "outside region should interpolate normally");
    }
}
