// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Low-level `f64` interpolation kernels shared by typed grid wrappers.
//!
//! The crate-internal kernels (`locate`, `locate_uniform`, `lerp`, `bilerp`,
//! `trilerp`) are used by the typed `Axis`-based `Grid1D`. The public
//! direction-aware functions (`validate_axis`, `linear_1d`, `bilinear`,
//! `bilinear_unit`, `trilinear`, `trilinear_unit`) are the canonical
//! interpolation kernels reused by `siderust::tables`.

#[cfg(not(feature = "std"))]
use qtty::Real;

use crate::grid::{AxisDirection, GridError, OutOfRange};

// ─── Axis-based internal kernels (Grid1D / Axis) ─────────────────────────────

/// Binary search: returns `(low_index, fraction)` for `x` in ascending `xs`.
#[must_use]
pub(crate) fn locate(xs: &[f64], x: f64) -> (usize, f64) {
    let last = xs.len() - 1;
    if !x.is_finite() {
        return if x.is_sign_positive() {
            (last - 1, 1.0)
        } else {
            (0, 0.0)
        };
    }
    if x <= xs[0] {
        return (0, 0.0);
    }
    if x >= xs[last] {
        return (last - 1, 1.0);
    }
    let upper = xs.partition_point(|v| *v <= x);
    let low = upper - 1;
    let denom = xs[upper] - xs[low];
    let fraction = if denom == 0.0 {
        0.0
    } else {
        (x - xs[low]) / denom
    };
    (low, fraction.clamp(0.0, 1.0))
}

/// Uniform locate: `O(1)`.
#[must_use]
pub(crate) fn locate_uniform(start: f64, step: f64, count: usize, x: f64) -> (usize, f64) {
    if count < 2 {
        return (0, 0.0);
    }
    let last = count - 1;
    let end = start + step * last as f64;
    if !x.is_finite() {
        return if x.is_sign_positive() {
            (last - 1, 1.0)
        } else {
            (0, 0.0)
        };
    }
    if x <= start {
        return (0, 0.0);
    }
    if x >= end {
        return (last - 1, 1.0);
    }
    let pos = (x - start) / step;
    let low = pos.floor() as usize;
    let clamped_low = low.min(last - 1);
    let fraction = (pos - clamped_low as f64).clamp(0.0, 1.0);
    (clamped_low, fraction)
}

/// Linear interpolation.
#[must_use]
pub(crate) fn lerp(y0: f64, y1: f64, t: f64) -> f64 {
    y0 + t * (y1 - y0)
}

/// Bilinear interpolation in axis-first (y-first) storage — kept for internal compatibility.
#[allow(dead_code)]
#[must_use]
pub(crate) fn bilerp(v00: f64, v01: f64, v10: f64, v11: f64, tx: f64, ty: f64) -> f64 {
    lerp(lerp(v00, v01, ty), lerp(v10, v11, ty), tx)
}

/// Trilinear interpolation in axis-first storage — kept for internal compatibility.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
#[must_use]
pub(crate) fn trilerp(
    v000: f64,
    v001: f64,
    v010: f64,
    v011: f64,
    v100: f64,
    v101: f64,
    v110: f64,
    v111: f64,
    tx: f64,
    ty: f64,
    tz: f64,
) -> f64 {
    let c00 = lerp(v000, v001, tz);
    let c01 = lerp(v010, v011, tz);
    let c10 = lerp(v100, v101, tz);
    let c11 = lerp(v110, v111, tz);
    bilerp(c00, c01, c10, c11, tx, ty)
}

// ─── Public direction-aware API ───────────────────────────────────────────────

/// Validates `xs` for use as an interpolation axis.
///
/// Returns the detected [`AxisDirection`] on success. The axis must be strictly
/// monotonic (ascending or descending), have at least two finite samples, and
/// contain no non-finite values.
///
/// # Errors
///
/// Returns [`GridError::TooFewSamples`], [`GridError::NonFinite`], or
/// [`GridError::NotMonotonic`] on invalid input.
///
/// # Examples
///
/// ```rust
/// use optica::grid::{algo, AxisDirection};
///
/// let dir = algo::validate_axis("x", &[1.0, 2.0, 3.0]).unwrap();
/// assert_eq!(dir, AxisDirection::Ascending);
///
/// let dir = algo::validate_axis("y", &[3.0, 2.0, 1.0]).unwrap();
/// assert_eq!(dir, AxisDirection::Descending);
/// ```
pub fn validate_axis(name: &'static str, xs: &[f64]) -> Result<AxisDirection, GridError> {
    if xs.len() < 2 {
        return Err(GridError::TooFewSamples {
            axis: name,
            len: xs.len(),
        });
    }
    for (index, &v) in xs.iter().enumerate() {
        if !v.is_finite() {
            return Err(GridError::NonFinite { axis: name, index });
        }
    }
    let ascending = xs[1] > xs[0];
    for i in 1..xs.len() {
        let monotone = if ascending {
            xs[i] > xs[i - 1]
        } else {
            xs[i] < xs[i - 1]
        };
        if !monotone {
            return Err(GridError::NotMonotonic {
                axis: name,
                at_index: i,
            });
        }
    }
    Ok(if ascending {
        AxisDirection::Ascending
    } else {
        AxisDirection::Descending
    })
}

/// Returns the `(lo, hi)` bounds of `xs` in value order (lo ≤ hi regardless of direction).
#[must_use]
pub(crate) fn axis_range(xs: &[f64], dir: AxisDirection) -> (f64, f64) {
    match dir {
        AxisDirection::Ascending => (xs[0], xs[xs.len() - 1]),
        AxisDirection::Descending => (xs[xs.len() - 1], xs[0]),
    }
}

/// Returns `(low_index, fraction)` for a value in a possibly-descending axis.
///
/// For descending axes the index still counts from position 0 to `len-2`, and
/// `fraction` is in `[0,1]`. The caller selects values at `[low_index]` and
/// `[low_index+1]` from the original (descending) storage.
#[must_use]
pub(crate) fn locate_dir(xs: &[f64], x: f64, dir: AxisDirection) -> (usize, f64) {
    match dir {
        AxisDirection::Ascending => locate(xs, x),
        AxisDirection::Descending => {
            let last = xs.len() - 1;
            if !x.is_finite() {
                return if x.is_sign_positive() {
                    (0, 0.0)
                } else {
                    (last - 1, 1.0)
                };
            }
            // xs is descending: xs[0] is highest, xs[last] is lowest
            if x >= xs[0] {
                return (0, 0.0);
            }
            if x <= xs[last] {
                return (last - 1, 1.0);
            }
            // Binary search for insertion point in descending slice
            // We want the last index where xs[i] > x
            let upper = xs.partition_point(|v| *v > x);
            let low = upper - 1;
            let denom = xs[low] - xs[upper]; // positive since descending
            let fraction = if denom == 0.0 {
                0.0
            } else {
                (xs[low] - x) / denom
            };
            (low, fraction.clamp(0.0, 1.0))
        }
    }
}

/// Applies `oor` for a scalar query `v` against `[lo, hi]`.
///
/// Returns `Ok(true)` when the caller should proceed with interpolation,
/// `Ok(false)` when the result should be zero, or an error.
pub(crate) fn check_oor(
    v: f64,
    lo: f64,
    hi: f64,
    oor: OutOfRange,
    axis: &'static str,
) -> Result<bool, GridError> {
    if v >= lo && v <= hi {
        return Ok(true);
    }
    match oor {
        OutOfRange::ClampToEndpoints => Ok(true),
        OutOfRange::Zero => Ok(false),
        OutOfRange::Error => Err(GridError::OutOfRange {
            axis,
            value: v,
            lo,
            hi,
        }),
    }
}

/// Interpolates a 1-D axis/value pair at `x`.
///
/// # Errors
///
/// Returns [`GridError::OutOfRange`] when `oor` is [`OutOfRange::Error`] and `x`
/// lies outside `xs`.
///
/// # Examples
///
/// ```rust
/// use optica::grid::{algo, AxisDirection, OutOfRange};
///
/// let xs = [0.0_f64, 1.0, 2.0];
/// let ys = [0.0_f64, 10.0, 20.0];
/// let v = algo::linear_1d(&xs, &ys, 0.5, OutOfRange::ClampToEndpoints, AxisDirection::Ascending).unwrap();
/// assert!((v - 5.0).abs() < 1e-12);
/// ```
pub fn linear_1d(
    xs: &[f64],
    ys: &[f64],
    x: f64,
    oor: OutOfRange,
    dir: AxisDirection,
) -> Result<f64, GridError> {
    let (lo, hi) = axis_range(xs, dir);
    if !check_oor(x, lo, hi, oor, "x")? {
        return Ok(0.0);
    }
    let (i, t) = locate_dir(xs, x, dir);
    Ok(lerp(ys[i], ys[i + 1], t))
}

/// Bilinear interpolation with per-row slice storage.
///
/// `table[iy][ix]` stores values in y-major row-major order. `table.len()` must equal
/// `ys.len()` and each `table[iy].len()` must equal `xs.len()`.
///
/// The corner convention is x-first:
/// - `f00 = table[iy0][ix0]`, `f10 = table[iy0][ix1]` (same row, next column)
/// - `f01 = table[iy1][ix0]`, `f11 = table[iy1][ix1]` (next row)
///
/// # Errors
///
/// Returns [`GridError::OutOfRange`] when either axis has `OutOfRange::Error` active.
///
/// # Examples
///
/// ```rust
/// use optica::grid::{algo, AxisDirection, OutOfRange};
///
/// let xs = [400.0_f64, 500.0];
/// let ys = [0.0_f64, 1.0];
/// let row0: &[f64] = &[1.0, 2.0]; // y=0: V(400)=1, V(500)=2
/// let row1: &[f64] = &[3.0, 4.0]; // y=1: V(400)=3, V(500)=4
/// let table: &[&[f64]] = &[row0, row1];
/// let v = algo::bilinear(
///     &xs, &ys, table, 450.0, 0.0,
///     OutOfRange::ClampToEndpoints, OutOfRange::ClampToEndpoints,
///     AxisDirection::Ascending, AxisDirection::Ascending,
/// ).unwrap();
/// assert!((v - 1.5).abs() < 1e-12);
/// ```
#[allow(clippy::too_many_arguments)]
pub fn bilinear(
    xs: &[f64],
    ys: &[f64],
    table: &[&[f64]],
    x: f64,
    y: f64,
    oor_x: OutOfRange,
    oor_y: OutOfRange,
    dir_x: AxisDirection,
    dir_y: AxisDirection,
) -> Result<f64, GridError> {
    let (x_lo, x_hi) = axis_range(xs, dir_x);
    let (y_lo, y_hi) = axis_range(ys, dir_y);
    if !check_oor(x, x_lo, x_hi, oor_x, "x")? {
        return Ok(0.0);
    }
    if !check_oor(y, y_lo, y_hi, oor_y, "y")? {
        return Ok(0.0);
    }
    let (ix, tx) = locate_dir(xs, x, dir_x);
    let (iy, ty) = locate_dir(ys, y, dir_y);
    let f00 = table[iy][ix];
    let f10 = table[iy][ix + 1];
    let f01 = table[iy + 1][ix];
    let f11 = table[iy + 1][ix + 1];
    Ok(bilinear_unit(f00, f10, f01, f11, tx, ty))
}

/// Bilinear kernel: interpolates within a unit cell `[0,1]×[0,1]`.
///
/// Corner convention (x-first):
/// - `f00` = value at `(x=0, y=0)`, `f10` = value at `(x=1, y=0)`
/// - `f01` = value at `(x=0, y=1)`, `f11` = value at `(x=1, y=1)`
///
/// # Examples
///
/// ```rust
/// use optica::grid::algo;
///
/// // Midpoint of a unit square with corners 0, 1, 2, 3
/// let v = algo::bilinear_unit(0.0, 1.0, 2.0, 3.0, 0.5, 0.5);
/// assert!((v - 1.5).abs() < 1e-12);
/// ```
#[must_use]
pub fn bilinear_unit(f00: f64, f10: f64, f01: f64, f11: f64, tx: f64, ty: f64) -> f64 {
    lerp(lerp(f00, f10, tx), lerp(f01, f11, tx), ty)
}

/// Trilinear interpolation with flat `(iz*ny+iy)*nx+ix` storage.
///
/// The `table` slice must have length `nx * ny * nz`. Storage is
/// z-outermost, y-middle, x-innermost: `table[(iz * ny + iy) * nx + ix]`.
///
/// # Errors
///
/// Returns [`GridError::OutOfRange`] when any axis has `OutOfRange::Error` active.
#[allow(clippy::too_many_arguments)]
pub fn trilinear(
    xs: &[f64],
    ys: &[f64],
    zs: &[f64],
    table: &[f64],
    x: f64,
    y: f64,
    z: f64,
    oor_x: OutOfRange,
    oor_y: OutOfRange,
    oor_z: OutOfRange,
    dir_x: AxisDirection,
    dir_y: AxisDirection,
    dir_z: AxisDirection,
) -> Result<f64, GridError> {
    let (x_lo, x_hi) = axis_range(xs, dir_x);
    let (y_lo, y_hi) = axis_range(ys, dir_y);
    let (z_lo, z_hi) = axis_range(zs, dir_z);
    if !check_oor(x, x_lo, x_hi, oor_x, "x")? {
        return Ok(0.0);
    }
    if !check_oor(y, y_lo, y_hi, oor_y, "y")? {
        return Ok(0.0);
    }
    if !check_oor(z, z_lo, z_hi, oor_z, "z")? {
        return Ok(0.0);
    }
    let (ix, tx) = locate_dir(xs, x, dir_x);
    let (iy, ty) = locate_dir(ys, y, dir_y);
    let (iz, tz) = locate_dir(zs, z, dir_z);
    let nx = xs.len();
    let ny = ys.len();
    let idx = |iz: usize, iy: usize, ix: usize| (iz * ny + iy) * nx + ix;
    let v000 = table[idx(iz, iy, ix)];
    let v100 = table[idx(iz, iy, ix + 1)];
    let v010 = table[idx(iz, iy + 1, ix)];
    let v110 = table[idx(iz, iy + 1, ix + 1)];
    let v001 = table[idx(iz + 1, iy, ix)];
    let v101 = table[idx(iz + 1, iy, ix + 1)];
    let v011 = table[idx(iz + 1, iy + 1, ix)];
    let v111 = table[idx(iz + 1, iy + 1, ix + 1)];
    Ok(trilinear_unit(
        v000, v100, v010, v110, v001, v101, v011, v111, tx, ty, tz,
    ))
}

/// Trilinear kernel: interpolates within a unit cube `[0,1]³`.
///
/// Corner convention (x-first, then y, then z):
/// - First index digit = z, second = y, third = x (bit 0 = x, bit 1 = y, bit 2 = z)
/// - `v000` = (x=0,y=0,z=0), `v100` = (x=1,y=0,z=0), `v010` = (x=0,y=1,z=0), …
///
/// # Examples
///
/// ```rust
/// use optica::grid::algo;
///
/// // All corners 0 except v111=8 → midpoint = 1
/// let v = algo::trilinear_unit(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.5, 0.5, 0.5);
/// assert!((v - 1.0).abs() < 1e-12);
/// ```
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn trilinear_unit(
    v000: f64,
    v100: f64,
    v010: f64,
    v110: f64,
    v001: f64,
    v101: f64,
    v011: f64,
    v111: f64,
    tx: f64,
    ty: f64,
    tz: f64,
) -> f64 {
    let c0 = bilinear_unit(v000, v100, v010, v110, tx, ty);
    let c1 = bilinear_unit(v001, v101, v011, v111, tx, ty);
    lerp(c0, c1, tz)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locate_clamps_to_bounds() {
        let xs = [1.0, 2.0, 4.0];
        assert_eq!(locate(&xs, -1.0), (0, 0.0));
        assert_eq!(locate(&xs, 10.0), (1, 1.0));
    }

    #[test]
    fn locate_dir_descending_midpoint() {
        let xs = [10.0, 5.0, 0.0];
        // query 7.5 → between xs[0]=10 and xs[1]=5, fraction = (10-7.5)/(10-5) = 0.5
        let (i, t) = locate_dir(&xs, 7.5, AxisDirection::Descending);
        assert_eq!(i, 0);
        assert!((t - 0.5).abs() < 1e-12);
    }

    #[test]
    fn locate_dir_descending_boundary() {
        let xs = [10.0, 5.0, 0.0];
        assert_eq!(locate_dir(&xs, 15.0, AxisDirection::Descending), (0, 0.0));
        assert_eq!(locate_dir(&xs, -5.0, AxisDirection::Descending), (1, 1.0));
    }

    #[test]
    fn validate_axis_detects_direction() {
        let dir = validate_axis("x", &[1.0, 2.0, 3.0]).unwrap();
        assert_eq!(dir, AxisDirection::Ascending);
        let dir = validate_axis("y", &[3.0, 2.0, 1.0]).unwrap();
        assert_eq!(dir, AxisDirection::Descending);
    }

    #[test]
    fn bilinear_unit_x_first_convention() {
        // f00=0, f10=1, f01=0, f11=1 → pure x interpolation regardless of y
        let v = bilinear_unit(0.0, 1.0, 0.0, 1.0, 0.5, 0.5);
        assert!((v - 0.5).abs() < 1e-12);
        // Discriminative: y-first would give same result here, so use asymmetric corners
        // f00=0, f10=2, f01=4, f11=8, tx=0.5, ty=0 → interpolates x at y=0: lerp(0,2,0.5)=1
        let v = bilinear_unit(0.0, 2.0, 4.0, 8.0, 0.5, 0.0);
        assert!((v - 1.0).abs() < 1e-12);
        // at ty=1: lerp(4,8,0.5)=6
        let v = bilinear_unit(0.0, 2.0, 4.0, 8.0, 0.5, 1.0);
        assert!((v - 6.0).abs() < 1e-12);
    }

    #[test]
    fn trilinear_unit_midpoint() {
        let v = trilerp(0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 0.5, 0.5, 0.5);
        assert_eq!(v, 7.0);
    }

    #[test]
    fn bilinear_row_slice_api() {
        let xs = [400.0, 500.0];
        let ys = [0.0, 1.0];
        let row0: &[f64] = &[1.0, 2.0];
        let row1: &[f64] = &[3.0, 4.0];
        let table: &[&[f64]] = &[row0, row1];
        let v = bilinear(
            &xs,
            &ys,
            table,
            450.0,
            0.0,
            OutOfRange::ClampToEndpoints,
            OutOfRange::ClampToEndpoints,
            AxisDirection::Ascending,
            AxisDirection::Ascending,
        )
        .unwrap();
        assert!((v - 1.5).abs() < 1e-12);
    }
}
