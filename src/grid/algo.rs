// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Low-level `f64` interpolation kernels shared by typed grid wrappers.

/// Binary search: returns `(low_index, fraction)` for `x` in `xs`.
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

    let upper = xs.partition_point(|value| *value <= x);
    let low = upper - 1;
    let denom = xs[upper] - xs[low];
    let fraction = if denom == 0.0 { 0.0 } else { (x - xs[low]) / denom };
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

/// Bilinear interpolation in row-major storage.
#[must_use]
pub(crate) fn bilerp(v00: f64, v01: f64, v10: f64, v11: f64, tx: f64, ty: f64) -> f64 {
    lerp(lerp(v00, v01, ty), lerp(v10, v11, ty), tx)
}

/// Trilinear interpolation in row-major storage.
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
    fn trilinear_interpolation_matches_corner_midpoint() {
        let value = trilerp(0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 0.5, 0.5, 0.5);
        assert_eq!(value, 7.0);
    }
}
