// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed rays and ray segments.

use affn::{CartesianDirection, Position, ReferenceCenter, ReferenceFrame};
use qtty::length::LengthUnit;
use qtty::Quantity;

/// A ray in 3-D space: an origin point and a normalized direction.
///
/// # Examples
///
/// ```rust
/// use affn::{CartesianDirection, Position, ReferenceCenter, ReferenceFrame};
/// use optica::ray::Ray;
/// use qtty::unit::Kilometer;
///
/// #[derive(Debug, Copy, Clone)]
/// struct Center;
/// impl ReferenceCenter for Center {
///     type Params = ();
///     fn center_name() -> &'static str { "Center" }
/// }
///
/// #[derive(Debug, Copy, Clone)]
/// struct Frame;
/// impl ReferenceFrame for Frame {
///     fn frame_name() -> &'static str { "Frame" }
/// }
///
/// let ray = Ray::new(
///     Position::<Center, Frame, Kilometer>::new(0.0, 0.0, 0.0),
///     CartesianDirection::<Frame>::new(0.0, 0.0, 1.0),
/// );
/// assert_eq!(ray.origin.as_array()[2].value(), 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct Ray<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> {
    /// Ray origin.
    pub origin: Position<C, F, U>,
    /// Ray propagation direction.
    pub direction: CartesianDirection<F>,
}

impl<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> Ray<C, F, U> {
    /// Creates a new ray.
    #[must_use]
    pub fn new(origin: Position<C, F, U>, direction: CartesianDirection<F>) -> Self {
        Self { origin, direction }
    }
}

/// A parametric ray interval `[t_min, t_max]` in the ray's length unit.
///
/// The constructor stores endpoints in ascending order.
///
/// # Examples
///
/// ```rust
/// use optica::ray::RaySegment;
/// use qtty::length::Kilometers;
/// use qtty::unit::Kilometer;
///
/// let segment = RaySegment::<Kilometer>::new(Kilometers::new(3.0), Kilometers::new(1.0));
/// assert_eq!(segment.length().value(), 2.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RaySegment<U: LengthUnit> {
    /// Lower endpoint of the segment.
    pub t_min: Quantity<U>,
    /// Upper endpoint of the segment.
    pub t_max: Quantity<U>,
}

impl<U: LengthUnit> RaySegment<U> {
    /// Creates a new ray segment.
    #[must_use]
    pub fn new(t_min: Quantity<U>, t_max: Quantity<U>) -> Self {
        let (t_min, t_max) = if t_min <= t_max {
            (t_min, t_max)
        } else {
            (t_max, t_min)
        };
        Self { t_min, t_max }
    }

    /// Returns the segment length.
    #[must_use]
    pub fn length(&self) -> Quantity<U> {
        self.t_max - self.t_min
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::length::Kilometers;
    use qtty::unit::Kilometer;

    #[test]
    fn segment_orders_endpoints() {
        let segment = RaySegment::<Kilometer>::new(Kilometers::new(3.0), Kilometers::new(1.0));
        assert_eq!(segment.t_min.value(), 1.0);
        assert_eq!(segment.t_max.value(), 3.0);
    }
}
