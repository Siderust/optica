// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed interpolation grids for optics workloads.
//!
//! The grids in this module expose typed query and value units while keeping the
//! interpolation kernels on contiguous `f64` storage for fast hot-loop use.
//!
//! ## References
//!
//! - Standard linear, bilinear, and trilinear interpolation formulas.

pub mod algo;
pub mod axis;
pub mod error;
pub mod grid1d;
pub mod grid2d;
pub mod grid3d;

pub use axis::Axis;
pub use error::GridError;
pub use grid1d::Grid1D;
pub use grid2d::{ConstantRegion, Grid2D};
pub use grid3d::Grid3D;

/// Axis monotonicity direction for direction-aware interpolation kernels.
///
/// Used by [`algo`] functions and [`Grid2D`]/[`Grid3D`] constructors that accept
/// both ascending and descending axis data.
///
/// # Examples
///
/// ```rust
/// use optica::grid::AxisDirection;
///
/// let dir = AxisDirection::Ascending;
/// assert_ne!(dir, AxisDirection::Descending);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxisDirection {
    /// Values increase with index.
    Ascending,
    /// Values decrease with index.
    Descending,
}

/// Policy for interpolation values requested outside the grid's range.
///
/// # Examples
///
/// ```rust
/// use optica::grid::OutOfRange;
///
/// assert_eq!(OutOfRange::default(), OutOfRange::ClampToEndpoints);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutOfRange {
    /// Clamp to the nearest endpoint value (default, fastest).
    #[default]
    ClampToEndpoints,
    /// Return zero outside range.
    Zero,
    /// Return an error.
    Error,
}
