// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Error types for typed interpolation grids.

/// Errors produced while constructing or querying interpolation grids.
///
/// # Examples
///
/// ```rust
/// use optica::grid::GridError;
///
/// let error = GridError::TooFewPoints { axis: 0, got: 1 };
/// assert!(error.to_string().contains("at least 2 points"));
/// ```
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum GridError {
    /// An axis does not contain enough samples for interpolation.
    #[error("axis {axis} must have at least 2 points (got {got})")]
    TooFewPoints {
        /// Zero-based axis index.
        axis: usize,
        /// Number of points supplied for the axis.
        got: usize,
    },
    /// An axis is not strictly increasing.
    #[error("axis {axis} is not strictly monotone increasing at index {index}")]
    NotMonotone {
        /// Zero-based axis index.
        axis: usize,
        /// Index where monotonicity fails.
        index: usize,
    },
    /// An axis contains a non-finite value.
    #[error("axis {axis} contains non-finite value at index {index}")]
    NonFinite {
        /// Zero-based axis index.
        axis: usize,
        /// Index of the offending value.
        index: usize,
    },
    /// The flattened value buffer has the wrong length.
    #[error("expected {expected} values, got {got}")]
    ShapeMismatch {
        /// Expected number of values from the axis shape.
        expected: usize,
        /// Actual number of values supplied.
        got: usize,
    },
    /// A uniform axis used a non-positive step.
    #[error("step must be positive (got {step})")]
    NonPositiveStep {
        /// Invalid step size.
        step: f64,
    },
    /// A query point lies outside an axis range while `OutOfRange::Error` is active.
    #[error("query value {value} is outside axis {axis} range [{min}, {max}]")]
    OutOfRange {
        /// Zero-based axis index.
        axis: usize,
        /// Queried raw coordinate value.
        value: f64,
        /// Minimum valid value on the axis.
        min: f64,
        /// Maximum valid value on the axis.
        max: f64,
    },
}
