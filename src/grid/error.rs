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
/// let error = GridError::TooFewSamples { axis: "x", len: 1 };
/// assert!(error.to_string().contains("too few"));
/// ```
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum GridError {
    /// The value buffer length does not match the product of the axis lengths.
    #[error("expected {expected} values, got {got}")]
    ShapeMismatch {
        /// Expected count derived from the axis dimensions.
        expected: usize,
        /// Actual count supplied.
        got: usize,
    },
    /// An axis has fewer than two samples (interpolation requires at least two).
    #[error("{axis} axis has too few samples ({len}); need ≥ 2")]
    TooFewSamples {
        /// Name of the undersampled axis.
        axis: &'static str,
        /// Number of samples present.
        len: usize,
    },
    /// An axis is not strictly monotonic (either ascending or descending).
    #[error("{axis} axis is not strictly monotonic at index {at_index}")]
    NotMonotonic {
        /// Axis name.
        axis: &'static str,
        /// Index where monotonicity fails.
        at_index: usize,
    },
    /// An axis contains a non-finite value.
    #[error("{axis} axis contains non-finite value at index {index}")]
    NonFinite {
        /// Axis name.
        axis: &'static str,
        /// Index of the offending value.
        index: usize,
    },
    /// A uniform axis used a non-positive step.
    #[error("step must be positive (got {step})")]
    NonPositiveStep {
        /// Invalid step size.
        step: f64,
    },
    /// The y-axis supplied to a y-descending constructor is not uniformly spaced.
    #[error("{axis} axis is not uniformly spaced: expected step {expected}, got {got}")]
    NonUniformStep {
        /// Axis name.
        axis: &'static str,
        /// Expected step (from the first pair).
        expected: f64,
        /// Observed deviating step.
        got: f64,
    },
    /// A query point lies outside an axis range while `OutOfRange::Error` is active.
    #[error("{axis} value {value} is outside axis range [{lo}, {hi}]")]
    OutOfRange {
        /// Axis name.
        axis: &'static str,
        /// Queried raw coordinate value.
        value: f64,
        /// Lower bound of the valid axis range.
        lo: f64,
        /// Upper bound of the valid axis range.
        hi: f64,
    },
}
