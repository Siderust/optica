// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Error types for sampled-spectrum operations.

/// Errors produced while constructing or querying a sampled spectrum.
///
/// # Examples
///
/// ```rust
/// use optica::spectrum::SpectrumError;
///
/// let error = SpectrumError::TooFewSamples(1);
/// assert!(error.to_string().contains("too few"));
/// ```
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SpectrumError {
    /// The `xs` and `ys` slices have different lengths.
    #[error("length mismatch: xs has {xs} elements, ys has {ys}")]
    LengthMismatch {
        /// Number of x samples supplied.
        xs: usize,
        /// Number of y samples supplied.
        ys: usize,
    },
    /// The spectrum has fewer than two samples (interpolation requires ≥ 2).
    #[error("spectrum has too few samples ({0}); need ≥ 2")]
    TooFewSamples(usize),
    /// The x-axis is not strictly monotonically increasing at the given index.
    #[error("x-axis is not strictly increasing at index {index}")]
    NotMonotonic {
        /// Index where monotonicity fails.
        index: usize,
    },
    /// Query point is outside the spectrum's x-axis range when `OutOfRange::Error` is active.
    #[error("query {x} is outside [{lo}, {hi}]")]
    OutOfRange {
        /// Query x value.
        x: f64,
        /// Axis lower bound.
        lo: f64,
        /// Axis upper bound.
        hi: f64,
    },
    /// Failed to parse a spectrum from text data.
    #[error("parse error: {0}")]
    Parse(String),
}
