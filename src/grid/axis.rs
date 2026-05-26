// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Lookup-axis definitions for interpolation tables.

use alloc::boxed::Box;

use crate::grid::algo::{locate, locate_uniform};
use crate::grid::error::GridError;

/// Lookup axis for a 1-D table, supporting uniform (`O(1)`) and non-uniform (`O(log n)`) grids.
///
/// # Examples
///
/// ```rust
/// use optica::grid::Axis;
///
/// let axis = Axis::uniform(400.0, 50.0, 3).unwrap();
/// assert_eq!(axis.locate(425.0), (0, 0.5));
/// ```
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Axis {
    /// Uniformly-spaced axis: `O(1)` lookup.
    Uniform {
        /// First sample coordinate.
        start: f64,
        /// Positive spacing between consecutive samples.
        step: f64,
        /// Number of samples on the axis.
        count: usize,
    },
    /// Non-uniform axis with binary search.
    NonUniform(Box<[f64]>),
}

impl Axis {
    /// Returns the number of samples on the axis.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::Axis;
    ///
    /// let axis = Axis::non_uniform([1.0, 2.0, 4.0]).unwrap();
    /// assert_eq!(axis.len(), 3);
    /// ```
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Uniform { count, .. } => *count,
            Self::NonUniform(xs) => xs.len(),
        }
    }

    /// Returns `true` when the axis has no samples.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the `(low_index, fraction)` pair for interpolation.
    ///
    /// `low_index` is clamped to `[0, len - 2]`. `fraction` is clamped to `[0, 1]`.
    #[must_use]
    pub fn locate(&self, x: f64) -> (usize, f64) {
        match self {
            Self::Uniform { start, step, count } => locate_uniform(*start, *step, *count, x),
            Self::NonUniform(xs) => locate(xs, x),
        }
    }

    /// Returns the axis endpoints as `(min, max)` (always with `min ≤ max`).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::Axis;
    ///
    /// let axis = Axis::uniform(400.0, 50.0, 3).unwrap();
    /// let (lo, hi) = axis.bounds();
    /// assert_eq!(lo, 400.0);
    /// assert_eq!(hi, 500.0);
    /// ```
    #[must_use]
    pub fn bounds(&self) -> (f64, f64) {
        match self {
            Self::Uniform { start, step, count } => {
                let end = start + step * (*count as f64 - 1.0);
                if *step >= 0.0 {
                    (*start, end)
                } else {
                    (end, *start)
                }
            }
            Self::NonUniform(xs) => {
                let lo = xs.iter().copied().fold(f64::INFINITY, f64::min);
                let hi = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                (lo, hi)
            }
        }
    }

    /// Validates the axis shape and monotonicity.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::Axis;
    ///
    /// let axis = Axis::non_uniform([400.0, 500.0, 700.0]).unwrap();
    /// assert!(axis.validate().is_ok());
    /// ```
    pub fn validate(&self) -> Result<(), GridError> {
        self.validate_for_axis("axis")
    }

    /// Constructs a validated uniform axis.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::Axis;
    ///
    /// let axis = Axis::uniform(0.0, 0.5, 4).unwrap();
    /// assert_eq!(axis.len(), 4);
    /// ```
    pub fn uniform(start: f64, step: f64, count: usize) -> Result<Self, GridError> {
        let axis = Self::Uniform { start, step, count };
        axis.validate()?;
        Ok(axis)
    }

    /// Constructs a validated non-uniform axis from sample coordinates.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use optica::grid::Axis;
    ///
    /// let axis = Axis::non_uniform([1.0, 3.0, 5.0]).unwrap();
    /// assert_eq!(axis.len(), 3);
    /// ```
    pub fn non_uniform(xs: impl Into<Box<[f64]>>) -> Result<Self, GridError> {
        let axis = Self::NonUniform(xs.into());
        axis.validate()?;
        Ok(axis)
    }

    pub(crate) fn validate_for_axis(&self, name: &'static str) -> Result<(), GridError> {
        match self {
            Self::Uniform { start, step, count } => {
                if *count < 2 {
                    return Err(GridError::TooFewSamples {
                        axis: name,
                        len: *count,
                    });
                }
                if !start.is_finite() {
                    return Err(GridError::NonFinite {
                        axis: name,
                        index: 0,
                    });
                }
                if !step.is_finite() {
                    return Err(GridError::NonFinite {
                        axis: name,
                        index: 1,
                    });
                }
                if *step <= 0.0 {
                    return Err(GridError::NonPositiveStep { step: *step });
                }
                Ok(())
            }
            Self::NonUniform(xs) => {
                if xs.len() < 2 {
                    return Err(GridError::TooFewSamples {
                        axis: name,
                        len: xs.len(),
                    });
                }
                for (index, value) in xs.iter().copied().enumerate() {
                    if !value.is_finite() {
                        return Err(GridError::NonFinite { axis: name, index });
                    }
                    if index > 0 && value <= xs[index - 1] {
                        return Err(GridError::NotMonotonic {
                            axis: name,
                            at_index: index,
                        });
                    }
                }
                Ok(())
            }
        }
    }

    #[must_use]
    pub(crate) fn contains(&self, x: f64) -> bool {
        let (min, max) = self.bounds();
        x.is_finite() && x >= min && x <= max
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_axis_locates_midpoint() {
        let axis = Axis::uniform(0.0, 2.0, 3).unwrap();
        assert_eq!(axis.locate(1.0), (0, 0.5));
    }

    #[test]
    fn non_uniform_axis_rejects_unsorted_values() {
        let error = Axis::non_uniform([1.0, 1.0, 2.0]).unwrap_err();
        assert!(matches!(error, GridError::NotMonotonic { .. }));
    }
}
