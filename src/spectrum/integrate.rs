// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Re-exports of untyped trapezoidal integration kernels.
//!
//! For typed integration over a [`super::SampledSpectrum`], use the methods
//! on that type directly.

pub use super::algo::{trapz, trapz_range, trapz_weighted};
