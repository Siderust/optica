// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # optica — participating-media and optics foundations
//!
//! `optica` provides fast, typed building blocks for optics and radiative-transfer
//! workloads: rays, optical coefficients, scattering phase functions, sampled
//! spectra, interpolation tables, and optical-depth integration kernels.
//!
//! The crate is intentionally domain-agnostic. It does **not** include astronomy,
//! planetary constants, ephemerides, or site-specific policies.
//!
//! ## References
//!
//! - Chandrasekhar, *Radiative Transfer*.
//! - Bodhaine et al. (1999), Rayleigh optical-depth approximation.

pub mod data;
pub mod grid;
pub mod medium;
pub mod phase;
pub mod prelude;
pub mod ray;
pub mod scatter;
pub mod spectrum;
pub mod transport;
