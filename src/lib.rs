// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

#![doc(html_root_url = "https://docs.rs/optica/0.3.0")]
#![cfg_attr(not(feature = "std"), no_std)]
//! # optica — participating-media and optics foundations
//!
//! `optica` provides fast, typed building blocks for optics and radiative-transfer
//! workloads: rays, optical coefficients, scattering phase functions, sampled
//! spectra, interpolation tables, and optical-depth integration kernels.
//!
//! The crate is intentionally domain-agnostic. It does **not** include astronomy,
//! planetary constants, ephemerides, or site-specific policies.
//!
//! ## Feature flags
//!
//! - `std` *(default)*: enables the standard library and implies `alloc`;
//!   forwards `qtty/std` and `affn/std`.
//! - `alloc`: enables heap-backed containers in `no_std` builds. All public types
//!   that use `Vec`/`Box`/`String` are gated behind `alloc` or `std`; forwards
//!   `qtty/alloc` and `affn/alloc`.
//! - `serde`: derives `Serialize`/`Deserialize` for the public data, error,
//!   and policy enums (implies `alloc`).
//!
//! Pure `core`-only builds use `--no-default-features`. Combine with
//! `--features alloc` and/or `--features serde` as needed.
//!
//! ## References
//!
//! - Chandrasekhar, *Radiative Transfer*.
//! - Bodhaine et al. (1999), Rayleigh optical-depth approximation.

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(feature = "alloc")]
pub mod data;
#[cfg(feature = "alloc")]
pub mod grid;
pub mod medium;
#[cfg(feature = "alloc")]
pub mod phase;
#[cfg(feature = "alloc")]
pub mod prelude;
pub mod ray;
pub mod scatter;
#[cfg(feature = "alloc")]
pub mod spectrum;
pub mod transport;

#[cfg(feature = "alloc")]
pub use grid::{AxisDirection, ConstantRegion};
