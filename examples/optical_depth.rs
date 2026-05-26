// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

#![allow(clippy::print_stdout)]
//! Integrate optical depth along a ray through a homogeneous medium.
//!
//! This example shows the midpoint-rule integration helper and the Beer–Lambert
//! transmittance computation. A homogeneous medium with known σ_a and σ_s is used
//! so that the result can be verified analytically: τ = σ_t × path_length.

use affn::{CartesianDirection, Position, ReferenceCenter, ReferenceFrame};
use optica::medium::HomogeneousMedium;
use optica::ray::{Ray, RaySegment};
use optica::transport::{integrate_optical_depth, transmittance, IntegrationOpts};
use qtty::length::{Kilometers, Nanometers};
use qtty::unit::Kilometer;

#[derive(Debug, Copy, Clone)]
struct EarthCenter;

impl ReferenceCenter for EarthCenter {
    type Params = ();
    fn center_name() -> &'static str {
        "EarthCenter"
    }
}

#[derive(Debug, Copy, Clone)]
struct Enu;

impl ReferenceFrame for Enu {
    fn frame_name() -> &'static str {
        "ENU"
    }
}

fn main() {
    // Medium: σ_a = 0.05 km⁻¹, σ_s = 0.10 km⁻¹ → σ_t = 0.15 km⁻¹
    let medium = HomogeneousMedium::<Kilometer>::new(0.05, 0.10);

    // Vertical ray from ground to 20 km altitude
    let ray = Ray::new(
        Position::<EarthCenter, Enu, Kilometer>::new(0.0, 0.0, 0.0),
        CartesianDirection::<Enu>::new(0.0, 0.0, 1.0),
    );
    let segment = RaySegment::new(Kilometers::new(0.0), Kilometers::new(20.0));

    // Analytical expectation: τ = 0.15 km⁻¹ × 20 km = 3.0
    let tau = integrate_optical_depth(
        &medium,
        &ray,
        segment,
        Nanometers::new(550.0),
        IntegrationOpts { n_steps: 100 },
    );
    let t = transmittance(tau);

    println!("Optical depth τ = {:.6}", tau.value());
    println!("Transmittance  T = exp(-τ) = {:.6}", t.value());
    println!("Expected τ = 3.000000, T = {:.6}", (-3.0_f64).exp());

    assert!(
        (tau.value() - 3.0).abs() < 1e-10,
        "unexpected τ: {}",
        tau.value()
    );
    println!("✓ passes analytical check");
}
