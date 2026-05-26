// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

#![allow(clippy::print_stdout)]
//! Load a two-column ASCII spectrum and query it.
//!
//! This example demonstrates parsing a whitespace-separated table of
//! wavelength/transmittance pairs and interpolating at an intermediate point.
//! The same pattern applies to Bessell filter curves, SVO passband files,
//! and atmospheric transmission tables stored in two-column ASCII format.

use optica::grid::OutOfRange;
use optica::spectrum::{loaders::ascii, Interpolation};
use qtty::unit::{Nanometer, Ratio};
use qtty::Quantity;

fn main() {
    let raw = "\
# Simplified ozone transmittance — illustrative values only
# wavelength_nm  transmittance
 290.0  0.00
 300.0  0.02
 320.0  0.30
 340.0  0.72
 360.0  0.95
 400.0  0.99
 700.0  1.00
";

    let spectrum = ascii::two_column::<Nanometer, Ratio>(
        raw,
        1.0, // x already in nanometres
        1.0, // y already dimensionless ratio
        Interpolation::Linear,
        OutOfRange::ClampToEndpoints,
        None,
    )
    .expect("valid ASCII data");

    println!("Loaded {} samples", spectrum.len());
    println!(
        "  x range: [{:.1}, {:.1}] nm",
        spectrum.xs_raw()[0],
        spectrum.xs_raw()[spectrum.len() - 1]
    );

    for lam_nm in [290.0_f64, 310.0, 330.0, 400.0, 550.0, 700.0] {
        let t = spectrum.interp_at(Quantity::<Nanometer>::new(lam_nm));
        println!("  T({lam_nm:.0} nm) = {:.4}", t.value());
    }
}
