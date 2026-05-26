// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

#![allow(clippy::print_stdout)]
//! Bilinear interpolation on a 2-D typed grid.
//!
//! This example builds a small wavelength × zenith-angle table and queries it
//! with `Grid2D::interp_at`. The same approach is used for phase tables, zodiacal
//! background tables, and any other two-parameter optical lookup.

use optica::grid::{Grid2D, OutOfRange};
use qtty::unit::{Nanometer, Radian, Ratio};
use qtty::Quantity;

fn main() {
    // Tiny transmittance table: T(wavelength, zenith_angle)
    // Wavelengths: 400, 550, 700 nm
    // Zenith angles: 0, 45, 90 degrees (stored as radians)
    let lam_nm = [400.0_f64, 550.0, 700.0];
    let za_rad = [
        0.0_f64,
        std::f64::consts::FRAC_PI_4,
        std::f64::consts::FRAC_PI_2,
    ];

    // Row-major values: row i (za=za_rad[i]) over all wavelengths
    //   row 0 (z=0°):  T = 0.90, 0.95, 0.97
    //   row 1 (z=45°): T = 0.78, 0.85, 0.90
    //   row 2 (z=90°): T = 0.00, 0.00, 0.00  (horizon, full extinction)
    #[rustfmt::skip]
    let values = [
        0.90_f64, 0.95, 0.97,
        0.78,     0.85, 0.90,
        0.00,     0.00, 0.00,
    ];

    let grid = Grid2D::<Nanometer, Radian, Ratio>::from_raw_row_major(
        &lam_nm,
        &za_rad,
        &values,
        OutOfRange::ClampToEndpoints,
    )
    .expect("valid grid");

    println!(
        "Grid: {}×{} (wavelength × zenith angle)",
        grid.nx(),
        grid.ny()
    );

    let queries = [
        (475.0_f64, 0.0_f64),
        (550.0, std::f64::consts::FRAC_PI_4),
        (625.0, std::f64::consts::FRAC_PI_4 / 2.0),
    ];
    for (lam, za) in queries {
        let t = grid.interp_at(Quantity::<Nanometer>::new(lam), Quantity::<Radian>::new(za));
        println!("  T({lam:.0} nm, {:.2} rad) = {:.4}", za, t.value());
    }
}
