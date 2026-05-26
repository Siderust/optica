// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Smoke tests for the `serde` feature.

#![cfg(feature = "serde")]

use optica::grid::{AxisDirection, OutOfRange};
use optica::scatter::MieParams;
use optica::spectrum::Interpolation;
use optica::transport::{IntegrationMethod, IntegrationOpts};

fn assert_json_roundtrip<T>(value: &T)
where
    T: serde::Serialize + for<'de> serde::Deserialize<'de> + PartialEq + core::fmt::Debug,
{
    let json = serde_json::to_string(value).expect("serialize");
    let back: T = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(value, &back);
}

#[test]
fn axis_direction_roundtrips() {
    assert_json_roundtrip(&AxisDirection::Ascending);
    assert_json_roundtrip(&AxisDirection::Descending);
}

#[test]
fn out_of_range_roundtrips() {
    assert_json_roundtrip(&OutOfRange::ClampToEndpoints);
    assert_json_roundtrip(&OutOfRange::Zero);
    assert_json_roundtrip(&OutOfRange::Error);
}

#[test]
fn interpolation_roundtrips() {
    assert_json_roundtrip(&Interpolation::Linear);
    assert_json_roundtrip(&Interpolation::Nearest);
    assert_json_roundtrip(&Interpolation::CubicSpline);
}

#[test]
fn integration_opts_roundtrips() {
    assert_json_roundtrip(&IntegrationOpts::new(32, IntegrationMethod::Simpson));
}

#[test]
fn mie_params_roundtrips() {
    let p = MieParams::try_new(0.12, 1.3, 550.0).unwrap();
    assert_json_roundtrip(&p);
}
