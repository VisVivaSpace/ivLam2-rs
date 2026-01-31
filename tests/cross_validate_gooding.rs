//! Cross-validation tests: Rust vercosine solver vs C Gooding solver.
//!
//! Both solvers are iterative with different formulations, so we use 1e-8
//! velocity agreement tolerance (documented: two iterative methods, different
//! formulations, both converged to ~1e-10 internally).
//!
//! The Gooding solver only handles prograde transfers (theta in [0, pi]),
//! so all tests here are prograde zero-rev unless noted.

#![cfg(feature = "gooding-ffi")]

mod common;

use common::gooding::gooding_lambert;
use lambert_solver::{solve_lambert, Direction};
use std::f64::consts::PI;

/// Tolerance for cross-validation: both iterative, different formulations
const CROSS_TOL: f64 = 1e-8;

fn vec_mag(v: &[f64; 3]) -> f64 {
    (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt()
}

fn compare_velocities(
    label: &str,
    rust_v1: &[f64; 3],
    rust_v2: &[f64; 3],
    good_v1: &[f64; 3],
    good_v2: &[f64; 3],
    tol: f64,
) {
    for i in 0..3 {
        assert!(
            (rust_v1[i] - good_v1[i]).abs() < tol,
            "{}: v1[{}] mismatch: rust={:.12e}, gooding={:.12e}, diff={:.2e}",
            label, i, rust_v1[i], good_v1[i], (rust_v1[i] - good_v1[i]).abs()
        );
        assert!(
            (rust_v2[i] - good_v2[i]).abs() < tol,
            "{}: v2[{}] mismatch: rust={:.12e}, gooding={:.12e}, diff={:.2e}",
            label, i, rust_v2[i], good_v2[i], (rust_v2[i] - good_v2[i]).abs()
        );
    }
}

fn cross_validate(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64, tol: f64) {
    let rust_sol = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0)
        .unwrap_or_else(|e| panic!("{}: Rust solver failed: {:?}", label, e));

    let good_sol = gooding_lambert(mu, &r1, &r2, 0, tof)
        .unwrap_or_else(|| panic!("{}: Gooding solver returned no solution", label));

    compare_velocities(label, &rust_sol.v1, &rust_sol.v2, &good_sol.v1, &good_sol.v2, tol);
}

// --- Transfer angle sweep ---

#[test]
fn test_cross_validate_10_deg() {
    let angle = 10.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    cross_validate("10°", r1, r2, angle, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_45_deg() {
    let angle = 45.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    cross_validate("45°", r1, r2, angle, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_90_deg() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    cross_validate("90°", r1, r2, PI / 2.0, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_120_deg() {
    let angle = 120.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    cross_validate("120°", r1, r2, angle, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_135_deg() {
    let angle = 135.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    cross_validate("135°", r1, r2, angle, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_170_deg() {
    let angle = 170.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    cross_validate("170°", r1, r2, 2.0, 1.0, CROSS_TOL);
}

// --- Different radii ---

#[test]
fn test_cross_validate_hohmann_like() {
    // r1=1, r2=2, 90° prograde
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 2.0, 0.0];
    cross_validate("Hohmann-like", r1, r2, 2.0, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_inner_transfer() {
    // r1=1, r2=0.7
    let angle = 90.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.7, 0.0];
    cross_validate("inner transfer", r1, r2, 1.5, 1.0, CROSS_TOL);
}

// --- Physical units ---

#[test]
fn test_cross_validate_leo_to_geo() {
    let mu = 398600.4418;
    let r1 = [6678.0, 0.0, 0.0];
    let r2 = [0.0, 42164.0, 0.0];
    let tof = 5.0 * 3600.0;
    cross_validate("LEO→GEO", r1, r2, tof, mu, CROSS_TOL);
}

// --- 3D transfers ---

#[test]
fn test_cross_validate_3d_transfer() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.8, 0.6]; // |r2|=1, out of plane
    cross_validate("3D", r1, r2, PI / 2.0, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_3d_different_radii() {
    let r1 = [1.0, 0.5, 0.0];
    let r2 = [0.0, 1.0, 0.8];
    cross_validate("3D diff-r", r1, r2, 2.0, 1.0, CROSS_TOL);
}

// --- Hyperbolic ---

#[test]
fn test_cross_validate_hyperbolic() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    cross_validate("hyperbolic", r1, r2, 0.1, 1.0, CROSS_TOL);
}
