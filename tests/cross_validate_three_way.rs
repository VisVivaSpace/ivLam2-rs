//! Three-way cross-validation: Rust vercosine vs C Gooding vs Fortran ivLam.
//!
//! All three solvers should agree to within 1e-10 velocity components.
//! ivLam assumes mu=1, so all tests use mu=1.

#![cfg(all(feature = "gooding-ffi", feature = "ivlam-ffi"))]

mod common;

use common::gooding::gooding_lambert;
use common::ivlam_ffi::ivlam_zero_rev;
use lambert_solver::{solve_lambert, Direction};
use std::f64::consts::PI;

/// Tolerance: three iterative methods, same problem, 1e-10 agreement
const THREE_WAY_TOL: f64 = 1e-10;
/// Looser tolerance for Gooding which uses a different formulation
const GOODING_TOL: f64 = 1e-8;

fn three_way_check(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64) {
    let mu = 1.0;

    // Rust solver
    let rust_sol = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0)
        .unwrap_or_else(|e| panic!("{}: Rust solver failed: {:?}", label, e));

    // ivLam (mu=1, direction=1 for prograde)
    let ivlam_sol = ivlam_zero_rev(&r1, &r2, tof, 1);
    assert_eq!(ivlam_sol.info, 0, "{}: ivLam returned info={}", label, ivlam_sol.info);

    // Gooding
    let good_sol = gooding_lambert(mu, &r1, &r2, 0, tof)
        .unwrap_or_else(|| panic!("{}: Gooding solver failed", label));

    // Compare Rust vs ivLam (both vercosine-based, should be very close)
    for i in 0..3 {
        assert!(
            (rust_sol.v1[i] - ivlam_sol.v1[i]).abs() < THREE_WAY_TOL,
            "{}: Rust vs ivLam v1[{}]: rust={:.12e}, ivlam={:.12e}, diff={:.2e}",
            label, i, rust_sol.v1[i], ivlam_sol.v1[i],
            (rust_sol.v1[i] - ivlam_sol.v1[i]).abs()
        );
        assert!(
            (rust_sol.v2[i] - ivlam_sol.v2[i]).abs() < THREE_WAY_TOL,
            "{}: Rust vs ivLam v2[{}]: rust={:.12e}, ivlam={:.12e}, diff={:.2e}",
            label, i, rust_sol.v2[i], ivlam_sol.v2[i],
            (rust_sol.v2[i] - ivlam_sol.v2[i]).abs()
        );
    }

    // Compare Rust vs Gooding (different formulations, looser tolerance)
    for i in 0..3 {
        assert!(
            (rust_sol.v1[i] - good_sol.v1[i]).abs() < GOODING_TOL,
            "{}: Rust vs Gooding v1[{}]: rust={:.12e}, good={:.12e}, diff={:.2e}",
            label, i, rust_sol.v1[i], good_sol.v1[i],
            (rust_sol.v1[i] - good_sol.v1[i]).abs()
        );
        assert!(
            (rust_sol.v2[i] - good_sol.v2[i]).abs() < GOODING_TOL,
            "{}: Rust vs Gooding v2[{}]: rust={:.12e}, good={:.12e}, diff={:.2e}",
            label, i, rust_sol.v2[i], good_sol.v2[i],
            (rust_sol.v2[i] - good_sol.v2[i]).abs()
        );
    }
}

#[test]
fn test_three_way_90_deg() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    three_way_check("90°", r1, r2, PI / 2.0);
}

#[test]
fn test_three_way_45_deg() {
    let angle = 45.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    three_way_check("45°", r1, r2, angle);
}

#[test]
fn test_three_way_135_deg() {
    let angle = 135.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    three_way_check("135°", r1, r2, angle);
}

#[test]
fn test_three_way_170_deg() {
    let angle = 170.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    three_way_check("170°", r1, r2, 2.0);
}

#[test]
fn test_three_way_different_radii() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 2.0, 0.0];
    three_way_check("diff-r", r1, r2, 2.0);
}

#[test]
fn test_three_way_3d() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.8, 0.6];
    three_way_check("3D", r1, r2, PI / 2.0);
}

#[test]
fn test_three_way_hyperbolic() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    three_way_check("hyperbolic", r1, r2, 0.1);
}

// ────────────────────────────────────────────────────────────────
// Angle sweep (equal radii, mu=1)
// ────────────────────────────────────────────────────────────────

#[test]
fn test_three_way_angle_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    for angle_deg in [5.0_f64, 15.0, 25.0, 30.0, 60.0, 75.0, 100.0, 110.0, 140.0, 150.0, 160.0] {
        let angle = angle_deg.to_radians();
        let r2 = [angle.cos(), angle.sin(), 0.0];
        three_way_check(&format!("angle_{}°", angle_deg), r1, r2, angle);
    }
}

// ────────────────────────────────────────────────────────────────
// Radius ratio sweep (90° transfer, mu=1)
// ────────────────────────────────────────────────────────────────

#[test]
fn test_three_way_radius_ratio_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    for ratio in [0.3_f64, 0.5, 0.7, 1.5, 2.0, 3.0, 5.0] {
        let r2 = [0.0, ratio, 0.0];
        let a = (1.0 + ratio) / 2.0;
        let tof = PI * (a * a * a).sqrt() * 0.5;
        three_way_check(&format!("ratio_{:.1}", ratio), r1, r2, tof);
    }
}

// ────────────────────────────────────────────────────────────────
// TOF variation (fixed 90° geometry, mu=1)
// ────────────────────────────────────────────────────────────────

#[test]
fn test_three_way_tof_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    for tof in [0.2_f64, 0.5, 1.0, 2.0, 5.0, 10.0] {
        three_way_check(&format!("tof_{:.1}", tof), r1, r2, tof);
    }
}

// ────────────────────────────────────────────────────────────────
// 3D transfers
// ────────────────────────────────────────────────────────────────

#[test]
fn test_three_way_3d_inclined() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.866, 0.5]; // 30° inclination, |r2|=1
    three_way_check("3D_30inc", r1, r2, PI / 2.0);
}

#[test]
fn test_three_way_3d_steep() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.5, 0.866]; // 60° inclination, |r2|=1
    three_way_check("3D_60inc", r1, r2, PI / 2.0);
}

#[test]
fn test_three_way_3d_both_off_plane() {
    let s: f64 = 1.0 / 3.0_f64.sqrt();
    let r1 = [s, s, s];
    let r2 = [-s, s, s];
    three_way_check("3D_both_off", r1, r2, 2.0);
}

// ────────────────────────────────────────────────────────────────
// Large eccentricity
// ────────────────────────────────────────────────────────────────

#[test]
fn test_three_way_high_eccentricity() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 7.0, 0.0];
    three_way_check("ecc_r7", r1, r2, 2.0);
}
