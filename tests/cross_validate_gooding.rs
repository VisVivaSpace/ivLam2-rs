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

// ────────────────────────────────────────────────────────────────
// Angle sweep (equal radii, r1=r2=1, prograde, TOF = circular orbit time)
// ────────────────────────────────────────────────────────────────

#[test]
fn test_cross_validate_angle_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    for angle_deg in [5.0_f64, 15.0, 25.0, 30.0, 60.0, 75.0, 100.0, 110.0, 140.0, 150.0, 160.0] {
        let angle = angle_deg.to_radians();
        let r2 = [angle.cos(), angle.sin(), 0.0];
        // TOF = circular orbit time for this angle (period = 2*pi, so t = angle)
        let tof = angle;
        cross_validate(
            &format!("angle_sweep_{}°", angle_deg),
            r1, r2, tof, 1.0, CROSS_TOL,
        );
    }
}

// ────────────────────────────────────────────────────────────────
// Radius ratio sweep (90° transfer, varying r2/r1)
// ────────────────────────────────────────────────────────────────

#[test]
fn test_cross_validate_radius_ratio_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    for ratio in [0.3_f64, 0.5, 0.7, 1.5, 2.0, 3.0, 5.0] {
        let r2 = [0.0, ratio, 0.0];
        // TOF: use a value that gives a reasonable elliptic transfer
        // For Hohmann-like semi-major axis a = (1 + ratio)/2, period ~ 2*pi*sqrt(a^3)
        // Use half the transfer time as a reasonable TOF
        let a = (1.0 + ratio) / 2.0;
        let tof = PI * (a * a * a).sqrt() * 0.5;
        cross_validate(
            &format!("ratio_{:.1}", ratio),
            r1, r2, tof, 1.0, CROSS_TOL,
        );
    }
}

// ────────────────────────────────────────────────────────────────
// TOF variation (fixed 90° geometry, r1=[1,0,0], r2=[0,1,0])
// ────────────────────────────────────────────────────────────────

#[test]
fn test_cross_validate_tof_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    // Spans hyperbolic (short TOF) through slow elliptic (long TOF)
    for tof in [0.2_f64, 0.5, 1.0, 2.0, 5.0, 10.0] {
        cross_validate(
            &format!("tof_{:.1}", tof),
            r1, r2, tof, 1.0, CROSS_TOL,
        );
    }
}

// ────────────────────────────────────────────────────────────────
// 3D transfers (out-of-plane)
// ────────────────────────────────────────────────────────────────

#[test]
fn test_cross_validate_3d_30_deg_inclination() {
    // 30° inclination: r2 tilted out of xy-plane
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.866, 0.5]; // |r2| = 1.0
    cross_validate("3D_30inc", r1, r2, PI / 2.0, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_3d_60_deg_inclination() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.5, 0.866]; // |r2| = 1.0
    cross_validate("3D_60inc", r1, r2, PI / 2.0, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_3d_45_deg_both_off_plane() {
    // Both vectors have z-components
    let s: f64 = 1.0 / 3.0_f64.sqrt(); // ≈ 0.577
    let r1 = [s, s, s]; // |r1| = 1.0
    let r2 = [-s, s, s]; // |r2| = 1.0
    cross_validate("3D_both_off", r1, r2, 2.0, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_3d_different_radii_inclined() {
    let r1 = [2.0, 0.0, 0.0];
    let r2 = [0.0, 1.5, 1.0]; // |r2| ≈ 1.803
    cross_validate("3D_diff_r_inc", r1, r2, 3.0, 1.0, CROSS_TOL);
}

// ────────────────────────────────────────────────────────────────
// Physical units
// ────────────────────────────────────────────────────────────────

#[test]
fn test_cross_validate_leo_to_meo() {
    let mu = 398600.4418; // km³/s²
    let r1 = [6678.0, 0.0, 0.0]; // LEO ~300 km alt
    let r2 = [0.0, 26560.0, 0.0]; // MEO (GPS orbit)
    let tof = 4.0 * 3600.0; // 4 hours
    cross_validate("LEO→MEO", r1, r2, tof, mu, CROSS_TOL);
}

#[test]
fn test_cross_validate_earth_mars_type() {
    // Sun-centered, AU-based (mu_sun ≈ 1.327e11 km³/s²)
    let mu = 1.327e11;
    let r_earth = 1.496e8; // km
    let r_mars = 2.279e8; // km
    let angle = 135.0_f64.to_radians();
    let r1 = [r_earth, 0.0, 0.0];
    let r2 = [r_mars * angle.cos(), r_mars * angle.sin(), 0.0];
    let tof = 200.0 * 86400.0; // 200 days in seconds
    cross_validate("Earth→Mars", r1, r2, tof, mu, CROSS_TOL);
}

#[test]
fn test_cross_validate_sun_centered_venus() {
    let mu = 1.327e11;
    let r_earth = 1.496e8;
    let r_venus = 1.082e8;
    let angle = 60.0_f64.to_radians();
    let r1 = [r_earth, 0.0, 0.0];
    let r2 = [r_venus * angle.cos(), r_venus * angle.sin(), 0.0];
    let tof = 120.0 * 86400.0; // 120 days
    cross_validate("Earth→Venus", r1, r2, tof, mu, CROSS_TOL);
}

// ────────────────────────────────────────────────────────────────
// Large eccentricity (r2/r1 > 5, short TOF → highly eccentric)
// ────────────────────────────────────────────────────────────────

#[test]
fn test_cross_validate_high_eccentricity_r_ratio_7() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 7.0, 0.0]; // r2/r1 = 7
    // Short TOF relative to circular: fast, eccentric transfer
    let tof = 2.0;
    cross_validate("ecc_r7", r1, r2, tof, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_high_eccentricity_r_ratio_10() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 10.0, 0.0]; // r2/r1 = 10
    let tof = 3.0;
    cross_validate("ecc_r10", r1, r2, tof, 1.0, CROSS_TOL);
}

#[test]
fn test_cross_validate_high_eccentricity_60_deg() {
    // Large radius ratio with non-90° angle
    let angle = 60.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [6.0 * angle.cos(), 6.0 * angle.sin(), 0.0]; // |r2| = 6
    let tof = 2.5;
    cross_validate("ecc_r6_60°", r1, r2, tof, 1.0, CROSS_TOL);
}
