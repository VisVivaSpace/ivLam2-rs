//! Round-trip integration tests for the Lambert solver.
//!
//! Solve Lambert → propagate (r1, v1) forward by TOF → verify arrival at r2.
//! Tolerance: 1e-8 (two iterative methods in series).

mod common;

use common::kepler::kepler_propagate;
use lambert_solver::{solve_lambert, Direction};
use std::f64::consts::PI;

fn vec_mag(v: &[f64; 3]) -> f64 {
    (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt()
}

fn vec_diff_mag(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Solve Lambert then propagate and check arrival.
fn round_trip_check(r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64, dir: Direction, tol: f64) {
    let sol = solve_lambert(&r1, &r2, tof, mu, dir, 0)
        .expect("Lambert solver failed");

    // Propagate r1, v1 forward by tof
    let (r2_prop, v2_prop) = kepler_propagate(&r1, &sol.v1, tof, mu);

    let pos_err = vec_diff_mag(&r2, &r2_prop);
    let r2_mag = vec_mag(&r2);
    assert!(
        pos_err < tol * r2_mag,
        "Position error: {:.2e} (relative: {:.2e}), tol={:.0e}",
        pos_err,
        pos_err / r2_mag,
        tol
    );

    // Also check v2 agreement
    let vel_err = vec_diff_mag(&sol.v2, &v2_prop);
    let v2_mag = vec_mag(&sol.v2);
    assert!(
        vel_err < tol * v2_mag,
        "Velocity error: {:.2e} (relative: {:.2e}), tol={:.0e}",
        vel_err,
        vel_err / v2_mag,
        tol
    );
}

#[test]
fn test_round_trip_circular_90_deg() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    round_trip_check(r1, r2, PI / 2.0, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn test_round_trip_circular_45_deg() {
    let angle = PI / 4.0;
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    round_trip_check(r1, r2, angle, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn test_round_trip_circular_150_deg() {
    let angle = 150.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    round_trip_check(r1, r2, angle, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn test_round_trip_elliptic_different_radii() {
    // r1 = 1, r2 = 1.5, 90° transfer
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.5, 0.0];
    round_trip_check(r1, r2, 2.0, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn test_round_trip_hyperbolic() {
    // Short TOF => hyperbolic
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    round_trip_check(r1, r2, 0.1, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn test_round_trip_retrograde() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    round_trip_check(r1, r2, 3.0 * PI / 2.0, 1.0, Direction::Retrograde, 1e-8);
}

#[test]
fn test_round_trip_3d() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.8, 0.6];
    round_trip_check(r1, r2, PI / 2.0, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn test_round_trip_physical_units() {
    // LEO to GEO-like transfer in km, s
    let mu = 398600.4418; // km^3/s^2
    let r1 = [6678.0, 0.0, 0.0]; // LEO ~300 km altitude
    let r2 = [0.0, 42164.0, 0.0]; // GEO
    let tof = 5.0 * 3600.0; // 5 hours
    round_trip_check(r1, r2, tof, mu, Direction::Prograde, 1e-8);
}
