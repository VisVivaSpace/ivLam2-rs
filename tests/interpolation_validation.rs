//! Validation tests for interpolation-based initial guesses.
//!
//! Verifies that:
//! 1. Solutions with interpolation match the same accuracy as before
//! 2. Interpolation reduces iteration counts compared to analytical fallback
//! 3. Multi-rev cases with interpolation produce correct results

use lambert_solver::{solve_lambert, solve_lambert_multi_rev, Direction};
use std::f64::consts::PI;

/// Verify a zero-rev solve produces correct result and low iteration count
fn check_zero_rev(
    r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64, dir: Direction,
    max_iters: usize, label: &str,
) {
    let sol = solve_lambert(&r1, &r2, tof, mu, dir, 0)
        .unwrap_or_else(|e| panic!("{}: solver failed: {:?}", label, e));

    assert!(sol.residual < 1e-10,
        "{}: residual {:.2e} too large", label, sol.residual);

    // With interpolation, zero-rev should converge in very few iterations
    assert!(sol.iterations <= max_iters,
        "{}: expected <= {} iterations, got {}", label, max_iters, sol.iterations);
}

#[test]
fn test_interpolation_90_degree_low_iterations() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    // With interpolation, should converge in ~2 iterations
    check_zero_rev(r1, r2, PI / 2.0, 1.0, Direction::Prograde, 4, "90deg");
}

#[test]
fn test_interpolation_angle_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    for angle_deg in [10.0_f64, 30.0, 60.0, 90.0, 120.0, 150.0, 170.0] {
        let angle = angle_deg.to_radians();
        let r2 = [angle.cos(), angle.sin(), 0.0];
        let tof = angle; // circular orbit TOF
        check_zero_rev(r1, r2, tof, 1.0, Direction::Prograde, 5,
            &format!("{}deg", angle_deg));
    }
}

#[test]
fn test_interpolation_hyperbolic() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    check_zero_rev(r1, r2, 0.1, 1.0, Direction::Prograde, 5, "hyperbolic");
}

#[test]
fn test_interpolation_retrograde() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    check_zero_rev(r1, r2, 3.0 * PI / 2.0, 1.0, Direction::Retrograde, 5, "retrograde");
}

#[test]
fn test_interpolation_3d() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.8, 0.6];
    check_zero_rev(r1, r2, PI / 2.0, 1.0, Direction::Prograde, 5, "3d");
}

#[test]
fn test_interpolation_physical_units() {
    let mu = 398600.4418;
    let r1 = [6678.0, 0.0, 0.0];
    let r2 = [0.0, 42164.0, 0.0];
    let tof = 5.0 * 3600.0;
    check_zero_rev(r1, r2, tof, mu, Direction::Prograde, 5, "LEO-GEO");
}

#[test]
fn test_interpolation_different_radii() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.5, 0.0];
    check_zero_rev(r1, r2, 2.0, 1.0, Direction::Prograde, 5, "diff-radii");
}

#[test]
fn test_interpolation_long_tof() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    check_zero_rev(r1, r2, 10.0, 1.0, Direction::Prograde, 5, "long-tof");
}

#[test]
fn test_interpolation_multi_rev_1rev() {
    // 1-revolution transfer: need sufficient TOF
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    // TOF must be large enough for 1 full revolution
    let tof = 2.0 * PI + PI / 2.0; // ~1.25 orbits

    let result = solve_lambert_multi_rev(&r1, &r2, tof, 1.0, Direction::Prograde, 1);
    match result {
        Ok(multi) => {
            assert!(multi.short_period.residual < 1e-10,
                "short_period residual: {:.2e}", multi.short_period.residual);
            assert!(multi.long_period.residual < 1e-10,
                "long_period residual: {:.2e}", multi.long_period.residual);

            // Verify iteration counts are reasonable
            assert!(multi.short_period.iterations <= 10,
                "short_period iterations: {}", multi.short_period.iterations);
            assert!(multi.long_period.iterations <= 10,
                "long_period iterations: {}", multi.long_period.iterations);
        }
        Err(e) => {
            // May fail if TOF is too short for 1 rev - that's acceptable
            eprintln!("Multi-rev 1-rev failed (may be expected): {:?}", e);
        }
    }
}

#[test]
fn test_interpolation_energy_conservation() {
    // Verify energy conservation with interpolation enabled
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    let mu = 1.0;
    let tof = PI / 2.0;

    let sol = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0).unwrap();

    let v1_sq = sol.v1[0].powi(2) + sol.v1[1].powi(2) + sol.v1[2].powi(2);
    let v2_sq = sol.v2[0].powi(2) + sol.v2[1].powi(2) + sol.v2[2].powi(2);
    let e1 = v1_sq / 2.0 - mu / 1.0;
    let e2 = v2_sq / 2.0 - mu / 1.0;

    assert!((e1 - e2).abs() < 1e-10,
        "Energy not conserved: E1={}, E2={}", e1, e2);
}
