//! Validation tests for the Lambert Jacobian (first-order sensitivities).
//!
//! Three validation strategies:
//! 1. Finite-difference check (central difference, h=1e-7)
//! 2. Comparison against Fortran ivLam reference values (feature-gated)
//! 3. Energy consistency check (dE₁/dy = dE₂/dy)

mod common;

use lambert_solver::{solve_lambert, solve_lambert_with_jacobian, Direction};
use std::f64::consts::PI;

/// Tolerance for finite-difference vs analytical Jacobian.
/// Central difference error is O(h²) ≈ 1e-14, but amplified by condition number.
/// We use h=1e-7, giving h²=1e-14. Allow 1e-6 relative tolerance.
const FD_REL_TOL: f64 = 1e-6;
const FD_ABS_TOL: f64 = 1e-8;

fn finite_difference_jacobian(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    h: f64,
) -> [[f64; 7]; 6] {
    let mut jac = [[0.0; 7]; 6];

    // Pack inputs into a 7-vector: [r1x, r1y, r1z, r2x, r2y, r2z, tof]
    let y0 = [r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], tof];

    for j in 0..7 {
        let mut y_plus = y0;
        let mut y_minus = y0;
        y_plus[j] += h;
        y_minus[j] -= h;

        let r1_p = [y_plus[0], y_plus[1], y_plus[2]];
        let r2_p = [y_plus[3], y_plus[4], y_plus[5]];
        let r1_m = [y_minus[0], y_minus[1], y_minus[2]];
        let r2_m = [y_minus[3], y_minus[4], y_minus[5]];

        let sol_p = solve_lambert(&r1_p, &r2_p, y_plus[6], mu, direction, 0)
            .expect("FD+ solve failed");
        let sol_m = solve_lambert(&r1_m, &r2_m, y_minus[6], mu, direction, 0)
            .expect("FD- solve failed");

        let z_plus = [sol_p.v1[0], sol_p.v1[1], sol_p.v1[2],
                       sol_p.v2[0], sol_p.v2[1], sol_p.v2[2]];
        let z_minus = [sol_m.v1[0], sol_m.v1[1], sol_m.v1[2],
                        sol_m.v2[0], sol_m.v2[1], sol_m.v2[2]];

        for i in 0..6 {
            jac[i][j] = (z_plus[i] - z_minus[i]) / (2.0 * h);
        }
    }

    jac
}

fn check_jacobian_against_fd(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64) {
    let (sol, sens) = solve_lambert_with_jacobian(&r1, &r2, tof, mu, Direction::Prograde, 0)
        .unwrap_or_else(|e| panic!("{}: solver failed: {:?}", label, e));

    let jac_analytical = sens.jacobian;
    let jac_fd = finite_difference_jacobian(&r1, &r2, tof, mu, Direction::Prograde, 1e-7);

    // Check each element
    let z = [sol.v1[0], sol.v1[1], sol.v1[2], sol.v2[0], sol.v2[1], sol.v2[2]];
    for i in 0..6 {
        for j in 0..7 {
            let a = jac_analytical[i][j];
            let fd = jac_fd[i][j];
            let diff = (a - fd).abs();
            let scale = a.abs().max(fd.abs()).max(z[i].abs()).max(1e-10);
            let rel = diff / scale;

            assert!(
                diff < FD_ABS_TOL || rel < FD_REL_TOL,
                "{}: Jacobian[{}][{}] mismatch: analytical={:.10e}, fd={:.10e}, diff={:.2e}, rel={:.2e}",
                label, i, j, a, fd, diff, rel
            );
        }
    }
}

// ────────────────────────────────────────────────────────────────
// Finite-difference validation tests
// ────────────────────────────────────────────────────────────────

#[test]
fn test_jacobian_fd_90_deg() {
    check_jacobian_against_fd("90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0, 1.0);
}

#[test]
fn test_jacobian_fd_45_deg() {
    let angle = 45.0_f64.to_radians();
    check_jacobian_against_fd(
        "45°",
        [1.0, 0.0, 0.0],
        [angle.cos(), angle.sin(), 0.0],
        angle,
        1.0,
    );
}

#[test]
fn test_jacobian_fd_135_deg() {
    let angle = 135.0_f64.to_radians();
    check_jacobian_against_fd(
        "135°",
        [1.0, 0.0, 0.0],
        [angle.cos(), angle.sin(), 0.0],
        2.0,
        1.0,
    );
}

#[test]
fn test_jacobian_fd_3d() {
    check_jacobian_against_fd("3D", [1.0, 0.0, 0.0], [0.0, 0.8, 0.6], PI / 2.0, 1.0);
}

#[test]
fn test_jacobian_fd_different_radii() {
    check_jacobian_against_fd("diff-r", [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], 2.0, 1.0);
}

#[test]
fn test_jacobian_fd_hyperbolic() {
    check_jacobian_against_fd("hyper", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1, 1.0);
}

#[test]
fn test_jacobian_fd_physical_units() {
    // LEO to GEO-like transfer with mu in km³/s²
    let mu = 398600.4418;
    let r1 = [6678.0, 0.0, 0.0]; // LEO
    let r2 = [0.0, 42164.0, 0.0]; // GEO
    let tof = 19000.0; // seconds
    check_jacobian_against_fd("phys", r1, r2, tof, mu);
}

// ────────────────────────────────────────────────────────────────
// Energy consistency check
// ────────────────────────────────────────────────────────────────

#[test]
fn test_jacobian_energy_consistency() {
    // For a Keplerian arc, orbital energy E = v²/2 - μ/r is the same at both endpoints.
    // Therefore dE₁/dy = dE₂/dy for all inputs y.
    //
    // E₁ = v₁·v₁/2 - μ/r₁
    // dE₁/dyⱼ = Σᵢ v₁ᵢ · ∂v₁ᵢ/∂yⱼ + (μ/r₁²)·∂r₁/∂yⱼ
    //
    // E₂ = v₂·v₂/2 - μ/r₂
    // dE₂/dyⱼ = Σᵢ v₂ᵢ · ∂v₂ᵢ/∂yⱼ + (μ/r₂²)·∂r₂/∂yⱼ

    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    let tof = PI / 2.0;
    let mu = 1.0;

    let (sol, sens) = solve_lambert_with_jacobian(&r1, &r2, tof, mu, Direction::Prograde, 0)
        .unwrap();

    let r1_mag = (r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]).sqrt();
    let r2_mag = (r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]).sqrt();

    // Unit vectors
    let r1_hat = [r1[0] / r1_mag, r1[1] / r1_mag, r1[2] / r1_mag];
    let r2_hat = [r2[0] / r2_mag, r2[1] / r2_mag, r2[2] / r2_mag];

    for j in 0..7 {
        // dE₁/dyⱼ
        let mut de1 = 0.0;
        for i in 0..3 {
            de1 += sol.v1[i] * sens.jacobian[i][j];
        }
        // Add ∂(-μ/r₁)/∂yⱼ = (μ/r₁²)·∂r₁/∂yⱼ
        // ∂r₁/∂r₁ⱼ = r̂₁ⱼ (only for j=0,1,2)
        if j < 3 {
            de1 += (mu / (r1_mag * r1_mag)) * r1_hat[j];
        }

        // dE₂/dyⱼ
        let mut de2 = 0.0;
        for i in 0..3 {
            de2 += sol.v2[i] * sens.jacobian[3 + i][j];
        }
        if j >= 3 && j < 6 {
            de2 += (mu / (r2_mag * r2_mag)) * r2_hat[j - 3];
        }

        assert!(
            (de1 - de2).abs() < 1e-8,
            "Energy derivative mismatch for j={}: dE₁/dy={:.10e}, dE₂/dy={:.10e}, diff={:.2e}",
            j, de1, de2, (de1 - de2).abs()
        );
    }
}

// ────────────────────────────────────────────────────────────────
// Fortran ivLam comparison (feature-gated)
// ────────────────────────────────────────────────────────────────

#[cfg(feature = "ivlam-ffi")]
mod ivlam_comparison {
    use super::*;
    use crate::common::ivlam_ffi::ivlam_with_derivs;

    /// Tolerance for Rust vs Fortran ivLam Jacobian.
    /// Both are analytical, should agree to near machine precision.
    const IVLAM_TOL: f64 = 1e-10;

    fn check_vs_ivlam(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64) {
        // Rust solver (mu=1 to match ivLam)
        let (_sol, sens) =
            solve_lambert_with_jacobian(&r1, &r2, tof, 1.0, Direction::Prograde, 0)
                .unwrap_or_else(|e| panic!("{}: Rust solver failed: {:?}", label, e));

        // Fortran ivLam
        let ivlam_result = ivlam_with_derivs(&r1, &r2, tof, 1, 0);
        assert_eq!(
            ivlam_result.info, 0,
            "{}: ivLam failed: info={}",
            label, ivlam_result.info
        );
        // Compare — ivlam_result.jacobian is already in our convention
        for i in 0..6 {
            for j in 0..7 {
                let rust = sens.jacobian[i][j];
                let fortran = ivlam_result.jacobian[i][j];
                let diff = (rust - fortran).abs();
                let scale = rust.abs().max(fortran.abs()).max(1e-15);

                assert!(
                    diff < IVLAM_TOL || diff / scale < IVLAM_TOL,
                    "{}: Jacobian[{}][{}] Rust={:.12e} vs ivLam={:.12e}, diff={:.2e}",
                    label, i, j, rust, fortran, diff
                );
            }
        }
    }

    #[test]
    fn test_vs_ivlam_90_deg() {
        check_vs_ivlam("90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0);
    }

    #[test]
    fn test_vs_ivlam_45_deg() {
        let angle = 45.0_f64.to_radians();
        check_vs_ivlam("45°", [1.0, 0.0, 0.0], [angle.cos(), angle.sin(), 0.0], angle);
    }

    #[test]
    fn test_vs_ivlam_3d() {
        check_vs_ivlam("3D", [1.0, 0.0, 0.0], [0.0, 0.8, 0.6], PI / 2.0);
    }
}
