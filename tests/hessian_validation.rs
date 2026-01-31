//! Validation tests for the Lambert Hessians (second-order sensitivities).
//!
//! Three validation strategies:
//! 1. Finite-difference Hessians (central difference of Jacobians)
//! 2. Symmetry check: H_i[j][l] == H_i[l][j]
//! 3. Comparison against Fortran ivLam reference values (feature-gated)

mod common;

use lambert_solver::{solve_lambert_with_hessian, solve_lambert_with_jacobian, Direction};
use std::f64::consts::PI;

/// Tolerance for finite-difference vs analytical Hessian.
/// We compute FD Hessians by central-differencing the analytical Jacobians with h=1e-5.
/// FD error is O(h²) = 1e-10. Allow 1e-4 relative, 1e-6 absolute.
const FD_REL_TOL: f64 = 1e-4;
const FD_ABS_TOL: f64 = 1e-6;

/// Tolerance for symmetry check (should be machine precision).
const SYMMETRY_TOL: f64 = 1e-13;

/// Compute finite-difference Hessian by central-differencing the analytical Jacobian.
/// H_FD[i][j][l] ≈ (J(y+h·e_l) - J(y-h·e_l))[i][j] / (2h)
fn finite_difference_hessian(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    h: f64,
) -> [[[f64; 7]; 7]; 6] {
    let mut hess = [[[0.0; 7]; 7]; 6];
    let y0 = [r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], tof];

    for l in 0..7 {
        let mut y_plus = y0;
        let mut y_minus = y0;
        y_plus[l] += h;
        y_minus[l] -= h;

        let r1_p = [y_plus[0], y_plus[1], y_plus[2]];
        let r2_p = [y_plus[3], y_plus[4], y_plus[5]];
        let r1_m = [y_minus[0], y_minus[1], y_minus[2]];
        let r2_m = [y_minus[3], y_minus[4], y_minus[5]];

        let (_, sens_p) =
            solve_lambert_with_jacobian(&r1_p, &r2_p, y_plus[6], mu, direction, 0)
                .expect("FD+ Jacobian solve failed");
        let (_, sens_m) =
            solve_lambert_with_jacobian(&r1_m, &r2_m, y_minus[6], mu, direction, 0)
                .expect("FD- Jacobian solve failed");

        let inv_2h = 1.0 / (2.0 * h);
        for i in 0..6 {
            for j in 0..7 {
                hess[i][j][l] = (sens_p.jacobian[i][j] - sens_m.jacobian[i][j]) * inv_2h;
            }
        }
    }
    hess
}

fn check_hessian_against_fd(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64) {
    let (sol, sens) =
        solve_lambert_with_hessian(&r1, &r2, tof, mu, Direction::Prograde, 0)
            .unwrap_or_else(|e| panic!("{}: solver failed: {:?}", label, e));

    let analytical = sens.hessians.expect("Hessians should be Some");
    let fd = finite_difference_hessian(&r1, &r2, tof, mu, Direction::Prograde, 1e-5);

    let z = [sol.v1[0], sol.v1[1], sol.v1[2], sol.v2[0], sol.v2[1], sol.v2[2]];

    for i in 0..6 {
        for j in 0..7 {
            for l in 0..7 {
                let a = analytical[i][j][l];
                let f = fd[i][j][l];
                let diff = (a - f).abs();
                let scale = a.abs().max(f.abs()).max(z[i].abs()).max(1e-10);
                let rel = diff / scale;

                assert!(
                    diff < FD_ABS_TOL || rel < FD_REL_TOL,
                    "{}: Hessian[{}][{}][{}] mismatch: analytical={:.8e}, fd={:.8e}, diff={:.2e}, rel={:.2e}",
                    label, i, j, l, a, f, diff, rel
                );
            }
        }
    }
}

fn check_hessian_symmetry(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64) {
    let (_, sens) =
        solve_lambert_with_hessian(&r1, &r2, tof, mu, Direction::Prograde, 0)
            .unwrap_or_else(|e| panic!("{}: solver failed: {:?}", label, e));

    let hess = sens.hessians.expect("Hessians should be Some");

    for i in 0..6 {
        for j in 0..7 {
            for l in 0..7 {
                let diff = (hess[i][j][l] - hess[i][l][j]).abs();
                assert!(
                    diff < SYMMETRY_TOL,
                    "{}: Hessian[{}][{}][{}]={:.10e} != Hessian[{}][{}][{}]={:.10e}, diff={:.2e}",
                    label, i, j, l, hess[i][j][l], i, l, j, hess[i][l][j], diff
                );
            }
        }
    }
}

// ────────────────────────────────────────────────────────────────
// Finite-difference validation tests
// ────────────────────────────────────────────────────────────────

#[test]
fn test_hessian_fd_90_deg() {
    check_hessian_against_fd("90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0, 1.0);
}

#[test]
fn test_hessian_fd_45_deg() {
    let angle = 45.0_f64.to_radians();
    check_hessian_against_fd(
        "45°",
        [1.0, 0.0, 0.0],
        [angle.cos(), angle.sin(), 0.0],
        angle,
        1.0,
    );
}

#[test]
fn test_hessian_fd_135_deg() {
    let angle = 135.0_f64.to_radians();
    check_hessian_against_fd(
        "135°",
        [1.0, 0.0, 0.0],
        [angle.cos(), angle.sin(), 0.0],
        2.0,
        1.0,
    );
}

#[test]
fn test_hessian_fd_3d() {
    check_hessian_against_fd("3D", [1.0, 0.0, 0.0], [0.0, 0.8, 0.6], PI / 2.0, 1.0);
}

#[test]
fn test_hessian_fd_different_radii() {
    check_hessian_against_fd("diff-r", [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], 2.0, 1.0);
}

#[test]
fn test_hessian_fd_hyperbolic() {
    check_hessian_against_fd("hyper", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1, 1.0);
}

// ────────────────────────────────────────────────────────────────
// Symmetry validation tests
// ────────────────────────────────────────────────────────────────

#[test]
fn test_hessian_symmetry_90_deg() {
    check_hessian_symmetry("90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0, 1.0);
}

#[test]
fn test_hessian_symmetry_45_deg() {
    let angle = 45.0_f64.to_radians();
    check_hessian_symmetry("45°", [1.0, 0.0, 0.0], [angle.cos(), angle.sin(), 0.0], angle, 1.0);
}

#[test]
fn test_hessian_symmetry_135_deg() {
    let angle = 135.0_f64.to_radians();
    check_hessian_symmetry("135°", [1.0, 0.0, 0.0], [angle.cos(), angle.sin(), 0.0], 2.0, 1.0);
}

#[test]
fn test_hessian_symmetry_3d() {
    check_hessian_symmetry("3D", [1.0, 0.0, 0.0], [0.0, 0.8, 0.6], PI / 2.0, 1.0);
}

#[test]
fn test_hessian_symmetry_different_radii() {
    check_hessian_symmetry("diff-r", [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], 2.0, 1.0);
}

#[test]
fn test_hessian_symmetry_hyperbolic() {
    check_hessian_symmetry("hyper", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1, 1.0);
}

// ────────────────────────────────────────────────────────────────
// Fortran ivLam comparison (feature-gated)
// ────────────────────────────────────────────────────────────────

#[cfg(feature = "ivlam-ffi")]
mod ivlam_hessian_comparison {
    use super::*;
    use crate::common::ivlam_ffi::ivlam_with_hessians;

    /// Tolerance for Rust vs Fortran ivLam Hessians.
    /// Both are analytical, should agree to near machine precision.
    const IVLAM_TOL: f64 = 1e-8;

    fn check_vs_ivlam(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64) {
        let (_, sens) =
            solve_lambert_with_hessian(&r1, &r2, tof, 1.0, Direction::Prograde, 0)
                .unwrap_or_else(|e| panic!("{}: Rust solver failed: {:?}", label, e));

        let ivlam_result = ivlam_with_hessians(&r1, &r2, tof, 1, 0);
        assert_eq!(
            ivlam_result.info, 0,
            "{}: ivLam failed: info={}",
            label, ivlam_result.info
        );

        let rust_hess = sens.hessians.expect("Hessians should be Some");

        for i in 0..6 {
            for j in 0..7 {
                for l in 0..7 {
                    let rust = rust_hess[i][j][l];
                    let fortran = ivlam_result.hessians[i][j][l];
                    let diff = (rust - fortran).abs();
                    let scale = rust.abs().max(fortran.abs()).max(1e-15);

                    assert!(
                        diff < IVLAM_TOL || diff / scale < IVLAM_TOL,
                        "{}: Hessian[{}][{}][{}] Rust={:.12e} vs ivLam={:.12e}, diff={:.2e}",
                        label, i, j, l, rust, fortran, diff
                    );
                }
            }
        }
    }

    #[test]
    fn test_hessian_vs_ivlam_90_deg() {
        check_vs_ivlam("90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0);
    }

    #[test]
    fn test_hessian_vs_ivlam_45_deg() {
        let angle = 45.0_f64.to_radians();
        check_vs_ivlam("45°", [1.0, 0.0, 0.0], [angle.cos(), angle.sin(), 0.0], angle);
    }

    #[test]
    fn test_hessian_vs_ivlam_3d() {
        check_vs_ivlam("3D", [1.0, 0.0, 0.0], [0.0, 0.8, 0.6], PI / 2.0);
    }

    #[test]
    fn test_hessian_vs_ivlam_diff_radii() {
        check_vs_ivlam("diff_r", [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], 2.0);
    }

    #[test]
    fn test_hessian_vs_ivlam_hyperbolic() {
        check_vs_ivlam("hyper", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1);
    }

    #[test]
    fn test_hessian_vs_ivlam_large_angle() {
        let angle = 160.0_f64.to_radians();
        check_vs_ivlam("160°", [1.0, 0.0, 0.0], [angle.cos(), angle.sin(), 0.0], 2.0);
    }
}
