//! Capture derivative reference values from Fortran ivLam.
//!
//! These tests call ivLam_NtildeWithDerivs to get first-order Jacobians
//! on canonical test cases. The values are printed and can be used as
//! reference for validating our Rust derivative implementation.

#![cfg(feature = "ivlam-ffi")]

mod common;

// Import from library crate to ensure build.rs link flags propagate
use lambert_solver::solve_lambert as _;
use common::ivlam_ffi::{ivlam_with_derivs, ivlam_zero_rev};
use std::f64::consts::PI;

#[test]
fn test_ivlam_derivatives_90_deg() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    let tof = PI / 2.0;

    let result = ivlam_with_derivs(&r1, &r2, tof, 1, 0);
    assert_eq!(result.info, 0, "ivLam deriv solver failed: info={}", result.info);

    let jac = &result.jacobian;

    // Verify velocities match zero-rev solver
    let basic = ivlam_zero_rev(&r1, &r2, tof, 1);
    for i in 0..3 {
        assert!((result.v1[i] - basic.v1[i]).abs() < 1e-14,
            "Deriv v1[{}] doesn't match basic: {} vs {}", i, result.v1[i], basic.v1[i]);
        assert!((result.v2[i] - basic.v2[i]).abs() < 1e-14,
            "Deriv v2[{}] doesn't match basic: {} vs {}", i, result.v2[i], basic.v2[i]);
    }

    // Print reference Jacobian for Phase 4 validation
    println!("=== ivLam Reference Jacobian for 90° transfer (mu=1) ===");
    println!("r1 = {:?}", r1);
    println!("r2 = {:?}", r2);
    println!("tof = {}", tof);
    println!("v1 = {:?}", result.v1);
    println!("v2 = {:?}", result.v2);
    println!("Jacobian dz/dy [6x7]:");
    for i in 0..6 {
        print!("  [");
        for j in 0..7 {
            if j > 0 { print!(", "); }
            print!("{:+.14e}", jac[i][j]);
        }
        println!("]");
    }

    // Jacobian should have non-zero entries
    let mut has_nonzero = false;
    for i in 0..6 {
        for j in 0..7 {
            if jac[i][j].abs() > 1e-15 {
                has_nonzero = true;
            }
        }
    }
    assert!(has_nonzero, "Jacobian is all zeros — derivative computation likely failed");
}

#[test]
fn test_ivlam_derivatives_45_deg() {
    let angle = 45.0_f64.to_radians();
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [angle.cos(), angle.sin(), 0.0];
    let tof = angle;

    let result = ivlam_with_derivs(&r1, &r2, tof, 1, 0);
    assert_eq!(result.info, 0, "ivLam deriv solver failed: info={}", result.info);

    let jac = &result.jacobian;

    println!("\n=== ivLam Reference Jacobian for 45° transfer (mu=1) ===");
    println!("r1 = {:?}", r1);
    println!("r2 = {:?}", r2);
    println!("tof = {}", tof);
    println!("v1 = {:?}", result.v1);
    println!("v2 = {:?}", result.v2);
    println!("Jacobian dz/dy [6x7]:");
    for i in 0..6 {
        print!("  [");
        for j in 0..7 {
            if j > 0 { print!(", "); }
            print!("{:+.14e}", jac[i][j]);
        }
        println!("]");
    }
}

#[test]
fn test_ivlam_derivatives_3d() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 0.8, 0.6];
    let tof = PI / 2.0;

    let result = ivlam_with_derivs(&r1, &r2, tof, 1, 0);
    assert_eq!(result.info, 0, "ivLam deriv solver failed: info={}", result.info);

    let jac = &result.jacobian;

    println!("\n=== ivLam Reference Jacobian for 3D transfer (mu=1) ===");
    println!("r1 = {:?}", r1);
    println!("r2 = {:?}", r2);
    println!("tof = {}", tof);
    println!("v1 = {:?}", result.v1);
    println!("v2 = {:?}", result.v2);
    println!("Jacobian dz/dy [6x7]:");
    for i in 0..6 {
        print!("  [");
        for j in 0..7 {
            if j > 0 { print!(", "); }
            print!("{:+.14e}", jac[i][j]);
        }
        println!("]");
    }
}
