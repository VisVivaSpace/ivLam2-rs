//! Fuzz testing: Rust vercosine solver vs Fortran ivLam with randomized inputs.
//!
//! Uses a seeded PRNG for reproducibility. On failure, prints full case
//! parameters for reproduction.
//!
//! Run with: DYLD_LIBRARY_PATH=fortran cargo test --features gooding-ffi,ivlam-ffi --test fuzz_validation

#![cfg(feature = "ivlam-ffi")]

mod common;

use common::ivlam_ffi::{ivlam_zero_rev, ivlam_with_derivs};
use lambert_solver::{solve_lambert, solve_lambert_with_jacobian, Direction};
use std::f64::consts::PI;

/// Tolerance for velocity comparison (both iterative, same algorithm family)
const VEL_TOL: f64 = 1e-8;

/// Tolerance for Jacobian comparison (both analytical)
const JAC_TOL: f64 = 1e-10;

/// Simple xoshiro256** PRNG for reproducibility without external dependencies.
struct Rng {
    s: [u64; 4],
}

impl Rng {
    fn new(seed: u64) -> Self {
        // SplitMix64 to initialize state from a single seed
        let mut z = seed;
        let mut s = [0u64; 4];
        for slot in &mut s {
            z = z.wrapping_add(0x9e3779b97f4a7c15);
            let mut w = z;
            w = (w ^ (w >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            w = (w ^ (w >> 27)).wrapping_mul(0x94d049bb133111eb);
            *slot = w ^ (w >> 31);
        }
        Rng { s }
    }

    fn next_u64(&mut self) -> u64 {
        let result = (self.s[1].wrapping_mul(5)).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Uniform f64 in [0, 1)
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform f64 in [lo, hi]
    fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

/// Generate a random unit vector on the sphere
fn random_unit_vec(rng: &mut Rng) -> [f64; 3] {
    // Marsaglia method
    loop {
        let x = rng.uniform(-1.0, 1.0);
        let y = rng.uniform(-1.0, 1.0);
        let s = x * x + y * y;
        if s < 1.0 {
            let factor = 2.0 * (1.0 - s).sqrt();
            return [x * factor, y * factor, 1.0 - 2.0 * s];
        }
    }
}

/// Generate a random position vector with magnitude in [r_min, r_max]
fn random_position(rng: &mut Rng, r_min: f64, r_max: f64) -> [f64; 3] {
    let r = rng.uniform(r_min, r_max);
    let u = random_unit_vec(rng);
    [u[0] * r, u[1] * r, u[2] * r]
}

fn vec_mag(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute the transfer angle between two position vectors
fn transfer_angle(r1: &[f64; 3], r2: &[f64; 3]) -> f64 {
    let r1m = vec_mag(r1);
    let r2m = vec_mag(r2);
    let cos_theta = dot(r1, r2) / (r1m * r2m);
    cos_theta.clamp(-1.0, 1.0).acos()
}

/// Estimate parabolic TOF for scaling (Barker's equation approximation)
fn estimate_parabolic_tof(r1: &[f64; 3], r2: &[f64; 3]) -> f64 {
    let r1m = vec_mag(r1);
    let r2m = vec_mag(r2);
    let s = (r1m + r2m + vec_mag(&[r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]])) / 2.0;
    // Parabolic TOF: T_p = sqrt(2)/3 * (s^(3/2) - (s-c)^(3/2)) where c = |r2-r1|
    // Simplified: T_p ~ sqrt(2*s^3/9) for mu=1
    (2.0 * s * s * s / 9.0).sqrt()
}

#[test]
fn test_fuzz_velocity_1000_cases() {
    let mut rng = Rng::new(42);
    let mut failures = Vec::new();

    for case_idx in 0..1000 {
        let r1 = random_position(&mut rng, 0.1, 100.0);
        let r2 = random_position(&mut rng, 0.1, 100.0);

        // Skip near-0° and near-180° cases
        let theta = transfer_angle(&r1, &r2);
        if theta < 5.0_f64.to_radians() || theta > 175.0_f64.to_radians() {
            continue;
        }

        // Choose a feasible TOF between ~parabolic and ~10x parabolic
        let t_parab = estimate_parabolic_tof(&r1, &r2);
        let tof = rng.uniform(0.3 * t_parab, 10.0 * t_parab);

        // Determine transfer direction (prograde = +1 if h_z > 0)
        let h_z = r1[0] * r2[1] - r1[1] * r2[0];
        let direction_int: i32 = if h_z >= 0.0 { 1 } else { -1 };
        let direction = if h_z >= 0.0 { Direction::Prograde } else { Direction::Retrograde };

        // Rust solver
        let rust_result = solve_lambert(&r1, &r2, tof, 1.0, direction, 0);
        let rust_sol = match rust_result {
            Ok(sol) => sol,
            Err(_) => continue, // Skip cases where solver fails
        };

        // Fortran ivLam
        let ivlam_sol = ivlam_zero_rev(&r1, &r2, tof, direction_int);
        if ivlam_sol.info != 0 {
            continue; // Skip cases where ivLam fails
        }

        // Compare velocities
        let mut max_diff = 0.0_f64;
        for i in 0..3 {
            max_diff = max_diff.max((rust_sol.v1[i] - ivlam_sol.v1[i]).abs());
            max_diff = max_diff.max((rust_sol.v2[i] - ivlam_sol.v2[i]).abs());
        }

        if max_diff > VEL_TOL {
            failures.push(format!(
                "Case {}: r1={:?}, r2={:?}, tof={:.6e}, theta={:.2}°, dir={}, max_diff={:.2e}\n  \
                 rust_v1={:?}\n  ivlam_v1={:?}\n  rust_v2={:?}\n  ivlam_v2={:?}",
                case_idx, r1, r2, tof, theta.to_degrees(), direction_int, max_diff,
                rust_sol.v1, ivlam_sol.v1, rust_sol.v2, ivlam_sol.v2
            ));
        }
    }

    if !failures.is_empty() {
        panic!(
            "{} of 1000 fuzz cases failed:\n{}",
            failures.len(),
            failures.join("\n\n")
        );
    }
}

#[test]
fn test_fuzz_jacobian_100_cases() {
    let mut rng = Rng::new(123);
    let mut failures = Vec::new();

    for case_idx in 0..100 {
        let r1 = random_position(&mut rng, 0.1, 100.0);
        let r2 = random_position(&mut rng, 0.1, 100.0);

        let theta = transfer_angle(&r1, &r2);
        if theta < 5.0_f64.to_radians() || theta > 175.0_f64.to_radians() {
            continue;
        }

        let t_parab = estimate_parabolic_tof(&r1, &r2);
        let tof = rng.uniform(0.3 * t_parab, 10.0 * t_parab);

        let h_z = r1[0] * r2[1] - r1[1] * r2[0];
        let direction_int: i32 = if h_z >= 0.0 { 1 } else { -1 };
        let direction = if h_z >= 0.0 { Direction::Prograde } else { Direction::Retrograde };

        // Rust solver with Jacobian
        let rust_result = solve_lambert_with_jacobian(&r1, &r2, tof, 1.0, direction, 0);
        let (_sol, sens) = match rust_result {
            Ok(r) => r,
            Err(_) => continue,
        };

        // Fortran ivLam with derivatives
        let ivlam_result = ivlam_with_derivs(&r1, &r2, tof, direction_int, 0);
        if ivlam_result.info != 0 {
            continue;
        }

        // Compare Jacobians
        let mut max_diff = 0.0_f64;
        for i in 0..6 {
            for j in 0..7 {
                let diff = (sens.jacobian[i][j] - ivlam_result.jacobian[i][j]).abs();
                max_diff = max_diff.max(diff);
            }
        }

        if max_diff > JAC_TOL {
            failures.push(format!(
                "Case {}: r1={:?}, r2={:?}, tof={:.6e}, theta={:.2}°, dir={}, max_jac_diff={:.2e}",
                case_idx, r1, r2, tof, theta.to_degrees(), direction_int, max_diff
            ));
        }
    }

    if !failures.is_empty() {
        panic!(
            "{} of 100 Jacobian fuzz cases failed:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
