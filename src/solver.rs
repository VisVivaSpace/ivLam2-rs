//! Lambert solver main module.
//!
//! This module contains the main solver routines that iterate on k to find
//! the solution to Lambert's equation using the vercosine formulation.

use crate::geometry::{Direction, Geometry, NpiRevInfo};
use crate::sensitivities::LambertSensitivities;
use crate::stumpff::compute_w_and_derivatives;
use crate::velocity::{compute_velocities, SolverState};

use std::f64::consts::PI;

/// Error types for the Lambert solver
#[derive(Debug, Clone, PartialEq)]
pub enum LambertError {
    /// The position vectors are identical (r1 = r2)
    IdenticalPositions,
    /// Time of flight is non-positive
    InvalidTimeOfFlight,
    /// Requested number of revolutions has no solution (TOF too short)
    NoSolutionForRevolutions(i32),
    /// Solver failed to converge
    ConvergenceFailed {
        iterations: usize,
        residual: f64,
    },
    /// Half-revolution case detected - velocity solution is undefined
    HalfRevolutionSingularity,
    /// Input validation failed
    InvalidInput(String),
}

/// Solution to the Lambert problem
#[derive(Debug, Clone)]
pub struct LambertSolution {
    /// Initial velocity vector
    pub v1: [f64; 3],
    /// Final velocity vector
    pub v2: [f64; 3],
    /// Number of complete revolutions (0 for direct transfer)
    pub n_rev: i32,
    /// Number of iterations required
    pub iterations: usize,
    /// Final residual (should be near zero)
    pub residual: f64,
    /// The solved k parameter
    pub k: f64,
    /// Information about near-singularity warnings
    pub warning: Option<String>,
}

/// Multi-revolution solution containing both short and long period solutions
#[derive(Debug, Clone)]
pub struct MultiRevSolution {
    /// Short-period solution (positive N tilde)
    pub short_period: LambertSolution,
    /// Long-period solution (negative N tilde)  
    pub long_period: LambertSolution,
}

// Solver parameters
const MAX_ITERATIONS: usize = 25;
const SQRT_2: f64 = std::f64::consts::SQRT_2;
const TWO_PI: f64 = 2.0 * PI;

/// Margin from k boundaries
const K_MARGIN: f64 = 1e-10;
const K_RANGE_LEFT_ZERO_REV: f64 = -SQRT_2 + K_MARGIN;
const K_RANGE_RIGHT_MULTI_REV: f64 = SQRT_2 - K_MARGIN;
const K_RANGE_LEFT_MULTI_REV: f64 = -SQRT_2 + K_MARGIN;

/// Maximum step size for multi-rev safeguard
const MAX_K_STEP_MULTI_REV: f64 = 0.5;

/// Solve the Lambert problem for zero or single-N revolutions.
///
/// # Arguments
/// * `r1` - Initial position vector [x, y, z]
/// * `r2` - Final position vector [x, y, z]
/// * `tof` - Time of flight (must be > 0)
/// * `mu` - Gravitational parameter
/// * `direction` - Transfer direction (Prograde or Retrograde)
/// * `n_rev` - Number of complete revolutions (0 for direct transfer)
///
/// # Returns
/// `Ok(LambertSolution)` on success, `Err(LambertError)` on failure
///
/// # Example
/// ```
/// use lambert_solver::{solve_lambert, Direction};
/// 
/// let r1 = [1.0, 0.0, 0.0];
/// let r2 = [0.0, 1.0, 0.0];
/// let tof = std::f64::consts::PI / 2.0;
/// let mu = 1.0;
///
/// let result = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0);
/// ```
pub fn solve_lambert(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    n_rev: i32,
) -> Result<LambertSolution, LambertError> {
    // Input validation
    if tof <= 0.0 {
        return Err(LambertError::InvalidTimeOfFlight);
    }
    
    // Check for identical positions
    let dr = [r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]];
    let dr_mag = (dr[0].powi(2) + dr[1].powi(2) + dr[2].powi(2)).sqrt();
    let r1_mag = (r1[0].powi(2) + r1[1].powi(2) + r1[2].powi(2)).sqrt();
    if dr_mag < 1e-14 * r1_mag {
        return Err(LambertError::IdenticalPositions);
    }
    
    // Compute geometry
    let geom = Geometry::new(r1, r2, tof, mu, direction);
    
    // Check for half-rev singularity
    if geom.n_pi_rev_info == NpiRevInfo::ExactHalfRev {
        return Err(LambertError::HalfRevolutionSingularity);
    }
    
    // For n_rev = 0, solve directly
    // For n_rev != 0, we need to check if solution exists
    let n_abs = n_rev.abs() as usize;
    let n_is_zero = n_rev == 0;
    let two_pi_n = TWO_PI * n_abs as f64;
    
    // Determine k bounds
    let (k_left, k_right) = if n_is_zero {
        let k_right = if geom.tau > 0.0 {
            (1.0 / geom.tau) * (1.0 - K_MARGIN)
        } else {
            f64::MAX
        };
        (K_RANGE_LEFT_ZERO_REV, k_right)
    } else {
        (K_RANGE_LEFT_MULTI_REV, K_RANGE_RIGHT_MULTI_REV)
    };

    // Get initial guess and clamp to valid range
    let k_initial = compute_initial_guess(&geom, n_is_zero, n_rev);
    let k_initial = k_initial.clamp(k_left + K_MARGIN, k_right - K_MARGIN);

    // Create initial solver state
    let mut state = SolverState::new(k_initial, geom.tau);
    
    // Main iteration loop
    let mut k_last = state.k_sol;
    let mut residual = f64::MAX;
    
    for iter in 0..MAX_ITERATIONS {
        state.iterations = iter + 1;
        k_last = state.k_sol;
        
        // Compute W and its derivatives
        state.dw = compute_w_and_derivatives(state.k_sol, n_is_zero, two_pi_n, 3);
        state.w = state.dw[0];
        
        // Compute the TOF function and its derivatives
        let (func, dfunc) = compute_tof_function(&state, &geom);
        residual = func.abs();
        
        // Check for convergence
        if residual < 1e-14 * geom.tof_by_s.max(1.0) {
            break;
        }
        
        // Compute correction using Newton-Raphson (with higher-order terms)
        let dk = compute_correction(func, &dfunc, n_is_zero);
        
        // Apply correction with safeguards
        let mut k_new = state.k_sol + dk;
        
        // Bound checking
        if k_new < k_left {
            k_new = 0.5 * (k_last + k_left);
        } else if k_new > k_right {
            k_new = 0.5 * (k_last + k_right);
        }
        
        state.k_sol = k_new;
        state.update_p(geom.tau);
        
        // Check for lack of improvement
        if (state.k_sol - k_last).abs() < f64::EPSILON * state.k_sol.abs().max(1.0) {
            break;
        }
    }
    
    // Check convergence
    if state.iterations >= MAX_ITERATIONS && residual > 1e-10 {
        return Err(LambertError::ConvergenceFailed {
            iterations: state.iterations,
            residual,
        });
    }
    
    // Compute velocities
    let (v1, v2) = compute_velocities(&state, &geom);
    
    // Build warning message if needed
    let warning = match geom.n_pi_rev_info {
        NpiRevInfo::NearHalfRev => Some("Near half-revolution transfer - velocity accuracy may be degraded".to_string()),
        NpiRevInfo::NearFullRev => Some("Near full-revolution transfer - close to degenerate case".to_string()),
        _ => None,
    };
    
    Ok(LambertSolution {
        v1,
        v2,
        n_rev,
        iterations: state.iterations,
        residual,
        k: state.k_sol,
        warning,
    })
}

/// Solve for both multi-revolution solutions (short and long period).
///
/// For |N| > 0, there are generally two solutions: a "short period" (higher energy)
/// and a "long period" (lower energy) transfer.
///
/// # Arguments
/// * `r1`, `r2`, `tof`, `mu`, `direction` - Same as `solve_lambert`
/// * `n_rev` - Number of complete revolutions (must be > 0)
///
/// # Returns
/// Both solutions if they exist, or error if TOF is too short for the requested N
pub fn solve_lambert_multi_rev(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    n_rev: i32,
) -> Result<MultiRevSolution, LambertError> {
    if n_rev == 0 {
        return Err(LambertError::InvalidInput(
            "n_rev must be non-zero for multi-rev solver".to_string(),
        ));
    }
    
    let n_abs = n_rev.abs();
    
    // Solve for positive N tilde (short period)
    let short = solve_lambert(r1, r2, tof, mu, direction, n_abs)?;
    
    // Solve for negative N tilde (long period)
    let long = solve_lambert(r1, r2, tof, mu, direction, -n_abs)?;
    
    Ok(MultiRevSolution {
        short_period: short,
        long_period: long,
    })
}

/// Solve Lambert's problem and compute the analytical Jacobian.
///
/// Returns both the velocity solution and the 6×7 Jacobian ∂[v₁,v₂]/∂[r₁,r₂,T].
/// The Jacobian is computed analytically via the implicit function theorem,
/// adding roughly one iteration's worth of computation to the base solve.
///
/// # Arguments
/// Same as [`solve_lambert`].
///
/// # Returns
/// `Ok((LambertSolution, LambertSensitivities))` on success
pub fn solve_lambert_with_jacobian(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    n_rev: i32,
) -> Result<(LambertSolution, LambertSensitivities), LambertError> {
    // Input validation
    if tof <= 0.0 {
        return Err(LambertError::InvalidTimeOfFlight);
    }

    let dr = [r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]];
    let dr_mag = (dr[0].powi(2) + dr[1].powi(2) + dr[2].powi(2)).sqrt();
    let r1_mag = (r1[0].powi(2) + r1[1].powi(2) + r1[2].powi(2)).sqrt();
    if dr_mag < 1e-14 * r1_mag {
        return Err(LambertError::IdenticalPositions);
    }

    let geom = Geometry::new(r1, r2, tof, mu, direction);

    if geom.n_pi_rev_info == NpiRevInfo::ExactHalfRev {
        return Err(LambertError::HalfRevolutionSingularity);
    }

    let n_abs = n_rev.abs() as usize;
    let n_is_zero = n_rev == 0;
    let two_pi_n = TWO_PI * n_abs as f64;

    let (k_left, k_right) = if n_is_zero {
        let k_right = if geom.tau > 0.0 {
            (1.0 / geom.tau) * (1.0 - K_MARGIN)
        } else {
            f64::MAX
        };
        (K_RANGE_LEFT_ZERO_REV, k_right)
    } else {
        (K_RANGE_LEFT_MULTI_REV, K_RANGE_RIGHT_MULTI_REV)
    };

    let k_initial = compute_initial_guess(&geom, n_is_zero, n_rev);
    let k_initial = k_initial.clamp(k_left + K_MARGIN, k_right - K_MARGIN);

    let mut state = SolverState::new(k_initial, geom.tau);
    let mut residual = f64::MAX;

    for iter in 0..MAX_ITERATIONS {
        state.iterations = iter + 1;
        let k_last = state.k_sol;

        state.dw = compute_w_and_derivatives(state.k_sol, n_is_zero, two_pi_n, 3);
        state.w = state.dw[0];

        let (func, dfunc) = compute_tof_function(&state, &geom);
        residual = func.abs();

        if residual < 1e-14 * geom.tof_by_s.max(1.0) {
            break;
        }

        let dk = compute_correction(func, &dfunc, n_is_zero);
        let mut k_new = state.k_sol + dk;

        if k_new < k_left {
            k_new = 0.5 * (k_last + k_left);
        } else if k_new > k_right {
            k_new = 0.5 * (k_last + k_right);
        }

        state.k_sol = k_new;
        state.update_p(geom.tau);

        if (state.k_sol - k_last).abs() < f64::EPSILON * state.k_sol.abs().max(1.0) {
            break;
        }
    }

    if state.iterations >= MAX_ITERATIONS && residual > 1e-10 {
        return Err(LambertError::ConvergenceFailed {
            iterations: state.iterations,
            residual,
        });
    }

    let (v1, v2) = compute_velocities(&state, &geom);

    let warning = match geom.n_pi_rev_info {
        NpiRevInfo::NearHalfRev => Some("Near half-revolution transfer - velocity accuracy may be degraded".to_string()),
        NpiRevInfo::NearFullRev => Some("Near full-revolution transfer - close to degenerate case".to_string()),
        _ => None,
    };

    // Compute Jacobian using the converged state
    let sens = LambertSensitivities::compute_first_order(&state, &geom);

    let solution = LambertSolution {
        v1,
        v2,
        n_rev,
        iterations: state.iterations,
        residual,
        k: state.k_sol,
        warning,
    };

    Ok((solution, sens))
}

/// Solve Lambert's problem and compute both the Jacobian and Hessians.
///
/// Returns the velocity solution, the 6×7 Jacobian, and the 6 symmetric 7×7 Hessian
/// matrices ∂²z_i/∂y_j∂y_l. This is the full second-order sensitivity computation.
///
/// # Arguments
/// Same as [`solve_lambert`].
///
/// # Returns
/// `Ok((LambertSolution, LambertSensitivities))` where `sensitivities.hessians` is `Some`.
pub fn solve_lambert_with_hessian(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    n_rev: i32,
) -> Result<(LambertSolution, LambertSensitivities), LambertError> {
    if tof <= 0.0 {
        return Err(LambertError::InvalidTimeOfFlight);
    }

    let dr = [r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2]];
    let dr_mag = (dr[0].powi(2) + dr[1].powi(2) + dr[2].powi(2)).sqrt();
    let r1_mag = (r1[0].powi(2) + r1[1].powi(2) + r1[2].powi(2)).sqrt();
    if dr_mag < 1e-14 * r1_mag {
        return Err(LambertError::IdenticalPositions);
    }

    let geom = Geometry::new(r1, r2, tof, mu, direction);

    if geom.n_pi_rev_info == NpiRevInfo::ExactHalfRev {
        return Err(LambertError::HalfRevolutionSingularity);
    }

    let n_abs = n_rev.abs() as usize;
    let n_is_zero = n_rev == 0;
    let two_pi_n = TWO_PI * n_abs as f64;

    let (k_left, k_right) = if n_is_zero {
        let k_right = if geom.tau > 0.0 {
            (1.0 / geom.tau) * (1.0 - K_MARGIN)
        } else {
            f64::MAX
        };
        (K_RANGE_LEFT_ZERO_REV, k_right)
    } else {
        (K_RANGE_LEFT_MULTI_REV, K_RANGE_RIGHT_MULTI_REV)
    };

    let k_initial = compute_initial_guess(&geom, n_is_zero, n_rev);
    let k_initial = k_initial.clamp(k_left + K_MARGIN, k_right - K_MARGIN);

    let mut state = SolverState::new(k_initial, geom.tau);
    let mut residual = f64::MAX;

    for iter in 0..MAX_ITERATIONS {
        state.iterations = iter + 1;
        let k_last = state.k_sol;

        state.dw = compute_w_and_derivatives(state.k_sol, n_is_zero, two_pi_n, 3);
        state.w = state.dw[0];

        let (func, dfunc) = compute_tof_function(&state, &geom);
        residual = func.abs();

        if residual < 1e-14 * geom.tof_by_s.max(1.0) {
            break;
        }

        let dk = compute_correction(func, &dfunc, n_is_zero);
        let mut k_new = state.k_sol + dk;

        if k_new < k_left {
            k_new = 0.5 * (k_last + k_left);
        } else if k_new > k_right {
            k_new = 0.5 * (k_last + k_right);
        }

        state.k_sol = k_new;
        state.update_p(geom.tau);

        if (state.k_sol - k_last).abs() < f64::EPSILON * state.k_sol.abs().max(1.0) {
            break;
        }
    }

    if state.iterations >= MAX_ITERATIONS && residual > 1e-10 {
        return Err(LambertError::ConvergenceFailed {
            iterations: state.iterations,
            residual,
        });
    }

    let (v1, v2) = compute_velocities(&state, &geom);

    let warning = match geom.n_pi_rev_info {
        NpiRevInfo::NearHalfRev => Some("Near half-revolution transfer - velocity accuracy may be degraded".to_string()),
        NpiRevInfo::NearFullRev => Some("Near full-revolution transfer - close to degenerate case".to_string()),
        _ => None,
    };

    let sens = LambertSensitivities::compute_with_hessians(&state, &geom);

    let solution = LambertSolution {
        v1,
        v2,
        n_rev,
        iterations: state.iterations,
        residual,
        k: state.k_sol,
        warning,
    };

    Ok((solution, sens))
}

/// Compute initial guess for k based on problem parameters.
///
/// When the `lightweight` feature is NOT enabled (default), this uses interpolation
/// tables from the ivLam coefficient file for accurate initial guesses (~1-2 iterations).
/// Falls back to analytical approximation if interpolation returns None or if
/// the `lightweight` feature is enabled.
fn compute_initial_guess(geom: &Geometry, n_is_zero: bool, n_rev: i32) -> f64 {
    #[cfg(not(feature = "lightweight"))]
    {
        if n_is_zero {
            if let Some(k) = crate::interpolation::interpolate_initial_k_zero_rev(
                geom.tau, geom.tof_by_s,
            ) {
                return k;
            }
        } else {
            if let Some(k) = crate::interpolation::interpolate_initial_k_multi_rev(
                geom.tau, geom.tof_by_s, n_rev,
            ) {
                return k;
            }
        }
    }

    if n_is_zero {
        // Zero-rev: start near the middle of the valid range
        // A better guess would use parabolic flight time
        let k_parab = SQRT_2;
        
        // Estimate whether we're elliptic or hyperbolic
        let t_parab = geom.parabolic_tof_by_s();
        
        if geom.tof_by_s < t_parab {
            // Hyperbolic: k > sqrt(2)
            // Simple linear extrapolation
            let ratio = t_parab / geom.tof_by_s;
            k_parab * ratio.min(5.0).max(1.0)
        } else {
            // Elliptic: k < sqrt(2)
            // Start somewhat below parabolic
            let ratio = geom.tof_by_s / t_parab;
            k_parab / ratio.sqrt().min(10.0).max(1.0)
        }
    } else {
        // Multi-rev: start in the middle of the elliptic range
        // A better implementation would interpolate from coefficient tables
        0.0
    }
}

/// Compute the TOF function F(k) = sqrt(p) * (tau + p*W) - T*/S
/// and its derivatives with respect to k.
fn compute_tof_function(state: &SolverState, geom: &Geometry) -> (f64, [f64; 4]) {
    let p = state.p;
    let sqrt_p = state.sqrt_p;
    let dw = &state.dw;
    let tau = geom.tau;
    let tau_sq = geom.tau_sq;
    let tau_cu = geom.tau_cu;
    
    // F = sqrt(p) * (tau + p*W) - T*/S
    let left_side = sqrt_p * (tau + p * dw[0]);
    
    let func = if geom.huge_tof_case {
        // Use log form for numerical stability with huge TOF
        left_side.ln() - geom.log_tof_by_s
    } else {
        left_side - geom.tof_by_s
    };
    
    // Derivatives
    let p_sq = p * p;
    let one_by_sqrt_p = 1.0 / sqrt_p;
    let one_by_p = one_by_sqrt_p * one_by_sqrt_p;
    let one_by_p_32 = one_by_sqrt_p * one_by_p;
    
    // dF/dk = (-3*p*tau*W + 2*p^2*dW - tau^2) / (2*sqrt(p))
    let df1 = (-3.0 * p * tau * dw[0] + 2.0 * p_sq * dw[1] - tau_sq) * one_by_sqrt_p * 0.5;
    
    // d²F/dk² = (3*p*tau²*W + 4*p³*d²W - 12*p²*tau*dW - tau³) / (4*p^(3/2))
    let df2 = (3.0 * p * tau_sq * dw[0] + 4.0 * p_sq * p * dw[2] 
               - 12.0 * p_sq * tau * dw[1] - tau_cu) * one_by_p_32 * 0.25;
    
    // d³F/dk³ (for higher-order convergence)
    let tau_4 = tau_sq * tau_sq;
    let df3 = (3.0 * p * tau_cu * dw[0] + 18.0 * p_sq * tau_sq * dw[1]
               - 36.0 * p_sq * p * tau * dw[2] + 8.0 * p_sq * p_sq * dw[3] 
               - 3.0 * tau_4) * one_by_p_32 * one_by_p * 0.125;
    
    // If using log form, need to transform derivatives
    let dfunc = if geom.huge_tof_case {
        let inv_f = 1.0 / left_side;
        let inv_f_sq = inv_f * inv_f;
        let t1 = df1 * inv_f;
        let t2 = df2 * inv_f - t1 * t1;
        let t3 = df3 * inv_f - 3.0 * df2 * df1 * inv_f_sq + 2.0 * t1.powi(3);
        [df1, t1, t2, t3]
    } else {
        [df1, df1, df2, df3]
    };
    
    (func, dfunc)
}

/// Compute the correction step using up to third-order Newton-Raphson
fn compute_correction(func: f64, dfunc: &[f64; 4], n_is_zero: bool) -> f64 {
    // First-order Newton
    let df1 = dfunc[1];
    if df1.abs() < 1e-30 {
        return 0.0;
    }
    
    let inv_df1 = -1.0 / df1;
    let dk1 = func * inv_df1;
    
    // Apply step limit for multi-rev
    let dk1_limited = if !n_is_zero {
        dk1.clamp(-MAX_K_STEP_MULTI_REV, MAX_K_STEP_MULTI_REV)
    } else {
        dk1
    };
    
    // Second-order correction
    let dk2 = 0.5 * dk1_limited * dk1_limited * inv_df1 * dfunc[2];
    
    // Third-order correction
    let dk1_sq = dk1_limited * dk1_limited;
    let dk3 = (dk1_sq * dk1_limited * dfunc[3] / 6.0 + dk1_sq * dk2 * dfunc[2]) * inv_df1;
    
    // Check for series convergence
    let total = dk1_limited + dk2 + dk3;
    
    // If corrections are diverging, just use first-order
    if dk2.abs() > dk1_limited.abs() || dk3.abs() > dk2.abs() {
        dk1_limited
    } else {
        total
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_solve_lambert_90_degree() {
        // 90-degree transfer in circular orbit
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let mu = 1.0;
        
        // For a circular orbit with r=1, v=1, period = 2*pi
        // Quarter orbit takes pi/2 time units
        let tof = PI / 2.0;
        
        let result = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0);
        
        match result {
            Ok(sol) => {
                // Velocity magnitude should be close to circular (v = 1)
                // Iterative solver with 1e-14 convergence criterion
                let v1_mag = (sol.v1[0].powi(2) + sol.v1[1].powi(2) + sol.v1[2].powi(2)).sqrt();
                assert!((v1_mag - 1.0).abs() < 1e-10, "v1_mag = {}", v1_mag);
                
                // Should converge in few iterations
                assert!(sol.iterations < 10, "iterations = {}", sol.iterations);
            }
            Err(e) => panic!("Solver failed: {:?}", e),
        }
    }
    
    #[test]
    fn test_solve_lambert_180_degree() {
        // 180-degree (half-orbit) transfer
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [-1.0, 0.0, 0.0];
        let mu = 1.0;
        let tof = PI;  // Half period for circular orbit
        
        let result = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0);
        
        // This should produce a warning about near-half-rev
        match result {
            Ok(sol) => {
                // The transfer plane is undefined for exact 180°,
                // but we can still get a solution in the xy-plane
                assert!(sol.warning.is_some() || sol.v1[2].abs() < 1e-10);
            }
            Err(LambertError::HalfRevolutionSingularity) => {
                // This is also acceptable - exact 180° case
            }
            Err(e) => panic!("Unexpected error: {:?}", e),
        }
    }
    
    #[test]
    fn test_solve_lambert_hyperbolic() {
        // Very short TOF should give hyperbolic transfer
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let mu = 1.0;
        let tof = 0.1;  // Much shorter than circular quarter orbit
        
        let result = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0);
        
        match result {
            Ok(sol) => {
                // Should converge
                assert!(sol.residual < 1e-10);
                
                // Velocity should be higher than circular (hyperbolic)
                let v1_mag = (sol.v1[0].powi(2) + sol.v1[1].powi(2) + sol.v1[2].powi(2)).sqrt();
                assert!(v1_mag > 1.0, "Expected hyperbolic v > 1, got {}", v1_mag);
            }
            Err(e) => panic!("Solver failed: {:?}", e),
        }
    }
    
    #[test]
    fn test_solve_lambert_retrograde() {
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let mu = 1.0;
        let tof = 3.0 * PI / 2.0;  // 3/4 orbit for retrograde
        
        let result = solve_lambert(&r1, &r2, tof, mu, Direction::Retrograde, 0);
        
        match result {
            Ok(sol) => {
                // Retrograde should have negative angular momentum
                // h = r1 x v1
                let h_z = r1[0] * sol.v1[1] - r1[1] * sol.v1[0];
                assert!(h_z < 0.0, "Expected negative h_z for retrograde, got {}", h_z);
            }
            Err(e) => panic!("Solver failed: {:?}", e),
        }
    }
    
    #[test]
    fn test_error_identical_positions() {
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [1.0, 0.0, 0.0];  // Same as r1
        let mu = 1.0;
        let tof = 1.0;
        
        let result = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0);
        assert!(matches!(result, Err(LambertError::IdenticalPositions)));
    }
    
    #[test]
    fn test_error_invalid_tof() {
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let mu = 1.0;

        let result = solve_lambert(&r1, &r2, -1.0, mu, Direction::Prograde, 0);
        assert!(matches!(result, Err(LambertError::InvalidTimeOfFlight)));

        let result = solve_lambert(&r1, &r2, 0.0, mu, Direction::Prograde, 0);
        assert!(matches!(result, Err(LambertError::InvalidTimeOfFlight)));
    }

    #[test]
    fn test_solve_lambert_transfer_angle_sweep() {
        // Sweep through various transfer angles (prograde, zero-rev)
        // All with r1=r2=1, mu=1 (circular orbit reference)
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];

        for angle_deg in [10.0_f64, 30.0, 45.0, 60.0, 90.0, 120.0, 135.0, 150.0, 170.0] {
            let angle = angle_deg.to_radians();
            let r2 = [angle.cos(), angle.sin(), 0.0];
            // Use circular orbit TOF for this angle: t = angle (since period = 2*pi)
            let tof = angle;

            let result = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0);
            match result {
                Ok(sol) => {
                    // Iterative solver convergence (1e-10 tolerance)
                    assert!(sol.residual < 1e-10,
                        "angle={}°: residual {} too large", angle_deg, sol.residual);

                    // For circular orbit, v_mag should be ~1
                    let v1_mag = (sol.v1[0].powi(2) + sol.v1[1].powi(2) + sol.v1[2].powi(2)).sqrt();
                    assert!((v1_mag - 1.0).abs() < 1e-6,
                        "angle={}°: v1_mag = {}, expected ~1", angle_deg, v1_mag);
                }
                Err(e) => panic!("angle={}°: solver failed: {:?}", angle_deg, e),
            }
        }
    }

    #[test]
    fn test_solve_lambert_energy_conservation() {
        // Energy at departure should equal energy at arrival (same orbit)
        // E = v^2/2 - mu/r
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let mu = 1.0;
        let tof = PI / 2.0;

        let sol = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0).unwrap();

        let r1_mag = 1.0;
        let r2_mag = 1.0;
        let v1_mag_sq = sol.v1[0].powi(2) + sol.v1[1].powi(2) + sol.v1[2].powi(2);
        let v2_mag_sq = sol.v2[0].powi(2) + sol.v2[1].powi(2) + sol.v2[2].powi(2);

        let e1 = v1_mag_sq / 2.0 - mu / r1_mag;
        let e2 = v2_mag_sq / 2.0 - mu / r2_mag;

        // Energy must be conserved (both on same Keplerian arc)
        // Iterative solver tolerance
        assert!((e1 - e2).abs() < 1e-10,
            "Energy not conserved: E1={}, E2={}", e1, e2);
    }

    #[test]
    fn test_solve_lambert_3d_transfer() {
        // Out-of-plane transfer to verify 3D works
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 0.8, 0.6]; // |r2| = 1.0, in y-z plane
        let mu = 1.0;
        let tof = PI / 2.0;

        let sol = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0).unwrap();

        assert!(sol.residual < 1e-10);

        // v1 should have a non-zero z-component for out-of-plane transfer
        assert!(sol.v1[2].abs() > 1e-6, "Expected non-zero v1_z for 3D transfer");
    }

    #[test]
    fn test_solve_lambert_different_radii() {
        // Hohmann-like: r1 = 1, r2 = 2, 180° transfer
        // For Hohmann: a = (r1+r2)/2 = 1.5, period = 2*pi*sqrt(a^3/mu)
        // TOF = half period = pi*sqrt(1.5^3) = pi*sqrt(3.375)
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [-2.0, 0.0, 0.0]; // 180° away, r2=2
        let mu = 1.0;
        let a = 1.5_f64;
        let tof = PI * (a * a * a).sqrt(); // half-period of transfer ellipse

        let result = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0);

        // Near-180° may give warning or singularity
        match result {
            Ok(sol) => {
                // For Hohmann transfer: v1 = sqrt(mu*(2/r1 - 1/a))
                let v1_hohmann = (mu * (2.0 / 1.0 - 1.0 / a)).sqrt();
                let v1_mag = (sol.v1[0].powi(2) + sol.v1[1].powi(2) + sol.v1[2].powi(2)).sqrt();
                // This is near 180° so may have degraded accuracy
                assert!((v1_mag - v1_hohmann).abs() < 0.01,
                    "v1_mag = {}, expected Hohmann v1 = {}", v1_mag, v1_hohmann);
            }
            Err(LambertError::HalfRevolutionSingularity) => {
                // Acceptable for exact 180°
            }
            Err(e) => panic!("Unexpected error: {:?}", e),
        }
    }
}
