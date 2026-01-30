//! Velocity computation from the solved Lambert parameter.
//!
//! Once the iteration variable k (or equivalently p = 1 - k*tau) has been found,
//! this module computes the terminal velocity vectors using the f and g functions.

use crate::geometry::Geometry;

/// State of the Lambert iteration, containing the solved parameter
/// and intermediate values needed for velocity computation.
#[derive(Debug, Clone)]
pub struct SolverState {
    /// The iteration variable k at the solution
    pub k_sol: f64,
    /// p = 1 - k*tau
    pub p: f64,
    /// sqrt(p)
    pub sqrt_p: f64,
    /// W function value at k_sol
    pub w: f64,
    /// Derivatives of W (for sensitivities)
    pub dw: [f64; 5],
    /// Number of iterations taken
    pub iterations: usize,
}

impl SolverState {
    /// Create a new solver state with initial guess
    pub fn new(k_initial: f64, tau: f64) -> Self {
        let p = 1.0 - k_initial * tau;
        Self {
            k_sol: k_initial,
            p,
            sqrt_p: p.sqrt(),
            w: 0.0,
            dw: [0.0; 5],
            iterations: 0,
        }
    }
    
    /// Update p and sqrt_p from current k_sol
    #[inline]
    pub fn update_p(&mut self, tau: f64) {
        self.p = 1.0 - self.k_sol * tau;
        self.sqrt_p = self.p.sqrt();
    }
}

/// Lagrange f coefficient
/// f = 1 - p*(r1+r2)/r1
#[inline]
pub fn compute_f(p: f64, r1_plus_r2: f64, inv_r1: f64) -> f64 {
    1.0 - p * r1_plus_r2 * inv_r1
}

/// Lagrange g coefficient  
/// g = S * tau * sqrt(p)
#[inline]
pub fn compute_g(s: f64, tau: f64, sqrt_p: f64) -> f64 {
    s * tau * sqrt_p
}

/// Lagrange g_dot coefficient
/// g_dot = 1 - p*(r1+r2)/r2
#[inline]
pub fn compute_g_dot(p: f64, r1_plus_r2: f64, inv_r2: f64) -> f64 {
    1.0 - p * r1_plus_r2 * inv_r2
}

/// Compute terminal velocity vectors from the solved state.
///
/// Using the Lagrange f and g coefficients:
/// - v1 = (r2 - f*r1) / g
/// - v2 = (g_dot*r2 - r1) / g
///
/// # Arguments
/// * `state` - The solved Lambert state
/// * `geom` - The problem geometry
///
/// # Returns
/// Tuple of (v1, v2) velocity vectors
pub fn compute_velocities(state: &SolverState, geom: &Geometry) -> ([f64; 3], [f64; 3]) {
    let f = compute_f(state.p, geom.r1_plus_r2, geom.inv_r1);
    let g = compute_g(geom.s, geom.tau, state.sqrt_p);
    let g_dot = compute_g_dot(state.p, geom.r1_plus_r2, geom.inv_r2);
    
    // Add small value to g for near-half-rev cases to avoid division by zero
    let inv_g = 1.0 / (g + geom.add_to_g);
    
    let v1 = [
        (geom.r2_vec[0] - f * geom.r1_vec[0]) * inv_g,
        (geom.r2_vec[1] - f * geom.r1_vec[1]) * inv_g,
        (geom.r2_vec[2] - f * geom.r1_vec[2]) * inv_g,
    ];
    
    let v2 = [
        (g_dot * geom.r2_vec[0] - geom.r1_vec[0]) * inv_g,
        (g_dot * geom.r2_vec[1] - geom.r1_vec[1]) * inv_g,
        (g_dot * geom.r2_vec[2] - geom.r1_vec[2]) * inv_g,
    ];
    
    (v1, v2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Direction;
    use std::f64::consts::PI;
    
    #[test]
    fn test_lagrange_coefficients() {
        // For a simple test case
        let p = 0.5;
        let r1_plus_r2 = 2.0;
        let inv_r1 = 1.0;
        let inv_r2 = 1.0;
        
        let f = compute_f(p, r1_plus_r2, inv_r1);
        let g_dot = compute_g_dot(p, r1_plus_r2, inv_r2);
        
        // f = 1 - 0.5 * 2.0 * 1.0 = 0
        assert!((f - 0.0).abs() < 1e-15);
        
        // g_dot = 1 - 0.5 * 2.0 * 1.0 = 0
        assert!((g_dot - 0.0).abs() < 1e-15);
    }
    
    #[test]
    fn test_velocity_computation_basic() {
        // Create a simple geometry
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let tof = PI / 2.0;
        let mu = 1.0;
        
        let geom = Geometry::new(&r1, &r2, tof, mu, Direction::Prograde);
        
        // Create a state with some reasonable values
        let mut state = SolverState::new(0.5, geom.tau);
        state.w = 0.5;
        
        let (v1, v2) = compute_velocities(&state, &geom);
        
        // Velocities should be finite
        for i in 0..3 {
            assert!(v1[i].is_finite());
            assert!(v2[i].is_finite());
        }
    }
}
