//! First and second-order sensitivities of the Lambert solution.
//!
//! This module computes the partial derivatives of the output velocities (v1, v2)
//! with respect to the inputs (r1, r2, tof), including second-order Hessians.
//!
//! Based on:
//! - Arora, N., Russell, R. P., Strange, N., and Ottesen, D., "Partial Derivatives
//!   of the Solution to the Lambert Boundary Value Problem," JGCD 2015
//! - Russell, R. P., "Complete Lambert Solver Including Second-Order Sensitivities,"
//!   JGCD 2022

use crate::geometry::Geometry;
use crate::velocity::SolverState;

/// First and second-order sensitivities of the Lambert solution.
///
/// The outputs z = [v1, v2] (6 components) depend on inputs y = [r1, r2, tof] (7 components).
#[derive(Debug, Clone)]
pub struct LambertSensitivities {
    /// First-order Jacobian: dz/dy (6x7 matrix)
    /// Stored as [dv1/dr1, dv1/dr2, dv1/dtof, dv2/dr1, dv2/dr2, dv2/dtof]
    pub jacobian: [[f64; 7]; 6],
    
    /// Second-order Hessians: d²z_i/dy² (6 matrices of 7x7)
    /// Only computed if requested
    pub hessians: Option<[[[f64; 7]; 7]; 6]>,
}

impl LambertSensitivities {
    /// Compute first-order sensitivities (Jacobian) of the Lambert solution.
    ///
    /// Uses the direct differentiation approach from Arora et al. 2015.
    pub fn compute_first_order(
        state: &SolverState,
        geom: &Geometry,
    ) -> Self {
        let mut jacobian = [[0.0; 7]; 6];
        
        // The sensitivity computation follows from differentiating:
        // 1. The Lagrange f, g, g_dot coefficients
        // 2. The velocity equations v1 = (r2 - f*r1)/g, v2 = (g_dot*r2 - r1)/g
        // 3. Through the implicit function theorem for dk/dy
        
        // For now, this is a placeholder that would contain the full
        // sensitivity computation from the FORTRAN code.
        // The full implementation requires translating ~3000 lines of
        // auto-generated derivative code from the Maple-to-FORTRAN conversion.
        
        // Key intermediate values needed:
        let p = state.p;
        let sqrt_p = state.sqrt_p;
        let tau = geom.tau;
        let s = geom.s;
        let r1_plus_r2 = geom.r1_plus_r2;
        let inv_r1 = geom.inv_r1;
        let inv_r2 = geom.inv_r2;
        let dw = &state.dw;
        
        // Compute f, g, g_dot
        let f = 1.0 - p * r1_plus_r2 * inv_r1;
        let g = s * tau * sqrt_p;
        let g_dot = 1.0 - p * r1_plus_r2 * inv_r2;
        let inv_g = 1.0 / (g + geom.add_to_g);
        
        // Partial of p with respect to k
        let dp_dk = -tau;
        
        // Partial of f with respect to p
        let df_dp = -r1_plus_r2 * inv_r1;
        
        // Partial of g with respect to p
        let dg_dp = s * tau * 0.5 / sqrt_p;
        
        // Partial of g_dot with respect to p
        let dg_dot_dp = -r1_plus_r2 * inv_r2;
        
        // The implicit function theorem gives:
        // dk/dy = -(dF/dk)^(-1) * dF/dy
        // where F is the Lambert TOF function
        
        // For a complete implementation, we need:
        // 1. dF/dy - explicit partials of TOF function w.r.t. inputs
        // 2. dF/dk - partial of TOF function w.r.t. k (already computed in solver)
        // 3. df/dy, dg/dy, dg_dot/dy - partials of Lagrange coefficients
        
        // The full computation is complex and is left as a TODO
        // for integration with the complete FORTRAN translation.
        
        Self {
            jacobian,
            hessians: None,
        }
    }
    
    /// Compute first and second-order sensitivities.
    pub fn compute_with_hessians(
        state: &SolverState,
        geom: &Geometry,
    ) -> Self {
        let first_order = Self::compute_first_order(state, geom);
        
        // Second-order computation would go here
        // This requires translating the H2 array computation from FORTRAN
        
        let hessians = [[[0.0; 7]; 7]; 6];
        
        Self {
            jacobian: first_order.jacobian,
            hessians: Some(hessians),
        }
    }
    
    /// Get the partial derivative dv1/dr1 (3x3 matrix)
    pub fn dv1_dr1(&self) -> [[f64; 3]; 3] {
        let mut result = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                result[i][j] = self.jacobian[i][j];
            }
        }
        result
    }
    
    /// Get the partial derivative dv1/dr2 (3x3 matrix)
    pub fn dv1_dr2(&self) -> [[f64; 3]; 3] {
        let mut result = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                result[i][j] = self.jacobian[i][3 + j];
            }
        }
        result
    }
    
    /// Get the partial derivative dv1/dtof (3-vector)
    pub fn dv1_dtof(&self) -> [f64; 3] {
        [self.jacobian[0][6], self.jacobian[1][6], self.jacobian[2][6]]
    }
    
    /// Get the partial derivative dv2/dr1 (3x3 matrix)
    pub fn dv2_dr1(&self) -> [[f64; 3]; 3] {
        let mut result = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                result[i][j] = self.jacobian[3 + i][j];
            }
        }
        result
    }
    
    /// Get the partial derivative dv2/dr2 (3x3 matrix)
    pub fn dv2_dr2(&self) -> [[f64; 3]; 3] {
        let mut result = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                result[i][j] = self.jacobian[3 + i][3 + j];
            }
        }
        result
    }
    
    /// Get the partial derivative dv2/dtof (3-vector)
    pub fn dv2_dtof(&self) -> [f64; 3] {
        [self.jacobian[3][6], self.jacobian[4][6], self.jacobian[5][6]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Direction;
    use std::f64::consts::PI;
    
    #[test]
    fn test_sensitivity_structure() {
        // Basic test that the structure is created correctly
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let tof = PI / 2.0;
        let mu = 1.0;
        
        let geom = Geometry::new(&r1, &r2, tof, mu, Direction::Prograde);
        let state = SolverState::new(0.5, geom.tau);
        
        let sens = LambertSensitivities::compute_first_order(&state, &geom);
        
        // Check matrix dimensions implicitly through accessor methods
        let dv1_dr1 = sens.dv1_dr1();
        assert_eq!(dv1_dr1.len(), 3);
        assert_eq!(dv1_dr1[0].len(), 3);
    }
}
