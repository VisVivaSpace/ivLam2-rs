//! Geometry computations for the vercosine Lambert formulation.
//!
//! This module handles the geometric preprocessing of the Lambert problem,
//! computing the key parameters tau, S, and related quantities.

use std::f64::consts::PI;

/// Transfer direction (prograde or retrograde)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Prograde (short way): 0 < theta < pi
    Prograde,
    /// Retrograde (long way): pi < theta < 2*pi
    Retrograde,
}

impl Direction {
    /// Returns +1.0 for prograde, -1.0 for retrograde
    #[inline]
    pub fn sign(self) -> f64 {
        match self {
            Direction::Prograde => 1.0,
            Direction::Retrograde => -1.0,
        }
    }
}

/// Geometric parameters computed from the boundary conditions.
///
/// These are computed once at the start and remain constant during iteration.
#[derive(Debug, Clone)]
pub struct Geometry {
    /// Initial position vector
    pub r1_vec: [f64; 3],
    /// Final position vector  
    pub r2_vec: [f64; 3],
    /// Magnitude of r1
    pub r1: f64,
    /// Magnitude of r2
    pub r2: f64,
    /// 1/r1
    pub inv_r1: f64,
    /// 1/r2
    pub inv_r2: f64,
    /// r1 + r2
    pub r1_plus_r2: f64,
    /// r1 * r2
    pub r1_times_r2: f64,
    /// cos(theta) where theta is the transfer angle
    pub cos_theta: f64,
    /// 1 + cos(theta)
    pub one_plus_cos_theta: f64,
    /// 1 - cos(theta)
    pub one_minus_cos_theta: f64,
    /// The geometry parameter S = sqrt((r1+r2)^3 / mu)
    pub s: f64,
    /// The normalized geometry parameter tau = d * sqrt(r1*r2*(1+cos_theta)) / (r1+r2)
    pub tau: f64,
    /// |tau|
    pub abs_tau: f64,
    /// tau^2
    pub tau_sq: f64,
    /// tau^3
    pub tau_cu: f64,
    /// Time of flight divided by S
    pub tof_by_s: f64,
    /// log(tof/S) - precomputed for huge TOF cases
    pub log_tof_by_s: f64,
    /// Whether this is a huge TOF case requiring log formulation
    pub huge_tof_case: bool,
    /// Small value added to g denominator for near-half-rev cases (normally 0)
    pub add_to_g: f64,
    /// Info about near-n*pi transfers
    pub n_pi_rev_info: NpiRevInfo,
}

/// Information about near-n*pi transfer cases
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NpiRevInfo {
    /// Normal case, no warnings
    Normal,
    /// Warning: close to half-rev, velocities may be degraded
    NearHalfRev,
    /// Error: exactly half-rev detected, velocities are bogus
    ExactHalfRev,
    /// Warning: close to full-rev (2*N*pi), not a singularity unless r1=r2
    NearFullRev,
}

// Threshold constants
const TOF_BY_S_HUGE_BOUNDARY: f64 = 1e4;
const NEAR_N_PI_REV_WARNING: f64 = 3.8e-5;
const ALTERNATE_TAU_THRESHOLD: f64 = 1e-8;
const SQRT_TINY: f64 = 1e-154; // sqrt of f64::MIN_POSITIVE approximately
const FOUR_RT_TINY: f64 = 1e-77;

impl Geometry {
    /// Compute geometry parameters from the problem inputs.
    ///
    /// # Arguments
    /// * `r1_vec` - Initial position vector [x, y, z]
    /// * `r2_vec` - Final position vector [x, y, z]
    /// * `tof` - Time of flight
    /// * `mu` - Gravitational parameter
    /// * `direction` - Transfer direction (prograde or retrograde)
    pub fn new(
        r1_vec: &[f64; 3],
        r2_vec: &[f64; 3],
        tof: f64,
        mu: f64,
        direction: Direction,
    ) -> Self {
        let r1 = (r1_vec[0].powi(2) + r1_vec[1].powi(2) + r1_vec[2].powi(2)).sqrt();
        let r2 = (r2_vec[0].powi(2) + r2_vec[1].powi(2) + r2_vec[2].powi(2)).sqrt();
        
        let inv_r1 = 1.0 / r1;
        let inv_r2 = 1.0 / r2;
        let r1_plus_r2 = r1 + r2;
        let r1_times_r2 = r1 * r2;
        
        // Compute cos(theta) from dot product
        let dot = r1_vec[0] * r2_vec[0] + r1_vec[1] * r2_vec[1] + r1_vec[2] * r2_vec[2];
        let cos_theta = dot * inv_r1 * inv_r2;
        let one_plus_cos_theta = cos_theta + 1.0;
        let one_minus_cos_theta = 1.0 - cos_theta;
        
        // Compute |tau| with precision-saving formula when close to pi
        let abs_tau = if one_plus_cos_theta < ALTERNATE_TAU_THRESHOLD {
            // Alternative formula for precision when theta is close to pi
            // |tau| = sin(theta) / sqrt((1-cos_theta) * r1 * r2) / (r1+r2)
            // Using: |r1 x r2| = r1 * r2 * sin(theta)
            let cross_x = r1_vec[1] * r2_vec[2] - r1_vec[2] * r2_vec[1];
            let cross_y = r1_vec[2] * r2_vec[0] - r1_vec[0] * r2_vec[2];
            let cross_z = r1_vec[0] * r2_vec[1] - r1_vec[1] * r2_vec[0];
            let sin_theta_r1r2 = (cross_x.powi(2) + cross_y.powi(2) + cross_z.powi(2)).sqrt();
            (1.0 / (one_minus_cos_theta * r1_times_r2)).sqrt() / r1_plus_r2 * sin_theta_r1r2
        } else {
            // Standard formula: |tau| = sqrt(r1*r2*(1+cos_theta)) / (r1+r2)
            (r1_times_r2 * one_plus_cos_theta).sqrt() / r1_plus_r2
        };
        
        let tau = direction.sign() * abs_tau;
        let tau_sq = tau * tau;
        let tau_cu = tau_sq * tau;
        
        // S = sqrt((r1+r2)^3 / mu)
        let s = r1_plus_r2 * (r1_plus_r2 / mu).sqrt();
        
        // Normalized time of flight
        let tof_by_s = tof / s;
        
        let huge_tof_case = tof_by_s > TOF_BY_S_HUGE_BOUNDARY;
        let log_tof_by_s = if huge_tof_case { tof_by_s.ln() } else { 0.0 };
        
        // Check for near-n*pi transfer warnings
        let (n_pi_rev_info, add_to_g) = if one_plus_cos_theta < NEAR_N_PI_REV_WARNING {
            if abs_tau < SQRT_TINY {
                // Exact half-rev - velocities will be bogus
                (NpiRevInfo::ExactHalfRev, FOUR_RT_TINY)
            } else {
                // Near half-rev warning
                (NpiRevInfo::NearHalfRev, 0.0)
            }
        } else if one_minus_cos_theta < NEAR_N_PI_REV_WARNING {
            // Near full-rev warning
            (NpiRevInfo::NearFullRev, 0.0)
        } else {
            (NpiRevInfo::Normal, 0.0)
        };
        
        Self {
            r1_vec: *r1_vec,
            r2_vec: *r2_vec,
            r1,
            r2,
            inv_r1,
            inv_r2,
            r1_plus_r2,
            r1_times_r2,
            cos_theta,
            one_plus_cos_theta,
            one_minus_cos_theta,
            s,
            tau,
            abs_tau,
            tau_sq,
            tau_cu,
            tof_by_s,
            log_tof_by_s,
            huge_tof_case,
            add_to_g,
            n_pi_rev_info,
        }
    }
    
    /// Compute the parabolic time of flight (T_p / S) for zero-rev case.
    /// This is used for the initial guess transformation.
    #[inline]
    pub fn parabolic_tof_by_s(&self) -> f64 {
        // For parabolic case, k = sqrt(2) and the TOF formula simplifies
        // T_p/S = sqrt(p) * (tau + p*W) where p = 1 - k*tau and W is the parabolic W value
        // For k = sqrt(2): p = 1 - sqrt(2)*tau
        let sqrt_2 = std::f64::consts::SQRT_2;
        let p = 1.0 - sqrt_2 * self.tau;
        let sqrt_p = p.sqrt();
        
        // W at k = sqrt(2) (parabola) can be computed from series
        // W_parab ≈ 0.471404... (see getD4W in FORTRAN)
        const W_PARAB: f64 = 0.47140452079103168;
        
        sqrt_p * (self.tau + p * W_PARAB)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_geometry_basic() {
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let tof = PI / 2.0;
        let mu = 1.0;
        
        let geom = Geometry::new(&r1, &r2, tof, mu, Direction::Prograde);
        
        // Both vectors have magnitude 1
        assert!((geom.r1 - 1.0).abs() < 1e-15);
        assert!((geom.r2 - 1.0).abs() < 1e-15);
        
        // cos(90°) = 0
        assert!(geom.cos_theta.abs() < 1e-15);
        
        // tau should be positive for prograde
        assert!(geom.tau > 0.0);
    }
    
    #[test]
    fn test_geometry_retrograde() {
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let tof = PI / 2.0;
        let mu = 1.0;
        
        let geom = Geometry::new(&r1, &r2, tof, mu, Direction::Retrograde);
        
        // tau should be negative for retrograde
        assert!(geom.tau < 0.0);
        
        // |tau| should be same as prograde
        let geom_pro = Geometry::new(&r1, &r2, tof, mu, Direction::Prograde);
        assert!((geom.abs_tau - geom_pro.abs_tau).abs() < 1e-15);
    }
}
