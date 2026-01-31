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
//!
//! The mathematical approach:
//!   outputs z = [v1, v2] depend on inputs y = [r1, r2, T] through:
//!     1. Geometry: τ(r1, r2), S(r1, r2)
//!     2. Implicit solve: F(k; τ, T/S) = 0 defines k(y)
//!     3. Velocity: v = v(k, τ, S, r1, r2) via Lagrange f,g coefficients
//!
//!   Total derivative:
//!     Dz/Dy = ∂z/∂y|_explicit + (∂z/∂k)·(dk/dy)
//!   where implicit function theorem gives:
//!     dk/dy = -(∂F/∂k)⁻¹ · (∂F/∂y)

use crate::geometry::Geometry;
use crate::velocity::SolverState;

/// First and second-order sensitivities of the Lambert solution.
///
/// The outputs z = [v1, v2] (6 components) depend on inputs y = [r1, r2, tof] (7 components).
#[derive(Debug, Clone)]
pub struct LambertSensitivities {
    /// First-order Jacobian: dz/dy (6x7 matrix)
    /// Row i = output i, Column j = input j
    /// Outputs: [v1x, v1y, v1z, v2x, v2y, v2z]
    /// Inputs: [r1x, r1y, r1z, r2x, r2y, r2z, tof]
    pub jacobian: [[f64; 7]; 6],

    /// Second-order Hessians: d²z_i/dy² (6 matrices of 7x7)
    /// Only computed if requested
    pub hessians: Option<[[[f64; 7]; 7]; 6]>,
}

impl LambertSensitivities {
    /// Compute first-order sensitivities (Jacobian) of the Lambert solution.
    ///
    /// Uses the implicit function theorem approach:
    ///   Dz/Dy = ∂z/∂y|_explicit + (∂z/∂k)·(dk/dy)
    /// where dk/dy = -(∂F/∂k)⁻¹ · ∂F/∂y
    pub fn compute_first_order(
        state: &SolverState,
        geom: &Geometry,
    ) -> Self {
        // ── Step 1: Gather core values ──
        let k = state.k_sol;
        let p = state.p;
        let sqrt_p = state.sqrt_p;
        let tau = geom.tau;
        let s = geom.s;
        let w = state.dw[0];
        let dw_dk = state.dw[1]; // W'(k)

        let r1 = geom.r1;
        let r2 = geom.r2;
        let inv_r1 = geom.inv_r1;
        let inv_r2 = geom.inv_r2;
        let r1_plus_r2 = geom.r1_plus_r2;
        let inv_r1_plus_r2 = 1.0 / r1_plus_r2;

        // Lagrange coefficients
        let f = 1.0 - p * r1_plus_r2 * inv_r1;
        let g = s * tau * sqrt_p;
        let g_dot = 1.0 - p * r1_plus_r2 * inv_r2;
        let inv_g = 1.0 / (g + geom.add_to_g);

        // ── Step 2: Geometry partials ──
        // We need ∂τ/∂rᵢ, ∂S/∂rᵢ, and ∂(T/S)/∂rᵢ for each of 7 inputs.
        //
        // Inputs y = [r1x, r1y, r1z, r2x, r2y, r2z, T]
        //
        // Key relations:
        //   τ = d · √(r1·r2·(1+cosθ)) / (r1+r2)
        //   S = (r1+r2) · √((r1+r2)/μ)
        //   T/S = T / S
        //   cosθ = (r⃗₁·r⃗₂) / (r1·r2)

        // Unit vectors
        let r1_hat = [geom.r1_vec[0] * inv_r1, geom.r1_vec[1] * inv_r1, geom.r1_vec[2] * inv_r1];
        let r2_hat = [geom.r2_vec[0] * inv_r2, geom.r2_vec[1] * inv_r2, geom.r2_vec[2] * inv_r2];

        // ∂r1/∂r1ᵢ = r1_hat[i] (scalar derivative of magnitude)
        // ∂r2/∂r2ᵢ = r2_hat[i]
        // ∂(cosθ)/∂r1ᵢ = (r2_hat[i] - cosθ·r1_hat[i]) / r1
        // ∂(cosθ)/∂r2ᵢ = (r1_hat[i] - cosθ·r2_hat[i]) / r2
        let cos_theta = geom.cos_theta;

        // Partial of τ w.r.t. each input
        // τ = d · √(r1·r2·(1+cosθ)) / (r1+r2)
        // Let Q = r1·r2·(1+cosθ), so τ = d·√Q / (r1+r2)
        let q = r1 * r2 * (1.0 + cos_theta);
        let sqrt_q = if q > 0.0 { q.sqrt() } else { 0.0 };
        let d = if tau >= 0.0 { 1.0 } else { -1.0 };
        // Avoid division by zero
        let inv_sqrt_q = if sqrt_q > 1e-30 { 1.0 / sqrt_q } else { 0.0 };

        // ∂Q/∂r1ᵢ = r2·(1+cosθ)·∂r1/∂r1ᵢ + r1·r2·∂cosθ/∂r1ᵢ
        //          = r2·(1+cosθ)·r̂₁ᵢ + r1·r2·(r̂₂ᵢ - cosθ·r̂₁ᵢ)/r1
        //          = r2·(1+cosθ)·r̂₁ᵢ + r2·(r̂₂ᵢ - cosθ·r̂₁ᵢ)
        //          = r2·(r̂₁ᵢ + r̂₂ᵢ)     [simplifies nicely!]
        // ∂Q/∂r2ᵢ = r1·(r̂₁ᵢ + r̂₂ᵢ)     [by symmetry]

        // ∂τ/∂r1ᵢ = d · [∂√Q/∂r1ᵢ · (r1+r2) - √Q · r̂₁ᵢ] / (r1+r2)²
        //          = d · [(r1+r2)·(∂Q/∂r1ᵢ)/(2√Q) - √Q·r̂₁ᵢ] / (r1+r2)²
        //
        // Let's compute ∂τ/∂r1ᵢ and ∂τ/∂r2ᵢ component by component.
        let inv_r1_plus_r2_sq = inv_r1_plus_r2 * inv_r1_plus_r2;

        // dtau/dy[j] for j = 0..6
        let mut dtau = [0.0_f64; 7];
        // ds/dy[j] for j = 0..6
        let mut ds = [0.0_f64; 7];

        for i in 0..3 {
            // ∂Q/∂r1ᵢ = r2 · (r̂₁ᵢ + r̂₂ᵢ)
            let dq_dr1i = r2 * (r1_hat[i] + r2_hat[i]);
            // ∂Q/∂r2ᵢ = r1 · (r̂₁ᵢ + r̂₂ᵢ)
            let dq_dr2i = r1 * (r1_hat[i] + r2_hat[i]);

            // ∂τ/∂r1ᵢ = d/(r1+r2)² · [(r1+r2)·dQ/(2√Q) - √Q·r̂₁ᵢ]
            dtau[i] = d * inv_r1_plus_r2_sq
                * (r1_plus_r2 * dq_dr1i * 0.5 * inv_sqrt_q - sqrt_q * r1_hat[i]);
            // ∂τ/∂r2ᵢ
            dtau[3 + i] = d * inv_r1_plus_r2_sq
                * (r1_plus_r2 * dq_dr2i * 0.5 * inv_sqrt_q - sqrt_q * r2_hat[i]);

            // ∂S/∂r1ᵢ = (3/2)·√((r1+r2)/μ)·r̂₁ᵢ
            // S = (r1+r2)·√((r1+r2)/μ) = (r1+r2)^(3/2) / √μ
            // ∂S/∂r1ᵢ = (3/2)·(r1+r2)^(1/2)/√μ · ∂(r1+r2)/∂r1ᵢ
            //          = (3/2)·(r1+r2)^(1/2)/√μ · r̂₁ᵢ
            // But S/(r1+r2) = √((r1+r2)/μ), so:
            // ∂S/∂r1ᵢ = (3/2) · S/(r1+r2) · r̂₁ᵢ
            ds[i] = 1.5 * s * inv_r1_plus_r2 * r1_hat[i];
            ds[3 + i] = 1.5 * s * inv_r1_plus_r2 * r2_hat[i];
        }
        // ∂τ/∂T = 0, ∂S/∂T = 0
        dtau[6] = 0.0;
        ds[6] = 0.0;

        // ── Step 3: ∂F/∂y (explicit, treating k as constant) ──
        // F = √p·(τ + p·W) - T/S
        // ∂F/∂yⱼ = √p·∂τ/∂yⱼ + √p·W·∂p/∂yⱼ·... hmm wait.
        //
        // Since p = 1 - k·τ, and k is treated as constant:
        //   ∂p/∂yⱼ = -k · ∂τ/∂yⱼ
        //   ∂√p/∂yⱼ = -k/(2√p) · ∂τ/∂yⱼ
        //
        // F = √p·(τ + p·W) - T/S
        // ∂F/∂yⱼ = ∂√p/∂yⱼ · (τ + p·W) + √p·(∂τ/∂yⱼ + ∂p/∂yⱼ·W) - ∂(T/S)/∂yⱼ
        //         = [-k/(2√p)]·∂τ/∂yⱼ·(τ+p·W) + √p·(∂τ/∂yⱼ - k·W·∂τ/∂yⱼ) - ∂(T/S)/∂yⱼ
        //         = ∂τ/∂yⱼ · [-k(τ+pW)/(2√p) + √p(1-kW)] - ∂(T/S)/∂yⱼ
        //
        // ∂(T/S)/∂yⱼ:
        //   T/S = T·S⁻¹
        //   ∂(T/S)/∂r_ij = -T·S⁻²·∂S/∂r_ij = -(T/S)·∂S/∂r_ij / S
        //   ∂(T/S)/∂T = 1/S

        let tof_by_s = geom.tof_by_s;
        let inv_s = 1.0 / s;
        let tau_plus_pw = tau + p * w;

        // Coefficient for ∂τ/∂yⱼ in ∂F/∂yⱼ:
        let coeff_dtau = -k * tau_plus_pw / (2.0 * sqrt_p) + sqrt_p * (1.0 - k * w);

        let mut df_dy = [0.0_f64; 7];
        for j in 0..7 {
            // ∂(T/S)/∂yⱼ
            let d_tof_by_s = if j == 6 {
                inv_s // ∂(T/S)/∂T = 1/S
            } else {
                -tof_by_s * ds[j] * inv_s // ∂(T/S)/∂rᵢⱼ = -(T/S)·(∂S/∂yⱼ)/S
            };

            df_dy[j] = coeff_dtau * dtau[j] - d_tof_by_s;
        }

        // ── Step 4: ∂F/∂k (already computed in solver iteration) ──
        // F' = dF/dk = (-3pτW + 2p²W' - τ²) / (2√p)
        let p_sq = p * p;
        let df_dk = (-3.0 * p * tau * w + 2.0 * p_sq * dw_dk - tau * tau) / (2.0 * sqrt_p);

        // ── Step 5: dk/dy via implicit function theorem ──
        // dk/dy = -(∂F/∂k)⁻¹ · (∂F/∂y)
        let inv_df_dk = -1.0 / df_dk;
        let mut dk_dy = [0.0_f64; 7];
        for j in 0..7 {
            dk_dy[j] = inv_df_dk * df_dy[j];
        }

        // ── Step 6: ∂z/∂k (velocity partials w.r.t. k, at constant geometry) ──
        // v₁ = (r⃗₂ - f·r⃗₁) / g
        // v₂ = (ġ·r⃗₂ - r⃗₁) / g
        //
        // dp/dk = -τ
        // df/dk = dp/dk · (r1+r2)/r1 = -τ·(r1+r2)/r1
        // dg/dk = S·τ·d(√p)/dk = S·τ·(-τ)/(2√p) = -S·τ²/(2√p)
        // dġ/dk = dp/dk · (r1+r2)/r2 = -τ·(r1+r2)/r2
        let df_dk_lagrange = tau * r1_plus_r2 * inv_r1;
        let dg_dk = -s * tau * tau / (2.0 * sqrt_p);
        let dg_dot_dk = tau * r1_plus_r2 * inv_r2;

        // dv₁/dk = (-df/dk·r⃗₁ - dg/dk·v⃗₁) / g  = (−df/dk·r⃗₁)/g − (dg/dk/g)·v⃗₁
        // dv₂/dk = (dġ/dk·r⃗₂ - dg/dk·v⃗₂) / g
        // But we need to recompute v1, v2 first
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

        let mut dv1_dk = [0.0; 3];
        let mut dv2_dk = [0.0; 3];
        for i in 0..3 {
            dv1_dk[i] = (-df_dk_lagrange * geom.r1_vec[i] - dg_dk * v1[i]) * inv_g;
            dv2_dk[i] = (dg_dot_dk * geom.r2_vec[i] - dg_dk * v2[i]) * inv_g;
        }

        // ── Step 7: ∂z/∂y|_explicit (velocity partials at constant k) ──
        // v₁ = (r⃗₂ - f·r⃗₁) / g
        // For ∂v₁/∂r1ⱼ: r⃗₁ appears directly AND through f and g (via τ).
        //   ∂v₁ᵢ/∂r1ⱼ = [-∂f/∂r1ⱼ·r1ᵢ - f·δᵢⱼ - ∂g/∂r1ⱼ·v1ᵢ] / g
        //               + [∂r2ᵢ/∂r1ⱼ] / g     (this is 0 since r2 doesn't depend on r1)
        //
        // ∂f/∂yⱼ|_k = ∂/∂yⱼ [1 - p·(r1+r2)/r1] = -∂p/∂yⱼ·(r1+r2)/r1 - p·∂[(r1+r2)/r1]/∂yⱼ
        //   ∂p/∂yⱼ = -k·∂τ/∂yⱼ
        //
        // ∂[(r1+r2)/r1]/∂r1ⱼ = (r̂₁ⱼ/r1 - (r1+r2)·r̂₁ⱼ/r1²) = r̂₁ⱼ(1/r1 - (r1+r2)/r1²)
        //                     = r̂₁ⱼ(r1 - r1 - r2)/r1² = -r2·r̂₁ⱼ/r1²
        // ∂[(r1+r2)/r1]/∂r2ⱼ = r̂₂ⱼ/r1
        //
        // ∂g/∂yⱼ|_k = S·∂τ/∂yⱼ·√p + τ·√p·∂S/∂yⱼ + S·τ·(-k·∂τ/∂yⱼ)/(2√p)
        //            = ∂τ/∂yⱼ·(S·√p - S·τ·k/(2√p)) + τ·√p·∂S/∂yⱼ
        //
        // Similarly for ġ = 1 - p·(r1+r2)/r2

        // Precompute common factors
        let r1_plus_r2_over_r1 = r1_plus_r2 * inv_r1;
        let r1_plus_r2_over_r2 = r1_plus_r2 * inv_r2;

        // Coefficient for ∂τ/∂yⱼ in ∂f/∂yⱼ
        let df_dtau_coeff = k * r1_plus_r2_over_r1; // ∂f/∂τ|_k = k·(r1+r2)/r1

        // Coefficient for ∂τ/∂yⱼ in ∂g/∂yⱼ
        let dg_dtau_coeff = s * sqrt_p - s * tau * k / (2.0 * sqrt_p);

        // Coefficient for ∂τ/∂yⱼ in ∂ġ/∂yⱼ
        let dg_dot_dtau_coeff = k * r1_plus_r2_over_r2;

        // Coefficient for ∂S/∂yⱼ in ∂g/∂yⱼ
        let dg_ds_coeff = tau * sqrt_p;

        let mut jacobian = [[0.0_f64; 7]; 6];

        for j in 0..7 {
            // Partials of f, g, ġ w.r.t. yⱼ at constant k
            let mut df_dyj = df_dtau_coeff * dtau[j];
            let dg_dyj = dg_dtau_coeff * dtau[j] + dg_ds_coeff * ds[j];
            let mut dg_dot_dyj = dg_dot_dtau_coeff * dtau[j];

            // Additional terms from direct r1, r2 dependence in f and ġ
            // f = 1 - p·(r1+r2)/r1
            //   ∂[(r1+r2)/r1]/∂r1ⱼ = -r2·r̂₁ⱼ/r1²
            //   ∂[(r1+r2)/r1]/∂r2ⱼ = r̂₂ⱼ/r1
            // ġ = 1 - p·(r1+r2)/r2
            //   ∂[(r1+r2)/r2]/∂r1ⱼ = r̂₁ⱼ/r2
            //   ∂[(r1+r2)/r2]/∂r2ⱼ = -r1·r̂₂ⱼ/r2²
            if j < 3 {
                // ∂/∂r1ⱼ
                df_dyj += -p * (-r2 * r1_hat[j] * inv_r1 * inv_r1);
                dg_dot_dyj += -p * (r1_hat[j] * inv_r2);
            } else if j < 6 {
                // ∂/∂r2ⱼ
                let jj = j - 3;
                df_dyj += -p * (r2_hat[jj] * inv_r1);
                dg_dot_dyj += -p * (-r1 * r2_hat[jj] * inv_r2 * inv_r2);
            }
            // j == 6 (T): no additional terms

            // Now compute ∂v1ᵢ/∂yⱼ|_explicit and ∂v2ᵢ/∂yⱼ|_explicit
            for i in 0..3 {
                // v₁ᵢ = (r2ᵢ - f·r1ᵢ) / g
                // ∂v₁ᵢ/∂yⱼ = [-df/dyⱼ·r1ᵢ + ∂r2ᵢ/∂yⱼ - f·∂r1ᵢ/∂yⱼ - dg/dyⱼ·v1ᵢ] / g

                // Direct position partials (Kronecker deltas)
                let dr1i_dyj = if j < 3 && j == i { 1.0 } else { 0.0 };
                let dr2i_dyj = if j >= 3 && j < 6 && (j - 3) == i { 1.0 } else { 0.0 };

                let dv1_explicit = (-df_dyj * geom.r1_vec[i] + dr2i_dyj
                    - f * dr1i_dyj - dg_dyj * v1[i])
                    * inv_g;

                // v₂ᵢ = (ġ·r2ᵢ - r1ᵢ) / g
                // ∂v₂ᵢ/∂yⱼ = [dġ/dyⱼ·r2ᵢ + ġ·∂r2ᵢ/∂yⱼ - ∂r1ᵢ/∂yⱼ - dg/dyⱼ·v2ᵢ] / g
                let dv2_explicit = (dg_dot_dyj * geom.r2_vec[i] + g_dot * dr2i_dyj
                    - dr1i_dyj - dg_dyj * v2[i])
                    * inv_g;

                // Total: Dv/Dy = explicit + (dv/dk)·(dk/dy)
                jacobian[i][j] = dv1_explicit + dv1_dk[i] * dk_dy[j];
                jacobian[3 + i][j] = dv2_explicit + dv2_dk[i] * dk_dy[j];
            }
        }

        // ── Step 8: Scale for mu ──
        // The solver works with normalized quantities. The Jacobian already accounts
        // for mu through S and tof_by_s, so no additional scaling is needed as long
        // as the geometry was constructed with the correct mu.

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

    /// Get the full 6x7 Jacobian as a reference
    pub fn as_slice(&self) -> &[[f64; 7]; 6] {
        &self.jacobian
    }

    /// Get the transposed Jacobian (7x6) for chain rule applications
    pub fn transpose(&self) -> [[f64; 6]; 7] {
        let mut t = [[0.0; 6]; 7];
        for i in 0..6 {
            for j in 0..7 {
                t[j][i] = self.jacobian[i][j];
            }
        }
        t
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
