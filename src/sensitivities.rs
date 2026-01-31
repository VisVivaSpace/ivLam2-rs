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

    /// Compute first and second-order sensitivities (Jacobian + Hessians).
    ///
    /// The Hessians are 6 symmetric 7×7 matrices: one per velocity component.
    /// H_i[j][l] = ∂²z_i / ∂y_j ∂y_l
    ///
    /// Uses the second-order implicit function theorem:
    ///   d²k/dy_j∂y_l = -(1/F_k)·[d²F/dy² + d²F/dkdy·dk/dy + dk/dy·d²F/dkdy + d²F/dk²·dk/dy⊗dk/dy]
    ///
    /// Then the total Hessian is assembled as:
    ///   H_i[j][l] = d²z_i/dy²  +  d²z_i/dkdy_j · dk/dy_l  +  dk/dy_j · d²z_i/dkdy_l
    ///             + d²z_i/dk² · dk/dy_j · dk/dy_l  +  dz_i/dk · d²k/dy_j∂y_l
    pub fn compute_with_hessians(
        state: &SolverState,
        geom: &Geometry,
    ) -> Self {
        // ── Recompute all first-order intermediates ──
        // (We need them for the second-order computation; compute_first_order
        //  keeps them local, so we rebuild rather than refactoring to avoid
        //  touching the working first-order code.)
        let k = state.k_sol;
        let p = state.p;
        let sqrt_p = state.sqrt_p;
        let tau = geom.tau;
        let s = geom.s;
        let w = state.dw[0];
        let dw_dk = state.dw[1];
        let d2w_dk2 = state.dw[2];

        let r1 = geom.r1;
        let r2 = geom.r2;
        let inv_r1 = geom.inv_r1;
        let inv_r2 = geom.inv_r2;
        let r1_plus_r2 = geom.r1_plus_r2;
        let inv_r1_plus_r2 = 1.0 / r1_plus_r2;

        // Lagrange coefficients
        let f_lag = 1.0 - p * r1_plus_r2 * inv_r1;
        let g_lag = s * tau * sqrt_p;
        let g_dot = 1.0 - p * r1_plus_r2 * inv_r2;
        let inv_g = 1.0 / (g_lag + geom.add_to_g);

        // Unit vectors
        let r1_hat = [geom.r1_vec[0] * inv_r1, geom.r1_vec[1] * inv_r1, geom.r1_vec[2] * inv_r1];
        let r2_hat = [geom.r2_vec[0] * inv_r2, geom.r2_vec[1] * inv_r2, geom.r2_vec[2] * inv_r2];
        let cos_theta = geom.cos_theta;

        let q = r1 * r2 * (1.0 + cos_theta);
        let sqrt_q = if q > 0.0 { q.sqrt() } else { 0.0 };
        let d = if tau >= 0.0 { 1.0 } else { -1.0 };
        let inv_sqrt_q = if sqrt_q > 1e-30 { 1.0 / sqrt_q } else { 0.0 };
        let inv_r1_plus_r2_sq = inv_r1_plus_r2 * inv_r1_plus_r2;

        // ── First-order geometry partials (same as compute_first_order) ──
        let mut dtau = [0.0_f64; 7];
        let mut ds = [0.0_f64; 7];
        let mut dq = [0.0_f64; 7]; // We also need dQ/dy for second derivatives

        for i in 0..3 {
            let dq_dr1i = r2 * (r1_hat[i] + r2_hat[i]);
            let dq_dr2i = r1 * (r1_hat[i] + r2_hat[i]);
            dq[i] = dq_dr1i;
            dq[3 + i] = dq_dr2i;

            dtau[i] = d * inv_r1_plus_r2_sq
                * (r1_plus_r2 * dq_dr1i * 0.5 * inv_sqrt_q - sqrt_q * r1_hat[i]);
            dtau[3 + i] = d * inv_r1_plus_r2_sq
                * (r1_plus_r2 * dq_dr2i * 0.5 * inv_sqrt_q - sqrt_q * r2_hat[i]);

            ds[i] = 1.5 * s * inv_r1_plus_r2 * r1_hat[i];
            ds[3 + i] = 1.5 * s * inv_r1_plus_r2 * r2_hat[i];
        }
        dtau[6] = 0.0;
        ds[6] = 0.0;
        dq[6] = 0.0;

        // First-order F partials
        let tof_by_s = geom.tof_by_s;
        let inv_s = 1.0 / s;
        let p_sq = p * p;
        let tau_plus_pw = tau + p * w;
        let coeff_dtau = -k * tau_plus_pw / (2.0 * sqrt_p) + sqrt_p * (1.0 - k * w);

        let mut df_dy = [0.0_f64; 7];
        let mut d_tof_by_s_dy = [0.0_f64; 7];
        for j in 0..7 {
            d_tof_by_s_dy[j] = if j == 6 {
                inv_s
            } else {
                -tof_by_s * ds[j] * inv_s
            };
            df_dy[j] = coeff_dtau * dtau[j] - d_tof_by_s_dy[j];
        }

        // dF/dk
        let df_dk = (-3.0 * p * tau * w + 2.0 * p_sq * dw_dk - tau * tau) / (2.0 * sqrt_p);
        let inv_df_dk = -1.0 / df_dk;

        // dk/dy (first-order IFT)
        let mut dk_dy = [0.0_f64; 7];
        for j in 0..7 {
            dk_dy[j] = inv_df_dk * df_dy[j];
        }

        // Velocity vectors
        let v1 = [
            (geom.r2_vec[0] - f_lag * geom.r1_vec[0]) * inv_g,
            (geom.r2_vec[1] - f_lag * geom.r1_vec[1]) * inv_g,
            (geom.r2_vec[2] - f_lag * geom.r1_vec[2]) * inv_g,
        ];
        let v2 = [
            (g_dot * geom.r2_vec[0] - geom.r1_vec[0]) * inv_g,
            (g_dot * geom.r2_vec[1] - geom.r1_vec[1]) * inv_g,
            (g_dot * geom.r2_vec[2] - geom.r1_vec[2]) * inv_g,
        ];

        // dz/dk (velocity partials w.r.t. k)
        let r1_plus_r2_over_r1 = r1_plus_r2 * inv_r1;
        let r1_plus_r2_over_r2 = r1_plus_r2 * inv_r2;
        let df_dk_lag = tau * r1_plus_r2_over_r1;  // note: df/dk of Lagrange f = -dp/dk·(r1+r2)/r1 (positive!)
        let dg_dk = -s * tau * tau / (2.0 * sqrt_p);
        let dg_dot_dk = tau * r1_plus_r2_over_r2;

        let mut dz_dk = [0.0_f64; 6];
        for i in 0..3 {
            dz_dk[i] = (-df_dk_lag * geom.r1_vec[i] - dg_dk * v1[i]) * inv_g;
            dz_dk[3 + i] = (dg_dot_dk * geom.r2_vec[i] - dg_dk * v2[i]) * inv_g;
        }

        // ── First-order Lagrange coefficient partials w.r.t. y (at constant k) ──
        let df_dtau_coeff = k * r1_plus_r2_over_r1;
        let dg_dtau_coeff = s * sqrt_p - s * tau * k / (2.0 * sqrt_p);
        let dg_dot_dtau_coeff = k * r1_plus_r2_over_r2;
        let dg_ds_coeff = tau * sqrt_p;

        let mut df_lag_dy = [0.0_f64; 7];
        let mut dg_lag_dy = [0.0_f64; 7];
        let mut dg_dot_dy = [0.0_f64; 7];

        for j in 0..7 {
            df_lag_dy[j] = df_dtau_coeff * dtau[j];
            dg_lag_dy[j] = dg_dtau_coeff * dtau[j] + dg_ds_coeff * ds[j];
            dg_dot_dy[j] = dg_dot_dtau_coeff * dtau[j];

            if j < 3 {
                df_lag_dy[j] += -p * (-r2 * r1_hat[j] * inv_r1 * inv_r1);
                dg_dot_dy[j] += -p * (r1_hat[j] * inv_r2);
            } else if j < 6 {
                let jj = j - 3;
                df_lag_dy[j] += -p * (r2_hat[jj] * inv_r1);
                dg_dot_dy[j] += -p * (-r1 * r2_hat[jj] * inv_r2 * inv_r2);
            }
        }

        // ── Compute first-order Jacobian (same formula as compute_first_order) ──
        let mut jacobian = [[0.0_f64; 7]; 6];
        // Also store dz/dy_explicit for second-order use
        let mut dz_dy_explicit = [[0.0_f64; 7]; 6];
        for j in 0..7 {
            for i in 0..3 {
                let dr1i_dyj = if j < 3 && j == i { 1.0 } else { 0.0 };
                let dr2i_dyj = if j >= 3 && j < 6 && (j - 3) == i { 1.0 } else { 0.0 };

                let dv1_explicit = (-df_lag_dy[j] * geom.r1_vec[i] + dr2i_dyj
                    - f_lag * dr1i_dyj - dg_lag_dy[j] * v1[i])
                    * inv_g;
                let dv2_explicit = (dg_dot_dy[j] * geom.r2_vec[i] + g_dot * dr2i_dyj
                    - dr1i_dyj - dg_lag_dy[j] * v2[i])
                    * inv_g;

                dz_dy_explicit[i][j] = dv1_explicit;
                dz_dy_explicit[3 + i][j] = dv2_explicit;
                jacobian[i][j] = dv1_explicit + dz_dk[i] * dk_dy[j];
                jacobian[3 + i][j] = dv2_explicit + dz_dk[3 + i] * dk_dy[j];
            }
        }

        // ════════════════════════════════════════════════════════════════
        // ── SECOND-ORDER COMPUTATION ──
        // ════════════════════════════════════════════════════════════════

        // ── 1. Second-order geometry: d²τ/dy_j·dy_l and d²S/dy_j·dy_l ──
        // These are 7×7 symmetric matrices, but only the position-position
        // blocks (0..6 × 0..6) are non-zero (TOF derivatives of τ, S are zero).
        let mut d2tau = [[0.0_f64; 7]; 7];
        let mut d2s = [[0.0_f64; 7]; 7];

        // We need d²Q/dy_j·dy_l first, then chain through τ = d·√Q/(r1+r2).
        // Also need d²(r1+r2)/dy_j·dy_l and related terms.
        //
        // Key second derivatives of scalars:
        //   d²r1/dr1_i·dr1_j = (δ_ij - r̂₁_i·r̂₁_j) / r1
        //   d²cosθ/dr1_i·dr1_j = see below
        //   d²Q/dy_j·dy_l = r1·r2·d²cosθ + cross terms

        // Precompute P_ij = (δ_ij - r̂_i·r̂_j)/r for each position vector
        // These are the second derivatives of the magnitude function.

        // ∂²Q/∂r1_i∂r1_j:
        //   Q = r1·r2·(1+cosθ)
        //   dQ/dr1_i = r2·(r̂₁_i + r̂₂_i)  [from first order]
        //
        //   d²Q/∂r1_i∂r1_j = r2 · d(r̂₁_i + r̂₂_i)/dr1_j
        //     dr̂₁_i/dr1_j = (δ_ij - r̂₁_i·r̂₁_j)/r1
        //     dr̂₂_i/dr1_j = 0  (r̂₂ doesn't depend on r1)
        //   So: d²Q/∂r1_i∂r1_j = r2·(δ_ij - r̂₁_i·r̂₁_j)/r1
        //
        // ∂²Q/∂r2_i∂r2_j = r1·(δ_ij - r̂₂_i·r̂₂_j)/r2  [by symmetry]
        //
        // ∂²Q/∂r1_i∂r2_j:
        //   d(r2·(r̂₁_i + r̂₂_i))/dr2_j = r̂₂_j·(r̂₁_i + r̂₂_i) + r2·(δ_ij - r̂₂_i·r̂₂_j)/r2
        //                                 = r̂₂_j·(r̂₁_i + r̂₂_i) + (δ_ij - r̂₂_i·r̂₂_j)
        //   Alternatively: d(r1·(r̂₁_i + r̂₂_i))/dr1_j by symmetry via Q.
        //   Let's compute from dQ/dr2_j = r1·(r̂₁_j + r̂₂_j), differentiate w.r.t. r1_i:
        //     d²Q/dr1_i·dr2_j = d[r1·(r̂₁_j + r̂₂_j)]/dr1_i
        //                     = r̂₁_i·(r̂₁_j + r̂₂_j) + r1·(δ_ij - r̂₁_i·r̂₁_j)/r1
        //                     = r̂₁_i·(r̂₁_j + r̂₂_j) + δ_ij - r̂₁_i·r̂₁_j
        //                     = r̂₁_i·r̂₂_j + δ_ij

        let mut d2q = [[0.0_f64; 7]; 7];
        for i in 0..3 {
            for j in 0..3 {
                let delta_ij = if i == j { 1.0 } else { 0.0 };
                // d²Q/∂r1_i∂r1_j
                d2q[i][j] = r2 * (delta_ij - r1_hat[i] * r1_hat[j]) / r1;
                // d²Q/∂r2_i∂r2_j
                d2q[3 + i][3 + j] = r1 * (delta_ij - r2_hat[i] * r2_hat[j]) / r2;
                // d²Q/∂r1_i∂r2_j
                d2q[i][3 + j] = r1_hat[i] * r2_hat[j] + delta_ij;
                // d²Q/∂r2_j∂r1_i (symmetric)
                d2q[3 + j][i] = d2q[i][3 + j];
            }
        }

        // Now compute d²τ/dy_j·dy_l from τ = d·√Q/(r1+r2).
        // Let R = r1+r2.
        // τ = d·√Q/R
        //
        // dτ/dy_j = d · [dQ_j/(2√Q·R) - √Q·dR_j/R²]
        //         = d · [(R·dQ_j - 2Q·dR_j) / (2√Q·R²)]  ... let's use the quotient rule directly.
        //
        // d²τ/dy_j·dy_l = d · d/dy_l [ dQ_j/(2√Q·R) - √Q·dR_j/R² ]
        //
        // Term A = dQ_j/(2√Q·R):
        //   d/dy_l[A] = d²Q_jl/(2√Q·R) + dQ_j·d/dy_l[1/(2√Q·R)]
        //             = d²Q_jl/(2√Q·R) - dQ_j·[dQ_l/(4Q^(3/2)·R) + dR_l/(2√Q·R²)]
        //
        // Term B = √Q·dR_j/R²:
        //   d/dy_l[B] = dQ_l/(2√Q)·dR_j/R² + √Q·d²R_jl/R² - 2√Q·dR_j·dR_l/R³
        //
        // dR/dr1_i = r̂₁_i, dR/dr2_i = r̂₂_i, dR/dT = 0
        // d²R/dr1_i·dr1_j = (δ_ij - r̂₁_i·r̂₁_j)/r1
        // d²R/dr2_i·dr2_j = (δ_ij - r̂₂_i·r̂₂_j)/r2
        // d²R/dr1_i·dr2_j = 0

        let mut dr = [0.0_f64; 7]; // dR/dy
        for i in 0..3 {
            dr[i] = r1_hat[i];
            dr[3 + i] = r2_hat[i];
        }

        let inv_2sqrt_q = 0.5 * inv_sqrt_q;
        let inv_q_32 = inv_sqrt_q * inv_sqrt_q * inv_sqrt_q; // 1/Q^(3/2)
        let inv_r = inv_r1_plus_r2;
        let inv_r_sq = inv_r1_plus_r2_sq;
        let inv_r_cu = inv_r_sq * inv_r;

        for j in 0..7 {
            for l in j..7 {
                // d²R/dy_j·dy_l
                let d2r_jl = if j < 3 && l < 3 {
                    let delta = if j == l { 1.0 } else { 0.0 };
                    (delta - r1_hat[j] * r1_hat[l]) / r1
                } else if j >= 3 && j < 6 && l >= 3 && l < 6 {
                    let jj = j - 3;
                    let ll = l - 3;
                    let delta = if jj == ll { 1.0 } else { 0.0 };
                    (delta - r2_hat[jj] * r2_hat[ll]) / r2
                } else {
                    0.0 // cross r1-r2 or TOF terms
                };

                // d²τ/dy_j·dy_l using product/quotient rule on τ = d·√Q/R
                //
                // Let u = √Q, v = 1/R
                // τ = d·u·v
                // d²τ = d·(d²u·v + du_j·dv_l + du_l·dv_j + u·d²v)
                //
                // du/dy = dQ/(2√Q)
                // d²u/dy_j·dy_l = d²Q_jl/(2√Q) - dQ_j·dQ_l/(4Q^(3/2))
                //
                // dv/dy = -dR/R²
                // d²v/dy_j·dy_l = -d²R/R² + 2·dR_j·dR_l/R³

                let du_j = dq[j] * inv_2sqrt_q;
                let du_l = dq[l] * inv_2sqrt_q;
                let d2u_jl = d2q[j][l] * inv_2sqrt_q
                    - dq[j] * dq[l] * 0.25 * inv_q_32;

                let dv_j = -dr[j] * inv_r_sq;
                let dv_l = -dr[l] * inv_r_sq;
                let d2v_jl = -d2r_jl * inv_r_sq + 2.0 * dr[j] * dr[l] * inv_r_cu;

                d2tau[j][l] = d * (d2u_jl / inv_r + du_j * dv_l + du_l * dv_j
                    + sqrt_q * d2v_jl);
                // Wait, let me redo this more carefully.
                // τ = d · u · (1/R) where u = √Q
                // So τ = d · √Q / R = d · √Q · inv_R
                //
                // d²(u/R)/dy_j dy_l = d²u_jl · (1/R) + du_j · d(1/R)/dy_l + du_l · d(1/R)/dy_j + u · d²(1/R)_jl
                //                    = d²u_jl / R + du_j · (-dR_l/R²) + du_l · (-dR_j/R²) + u · (-d²R_jl/R² + 2·dR_j·dR_l/R³)

                d2tau[j][l] = d * (
                    d2u_jl * inv_r
                    + du_j * (-dr[l] * inv_r_sq)
                    + du_l * (-dr[j] * inv_r_sq)
                    + sqrt_q * (-d2r_jl * inv_r_sq + 2.0 * dr[j] * dr[l] * inv_r_cu)
                );

                d2tau[l][j] = d2tau[j][l]; // symmetric
            }
        }

        // d²S/dy_j·dy_l
        // S = (r1+r2)^(3/2) / √μ  = R^(3/2) / √μ
        // dS/dy_j = (3/2)·R^(1/2)/√μ · dR_j = (3/2)·(S/R)·dR_j
        //
        // d²S/dy_j·dy_l = (3/2)/√μ · [R^(1/2)·d²R_jl + (1/2)·R^(-1/2)·dR_j·dR_l]
        //               = (3/2)·(S/R)·d²R_jl + (3/4)·S/R²·dR_j·dR_l
        let s_over_r = s * inv_r;
        let s_over_r_sq = s * inv_r_sq;

        for j in 0..7 {
            for l in j..7 {
                let d2r_jl = if j < 3 && l < 3 {
                    let delta = if j == l { 1.0 } else { 0.0 };
                    (delta - r1_hat[j] * r1_hat[l]) / r1
                } else if j >= 3 && j < 6 && l >= 3 && l < 6 {
                    let jj = j - 3;
                    let ll = l - 3;
                    let delta = if jj == ll { 1.0 } else { 0.0 };
                    (delta - r2_hat[jj] * r2_hat[ll]) / r2
                } else {
                    0.0
                };

                d2s[j][l] = 1.5 * s_over_r * d2r_jl + 0.75 * s_over_r_sq * dr[j] * dr[l];
                d2s[l][j] = d2s[j][l];
            }
        }

        // ── 2. Second-order F partials ──

        // d²F/dk²  (scalar, already computed in solver as df2)
        // From solver.rs line 464:
        // d²F/dk² = (3pτ²W + 4p³W'' - 12p²τW' - τ³) / (4p^(3/2))
        let tau_sq = tau * tau;
        let tau_cu = tau_sq * tau;
        let p_cu = p_sq * p;
        let inv_p_32 = 1.0 / (sqrt_p * p);
        let d2f_dk2 = (3.0 * p * tau_sq * w + 4.0 * p_cu * d2w_dk2
            - 12.0 * p_sq * tau * dw_dk - tau_cu) * 0.25 * inv_p_32;

        // d²F/dk·dy_j  (7-vector)
        // F_k = (-3pτW + 2p²W' - τ²) / (2√p)
        //
        // Differentiate F_k w.r.t. y_j at constant k.
        // Since p = 1 - k·τ, we have dp/dy_j = -k·dτ/dy_j.
        //
        // F_k = (-3pτW + 2p²W' - τ²) / (2√p)
        //
        // Let N = -3pτW + 2p²W' - τ²    and  D = 2√p
        // dF_k/dy_j = (dN/dy_j · D - N · dD/dy_j) / D²
        //
        // dN/dy_j = -3(dp·τ·W + p·dτ·W) + 2·2p·dp·W' - 2τ·dτ
        //         = -3(-k·dτ·τ·W + p·dτ·W) + 4p·(-k·dτ)·W' - 2τ·dτ
        //         = dτ·(3kτW - 3pW - 4pkW' - 2τ)
        //
        // dD/dy_j = 2·(-k/(2√p))·dτ = -k·dτ/√p
        //
        // So: dF_k/dy_j = [dτ·(3kτW - 3pW - 4pkW' - 2τ) · 2√p
        //                   - (-3pτW + 2p²W' - τ²) · (-k·dτ/√p)] / (4p)
        //
        // Simplify: factor out dτ:
        // dF_k/dy_j = dτ · [(3kτW - 3pW - 4pkW' - 2τ)/(2√p)
        //              + k·(-3pτW + 2p²W' - τ²)/(2p√p) · (1/2)]
        //
        // Actually let's just compute it directly, it's cleaner:
        let n_fk = -3.0 * p * tau * w + 2.0 * p_sq * dw_dk - tau_sq;

        let mut d2f_dkdy = [0.0_f64; 7];
        for j in 0..7 {
            let dp_dyj = -k * dtau[j];
            let dn_dyj = -3.0 * (dp_dyj * tau * w + p * dtau[j] * w)
                + 4.0 * p * dp_dyj * dw_dk
                - 2.0 * tau * dtau[j];
            let dd_dyj = dp_dyj / sqrt_p;  // d(2√p)/dy = 2·dp/(2√p) = dp/√p

            d2f_dkdy[j] = (dn_dyj * 2.0 * sqrt_p - n_fk * dd_dyj) / (4.0 * p);
        }

        // d²F/dy_j·dy_l  (7×7 symmetric matrix)
        // F = √p·(τ + p·W) - T/S
        // We already know dF/dy_j = coeff_dtau · dτ_j - d(T/S)_j
        //
        // where coeff_dtau = -k·(τ+pW)/(2√p) + √p·(1-kW)
        //
        // d²F/dy_j·dy_l = d(coeff_dtau)/dy_l · dτ_j + coeff_dtau · d²τ_jl - d²(T/S)_jl
        //
        // d(coeff_dtau)/dy_l: coeff_dtau depends on y only through τ (since p = 1-kτ).
        //   Let C = -k(τ+pW)/(2√p) + √p(1-kW)
        //   dC/dy_l = dC/dτ · dτ_l
        //
        //   dC/dτ = d/dτ[-k(τ+pW)/(2√p) + √p(1-kW)]
        //     with dp/dτ = -k, d(√p)/dτ = -k/(2√p)
        //
        //   First term: -k·d[(τ+pW)/(2√p)]/dτ
        //     = -k·[(1-kW)·2√p - (τ+pW)·(-k/√p)] / (4p)
        //     = -k·[2√p(1-kW) + k(τ+pW)/√p] / (4p)
        //
        //   Second term: d[√p(1-kW)]/dτ
        //     = (-k/(2√p))·(1-kW) + √p·0  = -k(1-kW)/(2√p)
        //     (W doesn't depend on τ at constant k)
        //
        //   dC/dτ = -k·[2√p(1-kW) + k(τ+pW)/√p]/(4p) - k(1-kW)/(2√p)
        //
        // Let me collect terms more carefully:
        let one_minus_kw = 1.0 - k * w;
        let inv_sqrt_p = 1.0 / sqrt_p;

        let dc_dtau = -k * (2.0 * sqrt_p * one_minus_kw + k * tau_plus_pw * inv_sqrt_p) / (4.0 * p)
            - k * one_minus_kw * inv_sqrt_p * 0.5;

        // d²(T/S)/dy_j·dy_l:
        // d(T/S)/dy_j = -(T/S)·dS_j/S for j<6, = 1/S for j=6
        //
        // d²(T/S)/dy_j·dy_l:
        //   For j<6, l<6: d/dy_l [-(T/S)·dS_j/S]
        //     = -(d(T/S)/dy_l)·dS_j/S - (T/S)·(d²S_jl/S - dS_j·dS_l/S²)
        //     = (T/S)·dS_l·dS_j/S² - (T/S)·d²S_jl/S + (T/S)·dS_j·dS_l/S²
        //     = (T/S)·[2·dS_j·dS_l/S² - d²S_jl/S]
        //   For j=6, l<6: d/dy_l [1/S] = -dS_l/S²
        //   For j<6, l=6: d/dy_j [1/S] · ... actually:
        //     d(T/S)/dy_j = -(T/S)·dS_j/S
        //     d²(T/S)/dy_j·dy_6 = d/dT[-(T/S)·dS_j/S] = -dS_j/S²  (since d(T/S)/dT = 1/S)
        //   For j=6, l=6: d/dT[1/S] = 0

        let mut d2_tof_by_s = [[0.0_f64; 7]; 7];
        let inv_s_sq = inv_s * inv_s;
        for j in 0..7 {
            for l in j..7 {
                if j == 6 && l == 6 {
                    d2_tof_by_s[j][l] = 0.0;
                } else if j == 6 || l == 6 {
                    // One is TOF, one is position
                    let pos = if j == 6 { l } else { j };
                    d2_tof_by_s[j][l] = -ds[pos] * inv_s_sq;
                } else {
                    d2_tof_by_s[j][l] = tof_by_s * (2.0 * ds[j] * ds[l] * inv_s_sq - d2s[j][l] * inv_s);
                }
                d2_tof_by_s[l][j] = d2_tof_by_s[j][l];
            }
        }

        let mut d2f_dydy = [[0.0_f64; 7]; 7];
        for j in 0..7 {
            for l in j..7 {
                d2f_dydy[j][l] = dc_dtau * dtau[l] * dtau[j]
                    + coeff_dtau * d2tau[j][l]
                    - d2_tof_by_s[j][l];
                d2f_dydy[l][j] = d2f_dydy[j][l];
            }
        }

        // ── 3. Second-order IFT: d²k/dy_j·dy_l ──
        // dk/dy = -(1/F_k)·F_y
        // d²k/dy_j·dy_l = -(1/F_k)·[d²F/dy_j·dy_l
        //   + d²F/dk·dy_j · dk/dy_l + dk/dy_j · d²F/dk·dy_l
        //   + d²F/dk² · dk/dy_j · dk/dy_l]
        let neg_inv_fk = inv_df_dk; // = -1/F_k
        let mut d2k_dy = [[0.0_f64; 7]; 7];
        for j in 0..7 {
            for l in j..7 {
                d2k_dy[j][l] = neg_inv_fk * (
                    d2f_dydy[j][l]
                    + d2f_dkdy[j] * dk_dy[l]
                    + dk_dy[j] * d2f_dkdy[l]
                    + d2f_dk2 * dk_dy[j] * dk_dy[l]
                );
                d2k_dy[l][j] = d2k_dy[j][l];
            }
        }

        // ── 4. Second-order velocity: d²z/dk², d²z/dkdy, d²z/dydy ──

        // d²(Lagrange coefficients)/dk²:
        // f = 1 - p·(r1+r2)/r1,  dp/dk = -τ,  d²p/dk² = 0
        // df/dk = -dp/dk · (r1+r2)/r1 = τ·(r1+r2)/r1
        // d²f/dk² = 0  (since dp/dk = -τ is constant w.r.t. k)
        let d2f_dk2_lag = 0.0;

        // g = S·τ·√p
        // dg/dk = S·τ·dp/(2√p) = -S·τ²/(2√p)
        // d²g/dk² = S·τ²·τ/(4p^(3/2)) · (-1) ... let me redo:
        //   dg/dk = S·τ·(-τ)/(2√p) = -Sτ²/(2√p)
        //   d²g/dk² = d/dk[-Sτ²/(2√p)]
        //           = -Sτ² · d(1/(2√p))/dk
        //           = -Sτ² · (-1) · dp/dk / (4p^(3/2))
        //           = -Sτ² · (-1) · (-τ) / (4p^(3/2))
        //           = -Sτ³ / (4p^(3/2))
        let d2g_dk2 = -s * tau_cu / (4.0 * sqrt_p * p);

        // ġ = 1 - p·(r1+r2)/r2
        // dġ/dk = -dp/dk · (r1+r2)/r2 = τ·(r1+r2)/r2
        // d²ġ/dk² = 0
        let d2g_dot_dk2 = 0.0;

        // d²v1_i/dk² = d/dk[(−df/dk·r1_i − dg/dk·v1_i) / g]
        // This requires d²f/dk², d²g/dk², and the chain through 1/g.
        // v1_i = (r2_i - f·r1_i) / g
        // dv1_i/dk = (-df/dk·r1_i - dg/dk·v1_i) / g
        //
        // d²v1_i/dk² = [(-d²f/dk²·r1_i - d²g/dk²·v1_i - dg/dk·dv1_i/dk) · g
        //               - (-df/dk·r1_i - dg/dk·v1_i) · dg/dk] / g²
        //            = [-d²f/dk²·r1_i - d²g/dk²·v1_i - 2·dg/dk·dv1_i/dk] / g
        //   (since the last term is -(dg/dk/g)·dv1_i/dk from the quotient rule)
        //   Actually let me be more careful:
        //   dv1_i/dk = N_i / g where N_i = -df/dk·r1_i - dg/dk·v1_i·g ... no.
        //
        // Let me use the approach: v1_i = (r2_i - f·r1_i)/g, so
        //   dv1/dk = (-df/dk·r1_i)/g - (dg/dk/g)·v1_i  ... = (-df/dk·r1_i - dg/dk·v1_i)/g
        //
        //   d²v1_i/dk² = [(-d²f/dk²·r1_i)/g - (d²g/dk²·v1_i + dg/dk·dv1_i/dk)/g
        //                 - (dg/dk/g)·dv1_i/dk]
        //              = (-d²f/dk²·r1_i - d²g/dk²·v1_i)/g - 2·(dg/dk/g)·dv1_i/dk
        //              = (-d²f/dk²·r1_i - d²g/dk²·v1_i - 2·dg/dk·dz_dk_i)/g
        //   Wait, that's wrong too. Let me redo carefully with dv1/dk = dz_dk[i]:
        //
        //   d(dv1_i/dk)/dk = d/dk[(-df_dk_lag·r1_i - dg_dk·v1_i) / g]
        //     = (-d2f_dk2_lag·r1_i - d2g_dk2·v1_i - dg_dk·dz_dk[i]) / g
        //       - (dg_dk / g) · dz_dk[i]
        //     = (-d2f_dk2_lag·r1_i - d2g_dk2·v1_i - 2·dg_dk·dz_dk[i]) / g

        let mut d2z_dk2 = [0.0_f64; 6];
        for i in 0..3 {
            d2z_dk2[i] = (-d2f_dk2_lag * geom.r1_vec[i] - d2g_dk2 * v1[i]
                - 2.0 * dg_dk * dz_dk[i]) * inv_g;
            d2z_dk2[3 + i] = (-d2g_dot_dk2 * geom.r2_vec[i] - d2g_dk2 * v2[i]
                - 2.0 * dg_dk * dz_dk[3 + i]) * inv_g;
        }

        // d²z_i/dk·dy_j  (6×7)
        // dz_i/dk = (-df_dk_lag·r1_i - dg_dk·v1_i) / g   for v1
        //
        // d/dy_j[dz_i/dk] = d/dy_j [(-df_dk_lag·r1_i - dg_dk·v1_i) / g]
        //
        // df_dk_lag = τ·(r1+r2)/r1, so d(df_dk_lag)/dy_j involves dtau, d(r1+r2)/r1
        // dg_dk = -Sτ²/(2√p), so d(dg_dk)/dy_j involves dτ, dS, dp
        //
        // Let me compute d(df_dk_lag)/dy_j and d(dg_dk)/dy_j separately.

        // d(df_dk_lag)/dy_j = d(τ·(r1+r2)/r1)/dy_j
        //   = dτ_j·(r1+r2)/r1 + τ·d((r1+r2)/r1)/dy_j
        //
        // d((r1+r2)/r1)/dr1_j = -r2·r̂₁_j/r1²
        // d((r1+r2)/r1)/dr2_j = r̂₂_j/r1

        let mut d_df_dk_lag_dy = [0.0_f64; 7];
        for j in 0..7 {
            d_df_dk_lag_dy[j] = dtau[j] * r1_plus_r2_over_r1;
            if j < 3 {
                d_df_dk_lag_dy[j] += tau * (-r2 * r1_hat[j] * inv_r1 * inv_r1);
            } else if j < 6 {
                let jj = j - 3;
                d_df_dk_lag_dy[j] += tau * (r2_hat[jj] * inv_r1);
            }
        }

        // d(dg_dk)/dy_j = d(-Sτ²/(2√p))/dy_j
        //   = -(dS_j·τ²/(2√p) + S·2τ·dτ_j/(2√p) + S·τ²·d(1/(2√p))/dy_j)
        //   = -(dS_j·τ²/(2√p) + S·τ·dτ_j/√p + S·τ²·k·dτ_j/(4p^(3/2)))
        //   [since d(1/(2√p))/dy_j = -dp_j/(4p^(3/2)) = k·dτ_j/(4p^(3/2))]

        let mut d_dg_dk_dy = [0.0_f64; 7];
        for j in 0..7 {
            d_dg_dk_dy[j] = -(ds[j] * tau_sq / (2.0 * sqrt_p)
                + s * tau * dtau[j] / sqrt_p
                + s * tau_sq * k * dtau[j] / (4.0 * sqrt_p * p));
        }

        // d(dg_dot_dk)/dy_j = d(τ·(r1+r2)/r2)/dy_j
        let mut d_dg_dot_dk_dy = [0.0_f64; 7];
        for j in 0..7 {
            d_dg_dot_dk_dy[j] = dtau[j] * r1_plus_r2_over_r2;
            if j < 3 {
                d_dg_dot_dk_dy[j] += tau * (r1_hat[j] * inv_r2);
            } else if j < 6 {
                let jj = j - 3;
                d_dg_dot_dk_dy[j] += tau * (-r1 * r2_hat[jj] * inv_r2 * inv_r2);
            }
        }

        // Now d²z_i/dk·dy_j:
        // dv1_i/dk = (-df_dk_lag·r1_i - dg_dk·v1_i) / g
        // d/dy_j = [(-d_df_dk_lag_j·r1_i - df_dk_lag·δ_{ij,r1} - d_dg_dk_j·v1_i - dg_dk·dv1_i/dy_j) / g
        //           - (dg_lag_dy_j / g) · dv1_i/dk]
        //        = [(-d_df_dk_lag_j·r1_i - df_dk_lag·δ_{ij,r1} - d_dg_dk_j·v1_i
        //            - dg_dk·dz_dy_explicit[i][j]) / g
        //           - (dg_lag_dy_j / g) · dz_dk[i]]
        //
        // Note: dv1_i/dy_j at constant k = dz_dy_explicit[i][j]
        // And the g in denominator gives the -dg/g·(dv/dk) term.

        let mut d2z_dkdy = [[0.0_f64; 7]; 6];
        for j in 0..7 {
            let dr1i_dyj = [
                if j == 0 { 1.0 } else { 0.0 },
                if j == 1 { 1.0 } else { 0.0 },
                if j == 2 { 1.0 } else { 0.0 },
            ];
            let dr2i_dyj = [
                if j == 3 { 1.0 } else { 0.0 },
                if j == 4 { 1.0 } else { 0.0 },
                if j == 5 { 1.0 } else { 0.0 },
            ];

            for i in 0..3 {
                // v1: dv1_i/dk = (-df_dk_lag·r1_i - dg_dk·v1_i) / g
                let num_dk_v1 = -d_df_dk_lag_dy[j] * geom.r1_vec[i]
                    - df_dk_lag * dr1i_dyj[i]
                    - d_dg_dk_dy[j] * v1[i]
                    - dg_dk * dz_dy_explicit[i][j];
                d2z_dkdy[i][j] = num_dk_v1 * inv_g - dg_lag_dy[j] * inv_g * dz_dk[i];

                // v2: dv2_i/dk = (dg_dot_dk·r2_i - dg_dk·v2_i) / g
                let num_dk_v2 = d_dg_dot_dk_dy[j] * geom.r2_vec[i]
                    + dg_dot_dk * dr2i_dyj[i]
                    - d_dg_dk_dy[j] * v2[i]
                    - dg_dk * dz_dy_explicit[3 + i][j];
                d2z_dkdy[3 + i][j] = num_dk_v2 * inv_g - dg_lag_dy[j] * inv_g * dz_dk[3 + i];
            }
        }

        // d²z_i/dy_j·dy_l  (6×7×7)
        // v1_i = (r2_i - f·r1_i) / g
        // dv1_i/dy_j = (-df_j·r1_i + δ_{r2i,j} - f·δ_{r1i,j} - dg_j·v1_i) / g
        //
        // d²v1_i/dy_j·dy_l = d/dy_l of the above, at constant k.
        //
        // Numerator N = -df_j·r1_i + δ_{r2i,j} - f·δ_{r1i,j} - dg_j·v1_i
        //
        // d/dy_l[N/g] = (dN/dy_l · g - N · dg_l) / g²
        //             = dN/dy_l / g - (N/g) · (dg_l/g)
        //             = dN/dy_l / g - dz_dy_explicit[i][j] · dg_l / g
        //
        // dN/dy_l = -d²f_jl·r1_i - df_j·δ_{r1i,l} + 0 - df_l·δ_{r1i,j}
        //           - d²g_jl·v1_i - dg_j·dv1_i/dy_l
        //
        // But dv1_i/dy_l = dz_dy_explicit[i][l] (at constant k)
        //
        // We need d²f/dy_j·dy_l, d²g/dy_j·dy_l, d²ġ/dy_j·dy_l at constant k.

        // d²f/dy_j·dy_l at constant k:
        // f = 1 - p·(r1+r2)/r1
        // df/dy_j = k·dτ_j·(r1+r2)/r1 - p·d((r1+r2)/r1)/dy_j
        // (df_dtau_coeff · dτ_j + geometric term)
        //
        // d²f/dy_j·dy_l = k·d²τ_jl·(r1+r2)/r1
        //   + k·dτ_j·d((r1+r2)/r1)/dy_l + k·dτ_l·d((r1+r2)/r1)/dy_j
        //   + k·dτ_l·(-k·dτ_j)·... wait, dp/dy_l = -k·dτ_l, so
        //   - dp_l·d((r1+r2)/r1)/dy_j - p·d²((r1+r2)/r1)/dy_j·dy_l
        //
        // Let A = (r1+r2)/r1. Then f = 1 - p·A.
        // df/dy_j = -dp_j·A - p·dA_j = k·dτ_j·A - p·dA_j
        // d²f/dy_j·dy_l = k·d²τ_jl·A + k·dτ_j·dA_l + k·dτ_l·dA_j + k²·dτ_j·dτ_l·... wait no.
        //   dp_j = -k·dτ_j, d²p_jl = -k·d²τ_jl
        //   d²f = -d²p·A - dp_j·dA_l - dp_l·dA_j - p·d²A_jl
        //       = k·d²τ_jl·A + k·dτ_j·dA_l + k·dτ_l·dA_j - p·d²A_jl

        // dA/dy_j where A = (r1+r2)/r1:
        //   dA/dr1_j = -r2·r̂₁_j/r1²
        //   dA/dr2_j = r̂₂_j/r1
        //   dA/dT = 0
        let mut da_dy = [0.0_f64; 7];
        for j in 0..3 {
            da_dy[j] = -r2 * r1_hat[j] * inv_r1 * inv_r1;
            da_dy[3 + j] = r2_hat[j] * inv_r1;
        }

        // d²A/dy_j·dy_l:
        //   d²((r1+r2)/r1)/dr1_j·dr1_l: A = 1 + r2/r1
        //     dA/dr1_j = -r2·r̂₁_j/r1²
        //     d²A/dr1_j·dr1_l = -r2 · d(r̂₁_j/r1²)/dr1_l
        //       = -r2 · [(δ_jl - r̂₁_j·r̂₁_l)/(r1·r1²) - 2·r̂₁_j·r̂₁_l/r1³]
        //       = -r2 · [(δ_jl - 3·r̂₁_j·r̂₁_l)/r1³]
        //       = r2·(3·r̂₁_j·r̂₁_l - δ_jl)/r1³
        //
        //   d²A/dr1_j·dr2_l: d/dr2_l[-r2·r̂₁_j/r1²]
        //     = -r̂₂_l·r̂₁_j/r1²   (only r2 magnitude depends on r2_l)
        //     Wait: d(r2)/dr2_l = r̂₂_l, and r̂₁_j doesn't depend on r2.
        //     So d²A/dr1_j·dr2_l = -r̂₂_l·r̂₁_j/r1²
        //
        //   d²A/dr2_j·dr2_l: A involves r̂₂_j/r1 for the dr2 part.
        //     dA/dr2_j = r̂₂_j/r1
        //     d²A/dr2_j·dr2_l = (δ_jl - r̂₂_j·r̂₂_l)/(r2·r1)
        //
        //   d²A/dr2_j·dr1_l: d/dr1_l[r̂₂_j/r1] = -r̂₁_l·r̂₂_j/r1²

        // Similarly for B = (r1+r2)/r2 (for ġ):
        let mut db_dy = [0.0_f64; 7];
        for j in 0..3 {
            db_dy[j] = r1_hat[j] * inv_r2;
            db_dy[3 + j] = -r1 * r2_hat[j] * inv_r2 * inv_r2;
        }

        // d²g/dy_j·dy_l at constant k:
        // g = S·τ·√p
        // dg/dy_j = dg_dtau_coeff · dτ_j + dg_ds_coeff · dS_j
        //
        // This is getting complex. Let me factor it cleanly:
        // g = S·τ·√p where p = 1-kτ
        //
        // dg/dy = S·√p·dτ + τ·√p·dS + S·τ·(-k·dτ)/(2√p)
        //       = (S·√p - Sτk/(2√p))·dτ + τ√p·dS
        //       = α·dτ + β·dS
        //
        // where α = S√p - Sτk/(2√p) = dg_dtau_coeff
        //       β = τ√p = dg_ds_coeff
        //
        // d²g/dy_j·dy_l = dα/dy_l·dτ_j + α·d²τ_jl + dβ/dy_l·dS_j + β·d²S_jl
        //
        // dα/dy_l = dS_l·(√p - τk/(2√p)) + S·d(√p - τk/(2√p))/dy_l
        //   d(√p)/dy_l = -k·dτ_l/(2√p)
        //   d(τk/(2√p))/dy_l = k·dτ_l/(2√p) + τk·k·dτ_l/(4p√p)  [chain through p]
        //     = k·dτ_l/(2√p) + τk²·dτ_l/(4p√p)
        //
        //   d(√p - τk/(2√p))/dy_l = -k·dτ_l/(2√p) - k·dτ_l/(2√p) - τk²·dτ_l/(4p√p)
        //     = -k·dτ_l/√p - τk²·dτ_l/(4p√p)
        //
        //   dα/dy_l = dS_l·(√p - τk/(2√p)) + S·(-k·dτ_l/√p - τk²·dτ_l/(4p√p))
        //
        // dβ/dy_l = dτ_l·√p + τ·(-k·dτ_l/(2√p)) = dτ_l·(√p - τk/(2√p)) = dτ_l·(α/S)
        //   Wait: β = τ√p, dβ/dy_l = dτ_l·√p + τ·(-k/(2√p))·dτ_l = dτ_l·(√p - τk/(2√p))
        //
        // This is correct. So:
        // For the Hessian assembly, I need d²z/dy_j·dy_l for each velocity component.
        // Rather than computing all intermediate d²f, d²g, d²ġ matrices, I'll compute
        // the Hessian directly in the assembly loop.
        //
        // For each (j, l) pair, compute:
        //   d²f_lag/dy_j·dy_l, d²g_lag/dy_j·dy_l, d²g_dot/dy_j·dy_l
        //   then d²v1_i/dy_j·dy_l and d²v2_i/dy_j·dy_l

        let mut hessians = [[[0.0_f64; 7]; 7]; 6];

        for j in 0..7 {
            for l in j..7 {
                // ── d²f_lag/dy_j·dy_l ──
                // f = 1 - p·A where A = (r1+r2)/r1, dp = -k·dτ
                // d²f = k·d²τ·A + k·dτ_j·dA_l + k·dτ_l·dA_j - p·d²A_jl
                let d2a_jl = compute_d2a_r1_r2(j, l, &r1_hat, &r2_hat, r1, r2, inv_r1, inv_r2);
                let d2f_lag_jl = k * d2tau[j][l] * r1_plus_r2_over_r1
                    + k * dtau[j] * da_dy[l]
                    + k * dtau[l] * da_dy[j]
                    - p * d2a_jl;

                // ── d²g_lag/dy_j·dy_l ──
                // g = S·τ·√p = α·dτ + β·dS integration... we use:
                // d²g = dα_l·dτ_j + α·d²τ_jl + dβ_l·dS_j + β·d²S_jl
                //       (+ symmetric partner from dα_j·dτ_l etc. — but the formula
                //        already handles it because dα and dβ depend on dτ and dS)
                //
                // Actually, g = S·τ·√p, so using product rule directly for 3 factors:
                // dg = dS·τ·√p + S·dτ·√p + S·τ·d(√p)
                // d²g = d²S·τ·√p + dS_j·dτ_l·√p + dS_j·τ·d√p_l
                //     + dS_l·dτ_j·√p + S·d²τ·√p + S·dτ_j·d√p_l
                //     + dS_l·τ·d√p_j + S·dτ_l·d√p_j + S·τ·d²√p
                //
                // d√p/dy = -k·dτ/(2√p)
                // d²√p/dy_j·dy_l = -k·d²τ/(2√p) + (-k·dτ_j)·d[1/(2√p)]/dy_l
                //   = -k·d²τ/(2√p) - k²·dτ_j·dτ_l/(4p√p)
                let dsqrtp_j = -k * dtau[j] / (2.0 * sqrt_p);
                let dsqrtp_l = -k * dtau[l] / (2.0 * sqrt_p);
                let d2sqrtp_jl = -k * d2tau[j][l] / (2.0 * sqrt_p)
                    - k * k * dtau[j] * dtau[l] / (4.0 * p * sqrt_p);

                let d2g_lag_jl = d2s[j][l] * tau * sqrt_p
                    + ds[j] * dtau[l] * sqrt_p + ds[j] * tau * dsqrtp_l
                    + ds[l] * dtau[j] * sqrt_p + s * d2tau[j][l] * sqrt_p + s * dtau[j] * dsqrtp_l
                    + ds[l] * tau * dsqrtp_j + s * dtau[l] * dsqrtp_j + s * tau * d2sqrtp_jl;

                // ── d²g_dot/dy_j·dy_l ──
                // ġ = 1 - p·B where B = (r1+r2)/r2
                // Same structure as f but with B instead of A
                let d2b_jl = compute_d2b_r1_r2(j, l, &r1_hat, &r2_hat, r1, r2, inv_r1, inv_r2);
                let d2g_dot_jl = k * d2tau[j][l] * r1_plus_r2_over_r2
                    + k * dtau[j] * db_dy[l]
                    + k * dtau[l] * db_dy[j]
                    - p * d2b_jl;

                // ── d²v/dy_j·dy_l at constant k ──
                for i in 0..3 {
                    let dr1i_dyj = if j < 3 && j == i { 1.0 } else { 0.0 };
                    let dr1i_dyl = if l < 3 && l == i { 1.0 } else { 0.0 };
                    let dr2i_dyj = if j >= 3 && j < 6 && (j - 3) == i { 1.0 } else { 0.0 };
                    let dr2i_dyl = if l >= 3 && l < 6 && (l - 3) == i { 1.0 } else { 0.0 };

                    // v1_i = (r2_i - f·r1_i) / g
                    // N1 = -df_j·r1_i + δ_{r2,j} - f·δ_{r1,j} - dg_j·v1_i
                    // dv1/dy_j = N1 / g
                    //
                    // dN1/dy_l = -d²f_jl·r1_i - df_j·δ_{r1,l} - df_l·δ_{r1,j}
                    //            - d²g_jl·v1_i - dg_j·dv1_explicit_l
                    //
                    // d²v1/dy_j·dy_l = (dN1_l / g) - (N1/g)·(dg_l/g)
                    //                = (dN1_l - dz_dy_explicit[i][j]·dg_l) / g

                    let dn1_dyl = -d2f_lag_jl * geom.r1_vec[i]
                        - df_lag_dy[j] * dr1i_dyl
                        - df_lag_dy[l] * dr1i_dyj
                        - d2g_lag_jl * v1[i]
                        - dg_lag_dy[j] * dz_dy_explicit[i][l];

                    let d2v1_dydy = (dn1_dyl - dz_dy_explicit[i][j] * dg_lag_dy[l]) * inv_g;

                    // v2_i = (ġ·r2_i - r1_i) / g
                    // N2 = dġ_j·r2_i + ġ·δ_{r2,j} - δ_{r1,j} - dg_j·v2_i
                    let dn2_dyl = d2g_dot_jl * geom.r2_vec[i]
                        + dg_dot_dy[j] * dr2i_dyl
                        + dg_dot_dy[l] * dr2i_dyj
                        - d2g_lag_jl * v2[i]
                        - dg_lag_dy[j] * dz_dy_explicit[3 + i][l];

                    let d2v2_dydy = (dn2_dyl - dz_dy_explicit[3 + i][j] * dg_lag_dy[l]) * inv_g;

                    // ── Full Hessian assembly ──
                    // H_i[j][l] = d²z/dy²  +  d²z/dkdy_j·dk/dy_l  +  dk/dy_j·d²z/dkdy_l
                    //            + d²z/dk²·dk/dy_j·dk/dy_l  +  dz/dk·d²k/dy_j·dy_l
                    let term1_v1 = d2v1_dydy;
                    let term2_v1 = d2z_dkdy[i][j] * dk_dy[l];
                    let term3_v1 = dk_dy[j] * d2z_dkdy[i][l];
                    let term4_v1 = d2z_dk2[i] * dk_dy[j] * dk_dy[l];
                    let term5_v1 = dz_dk[i] * d2k_dy[j][l];

                    hessians[i][j][l] = term1_v1 + term2_v1 + term3_v1 + term4_v1 + term5_v1;
                    hessians[i][l][j] = hessians[i][j][l]; // symmetric

                    let term1_v2 = d2v2_dydy;
                    let term2_v2 = d2z_dkdy[3 + i][j] * dk_dy[l];
                    let term3_v2 = dk_dy[j] * d2z_dkdy[3 + i][l];
                    let term4_v2 = d2z_dk2[3 + i] * dk_dy[j] * dk_dy[l];
                    let term5_v2 = dz_dk[3 + i] * d2k_dy[j][l];

                    hessians[3 + i][j][l] = term1_v2 + term2_v2 + term3_v2 + term4_v2 + term5_v2;
                    hessians[3 + i][l][j] = hessians[3 + i][j][l]; // symmetric
                }
            }
        }

        Self {
            jacobian,
            hessians: Some(hessians),
        }
    }

}

/// Compute d²A/dy_j·dy_l where A = (r1+r2)/r1.
///
/// A = 1 + r2/r1
/// dA/dr1_j = -r2·r̂₁_j/r1²
/// dA/dr2_j = r̂₂_j/r1
///
/// Second derivatives:
///   d²A/dr1_j·dr1_l = r2·(3·r̂₁_j·r̂₁_l - δ_jl)/r1³
///   d²A/dr1_j·dr2_l = -r̂₂_l·r̂₁_j/r1²
///   d²A/dr2_j·dr2_l = (δ_jl - r̂₂_j·r̂₂_l)/(r2·r1)
fn compute_d2a_r1_r2(
    j: usize, l: usize,
    r1_hat: &[f64; 3], r2_hat: &[f64; 3],
    r1: f64, r2: f64, inv_r1: f64, _inv_r2: f64,
) -> f64 {
    if j < 3 && l < 3 {
        // d²A/dr1_j·dr1_l
        let delta = if j == l { 1.0 } else { 0.0 };
        r2 * (3.0 * r1_hat[j] * r1_hat[l] - delta) * inv_r1 * inv_r1 * inv_r1
    } else if j < 3 && l >= 3 && l < 6 {
        // d²A/dr1_j·dr2_l
        -r2_hat[l - 3] * r1_hat[j] * inv_r1 * inv_r1
    } else if j >= 3 && j < 6 && l < 3 {
        // d²A/dr2_j·dr1_l (symmetric)
        -r2_hat[j - 3] * r1_hat[l] * inv_r1 * inv_r1
    } else if j >= 3 && j < 6 && l >= 3 && l < 6 {
        // d²A/dr2_j·dr2_l
        let jj = j - 3;
        let ll = l - 3;
        let delta = if jj == ll { 1.0 } else { 0.0 };
        (delta - r2_hat[jj] * r2_hat[ll]) / (r2 * r1)
    } else {
        0.0 // TOF terms
    }
}

/// Compute d²B/dy_j·dy_l where B = (r1+r2)/r2.
///
/// B = r1/r2 + 1
/// dB/dr1_j = r̂₁_j/r2
/// dB/dr2_j = -r1·r̂₂_j/r2²
///
/// Second derivatives:
///   d²B/dr1_j·dr1_l = (δ_jl - r̂₁_j·r̂₁_l)/(r1·r2)
///   d²B/dr1_j·dr2_l = -r̂₁_j·r̂₂_l/r2²
///   d²B/dr2_j·dr2_l = r1·(3·r̂₂_j·r̂₂_l - δ_jl)/r2³
fn compute_d2b_r1_r2(
    j: usize, l: usize,
    r1_hat: &[f64; 3], r2_hat: &[f64; 3],
    r1: f64, r2: f64, _inv_r1: f64, inv_r2: f64,
) -> f64 {
    if j < 3 && l < 3 {
        // d²B/dr1_j·dr1_l
        let delta = if j == l { 1.0 } else { 0.0 };
        (delta - r1_hat[j] * r1_hat[l]) / (r1 * r2)
    } else if j < 3 && l >= 3 && l < 6 {
        // d²B/dr1_j·dr2_l
        -r1_hat[j] * r2_hat[l - 3] * inv_r2 * inv_r2
    } else if j >= 3 && j < 6 && l < 3 {
        // d²B/dr2_j·dr1_l (symmetric)
        -r1_hat[l] * r2_hat[j - 3] * inv_r2 * inv_r2
    } else if j >= 3 && j < 6 && l >= 3 && l < 6 {
        // d²B/dr2_j·dr2_l
        let jj = j - 3;
        let ll = l - 3;
        let delta = if jj == ll { 1.0 } else { 0.0 };
        r1 * (3.0 * r2_hat[jj] * r2_hat[ll] - delta) * inv_r2 * inv_r2 * inv_r2
    } else {
        0.0 // TOF terms
    }
}

impl LambertSensitivities {
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
