//! Interpolation-based initial guess for the Lambert solver.
//!
//! Uses pre-computed polynomial coefficient data (generated from Fortran `.bin` file)
//! to provide accurate initial guesses for k, reducing Newton-Raphson iterations
//! from ~5-10 down to ~1-2.

use crate::generated_coefficients::*;

// ---- Constants from Fortran gammaBetaMod ----

const BETA_LOW_BOUND: f64 = -15.77152206023506;
const MINUS_1_BY_BETA_LOW_BOUND: f64 = -1.0 / BETA_LOW_BOUND;
const GAMMA_BAR_LOW_BOUND_NZERO: f64 = -6.90775528; // ln(1e-3)
const GAMMA_MAX_ADD_NEW: f64 = 13.30314718055995;
const GAMMA_LOW_BOUND: f64 = -16.11809565095832;

// Constants for GetGammaMaxNzero
const GAMMA_MAX_BETA_BOUND1: f64 = 0.0;
const GAMMA_MAX_BETA_BOUND2: f64 = -11.82864154517630;
const GAMMA_MAX_BETA_FLAT: f64 = 1.0;
const GAMMA_M_POS_NZERO: f64 = 9.0 / 16.0; // 0.5625
const GAMMA_B_POS_NZERO: f64 = 15.0;
const GAMMA_M_NEG_NZERO: f64 = 14.0 / 11.82864154517630; // 14 / (-GAMMA_MAX_BETA_BOUND2)
const GAMMA_B_NEG_NZERO: f64 = 15.0;

// Multi-rev gamma transform constants
const MREV_LCOF: f64 = 1.0 / (GAMMA_MAX_ADD_NEW - GAMMA_LOW_BOUND);
const MREV_BCOF: f64 = GAMMA_LOW_BOUND / (GAMMA_LOW_BOUND - GAMMA_MAX_ADD_NEW);
const MIN_TB_MARGIN: f64 = -0.02;

// Nudge epsilon for bounds checking
const LIMS_NUDGE_EPS: f64 = 1000.0 * f64::EPSILON;

const SQRT_2: f64 = std::f64::consts::SQRT_2;

/// Interpolate initial k for zero-revolution transfer.
/// Returns None if inputs are out of the interpolation domain.
pub fn interpolate_initial_k_zero_rev(tau: f64, tof_by_s: f64) -> Option<f64> {
    // Step 1: Transform tau -> x via getXfromTau
    let beta = get_x_from_tau_beta(tau);
    let x = beta * MINUS_1_BY_BETA_LOW_BOUND;

    // Step 2: Bounds check on x
    let x_range = ZREV_XHI[0] - ZREV_XLOW[0];
    let x_lo = ZREV_XLOW[0] + LIMS_NUDGE_EPS * x_range;
    let x_hi = ZREV_XHI[0] - LIMS_NUDGE_EPS * x_range;
    if x < x_lo || x > x_hi {
        return None;
    }

    // Step 3: Compute TpbyS (parabolic TOF/S)
    let tp_by_s = (1.0 - tau * SQRT_2).sqrt() * (tau + SQRT_2) / 3.0;

    // Step 4: Transform tof_by_s -> y via getYfromGammaTpNzero
    let gamma_bar_max = get_gamma_max_nzero(beta);
    let gamma_bar = (tof_by_s / tp_by_s).ln();
    let y = (gamma_bar - GAMMA_BAR_LOW_BOUND_NZERO) / (gamma_bar_max - GAMMA_BAR_LOW_BOUND_NZERO);

    // Step 5: Bounds check on y
    let y_range = ZREV_XHI[1] - ZREV_XLOW[1];
    let y_lo = ZREV_XLOW[1] + LIMS_NUDGE_EPS * y_range;
    let y_hi = ZREV_XHI[1] - LIMS_NUDGE_EPS * y_range;
    if y < y_lo {
        return None; // TOF too small
    }
    // Clamp y to upper bound when tof_by_s exceeds the interpolation domain.
    // This is safe because the initial guess only needs to be "close enough" for
    // Newton-Raphson to converge — the iteration corrects any error from clamping.
    let y = y.min(y_hi);

    // Step 6: Bin lookup
    let x_bin = get_x_bin_zero_rev(x);
    let y_bin = ((y - ZREV_XLOW[1]) * ZREV_ONE_BY_DELT_X[1]) as usize;

    // Clamp bins to valid range
    let x_bin = x_bin.min(ZREV_NX_BINS - 1);
    let y_bin = y_bin.min(ZREV_NY_BINS - 1);

    // Step 7: Patch lookup (column-major: x varies fastest)
    let q_fortran = ZREV_POINT[x_bin + y_bin * ZREV_NX_BINS] as usize;
    let q = q_fortran - 1; // Convert from 1-based to 0-based

    // Step 8: Normalize coordinates
    let patch = &ZREV_DATA[q];
    let xcorn = [patch[0], patch[1]];
    let inv_del = [patch[2], patch[3]];
    let xbar = [
        (x - xcorn[0]) * inv_del[0],
        (y - xcorn[1]) * inv_del[1],
    ];

    // Step 9: Evaluate polynomial
    let k = square_poly_eval_8(&xbar, &patch[4..]);
    Some(k)
}

/// Interpolate initial k for multi-revolution transfer.
/// n_rev is the signed revolution count (positive = short period, negative = long period).
/// Returns None if inputs are out of domain or no solution exists.
pub fn interpolate_initial_k_multi_rev(tau: f64, tof_by_s: f64, n_rev: i32) -> Option<f64> {
    let n_abs = n_rev.unsigned_abs() as usize;
    if n_abs == 0 {
        return None;
    }

    // Step 1: Transform tau -> x
    let beta = get_x_from_tau_beta(tau);
    let x = beta * MINUS_1_BY_BETA_LOW_BOUND;

    // Check x bounds against mrev lims
    let x_range = MREV_LIMS[1][0] - MREV_LIMS[0][0];
    let x_lo = MREV_LIMS[0][0] + LIMS_NUDGE_EPS * x_range;
    let x_hi = MREV_LIMS[1][0] - LIMS_NUDGE_EPS * x_range;
    if x < x_lo || x > x_hi {
        return None;
    }

    // Step 2: Get x bin (dimension 1)
    let x_pntd = (x - MREV_XLOW[0]) * MREV_ONE_BY_DELT_X[0];
    let x_pnti = x_pntd as usize;
    let x_pnti = x_pnti.min(MREV_NUM_EACH_DIR[0] - 1);

    // Step 3: Compute z (N-direction)
    // NmaxTree = exp(lims_hi[2]) rounded - 1
    let n_max_tree = (MREV_LIMS[1][2].exp() + 0.01) as usize - 1;

    let (z, n_huge) = if n_abs > n_max_tree {
        // Use clamped N
        ((n_max_tree as f64).ln(), true)
    } else {
        ((n_abs as f64).ln(), false)
    };

    // Get z bin (dimension 3)
    let z_pntd = (z - MREV_XLOW[2]) * MREV_ONE_BY_DELT_X[2];
    let z_pnti = z_pntd as usize;
    let z_pnti = z_pnti.min(MREV_NUM_EACH_DIR[2] - 1);

    // Step 4: Compute T/S_bottom
    let tof_by_s_bot = if n_huge {
        // Analytical kbottom at N->infinity
        let abs_tau = tau.abs();
        let k_sol_ninf = if abs_tau < 1e-3 {
            tau + 0.5 * tau * tau * tau
        } else {
            (1.0 - (1.0 - 2.0 * tau * tau).sqrt()) / tau
        };
        // Compute T/S from k for multi-rev (quick version)
        get_tof_from_k_mrev_quick(k_sol_ninf, tau, n_abs as f64 * 2.0 * std::f64::consts::PI)
    } else {
        // Interpolate kbot from coefficient data
        eval_oct_kbot(x_pntd, z_pntd, x_pnti, z_pnti)
    };

    // Step 5: gamma = tof_by_s - T/S_bottom
    let gamma = tof_by_s - tof_by_s_bot;

    // Step 6: Transform gamma -> y via getYfromGamma
    let (y, no_solution) = get_y_from_gamma(gamma, n_rev, n_abs);
    if no_solution {
        return None;
    }

    // Step 7: Bounds check on y
    let y_range = MREV_LIMS[1][1] - MREV_LIMS[0][1];
    let y_lo = MREV_LIMS[0][1] + LIMS_NUDGE_EPS * y_range;
    let y_hi = MREV_LIMS[1][1] - LIMS_NUDGE_EPS * y_range;
    // Clamp y to interpolation bounds. Safe because the initial guess only
    // needs to be "close enough" for Newton-Raphson — iteration corrects any
    // error from boundary clamping.
    let y = y.clamp(y_lo, y_hi);

    // Step 8: Get y bin (dimension 2)
    let y_pntd = (y - MREV_XLOW[1]) * MREV_ONE_BY_DELT_X[1];
    let y_pnti = y_pntd as usize;
    let y_pnti = y_pnti.min(MREV_NUM_EACH_DIR[1] - 1);

    // Step 9: 3D patch lookup (column-major: dim1 varies fastest)
    let flat_idx = x_pnti + y_pnti * MREV_NUM_EACH_DIR[0] + z_pnti * MREV_NUM_EACH_DIR[0] * MREV_NUM_EACH_DIR[1];
    let q_fortran = MREV_POINT[flat_idx] as usize;
    let q = q_fortran - 1;

    // Step 10: Normalize coordinates
    let patch = &MREV_DATA[q];
    let xcorn = [patch[0], patch[1], patch[2]];
    let inv_del = [patch[3], patch[4], patch[5]];
    let xbar = [
        (x - xcorn[0]) * inv_del[0],
        (y - xcorn[1]) * inv_del[1],
        (z - xcorn[2]) * inv_del[2],
    ];

    // Step 11: Evaluate polynomial
    let k = cube_poly_eval_5(&xbar, &patch[6..]);
    Some(k)
}

// ---- Internal helper functions ----

/// Compute beta from tau (the intermediate in tau->x transform)
#[inline]
fn get_x_from_tau_beta(tau: f64) -> f64 {
    if tau >= 0.0 {
        -0.5 * (0.5 * (SQRT_2 - 2.0 * tau).powi(2)).ln()
    } else {
        0.5 * (0.5 * (SQRT_2 + 2.0 * tau).powi(2)).ln()
    }
}

/// Compute the gamma_max for zero-rev as a function of beta
#[inline]
fn get_gamma_max_nzero(beta: f64) -> f64 {
    if beta > GAMMA_MAX_BETA_BOUND1 {
        GAMMA_M_POS_NZERO * beta + GAMMA_B_POS_NZERO
    } else if beta < GAMMA_MAX_BETA_BOUND2 {
        GAMMA_MAX_BETA_FLAT
    } else {
        GAMMA_M_NEG_NZERO * beta + GAMMA_B_NEG_NZERO
    }
}

/// Custom x-bin lookup for zero-rev (3-zone piecewise-uniform)
#[inline]
fn get_x_bin_zero_rev(x: f64) -> usize {
    if x > ZREV_DZONE_BOUND_CUSTOM[1] {
        // Zone 3: large bins, offset by addChunk[1]
        ZREV_ADD_CHUNK[1] as usize + ((x - ZREV_XLOW[0]) / ZREV_LARGE_DX) as usize
    } else if x > ZREV_DZONE_BOUND_CUSTOM[0] {
        // Zone 2: small bins near tau=0, offset by addChunk[0]
        ZREV_ADD_CHUNK[0] as usize + ((x - ZREV_DZONE_BOUND_CUSTOM[0]) / ZREV_SMALL_DX) as usize
    } else {
        // Zone 1: large bins from Xlow
        ((x - ZREV_XLOW[0]) / ZREV_LARGE_DX) as usize
    }
}

/// Transform gamma -> y for multi-rev
#[inline]
fn get_y_from_gamma(gamma: f64, n_signed: i32, n_abs: usize) -> (f64, bool) {
    if gamma <= 0.0 {
        return (0.0, true); // No solution
    }

    let gamma_log = (gamma / n_abs as f64).ln();
    let g1 = MREV_LCOF * gamma_log + MREV_BCOF;

    if g1 < MIN_TB_MARGIN {
        return (0.0, true); // No solution (too close to bottom)
    }

    let g1 = g1.max(1e-16);

    let y = if n_signed < 0 { -g1 } else { g1 };
    (y, false)
}

/// Evaluate T/S at the bottom of the multi-rev time curve
#[inline]
fn eval_oct_kbot(x_pntd: f64, z_pntd: f64, x_pnti: usize, z_pnti: usize) -> f64 {
    let xbar = [
        x_pntd - x_pnti as f64,
        z_pntd - z_pnti as f64,
    ];

    // Index into kbot data: (xi + zi * n1) in our generated layout
    let idx = x_pnti + z_pnti * MREV_NUM_EACH_DIR[0];
    let coefs = &MREV_KBOT_DATA[idx];

    square_poly_eval_7(&xbar, coefs)
}

/// Quick T/S computation from k for multi-rev (used for N > NmaxTree)
fn get_tof_from_k_mrev_quick(k: f64, tau: f64, two_pi_n: f64) -> f64 {
    let p = 1.0 - k * tau;
    let sqrt_p = p.sqrt();

    // Compute W(k) for the given k and n_rev
    // This is a simplified version - we use the full stumpff computation
    let w = crate::stumpff::compute_w_and_derivatives(k, false, two_pi_n, 0);

    sqrt_p * (tau + p * w[0])
}

// ---- Polynomial evaluators ----
// Direct ports from Fortran squarePolyEval8run, squarePolyEval7run, cubePolyEval5run.
// Coefficient ordering matches the Fortran: grouped by increasing y-power (then z-power),
// with increasing x-power within each group, subject to total degree constraint.

/// 2D order-8 polynomial, 45 coefficients
/// Monomials x1^a * x2^b where a+b <= 8
fn square_poly_eval_8(x: &[f64; 2], c: &[f64]) -> f64 {
    let x1 = x[0];
    let x2 = x[1];

    // Powers of x1
    let x1_2 = x1 * x1;
    let x1_3 = x1_2 * x1;
    let x1_4 = x1_3 * x1;
    let x1_5 = x1_4 * x1;
    let x1_6 = x1_5 * x1;
    let x1_7 = x1_6 * x1;
    let x1_8 = x1_7 * x1;

    // Powers of x2
    let x2_2 = x2 * x2;
    let x2_3 = x2_2 * x2;
    let x2_4 = x2_3 * x2;
    let x2_5 = x2_4 * x2;
    let x2_6 = x2_5 * x2;
    let x2_7 = x2_6 * x2;
    let x2_8 = x2_7 * x2;

    // Coefficient mapping (0-based):
    // c[0]=1, c[1]=x1, c[2]=x1^2, c[3]=x1^3, c[4]=x1^4, c[5]=x1^5, c[6]=x1^6, c[7]=x1^7, c[8]=x1^8
    // c[9]=x2, c[10]=x1*x2, c[11]=x1^2*x2, c[12]=x1^3*x2, c[13]=x1^4*x2, c[14]=x1^5*x2, c[15]=x1^6*x2, c[16]=x1^7*x2
    // c[17]=x2^2, c[18]=x1*x2^2, c[19]=x1^2*x2^2, c[20]=x1^3*x2^2, c[21]=x1^4*x2^2, c[22]=x1^5*x2^2, c[23]=x1^6*x2^2
    // c[24]=x2^3, c[25]=x1*x2^3, c[26]=x1^2*x2^3, c[27]=x1^3*x2^3, c[28]=x1^4*x2^3, c[29]=x1^5*x2^3
    // c[30]=x2^4, c[31]=x1*x2^4, c[32]=x1^2*x2^4, c[33]=x1^3*x2^4, c[34]=x1^4*x2^4
    // c[35]=x2^5, c[36]=x1*x2^5, c[37]=x1^2*x2^5, c[38]=x1^3*x2^5
    // c[39]=x2^6, c[40]=x1*x2^6, c[41]=x1^2*x2^6
    // c[42]=x2^7, c[43]=x1*x2^7
    // c[44]=x2^8

    c[0]
    + c[1]*x1 + c[2]*x1_2 + c[3]*x1_3 + c[4]*x1_4 + c[5]*x1_5 + c[6]*x1_6 + c[7]*x1_7 + c[8]*x1_8
    + (c[9] + c[10]*x1 + c[11]*x1_2 + c[12]*x1_3 + c[13]*x1_4 + c[14]*x1_5 + c[15]*x1_6 + c[16]*x1_7) * x2
    + (c[17] + c[18]*x1 + c[19]*x1_2 + c[20]*x1_3 + c[21]*x1_4 + c[22]*x1_5 + c[23]*x1_6) * x2_2
    + (c[24] + c[25]*x1 + c[26]*x1_2 + c[27]*x1_3 + c[28]*x1_4 + c[29]*x1_5) * x2_3
    + (c[30] + c[31]*x1 + c[32]*x1_2 + c[33]*x1_3 + c[34]*x1_4) * x2_4
    + (c[35] + c[36]*x1 + c[37]*x1_2 + c[38]*x1_3) * x2_5
    + (c[39] + c[40]*x1 + c[41]*x1_2) * x2_6
    + (c[42] + c[43]*x1) * x2_7
    + c[44] * x2_8
}

/// 2D order-7 polynomial, 36 coefficients (used for kbot evaluation)
fn square_poly_eval_7(x: &[f64; 2], c: &[f64]) -> f64 {
    let x1 = x[0];
    let x2 = x[1];

    let x1_2 = x1 * x1;
    let x1_3 = x1_2 * x1;
    let x1_4 = x1_3 * x1;
    let x1_5 = x1_4 * x1;
    let x1_6 = x1_5 * x1;
    let x1_7 = x1_6 * x1;

    let x2_2 = x2 * x2;
    let x2_3 = x2_2 * x2;
    let x2_4 = x2_3 * x2;
    let x2_5 = x2_4 * x2;
    let x2_6 = x2_5 * x2;
    let x2_7 = x2_6 * x2;

    c[0]
    + c[1]*x1 + c[2]*x1_2 + c[3]*x1_3 + c[4]*x1_4 + c[5]*x1_5 + c[6]*x1_6 + c[7]*x1_7
    + (c[8] + c[9]*x1 + c[10]*x1_2 + c[11]*x1_3 + c[12]*x1_4 + c[13]*x1_5 + c[14]*x1_6) * x2
    + (c[15] + c[16]*x1 + c[17]*x1_2 + c[18]*x1_3 + c[19]*x1_4 + c[20]*x1_5) * x2_2
    + (c[21] + c[22]*x1 + c[23]*x1_2 + c[24]*x1_3 + c[25]*x1_4) * x2_3
    + (c[26] + c[27]*x1 + c[28]*x1_2 + c[29]*x1_3) * x2_4
    + (c[30] + c[31]*x1 + c[32]*x1_2) * x2_5
    + (c[33] + c[34]*x1) * x2_6
    + c[35] * x2_7
}

/// 3D order-5 polynomial, 56 coefficients
fn cube_poly_eval_5(x: &[f64; 3], c: &[f64]) -> f64 {
    let x1 = x[0];
    let x2 = x[1];
    let x3 = x[2];

    let x1_2 = x1 * x1;
    let x1_3 = x1_2 * x1;
    let x1_4 = x1_3 * x1;
    let x1_5 = x1_4 * x1;

    let x2_2 = x2 * x2;
    let x2_3 = x2_2 * x2;
    let x2_4 = x2_3 * x2;
    let x2_5 = x2_4 * x2;

    let x3_2 = x3 * x3;
    let x3_3 = x3_2 * x3;
    let x3_4 = x3_3 * x3;
    let x3_5 = x3_4 * x3;

    // Pure x1 terms (a+0+0 <= 5)
    c[0] + c[1]*x1 + c[2]*x1_2 + c[3]*x1_3 + c[4]*x1_4 + c[5]*x1_5
    // x2 terms (a+b+0 <= 5, b >= 1)
    + (c[6] + c[7]*x1 + c[8]*x1_2 + c[9]*x1_3 + c[10]*x1_4) * x2
    + (c[11] + c[12]*x1 + c[13]*x1_2 + c[14]*x1_3) * x2_2
    + (c[15] + c[16]*x1 + c[17]*x1_2) * x2_3
    + (c[18] + c[19]*x1) * x2_4
    + c[20] * x2_5
    // x3 terms (a+b+c <= 5, c >= 1)
    + (c[21] + c[22]*x1 + c[23]*x1_2 + c[24]*x1_3 + c[25]*x1_4) * x3
    + (c[26] + c[27]*x1 + c[28]*x1_2 + c[29]*x1_3) * x2 * x3
    + (c[30] + c[31]*x1 + c[32]*x1_2) * x2_2 * x3
    + (c[33] + c[34]*x1) * x2_3 * x3
    + c[35] * x2_4 * x3
    // x3^2 terms
    + (c[36] + c[37]*x1 + c[38]*x1_2 + c[39]*x1_3) * x3_2
    + (c[40] + c[41]*x1 + c[42]*x1_2) * x2 * x3_2
    + (c[43] + c[44]*x1) * x2_2 * x3_2
    + c[45] * x2_3 * x3_2
    // x3^3 terms
    + (c[46] + c[47]*x1 + c[48]*x1_2) * x3_3
    + (c[49] + c[50]*x1) * x2 * x3_3
    + c[51] * x2_2 * x3_3
    // x3^4 terms
    + (c[52] + c[53]*x1) * x3_4
    + c[54] * x2 * x3_4
    // x3^5
    + c[55] * x3_5
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_x_from_tau_symmetry() {
        // beta should be antisymmetric: beta(-tau) = -beta(tau)
        let tau = 0.3;
        let beta_pos = get_x_from_tau_beta(tau);
        let beta_neg = get_x_from_tau_beta(-tau);
        assert!((beta_pos + beta_neg).abs() < 1e-14,
            "beta_pos={}, beta_neg={}", beta_pos, beta_neg);
    }

    #[test]
    fn test_get_x_from_tau_zero() {
        // tau = 0 should give beta = -0.5*ln(0.5*2) = -0.5*ln(1) = 0
        let beta = get_x_from_tau_beta(0.0);
        assert!(beta.abs() < 1e-14, "beta at tau=0 = {}", beta);
    }

    #[test]
    fn test_interpolate_zero_rev_returns_some() {
        // A standard 90-degree transfer should be in domain
        // tau ~ 0.5 for equal radius 90-deg, tof_by_s varies
        let tau = 0.3536; // roughly sqrt(2)/4 for 90-deg equal radius
        let tof_by_s = 0.5; // reasonable TOF
        let result = interpolate_initial_k_zero_rev(tau, tof_by_s);
        assert!(result.is_some(), "Expected Some for standard case, got None");
    }

    #[test]
    fn test_square_poly_eval_8_constant() {
        // All zeros except c[0] = 5.0 => result = 5.0
        let mut c = [0.0f64; 45];
        c[0] = 5.0;
        let x = [0.5, 0.5];
        let result = square_poly_eval_8(&x, &c);
        assert!((result - 5.0).abs() < 1e-14);
    }

    #[test]
    fn test_square_poly_eval_8_linear() {
        // c[1] = 2.0 (x1 term), c[9] = 3.0 (x2 term)
        let mut c = [0.0f64; 45];
        c[1] = 2.0;
        c[9] = 3.0;
        let x = [0.4, 0.6];
        let result = square_poly_eval_8(&x, &c);
        let expected = 2.0 * 0.4 + 3.0 * 0.6;
        assert!((result - expected).abs() < 1e-14);
    }
}
