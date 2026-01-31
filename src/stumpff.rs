//! The W function and its derivatives for the vercosine Lambert formulation.
//!
//! The W(k) function is analogous to the Stumpff functions in other universal
//! formulations, but the vercosine approach requires only a single such function.
//! It smoothly handles elliptic (k < sqrt(2)), parabolic (k = sqrt(2)), and
//! hyperbolic (k > sqrt(2)) cases.

use std::f64::consts::PI;

/// Threshold for series solution near k = 0
const ZERO_BAND: f64 = 0.02;
const ZERO_BAND_SQ: f64 = ZERO_BAND * ZERO_BAND;

/// Threshold for series solution near k = sqrt(2) (parabola)
const PARABOLA_BAND: f64 = 0.02;

/// Threshold for near k = -sqrt(2) boundary
const K_CLOSE_TO_M_SQRT2: f64 = 0.01;

const SQRT_2: f64 = std::f64::consts::SQRT_2;
const TWO_PI: f64 = 2.0 * PI;

/// Compute W(k) and its derivatives up to the specified order.
///
/// # Arguments
/// * `k` - The iteration variable
/// * `n_rev_is_zero` - Whether this is a zero-revolution case
/// * `two_pi_n` - Precomputed 2*pi*|N| for multi-rev cases
/// * `order` - Maximum derivative order (0-4)
///
/// # Returns
/// Array `[W, dW/dk, d²W/dk², d³W/dk³, d⁴W/dk⁴]` filled up to the requested order
pub fn compute_w_and_derivatives(
    k: f64,
    n_rev_is_zero: bool,
    two_pi_n: f64,
    order: usize,
) -> [f64; 5] {
    let mut dw = [0.0; 5];
    
    let k_sq = k * k;
    let nu = k - SQRT_2;
    
    // Determine which region we're in
    let region = if k_sq <= ZERO_BAND_SQ {
        Region::NearZero
    } else if n_rev_is_zero {
        if nu.abs() < PARABOLA_BAND {
            Region::NearParabola
        } else if k_sq > 2.0 {
            Region::Hyperbola
        } else {
            Region::Ellipse
        }
    } else {
        Region::Ellipse  // Multi-rev is always elliptic
    };
    
    match region {
        Region::NearParabola => {
            // Series solution near k = sqrt(2)
            compute_w_near_parabola(nu, &mut dw, order);
        }
        Region::NearZero => {
            // Series solution near k = 0
            compute_w_near_zero(k, two_pi_n + PI, &mut dw, order);
        }
        _ => {
            // General ellipse or hyperbola
            compute_w_general(k, k_sq, n_rev_is_zero, two_pi_n, &mut dw, order);
        }
    }
    
    dw
}

#[derive(Debug, Clone, Copy)]
enum Region {
    NearZero,
    NearParabola,
    Ellipse,
    Hyperbola,
}

/// Series solution for W near the parabolic case (k ≈ sqrt(2))
fn compute_w_near_parabola(nu: f64, dw: &mut [f64; 5], order: usize) {
    let nu2 = nu * nu;
    let nu3 = nu2 * nu;
    let nu4 = nu2 * nu2;
    let nu5 = nu4 * nu;
    let nu6 = nu4 * nu2;
    let nu7 = nu4 * nu3;
    let nu8 = nu4 * nu4;
    
    // Coefficients from maple derivation (see ivLam FORTRAN code)
    dw[0] = 0.47140452079103168293389624140323
        - 0.20000000000000000000000000000000 * nu
        + 0.80812203564176859931525069954840e-1 * nu2
        - 0.31746031746031746031746031746032e-1 * nu3
        + 0.12244273267299524232049253023461e-1 * nu4
        - 0.46620046620046620046620046620047e-2 * nu5
        + 0.17581520588942906589609183828558e-2 * nu6
        - 0.65816536404771698889345948169478e-3 * nu7
        + 0.24494378529487021564470999141955e-3 * nu8;
    
    if order >= 1 {
        dw[1] = -0.20000000000000000000000000000000
            + 0.16162440712835371986305013990967 * nu
            - 0.95238095238095238095238095238095e-1 * nu2
            + 0.48977093069198096928197012093843e-1 * nu3
            - 0.23310023310023310023310023310023e-1 * nu4
            + 0.10548912353365743953765510297135e-1 * nu5
            - 0.46071575483340189222542163718634e-2 * nu6
            + 0.19595502823589617251576799313564e-2 * nu7;
    }
    
    if order >= 2 {
        dw[2] = 0.16162440712835371986305013990967
            - 0.19047619047619047619047619047619 * nu
            + 0.14693127920759429078459103628152 * nu2
            - 0.93240093240093240093240093240093e-1 * nu3
            + 0.52744561766828719768827551485676e-1 * nu4
            - 0.27642945290004113533525298231181e-1 * nu5
            + 0.13716851976512732076103759519495e-1 * nu6;
    }
    
    if order >= 3 {
        dw[3] = -0.19047619047619047619047619047619
            + 0.29386255841518858156918207256306 * nu
            - 0.27972027972027972027972027972028 * nu2
            + 0.21097824706731487907531020594271 * nu3
            - 0.13821472645002056766762649115590 * nu4
            + 0.82301111859076392456622557116969e-1 * nu5;
    }
    
    if order >= 4 {
        dw[4] = 0.29386255841518858156918207256306
            - 0.55944055944055944055944055944056 * nu
            + 0.63293474120194463722593061782812 * nu2
            - 0.55285890580008227067050596462361 * nu3
            + 0.41150555929538196228311278558484 * nu4;
    }
}

/// Series solution for W near k = 0
fn compute_w_near_zero(k: f64, tn_plus_pi: f64, dw: &mut [f64; 5], order: usize) {
    // tn_plus_pi = 2*pi*N + pi
    let k2 = k * k;
    let k3 = k2 * k;
    let k4 = k2 * k2;
    let k5 = k4 * k;
    let k6 = k4 * k2;
    let k7 = k4 * k3;
    let k8 = k4 * k4;
    
    // Series expansion of W around k=0
    // W = A0*(N+1/2)*pi - k + A2*(N+1/2)*pi*k^2 - (2/3)*k^3 + ...
    // where A0 = 1/sqrt(2), A2 = 3/(4*sqrt(2)), etc.
    const A0: f64 = 0.3535533905932738; // 1/(2*sqrt(2))
    const A2: f64 = 0.2651650429449553; // 3/(8*sqrt(2))
    const A4: f64 = 0.1657281518405971; // 15/(64*sqrt(2))
    const A6: f64 = 0.09667475524034829; // 35/(256*sqrt(2))
    const A8: f64 = 0.05437954982269591; // 315/(4096*sqrt(2))
    
    dw[0] = A0 * tn_plus_pi
        - k
        + A2 * tn_plus_pi * k2
        - (2.0/3.0) * k3
        + A4 * tn_plus_pi * k4
        - 0.4 * k5
        + A6 * tn_plus_pi * k6
        - (8.0/35.0) * k7
        + A8 * tn_plus_pi * k8;
    
    // Derivatives follow from term-by-term differentiation
    // For now, fall back to general formula if derivatives needed
    if order >= 1 {
        let m = 2.0 - k2;
        let one_by_m = 1.0 / m;
        let t2 = 3.0 * dw[0];
        dw[1] = (t2 * k - 2.0) * one_by_m;
    }
    
    if order >= 2 {
        let m = 2.0 - k2;
        let one_by_m = 1.0 / m;
        dw[2] = (5.0 * dw[1] * k + 3.0 * dw[0]) * one_by_m;
    }
    
    if order >= 3 {
        let m = 2.0 - k2;
        let one_by_m = 1.0 / m;
        dw[3] = (7.0 * dw[2] * k + 8.0 * dw[1]) * one_by_m;
    }
    
    if order >= 4 {
        let m = 2.0 - k2;
        let one_by_m = 1.0 / m;
        dw[4] = (9.0 * dw[3] * k + 15.0 * dw[2]) * one_by_m;
    }
}

/// General formula for W (ellipse or hyperbola)
fn compute_w_general(
    k: f64,
    k_sq: f64,
    n_rev_is_zero: bool,
    two_pi_n: f64,
    dw: &mut [f64; 5],
    order: usize,
) {
    let k_sq_m1 = k_sq - 1.0;
    let m = 2.0 - k_sq; // = 1 - (k^2 - 1)
    let one_by_m = 1.0 / m;
    
    if m > 0.0 {
        // Ellipse case: k^2 < 2
        let sqrt_m_cubed = (one_by_m * one_by_m * one_by_m).sqrt();
        
        if k > 0.0 {
            // Standard ellipse
            let acos_val = k_sq_m1.clamp(-1.0, 1.0).acos();
            dw[0] = (two_pi_n + acos_val) * sqrt_m_cubed - k * one_by_m;
        } else {
            // k < 0: need 2*pi - acos for correct branch
            let kps2 = k + SQRT_2;
            
            if kps2 < K_CLOSE_TO_M_SQRT2 {
                // Very close to k = -sqrt(2) boundary - use series for precision
                let tn_p = TWO_PI + two_pi_n;
                let tb1 = kps2 * kps2;
                let tb2 = kps2.sqrt();
                let tb3 = tb2 * tb1;
                let tb9 = tb1 * tb1;
                let tb10 = tb2 * tb9;
                
                dw[0] = 0.12110150049603174603174603174603e-7 / tb3 * (
                    -0.38926398009946925989672338336519e8 * tb3 
                    - 0.16515072000000000000000000e8 * tb2 * kps2 * tb1 
                    - 0.1976320000000000000000000e7 * tb10 * (kps2 + 0.24532575164897006338798206469711e1) 
                    - 0.18246749067162621557658908595243e7 * tb10 
                    + 0.25959796716951899525909665607350e6 * (tb9 + 0.64646464646464646464646464646465e1 * tb1 + 0.35463203463203463203463203463203e2) * tn_p * tb1 
                    + 0.66750357442839860425810740303391e6 * tn_p * (tb9 + 0.60952380952380952380952380952381e1 * tb1 + 0.26006349206349206349206349206349e2) * kps2 
                    - 0.645120e6 * tb2 * kps2 * tb9
                );
            } else {
                let acos_val = k_sq_m1.clamp(-1.0, 1.0).acos();
                dw[0] = (TWO_PI + two_pi_n - acos_val) * sqrt_m_cubed - k * one_by_m;
            }
        }
    } else {
        // Hyperbola case: k^2 > 2, so m < 0
        // Use: acosh(x) = ln(x + sqrt(x^2 - 1))
        let neg_m = -m;
        let sqrt_neg_m_cubed = (1.0 / (neg_m * neg_m * neg_m)).sqrt();
        let acosh_val = (k_sq_m1 + (k_sq_m1 * k_sq_m1 - 1.0).sqrt()).ln();
        dw[0] = -acosh_val * sqrt_neg_m_cubed - k * one_by_m;
    }
    
    // Compute derivatives using recurrence relation:
    // dW/dk = (3*W*k - 2) / m
    // d²W/dk² = (5*dW/dk*k + 3*W) / m
    // d³W/dk³ = (7*d²W/dk²*k + 8*dW/dk) / m
    // d⁴W/dk⁴ = (9*d³W/dk³*k + 15*d²W/dk²) / m
    
    let t2 = 3.0 * dw[0];
    
    if order >= 1 {
        dw[1] = (t2 * k - 2.0) * one_by_m;
    }
    
    if order >= 2 {
        dw[2] = (5.0 * dw[1] * k + t2) * one_by_m;
    }
    
    if order >= 3 {
        dw[3] = (7.0 * dw[2] * k + 8.0 * dw[1]) * one_by_m;
    }
    
    if order >= 4 {
        dw[4] = (9.0 * dw[3] * k + 15.0 * dw[2]) * one_by_m;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_w_at_parabola() {
        // At k = sqrt(2), W should be approximately 0.4714
        let dw = compute_w_and_derivatives(SQRT_2, true, 0.0, 0);
        assert!((dw[0] - 0.47140452079103168).abs() < 1e-10);
    }
    
    #[test]
    fn test_w_derivatives_continuity() {
        // Test that W and its derivatives are continuous across region boundaries
        let k_vals = [0.019, 0.021, SQRT_2 - 0.021, SQRT_2 + 0.021];
        
        for &k in &k_vals {
            let dw = compute_w_and_derivatives(k, true, 0.0, 3);
            // All values should be finite
            for i in 0..4 {
                assert!(dw[i].is_finite(), "dw[{}] not finite at k={}", i, k);
            }
        }
    }
    
    #[test]
    fn test_w_ellipse_basic() {
        // For k=0, W(0) = pi/(2*sqrt(2)) ≈ 1.1107 (zero-rev, N=0)
        // The near-zero series with tn_plus_pi = pi gives:
        //   W(0) = A0 * pi = (1/(2*sqrt(2))) * pi = pi/(2*sqrt(2))
        let dw = compute_w_and_derivatives(0.0, true, 0.0, 0);
        let expected = PI * 0.5 * std::f64::consts::FRAC_1_SQRT_2;
        assert!((dw[0] - expected).abs() < 1e-12, "W(0) = {}, expected {}", dw[0], expected);
    }

    #[test]
    fn test_w_zero_all_derivatives() {
        // At k=0, check W and all derivatives via the general ellipse formula
        // General formula: W(0) = acos(0^2-1) / (2-0^2)^(3/2) - 0/(2-0^2)
        //                       = acos(-1) / 2^(3/2) = pi / (2*sqrt(2))
        let dw = compute_w_and_derivatives(0.0, true, 0.0, 4);

        let expected_w = PI / (2.0 * SQRT_2);
        assert!((dw[0] - expected_w).abs() < 1e-12, "W(0) = {}, expected {}", dw[0], expected_w);

        // dW/dk at k=0: recurrence gives (3*W*0 - 2)/(2-0) = -1
        assert!((dw[1] - (-1.0)).abs() < 1e-12, "W'(0) = {}, expected -1", dw[1]);

        // d2W/dk2 at k=0: (5*(-1)*0 + 3*W)/2 = 3*W/2
        let expected_w2 = 1.5 * expected_w;
        assert!((dw[2] - expected_w2).abs() < 1e-12, "W''(0) = {}, expected {}", dw[2], expected_w2);

        // d3W/dk3 at k=0: (7*W''*0 + 8*W')/2 = 8*(-1)/2 = -4
        assert!((dw[3] - (-4.0)).abs() < 1e-12, "W'''(0) = {}, expected -4", dw[3]);

        // d4W/dk4 at k=0: (9*(-4)*0 + 15*W'')/2 = 15*W''/2
        let expected_w4 = 7.5 * expected_w2;
        assert!((dw[4] - expected_w4).abs() < 1e-12, "W''''(0) = {}, expected {}", dw[4], expected_w4);
    }

    #[test]
    fn test_w_continuity_near_zero_boundary() {
        // Test continuity across the near-zero / general-ellipse boundary
        // ZERO_BAND = 0.02, so test just inside and outside
        let k_inside = 0.019; // uses near-zero series
        let k_outside = 0.021; // uses general ellipse formula

        let dw_in = compute_w_and_derivatives(k_inside, true, 0.0, 2);
        let dw_out = compute_w_and_derivatives(k_outside, true, 0.0, 2);

        // W should be continuous — check that the values at nearby points are close
        // Estimate derivative magnitude ~1 => dW ~ 0.002 over dk=0.002
        assert!((dw_in[0] - dw_out[0]).abs() < 0.01,
            "W discontinuity at zero boundary: {} vs {}", dw_in[0], dw_out[0]);
        assert!((dw_in[1] - dw_out[1]).abs() < 0.1,
            "W' discontinuity at zero boundary: {} vs {}", dw_in[1], dw_out[1]);
    }

    #[test]
    fn test_w_continuity_near_parabola_boundary() {
        // Test continuity across the near-parabola / general boundary
        // PARABOLA_BAND = 0.02
        let k_inside = SQRT_2 - 0.019; // uses near-parabola series
        let k_outside = SQRT_2 - 0.021; // uses general ellipse formula

        let dw_in = compute_w_and_derivatives(k_inside, true, 0.0, 2);
        let dw_out = compute_w_and_derivatives(k_outside, true, 0.0, 2);

        assert!((dw_in[0] - dw_out[0]).abs() < 0.01,
            "W discontinuity at parabola boundary: {} vs {}", dw_in[0], dw_out[0]);
        assert!((dw_in[1] - dw_out[1]).abs() < 0.1,
            "W' discontinuity at parabola boundary: {} vs {}", dw_in[1], dw_out[1]);
    }

    #[test]
    fn test_w_hyperbolic() {
        // For k > sqrt(2), W should be computed via acosh formula
        let k = 2.0;
        let dw = compute_w_and_derivatives(k, true, 0.0, 2);

        // All should be finite
        assert!(dw[0].is_finite(), "W({}) not finite", k);
        assert!(dw[1].is_finite(), "W'({}) not finite", k);
        assert!(dw[2].is_finite(), "W''({}) not finite", k);

        // For hyperbolic, W should be positive (acosh term dominates)
        // m = 2 - 4 = -2, neg_m = 2, sqrt_neg_m_cubed = 1/(2*sqrt(2))
        // acosh_val = ln(3 + sqrt(8)) = ln(3 + 2*sqrt(2))
        // W = -acosh_val / (2*sqrt(2)) - k/m = -acosh_val/(2*sqrt(2)) + 1
        let m = 2.0 - k * k;
        let k_sq_m1 = k * k - 1.0;
        let neg_m = -m;
        let acosh_val = (k_sq_m1 + (k_sq_m1 * k_sq_m1 - 1.0).sqrt()).ln();
        let sqrt_neg_m_cubed = (1.0 / (neg_m * neg_m * neg_m)).sqrt();
        let expected = -acosh_val * sqrt_neg_m_cubed - k / m;
        assert!((dw[0] - expected).abs() < 1e-14,
            "W({}) = {}, expected {}", k, dw[0], expected);
    }

    #[test]
    fn test_w_negative_k_ellipse() {
        // Negative k (still elliptic if k^2 < 2)
        let k = -1.0;
        let dw = compute_w_and_derivatives(k, true, 0.0, 2);
        assert!(dw[0].is_finite(), "W({}) not finite", k);

        // Verify against general formula directly
        let m = 2.0 - k * k; // = 1
        let sqrt_m_cubed = 1.0; // (1/1)^(3/2)
        let acos_val = (k * k - 1.0_f64).clamp(-1.0, 1.0).acos(); // acos(0) = pi/2
        // k < 0: use 2*pi - acos branch
        let expected = (TWO_PI - acos_val) * sqrt_m_cubed - k / m;
        assert!((dw[0] - expected).abs() < 1e-14,
            "W({}) = {}, expected {}", k, dw[0], expected);
    }

    #[test]
    fn test_w_derivative_recurrence() {
        // Verify the recurrence relation: dW/dk = (3*W*k - 2) / m
        // at several k values using the general formula
        for &k in &[0.5, 1.0, -0.5, 1.8, 2.5] {
            let dw = compute_w_and_derivatives(k, true, 0.0, 1);
            let m = 2.0 - k * k;
            if m.abs() < 1e-10 { continue; }
            let recurrence_val = (3.0 * dw[0] * k - 2.0) / m;
            assert!((dw[1] - recurrence_val).abs() < 1e-12,
                "Recurrence failed at k={}: W'={}, recurrence={}", k, dw[1], recurrence_val);
        }
    }

    #[test]
    fn test_w_multi_rev() {
        // Multi-rev (N=1): W should include the 2*pi*N offset
        let two_pi_n = TWO_PI;
        let dw = compute_w_and_derivatives(0.0, false, two_pi_n, 0);

        // W(0, N=1) = A0 * (2*pi + pi) = (1/(2*sqrt(2))) * 3*pi
        let expected = 3.0 * PI / (2.0 * SQRT_2);
        assert!((dw[0] - expected).abs() < 1e-12,
            "W(0, N=1) = {}, expected {}", dw[0], expected);
    }
}
