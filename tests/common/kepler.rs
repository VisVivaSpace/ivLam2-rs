//! Simple Kepler propagator for test-only use.
//!
//! Propagates a state (r, v) forward by time dt under two-body dynamics.
//! Uses the universal variable formulation with Stumpff functions.

/// Propagate (r, v) forward by dt under two-body dynamics with gravitational parameter mu.
/// Returns (r_final, v_final).
pub fn kepler_propagate(
    r0: &[f64; 3],
    v0: &[f64; 3],
    dt: f64,
    mu: f64,
) -> ([f64; 3], [f64; 3]) {
    let r0_mag = (r0[0].powi(2) + r0[1].powi(2) + r0[2].powi(2)).sqrt();
    let v0_mag_sq = v0[0].powi(2) + v0[1].powi(2) + v0[2].powi(2);
    let rdotv = r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2];

    // Orbital energy => semi-major axis
    let energy = v0_mag_sq / 2.0 - mu / r0_mag;
    let alpha = -2.0 * energy / mu; // = 1/a

    // Initial guess for universal variable chi
    let mut chi = if alpha > 1e-12 {
        // Elliptic
        mu.sqrt() * dt * alpha
    } else if alpha < -1e-12 {
        // Hyperbolic
        let a = 1.0 / alpha;
        let sign_dt = if dt >= 0.0 { 1.0 } else { -1.0 };
        sign_dt * (-a).sqrt() * ((-2.0 * mu * alpha * dt * dt)
            / (rdotv + sign_dt * (-mu * a).sqrt() * (1.0 - r0_mag * alpha)))
            .ln()
    } else {
        // Near-parabolic
        mu.sqrt() * dt / r0_mag
    };

    // Newton-Raphson iteration on Kepler's equation in universal variables
    let tol = 1e-14 * dt.abs().max(1.0);
    for _ in 0..50 {
        let chi2 = chi * chi;
        let psi = alpha * chi2;

        let (c2, c3) = stumpff_c2c3(psi);

        let r = chi2 * c2 + rdotv / mu.sqrt() * chi * (1.0 - psi * c3) + r0_mag * (1.0 - psi * c2);

        let f_val = r0_mag * chi * (1.0 - psi * c3)
            + rdotv / mu.sqrt() * chi2 * c2
            + chi2 * chi * c3
            - mu.sqrt() * dt;

        let df_val = r;

        let delta = f_val / df_val;
        chi -= delta;

        if delta.abs() < tol {
            break;
        }
    }

    // Compute f, g, f_dot, g_dot
    let chi2 = chi * chi;
    let psi = alpha * chi2;
    let (c2, c3) = stumpff_c2c3(psi);

    let r_mag = chi2 * c2 + rdotv / mu.sqrt() * chi * (1.0 - psi * c3) + r0_mag * (1.0 - psi * c2);

    let f = 1.0 - chi2 / r0_mag * c2;
    let g = dt - chi2 * chi / mu.sqrt() * c3;
    let g_dot = 1.0 - chi2 / r_mag * c2;
    let f_dot = mu.sqrt() / (r_mag * r0_mag) * chi * (psi * c3 - 1.0);

    let r_final = [
        f * r0[0] + g * v0[0],
        f * r0[1] + g * v0[1],
        f * r0[2] + g * v0[2],
    ];

    let v_final = [
        f_dot * r0[0] + g_dot * v0[0],
        f_dot * r0[1] + g_dot * v0[1],
        f_dot * r0[2] + g_dot * v0[2],
    ];

    (r_final, v_final)
}

/// Stumpff functions c2(psi) and c3(psi).
fn stumpff_c2c3(psi: f64) -> (f64, f64) {
    if psi > 1e-6 {
        let sqrt_psi = psi.sqrt();
        let c2 = (1.0 - sqrt_psi.cos()) / psi;
        let c3 = (sqrt_psi - sqrt_psi.sin()) / (psi * sqrt_psi);
        (c2, c3)
    } else if psi < -1e-6 {
        let sqrt_neg_psi = (-psi).sqrt();
        let c2 = (1.0 - sqrt_neg_psi.cosh()) / psi;
        let c3 = (sqrt_neg_psi.sinh() - sqrt_neg_psi) / ((-psi) * sqrt_neg_psi);
        (c2, c3)
    } else {
        // Series expansion near psi = 0
        let c2 = 1.0 / 2.0 - psi / 24.0 + psi * psi / 720.0;
        let c3 = 1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0;
        (c2, c3)
    }
}
