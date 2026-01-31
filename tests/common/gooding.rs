//! Safe Rust wrapper around the C Gooding Lambert solver.
//!
//! Only available when the `gooding-ffi` feature is enabled.

extern "C" {
    fn lambert(
        gm: f64,
        r1: *const f64,
        r2: *const f64,
        nrev: i32,
        dt: f64,
        v1: *mut f64,
        v2: *mut f64,
    ) -> i32;
}

/// Result from the Gooding Lambert solver.
#[derive(Debug, Clone)]
pub struct GoodingResult {
    pub v1: [f64; 3],
    pub v2: [f64; 3],
}

/// Solve Lambert's problem using the Gooding solver.
///
/// Note: The Gooding solver only handles prograde transfers (theta in [0, pi]).
/// nrev = 0 for zero-rev, nrev < 0 for short-period multi-rev, nrev > 0 for long-period.
///
/// Returns None if the solver fails (no solution or convergence failure).
pub fn gooding_lambert(
    mu: f64,
    r1: &[f64; 3],
    r2: &[f64; 3],
    nrev: i32,
    tof: f64,
) -> Option<GoodingResult> {
    let mut v1 = [0.0_f64; 3];
    let mut v2 = [0.0_f64; 3];

    let code = unsafe {
        lambert(
            mu,
            r1.as_ptr(),
            r2.as_ptr(),
            nrev,
            tof,
            v1.as_mut_ptr(),
            v2.as_mut_ptr(),
        )
    };

    if code > 0 {
        Some(GoodingResult { v1, v2 })
    } else {
        None
    }
}
