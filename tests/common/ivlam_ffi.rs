//! Safe Rust wrapper around the Fortran ivLam2 solver.
//!
//! Only available when the `ivlam-ffi` feature is enabled.
//! NOTE: ivLam assumes mu=1. Callers must scale inputs/outputs.

use std::ffi::c_int;
use std::sync::{Mutex, Once};

extern "C" {
    fn ivlam_init_c(nmax: c_int, path: *const u8, pathlen: c_int, info: *mut c_int);

    fn ivlam_zero_rev_c(
        r1: *const f64,
        r2: *const f64,
        tof: *const f64,
        direction: c_int,
        v1: *mut f64,
        v2: *mut f64,
        info: *mut c_int,
        halfrev: *mut c_int,
    );

    fn ivlam_ntilde_with_derivs_c(
        r1: *const f64,
        r2: *const f64,
        tof: *const f64,
        direction: c_int,
        ntilde: c_int,
        v1: *mut f64,
        v2: *mut f64,
        info: *mut c_int,
        halfrev: *mut c_int,
        include_second: c_int,
        dzdy_t: *mut f64,
        d2zdy_t: *mut f64,
    );
}

static INIT: Once = Once::new();
/// Fortran ivLam uses global state and is NOT thread-safe.
/// All calls must be serialized via this mutex.
static FORTRAN_LOCK: Mutex<()> = Mutex::new(());

fn ensure_initialized() {
    INIT.call_once(|| {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/fortran/ivLamTree_20210202_160219_i2d8.bin", manifest_dir);
        let path_bytes = path.as_bytes();
        let mut info: c_int = 0;
        unsafe {
            ivlam_init_c(-1, path_bytes.as_ptr(), path_bytes.len() as c_int, &mut info);
        }
        assert_eq!(info, 0, "ivLam initialization failed with info={}", info);
    });
}

/// Result from the ivLam solver.
#[derive(Debug, Clone)]
pub struct IvlamResult {
    pub v1: [f64; 3],
    pub v2: [f64; 3],
    pub info: i32,
    pub halfrev: i32,
}

/// Solve zero-rev Lambert with ivLam (mu=1).
pub fn ivlam_zero_rev(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    direction: i32,
) -> IvlamResult {
    ensure_initialized();
    let _lock = FORTRAN_LOCK.lock().unwrap();

    let mut v1 = [0.0_f64; 3];
    let mut v2 = [0.0_f64; 3];
    let mut info: c_int = 0;
    let mut halfrev: c_int = 0;

    unsafe {
        ivlam_zero_rev_c(
            r1.as_ptr(),
            r2.as_ptr(),
            &tof,
            direction as c_int,
            v1.as_mut_ptr(),
            v2.as_mut_ptr(),
            &mut info,
            &mut halfrev,
        );
    }

    IvlamResult {
        v1,
        v2,
        info: info as i32,
        halfrev: halfrev as i32,
    }
}

/// Result from ivLam with first-order derivatives.
#[derive(Debug, Clone)]
pub struct IvlamDerivResult {
    pub v1: [f64; 3],
    pub v2: [f64; 3],
    pub info: i32,
    pub halfrev: i32,
    /// Jacobian dz/dy transposed: dzdyT[j][i] = ∂z_i/∂y_j
    /// where z = [v1, v2] (6 components), y = [r1, r2, tof] (7 components)
    /// Fortran layout: column-major (7, 6)
    pub dzdy_t: [[f64; 6]; 7],
}

/// Solve Lambert with derivatives using ivLam (mu=1).
pub fn ivlam_with_derivs(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    direction: i32,
    ntilde: i32,
) -> IvlamDerivResult {
    ensure_initialized();
    let _lock = FORTRAN_LOCK.lock().unwrap();

    let mut v1 = [0.0_f64; 3];
    let mut v2 = [0.0_f64; 3];
    let mut info: c_int = 0;
    let mut halfrev: c_int = 0;
    let mut dzdy_t = [[0.0_f64; 6]; 7];
    let mut d2zdy_t = [[[0.0_f64; 6]; 7]; 7];

    unsafe {
        ivlam_ntilde_with_derivs_c(
            r1.as_ptr(),
            r2.as_ptr(),
            &tof,
            direction as c_int,
            ntilde as c_int,
            v1.as_mut_ptr(),
            v2.as_mut_ptr(),
            &mut info,
            &mut halfrev,
            0, // first-order only
            dzdy_t.as_mut_ptr() as *mut f64,
            d2zdy_t.as_mut_ptr() as *mut f64,
        );
    }

    IvlamDerivResult {
        v1,
        v2,
        info: info as i32,
        halfrev: halfrev as i32,
        dzdy_t,
    }
}

/// Convert ivLam's transposed Jacobian dzdyT(7,6) to our convention jacobian[6][7].
/// Fortran dzdyT(j,i) = ∂z_i/∂y_j, stored column-major.
/// Our jacobian[i][j] = ∂z_i/∂y_j.
pub fn ivlam_dzdy_to_jacobian(dzdy_t: &[[f64; 6]; 7]) -> [[f64; 7]; 6] {
    let mut jac = [[0.0; 7]; 6];
    for i in 0..6 {
        for j in 0..7 {
            jac[i][j] = dzdy_t[j][i];
        }
    }
    jac
}
