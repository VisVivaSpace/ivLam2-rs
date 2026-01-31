//! # Lambert Solver
//!
//! A Rust implementation of the vercosine Lambert solver based on the ivLam2 algorithm
//! by Ryan P. Russell. Supports all conic types, multi-revolution transfers, and
//! analytical first-order (Jacobian) and second-order (Hessian) sensitivities.
//!
//! For detailed algorithm documentation with mathematical derivations, see the
//! [`docs/`](https://github.com/ivansche/ivLam2-rs/tree/main/docs) directory:
//! - [`docs/algorithm.md`](https://github.com/ivansche/ivLam2-rs/blob/main/docs/algorithm.md) — Algorithm walkthrough
//! - [`docs/sensitivities.md`](https://github.com/ivansche/ivLam2-rs/blob/main/docs/sensitivities.md) — Sensitivity derivations
//!
//! ## References
//!
//! 1. Russell, Ryan P., "On the Solution to Every Lambert Problem,"
//!    Celestial Mechanics and Dynamical Astronomy, Vol. 131, Article 50, 2019, pp. 1–33,
//!    <https://dx.doi.org/10.1007/s10569-019-9927-z>
//!
//! 2. Russell, Ryan P., "Complete Lambert Solver Including Second-Order Sensitivities,"
//!    Journal of Guidance, Control, and Dynamics, vol. 45, no. 2, pp. 196–212, 2022,
//!    <https://doi.org/10.2514/1.G006089>
//!
//! 3. Arora, N., Russell, R. P., Strange, N., and Ottesen, D., "Partial Derivatives
//!    of the Solution to the Lambert Boundary Value Problem," Journal of
//!    Guidance, Control, and Dynamics, Vol. 38, No. 9, 2015, pp. 1563–1572.
//!
//! ## Examples
//!
//! ### Basic solve
//!
//! ```rust
//! use lambert_solver::{solve_lambert, solve_lambert_with_jacobian, Direction};
//!
//! // Define the problem: two position vectors and time of flight
//! let r1 = [1.0, 0.0, 0.0];        // Initial position (DU)
//! let r2 = [0.0, 1.0, 0.0];        // Final position (DU)
//! let tof = std::f64::consts::PI;  // Time of flight (TU), ~half orbit for circular
//! let mu = 1.0;                     // Gravitational parameter
//!
//! // Solve for zero-revolution transfer
//! let result = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0);
//!
//! match result {
//!     Ok(solution) => {
//!         println!("v1 = {:?}", solution.v1);
//!         println!("v2 = {:?}", solution.v2);
//!     }
//!     Err(e) => eprintln!("Error: {:?}", e),
//! }
//!
//! // Solve with analytical Jacobian
//! let (sol, sens) = solve_lambert_with_jacobian(
//!     &r1, &r2, tof, mu, Direction::Prograde, 0,
//! ).unwrap();
//! let dv1_dr1 = sens.dv1_dr1(); // 3x3 partial derivative matrix
//! ```
//!
//! ### Second-order sensitivities (Hessian)
//!
//! ```rust
//! use lambert_solver::{solve_lambert_with_hessian, Direction};
//!
//! let r1 = [1.0, 0.0, 0.0];
//! let r2 = [0.0, 1.0, 0.0];
//! let tof = std::f64::consts::PI / 2.0;
//! let mu = 1.0;
//!
//! let (sol, sens) = solve_lambert_with_hessian(
//!     &r1, &r2, tof, mu, Direction::Prograde, 0,
//! ).unwrap();
//!
//! // First-order: Jacobian is always populated
//! let dv1_dr1 = sens.dv1_dr1(); // 3x3 matrix
//!
//! // Second-order: Hessians are populated by solve_lambert_with_hessian
//! let hessians = sens.hessians.unwrap(); // [[[f64; 7]; 7]; 6]
//! let d2_v1x = hessians[0];             // 7x7 Hessian of v1.x
//! ```
//!
//! ### Multi-revolution transfer
//!
//! ```rust
//! use lambert_solver::{solve_lambert_multi_rev, Direction};
//!
//! let r1 = [1.0, 0.0, 0.0];
//! let r2 = [-1.0, 0.1, 0.0];
//! let tof = 20.0; // long enough for 1-rev solutions
//! let mu = 1.0;
//!
//! let result = solve_lambert_multi_rev(&r1, &r2, tof, mu, Direction::Prograde, 1);
//! if let Ok(multi) = result {
//!     println!("Short-period v1: {:?}", multi.short_period.v1);
//!     println!("Long-period  v1: {:?}", multi.long_period.v1);
//! }
//! ```
//!
//! ### Retrograde direction
//!
//! ```rust
//! use lambert_solver::{solve_lambert, Direction};
//!
//! let r1 = [1.0, 0.0, 0.0];
//! let r2 = [0.0, 1.0, 0.0];
//! let tof = std::f64::consts::PI / 2.0;
//! let mu = 1.0;
//!
//! // Retrograde (long way) uses transfer angle > 180°
//! let sol = solve_lambert(&r1, &r2, tof, mu, Direction::Retrograde, 0).unwrap();
//! // Retrograde takes the longer path around; different velocity than prograde
//! let pro = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0).unwrap();
//! assert!((sol.v1[0] - pro.v1[0]).abs() > 1e-6);
//! ```
//!
//! ### Physical units (Earth orbit, km and seconds)
//!
//! ```rust
//! use lambert_solver::{solve_lambert, Direction};
//!
//! let mu_earth = 398600.4418; // km^3/s^2
//!
//! // LEO departure, GEO arrival
//! let r1 = [6678.0, 0.0, 0.0];       // LEO radius (km)
//! let r2 = [0.0, 42164.0, 0.0];      // GEO radius (km)
//! let tof = 5.0 * 3600.0;            // 5-hour transfer (seconds)
//!
//! let sol = solve_lambert(&r1, &r2, tof, mu_earth, Direction::Prograde, 0).unwrap();
//! // Velocities are in km/s
//! let speed = (sol.v1[0].powi(2) + sol.v1[1].powi(2) + sol.v1[2].powi(2)).sqrt();
//! assert!(speed > 5.0 && speed < 15.0); // reasonable LEO departure speed
//! ```
//!
//! ### Error handling
//!
//! ```rust
//! use lambert_solver::{solve_lambert, LambertError, Direction};
//!
//! // Invalid mu (negative)
//! let result = solve_lambert(
//!     &[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0], 1.0, -1.0, Direction::Prograde, 0,
//! );
//! assert!(matches!(result, Err(LambertError::InvalidInput(_))));
//!
//! // Invalid mu (NaN)
//! let result = solve_lambert(
//!     &[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0], 1.0, f64::NAN, Direction::Prograde, 0,
//! );
//! assert!(matches!(result, Err(LambertError::InvalidInput(_))));
//! ```

mod geometry;
mod stumpff;
mod solver;
mod velocity;
mod sensitivities;
#[cfg(not(feature = "lightweight"))]
#[allow(dead_code, clippy::excessive_precision, clippy::unreadable_literal)]
mod generated_coefficients;
#[cfg(not(feature = "lightweight"))]
#[allow(clippy::excessive_precision)]
mod interpolation;

pub use geometry::Direction;
pub use solver::{solve_lambert, solve_lambert_multi_rev, solve_lambert_with_jacobian,
                 solve_lambert_with_hessian, LambertSolution, MultiRevSolution, LambertError};
pub use sensitivities::LambertSensitivities;
