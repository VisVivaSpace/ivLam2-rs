//! # Lambert Solver
//!
//! A Rust implementation of the vercosine Lambert solver based on the ivLam2 algorithm
//! by Ryan P. Russell.
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
//! ## Example
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

mod geometry;
mod stumpff;
mod solver;
mod velocity;
mod sensitivities;

pub use geometry::Direction;
pub use solver::{solve_lambert, solve_lambert_multi_rev, solve_lambert_with_jacobian,
                 LambertSolution, LambertError};
pub use sensitivities::LambertSensitivities;
