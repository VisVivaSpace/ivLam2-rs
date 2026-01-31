# Lambert Solver

A Rust implementation of the vercosine Lambert solver based on the ivLam2 algorithm by Ryan P. Russell.

## Overview

This crate solves **Lambert's problem**, the classic two-point boundary value problem in orbital mechanics: given two position vectors **r₁** and **r₂**, a time of flight *T*, and the gravitational parameter *μ*, find the initial and final velocity vectors **v₁** and **v₂**.

The implementation uses the **vercosine formulation** which:
- Is singularity-free except for the true singularity (r₁ = r₂)
- Handles all conic types (elliptic, parabolic, hyperbolic) with a single unified equation
- Supports multi-revolution transfers (N complete orbits)
- Includes analytical first- and second-order sensitivities (Jacobian and Hessian)

## References

1. Russell, Ryan P., **"On the Solution to Every Lambert Problem,"**
   Celestial Mechanics and Dynamical Astronomy, Vol. 131, Article 50, 2019, pp. 1–33,
   https://dx.doi.org/10.1007/s10569-019-9927-z

2. Russell, Ryan P., **"Complete Lambert Solver Including Second-Order Sensitivities,"**
   Journal of Guidance, Control, and Dynamics, vol. 45, no. 2, pp. 196–212, 2022,
   https://doi.org/10.2514/1.G006089

3. Arora, N., Russell, R. P., Strange, N., and Ottesen, D., **"Partial Derivatives
   of the Solution to the Lambert Boundary Value Problem,"** Journal of
   Guidance, Control, and Dynamics, Vol. 38, No. 9, 2015, pp. 1563–1572.

## Documentation

Detailed mathematical documentation is available in the [`docs/`](docs/) directory:

- **[Algorithm Description](docs/algorithm.md)** — Complete walkthrough of the vercosine Lambert solver with LaTeX equations
- **[Sensitivities](docs/sensitivities.md)** — First-order (Jacobian) and second-order (Hessian) sensitivity derivations
- **[LLM Context](docs/llm-context.md)** — Compact technical reference for AI assistants

## Usage

```rust
use lambert_solver::{solve_lambert, Direction};
use std::f64::consts::PI;

fn main() {
    // Define the problem in canonical units (DU, TU where mu = 1)
    let r1 = [1.0, 0.0, 0.0];        // Initial position
    let r2 = [0.0, 1.0, 0.0];        // Final position
    let tof = PI / 2.0;              // Time of flight (quarter orbit)
    let mu = 1.0;                     // Gravitational parameter

    // Solve for zero-revolution prograde transfer
    let result = solve_lambert(&r1, &r2, tof, mu, Direction::Prograde, 0);

    match result {
        Ok(solution) => {
            println!("Initial velocity: {:?}", solution.v1);
            println!("Final velocity: {:?}", solution.v2);
            println!("Iterations: {}", solution.iterations);
        }
        Err(e) => eprintln!("Error: {:?}", e),
    }
}
```

## API

### Main Functions

#### `solve_lambert`

```rust
pub fn solve_lambert(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    n_rev: i32,
) -> Result<LambertSolution, LambertError>
```

Solves Lambert's problem for a single transfer.

- `r1`, `r2`: Position vectors
- `tof`: Time of flight (must be positive)
- `mu`: Gravitational parameter (must be finite and positive)
- `direction`: `Direction::Prograde` (short way) or `Direction::Retrograde` (long way)
- `n_rev`: Number of complete revolutions (0 for direct transfer)

#### `solve_lambert_with_jacobian`

```rust
pub fn solve_lambert_with_jacobian(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    n_rev: i32,
) -> Result<(LambertSolution, LambertSensitivities), LambertError>
```

Solves Lambert's problem and computes the analytical 6×7 Jacobian ∂[v₁,v₂]/∂[r₁,r₂,T] via the implicit function theorem. Adds roughly one iteration's worth of computation to the base solve.

```rust
use lambert_solver::{solve_lambert_with_jacobian, Direction};

let r1 = [1.0, 0.0, 0.0];
let r2 = [0.0, 1.0, 0.0];
let tof = std::f64::consts::PI / 2.0;

let (sol, sens) = solve_lambert_with_jacobian(&r1, &r2, tof, 1.0, Direction::Prograde, 0)
    .unwrap();

let dv1_dr1 = sens.dv1_dr1();   // 3×3 matrix
let dv1_dtof = sens.dv1_dtof();  // 3-vector
let full_jac = &sens.jacobian;   // 6×7 [[f64; 7]; 6]
```

#### `solve_lambert_with_hessian`

```rust
pub fn solve_lambert_with_hessian(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    n_rev: i32,
) -> Result<(LambertSolution, LambertSensitivities), LambertError>
```

Solves Lambert's problem and computes both the 6×7 Jacobian and the 6×7×7 Hessian (second-order sensitivities d²[v₁,v₂]/d[r₁,r₂,T]²). The Hessian is stored in `sens.hessians` as `Option<[[[f64; 7]; 7]; 6]>`.

```rust
use lambert_solver::{solve_lambert_with_hessian, Direction};

let r1 = [1.0, 0.0, 0.0];
let r2 = [0.0, 1.0, 0.0];
let tof = std::f64::consts::PI / 2.0;

let (sol, sens) = solve_lambert_with_hessian(&r1, &r2, tof, 1.0, Direction::Prograde, 0)
    .unwrap();

let dv1_dr1 = sens.dv1_dr1();              // 3×3 Jacobian block
let hessians = sens.hessians.unwrap();      // 6×7×7
let d2_v1x = hessians[0];                   // 7×7 Hessian of v1.x
```

#### `solve_lambert_multi_rev`

```rust
pub fn solve_lambert_multi_rev(
    r1: &[f64; 3],
    r2: &[f64; 3],
    tof: f64,
    mu: f64,
    direction: Direction,
    n_rev: i32,
) -> Result<MultiRevSolution, LambertError>
```

Solves for both solutions when `|n_rev| > 0` (short-period and long-period).

```rust
use lambert_solver::{solve_lambert_multi_rev, Direction};

let r1 = [1.0, 0.0, 0.0];
let r2 = [-1.0, 0.1, 0.0];
let tof = 20.0;

let multi = solve_lambert_multi_rev(&r1, &r2, tof, 1.0, Direction::Prograde, 1).unwrap();
println!("Short-period v1: {:?}", multi.short_period.v1);
println!("Long-period  v1: {:?}", multi.long_period.v1);
```

### Physical Units Example (Earth, km/s)

```rust
use lambert_solver::{solve_lambert, Direction};

let mu_earth = 398600.4418; // km^3/s^2

// LEO departure, GEO arrival
let r1 = [6678.0, 0.0, 0.0];       // LEO radius (km)
let r2 = [0.0, 42164.0, 0.0];      // GEO radius (km)
let tof = 5.0 * 3600.0;            // 5-hour transfer (seconds)

let sol = solve_lambert(&r1, &r2, tof, mu_earth, Direction::Prograde, 0).unwrap();
// Velocities are in km/s
let speed = (sol.v1[0].powi(2) + sol.v1[1].powi(2) + sol.v1[2].powi(2)).sqrt();
println!("Departure speed: {:.3} km/s", speed);
```

### Types

#### `LambertSolution`

```rust
pub struct LambertSolution {
    pub v1: [f64; 3],          // Initial velocity vector
    pub v2: [f64; 3],          // Final velocity vector
    pub n_rev: i32,            // Number of revolutions
    pub iterations: usize,      // Solver iterations
    pub residual: f64,         // Final residual
    pub k: f64,                // Solved parameter
    pub warning: Option<String>, // Near-singularity warnings
}
```

#### `LambertSensitivities`

```rust
pub struct LambertSensitivities {
    pub jacobian: [[f64; 7]; 6],              // dz/dy: rows=outputs, cols=inputs
    pub hessians: Option<[[[f64; 7]; 7]; 6]>, // d²z/dy²: populated by solve_lambert_with_hessian
}
```

Accessor methods: `dv1_dr1()`, `dv1_dr2()`, `dv1_dtof()`, `dv2_dr1()`, `dv2_dr2()`, `dv2_dtof()`, `transpose()`.

#### `Direction`

```rust
pub enum Direction {
    Prograde,   // Short way (0 < θ < π)
    Retrograde, // Long way (π < θ < 2π)
}
```

#### `LambertError`

```rust
pub enum LambertError {
    IdenticalPositions,          // r1 = r2
    InvalidTimeOfFlight,         // TOF ≤ 0
    NoSolutionForRevolutions(i32), // TOF too short for N revs
    ConvergenceFailed { iterations: usize, residual: f64 },
    HalfRevolutionSingularity,   // Exact 180° transfer
    InvalidInput(String),
}
```

## The Vercosine Formulation

The core Lambert equation in this formulation is:

$$T = S\sqrt{1-k\tau}\left[\tau + (1-k\tau)W(k)\right]$$

Where:
- $T$ is the time of flight
- $S = \sqrt{(r_1+r_2)^3/\mu}$ is a geometry parameter
- $\tau = d\sqrt{r_1 r_2(1+\cos\theta)}/(r_1+r_2)$ is the normalized geometry variable
- $d = \pm 1$ for prograde/retrograde direction
- $k$ is the iteration variable to be found
- $W(k)$ is a Stumpff-like function

The variable $k$ determines the conic type:
- $k < \sqrt{2}$: Ellipse
- $k = \sqrt{2}$: Parabola
- $k > \sqrt{2}$: Hyperbola

## Features

First-order sensitivities (Jacobian) are always available via `solve_lambert_with_jacobian`.

Optional Cargo features for testing only:

- **`gooding-ffi`**: Cross-validate against C Gooding Lambert solver
- **`ivlam-ffi`**: Cross-validate against Fortran ivLam2 reference implementation

## Implementation Notes

### Numerical Precision

The solver uses several precision-preserving techniques:
- Series expansions near $k = 0$ and $k = \sqrt{2}$ (parabola)
- Alternative tau computation for near-180° transfers
- Log-form TOF function for extremely large times of flight
- Higher-order Newton-Raphson iteration

### Multi-Revolution Cases

For $|N| > 0$ revolutions, there are generally two solutions:
- **Short-period** (positive $\tilde{N}$): Higher energy, shorter semi-major axis
- **Long-period** (negative $\tilde{N}$): Lower energy, larger semi-major axis

The minimum time of flight for N revolutions can be computed from the geometry.

### Known Limitations

Half-revolution transfers (θ = π exactly) are a true singularity — the transfer plane is undefined, and the solver returns `HalfRevolutionSingularity`.

## License

GPL-3.0, following the original ivLam2 code.

## Acknowledgments

This implementation is based on the ivLam2 FORTRAN code by Ryan P. Russell at
The University of Texas at Austin.

- Official page: https://sites.utexas.edu/russell/publications/code/lambert/
- Code archive: https://doi.org/10.5281/zenodo.3479923
