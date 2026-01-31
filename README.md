# Lambert Solver

A Rust implementation of the vercosine Lambert solver based on the ivLam2 algorithm by Ryan P. Russell.

## Overview

This crate solves **Lambert's problem**, the classic two-point boundary value problem in orbital mechanics: given two position vectors **r₁** and **r₂**, a time of flight *T*, and the gravitational parameter *μ*, find the initial and final velocity vectors **v₁** and **v₂**.

The implementation uses the **vercosine formulation** which:
- Is singularity-free except for the true singularity (r₁ = r₂)
- Handles all conic types (elliptic, parabolic, hyperbolic) with a single unified equation
- Supports multi-revolution transfers (N complete orbits)
- Includes first and second-order sensitivities (with the `sensitivities` feature)

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
- `mu`: Gravitational parameter
- `direction`: `Direction::Prograde` (short way) or `Direction::Retrograde` (long way)
- `n_rev`: Number of complete revolutions (0 for direct transfer)

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

- **`sensitivities`**: Enable first and second-order partial derivatives of the solution

```toml
[dependencies]
lambert_solver = { version = "0.1", features = ["sensitivities"] }
```

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

1. **Initial guess**: The current implementation uses simple analytical approximations. 
   For production use, the interpolation tables from the original ivLam2 coefficient 
   file would provide better initial guesses.

2. **Sensitivities**: The sensitivity module is a placeholder. Full implementation
   requires translating ~3000 lines of auto-generated derivative code.

3. **No network access**: Cannot download the coefficient file for optimal initial guesses.

## License

GPL-3.0, following the original ivLam2 code.

## Acknowledgments

This implementation is based on the ivLam2 FORTRAN code by Ryan P. Russell at
The University of Texas at Austin.

- Official page: https://sites.utexas.edu/russell/publications/code/lambert/
- Code archive: https://doi.org/10.5281/zenodo.3479923
