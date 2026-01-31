# LLM Context: ivLam2-rs Technical Reference

Compact reference for LLMs and AI assistants working with this codebase. For detailed derivations, see [algorithm.md](algorithm.md) and [sensitivities.md](sensitivities.md).

---

## Architecture

Module pipeline (executed in order for each Lambert solve):

```
geometry.rs → stumpff.rs → solver.rs → velocity.rs → sensitivities.rs
```

| Module | Responsibility | Key function |
|--------|---------------|--------------|
| `geometry.rs` | Compute τ, S, cos(θ) from r₁, r₂ | `Geometry::new()` |
| `stumpff.rs` | W(k) and derivatives W', W'', W''', W'''' | `compute_w_and_derivatives()` |
| `solver.rs` | Newton iteration on k; main API | `solve_lambert()`, `solve_lambert_with_jacobian()`, `solve_lambert_with_hessian()` |
| `velocity.rs` | Lagrange f,g coefficients → v₁, v₂ | `compute_velocities()` |
| `sensitivities.rs` | Jacobian (6×7) and Hessian (6×7×7) via IFT | `LambertSensitivities::compute_first_order()`, `compute_with_hessians()` |
| `interpolation.rs` | Polynomial initial guess from coefficient tables | `interpolate_zero_rev()`, `interpolate_multi_rev()` |
| `generated_coefficients.rs` | Const arrays of polynomial coefficients | Data only, no logic |

---

## Key Variables

| Variable | Symbol | Meaning | Range | Defined in |
|----------|--------|---------|-------|------------|
| `k` / `k_sol` | $k$ | Iteration variable / conic parameter | $(-\sqrt{2}, +\infty)$ | `solver.rs` |
| `p` | $p = 1 - k\tau$ | Derived from k and τ | $> 0$ | `solver.rs` |
| `tau` | $\tau$ | Normalized geometry parameter | $[-1/\sqrt{2}, 1/\sqrt{2}]$ | `geometry.rs` |
| `s` | $S$ | Time-scaling: $\sqrt{(r_1+r_2)^3/\mu}$ | $> 0$ | `geometry.rs` |
| `tof_by_s` | $T_*/S$ | Normalized TOF | $> 0$ | `geometry.rs` |
| `w` / `dw[0]` | $W(k)$ | Vercosine Stumpff-like function | varies | `stumpff.rs` |
| `dw[1..4]` | $W', W'', W''', W''''$ | Derivatives of W | varies | `stumpff.rs` |
| `f`, `g`, `g_dot` | $f, g, \dot{g}$ | Lagrange coefficients | — | `velocity.rs` |
| `cos_theta` | $\cos\theta$ | Transfer angle cosine | $[-1, 1]$ | `geometry.rs` |

**Conic type from k:** $k < \sqrt{2}$ → ellipse, $k = \sqrt{2}$ → parabola, $k > \sqrt{2}$ → hyperbola.

---

## Equation Quick Reference

Core Lambert equation (Russell 2022, Eq. 3):

$$T/S = \sqrt{p}\,[\tau + p\,W(k)]$$

Geometry (Eqs. 1–2):

$$S = \sqrt{(r_1+r_2)^3/\mu}, \qquad \tau = d\sqrt{r_1 r_2(1+\cos\theta)}/(r_1+r_2)$$

Residual: $F(k) = \sqrt{p}[\tau + pW(k)] - T_*/S = 0$

TOF derivatives (Eqs. 8–10):

$$F' = \frac{-3p\tau W + 2p^2 W' - \tau^2}{2\sqrt{p}}$$

$$F'' = \frac{3p\tau^2 W + 4p^3 W'' - 12p^2\tau W' - \tau^3}{4p^{3/2}}$$

W recurrences ($m = 2 - k^2$):

$$W' = (3Wk-2)/m, \quad W'' = (5W'k+3W)/m, \quad W''' = (7W''k+8W')/m$$

Lagrange coefficients:

$$f = 1 - p(r_1+r_2)/r_1, \quad g = S\tau\sqrt{p}, \quad \dot{g} = 1 - p(r_1+r_2)/r_2$$

Velocities: $\mathbf{v}_1 = (\mathbf{r}_2 - f\mathbf{r}_1)/g$, $\mathbf{v}_2 = (\dot{g}\mathbf{r}_2 - \mathbf{r}_1)/g$

IFT first-order: $dk/d\mathbf{y} = -(F')^{-1} \cdot \partial F/\partial\mathbf{y}$

IFT total derivative: $d\mathbf{v}/d\mathbf{y} = \partial\mathbf{v}/\partial\mathbf{y} + (\partial\mathbf{v}/\partial k)(dk/d\mathbf{y})$

---

## Code Mapping

| Equation / concept | Function | File |
|-------------------|----------|------|
| Geometry (τ, S) | `Geometry::new()` | `geometry.rs` |
| W(k) and derivatives | `compute_w_and_derivatives()` | `stumpff.rs` |
| F(k), F', F'', F''' | `compute_tof_function()` | `solver.rs` |
| Newton iteration | `iterate_to_convergence()` | `solver.rs` |
| Initial guess (analytical) | `compute_initial_guess()` | `solver.rs` |
| Initial guess (interpolation) | `interpolate_zero_rev()` / `interpolate_multi_rev()` | `interpolation.rs` |
| Lagrange f, g, ġ | `compute_velocities()` | `velocity.rs` |
| Jacobian (dk/dy, dv/dy) | `LambertSensitivities::compute_first_order()` | `sensitivities.rs` |
| Hessian (d²k/dy², d²v/dy²) | `LambertSensitivities::compute_with_hessians()` | `sensitivities.rs` |
| Direction enum | `Direction::Prograde` / `Retrograde` | `geometry.rs` |
| Error types | `LambertError` | `solver.rs` |

---

## Common Pitfalls

### Numerical Issues

1. **Cancellation near k=0 and k=√2**: W(k) formulas lose precision. The code switches to Taylor series (trigger: |k| < 0.02 or |k−√2| < 0.02). If modifying `stumpff.rs`, preserve these branches.

2. **Half-revolution singularity (τ→0)**: When θ≈π, τ→0 and g→0. The alternative τ formula (using cross product) is triggered when 1+cos(θ) < 1e-8. The velocity solution is physically undefined in this case.

3. **Large TOF overflow**: When T/S > 1e4, the solver switches to log-form TOF to avoid overflow. Check `huge_tof_case` flag in Geometry.

4. **Multi-rev boundary**: Near k=−√2, special series are needed. The minimum TOF for N revolutions is computed separately.

### Sign Conventions

- `d = +1` for prograde (0 < θ < π), `d = -1` for retrograde (π < θ < 2π)
- τ has the same sign as d
- `n_rev` is always non-negative in the public API; internally `n_tilde` can be ±|N|

### Sensitivity-Specific Issues

- The Jacobian reuses `dw[0]` (W) and `dw[1]` (W') from the solver state — these must be from the final converged iteration
- The Hessian needs W'' through W'''' — it calls `compute_w_and_derivatives` with `order=4`
- Second-order geometry partials (d²τ/dy², d²S/dy²) are computed fresh in the Hessian path since they're not needed for the core solve

---

## Testing Patterns

### Tolerance Tiers

| Category | Tolerance | Example |
|----------|-----------|---------|
| Algebraic identities | 1e-14 | Lagrange f·ġ − ḟ·g = 1 |
| Converged solver | 1e-12 | v₁, v₂ vs known solutions |
| FD Jacobian validation | 1e-6 | Analytical vs central FD (h=1e-7) |
| FD Hessian validation | 1e-4 | Analytical vs central FD of Jacobians |
| Cross-validation (Fortran) | 1e-10 | Rust vs Fortran ivLam2 |

### Feature-Gated Tests

- Default: core solver + interpolation + Jacobian + Hessian tests
- `gooding-ffi`: adds C Gooding cross-validation (needs `cc` build)
- `ivlam-ffi`: adds Fortran ivLam cross-validation (needs `DYLD_LIBRARY_PATH=fortran`)
- `lightweight`: disables interpolation (analytical guess only)

---

## Public API Summary

```rust
// Core solve — velocities only
solve_lambert(r1, r2, tof, mu, direction, n_rev) -> Result<LambertSolution, LambertError>

// Both multi-rev solutions
solve_lambert_multi_rev(r1, r2, tof, mu, direction, n_rev) -> Result<MultiRevSolution, LambertError>

// Solve + Jacobian (6×7)
solve_lambert_with_jacobian(r1, r2, tof, mu, direction, n_rev) -> Result<(LambertSolution, LambertSensitivities), LambertError>

// Solve + Jacobian + Hessian (6×7×7)
solve_lambert_with_hessian(r1, r2, tof, mu, direction, n_rev) -> Result<(LambertSolution, LambertSensitivities), LambertError>
```

All inputs are `&[f64; 3]` for positions, `f64` for tof/mu, `Direction` enum, `i32` for n_rev.
