# Vercosine Lambert Solver: Algorithm Description

This document provides a detailed mathematical description of the vercosine Lambert solver implemented in this crate. The algorithm is based on Russell's formulation (Russell 2019, 2022).

**Contents:**

1. [Introduction](#introduction)
2. [Problem Statement](#problem-statement)
3. [The Vercosine Formulation](#the-vercosine-formulation)
4. [The W Function](#the-w-function)
5. [Root-Finding Iteration](#root-finding-iteration)
6. [Velocity Recovery](#velocity-recovery)
7. [Interpolation-Based Initial Guess](#interpolation-based-initial-guess)
8. [Multi-Revolution Solutions](#multi-revolution-solutions)
9. [Numerical Considerations](#numerical-considerations)
10. [References](#references)

---

## Introduction

**Lambert's problem** is a fundamental two-point boundary value problem in orbital mechanics. It appears throughout astrodynamics: trajectory design, orbit determination, rendezvous planning, and interplanetary mission optimization all require solving Lambert's problem — often millions of times in a single search.

The **vercosine formulation** (Russell 2019) reformulates the Lambert equation using a single iteration variable $k$ that smoothly spans all conic types. Unlike classical approaches that switch between different equations for elliptic, parabolic, and hyperbolic orbits, this formulation uses a single unified equation with a Stumpff-like function $W(k)$ that is analytic across all regimes.

Key advantages:
- **Unified equation** — one formula for ellipses, parabolas, and hyperbolas
- **Smooth across conic boundaries** — no branching or switching logic in the core iteration
- **Fast convergence** — Newton-Raphson with higher-order corrections typically converges in 2–5 iterations
- **Analytical sensitivities** — first- and second-order derivatives available via the implicit function theorem (see [sensitivities.md](sensitivities.md))

---

## Problem Statement

**Given:**
- Two position vectors $\mathbf{r}_1, \mathbf{r}_2 \in \mathbb{R}^3$
- Time of flight $T_* > 0$
- Gravitational parameter $\mu > 0$
- Transfer direction $d \in \{-1, +1\}$ (prograde or retrograde)
- Number of complete revolutions $N \geq 0$

**Find:**
- Initial velocity $\mathbf{v}_1$ and final velocity $\mathbf{v}_2$ such that a Keplerian orbit departing $\mathbf{r}_1$ with velocity $\mathbf{v}_1$ arrives at $\mathbf{r}_2$ with velocity $\mathbf{v}_2$ after time $T_*$ and $N$ complete revolutions.

The transfer angle $\theta$ is determined from the geometry:

$$\cos\theta = \frac{\mathbf{r}_1 \cdot \mathbf{r}_2}{r_1 \, r_2}$$

where $r_1 = |\mathbf{r}_1|$, $r_2 = |\mathbf{r}_2|$. The direction parameter selects:
- $d = +1$: Prograde (short way), $0 < \theta < \pi$
- $d = -1$: Retrograde (long way), $\pi < \theta < 2\pi$

---

## The Vercosine Formulation

### Geometry Parameters

The formulation uses two key geometric quantities computed from the input positions.

**Time-scaling parameter** $S$ (Russell 2022, Eq. 1):

$$S = \sqrt{\frac{(r_1 + r_2)^3}{\mu}}$$

**Normalized geometry parameter** $\tau$ (Russell 2022, Eq. 2):

$$\tau = d \cdot \frac{\sqrt{r_1 \, r_2 \, (1 + \cos\theta)}}{r_1 + r_2}$$

The parameter $\tau$ encodes the geometry of the transfer. Its properties:
- $|\tau| \in [0, 1/\sqrt{2}]$
- $\tau \to 0$ as $\theta \to \pi$ (half-revolution singularity — transfer plane is undefined)
- $|\tau| \to 1/\sqrt{2}$ as $r_1 \to r_2$ and $\theta \to 0$

### The Lambert Equation

The normalized time of flight in the vercosine formulation is (Russell 2022, Eq. 3):

$$\frac{T}{S} = \sqrt{p} \left[\tau + p \, W(k)\right]$$

where:
- $k$ is the **iteration variable** to be found by root-solving
- $p = 1 - k\tau$ is a derived quantity
- $W(k)$ is the **vercosine W function** (a Stumpff-like function, see [below](#the-w-function))
- $T$ is the actual time of flight, $S$ is the scaling parameter

The variable $k$ determines the **conic type** of the transfer orbit:

| $k$ range | Conic type |
|-----------|------------|
| $k < \sqrt{2}$ | Ellipse |
| $k = \sqrt{2}$ | Parabola |
| $k > \sqrt{2}$ | Hyperbola |

### The Root Equation

We define the residual function:

$$F(k) = \sqrt{p}\left[\tau + p \, W(k)\right] - \frac{T_*}{S} = 0$$

The solver finds $k^*$ such that $F(k^*) = 0$, then recovers the velocities from $k^*$.

### Valid Domain

**Zero-revolution** ($N = 0$):
- $k \in (-\sqrt{2}, \; 1/\tau)$ for $\tau > 0$
- $k \in (-\sqrt{2}, \; +\infty)$ for $\tau \leq 0$

**Multi-revolution** ($|N| > 0$):
- $k \in (-\sqrt{2}, \; +\sqrt{2})$ — always elliptic
- A minimum TOF exists; no solution if $T_* < T_{\min}(N)$

---

## The W Function

The $W(k)$ function is analogous to classical Stumpff functions but unifies all conic types into a single framework. It maps the iteration variable $k$ to the angular/hyperbolic component of the time of flight.

### Definition by Conic Type

**Elliptic** ($k^2 < 2$):

$$W(k) = \frac{2\pi N + \cos^{-1}(k^2 - 1)}{(2 - k^2)^{3/2}} - \frac{k}{2 - k^2}$$

For $k < 0$, use $2\pi - \cos^{-1}(k^2 - 1)$ instead of $\cos^{-1}(k^2 - 1)$.

**Parabolic** ($k = \sqrt{2}$):

$$W(\sqrt{2}) = \frac{\sqrt{2}}{3} \approx 0.47140452079103168$$

**Hyperbolic** ($k^2 > 2$):

$$W(k) = -\frac{\cosh^{-1}(k^2 - 1)}{(k^2 - 2)^{3/2}} - \frac{k}{2 - k^2}$$

### Derivatives of W

The derivatives satisfy a recurrence relation. Let $m = 2 - k^2$:

$$W' = \frac{3Wk - 2}{m}$$

$$W'' = \frac{5W'k + 3W}{m}$$

$$W''' = \frac{7W''k + 8W'}{m}$$

$$W'''' = \frac{9W'''k + 15W''}{m}$$

These derivatives are needed for:
- Newton-Raphson iteration (up to $W'''$ for third-order corrections)
- Sensitivity computation ($W'$ and $W''$)
- Hessian computation ($W'$, $W''$, $W'''$, $W''''$)

### Numerical Regions

Direct evaluation of the $W$ formulas loses precision near certain points. The implementation uses series expansions in these regions:

| Region | Trigger | Technique | Source module |
|--------|---------|-----------|---------------|
| $k \approx 0$ | $\|k\| < 0.02$ | Taylor series in $k$ | `stumpff.rs` |
| $k \approx \sqrt{2}$ | $\|k - \sqrt{2}\| < 0.02$ | Taylor series in $\nu = k - \sqrt{2}$ | `stumpff.rs` |
| $k \approx -\sqrt{2}$ | $k + \sqrt{2} < 0.01$ | Special series (multi-rev boundary) | `stumpff.rs` |
| General ellipse | otherwise | $\cos^{-1}$ formula with branch handling | `stumpff.rs` |
| Hyperbola | $k^2 > 2$ | $\log$ form of $\cosh^{-1}$ | `stumpff.rs` |

---

## Root-Finding Iteration

### Newton-Raphson with Higher-Order Corrections

The solver uses up to third-order corrections for fast convergence.

**Derivatives of** $F$ with respect to $k$ (Russell 2022, Eqs. 8–10):

$$F' = \frac{-3p\tau W + 2p^2 W' - \tau^2}{2\sqrt{p}}$$

$$F'' = \frac{3p\tau^2 W + 4p^3 W'' - 12p^2\tau W' - \tau^3}{4p^{3/2}}$$

$$F''' = \frac{3p\tau^3 W + 18p^2\tau^2 W' - 36p^3\tau W'' + 8p^4 W''' - 3\tau^4}{8p^{5/2}}$$

**Correction steps:**

First-order (Newton):

$$\delta k_1 = -\frac{F}{F'}$$

Second-order (Halley-like):

$$\delta k_2 = \frac{(\delta k_1)^2 \, F''}{2 F'}$$

Third-order:

$$\delta k_3 = \frac{(\delta k_1)^3 \, F''' / 6 + (\delta k_1)^2 \, \delta k_2 \, F''}{F'}$$

**Total correction:**

$$\Delta k = \delta k_1 + \delta k_2 + \delta k_3$$

If the series diverges ($|\delta k_2| > |\delta k_1|$ or $|\delta k_3| > |\delta k_2|$), the solver falls back to first-order Newton only.

### Convergence Criteria

The iteration terminates when:

$$|F(k)| < 10^{-14} \cdot \max(T_*/S, \; 1)$$

Maximum iterations: 25. If convergence is not achieved, a `ConvergenceFailed` error is returned.

### Safeguards

1. **Boundary clamping**: $k$ is kept within valid bounds with a margin $\varepsilon \approx 10^{-10}$
2. **Step limiting**: For multi-rev cases, $|\Delta k| \leq 0.5$ to avoid jumping to the wrong solution branch
3. **Bisection fallback**: If the correction would exit the valid domain, bisect toward the boundary
4. **Log-form TOF**: When $T_*/S > 10^4$, the solver works with $\log(T/S)$ instead of $T/S$ to prevent overflow

### Typical Convergence Behavior

| Case | Typical iterations | Notes |
|------|-------------------|-------|
| Zero-rev, good initial guess | 2–3 | Quadratic convergence |
| Zero-rev, poor initial guess | 4–6 | May need safeguards |
| Multi-rev | 3–5 | Step limiting active |
| Near boundaries | 5–10 | Reduced step sizes |
| Huge TOF | 3–5 | Log form helps |

---

## Velocity Recovery

Once $k^*$ is found, the velocities are computed from the **Lagrange $f$ and $g$ coefficients**.

### Lagrange Coefficients

$$f = 1 - p\,\frac{r_1 + r_2}{r_1}$$

$$g = S \, \tau \, \sqrt{p}$$

$$\dot{g} = 1 - p\,\frac{r_1 + r_2}{r_2}$$

where $p = 1 - k^*\tau$.

### Velocity Vectors

$$\mathbf{v}_1 = \frac{\mathbf{r}_2 - f\,\mathbf{r}_1}{g}$$

$$\mathbf{v}_2 = \frac{\dot{g}\,\mathbf{r}_2 - \mathbf{r}_1}{g}$$

**Half-revolution singularity**: When $\tau \approx 0$ (transfer angle $\approx \pi$), $g \to 0$ causing division by zero. The transfer plane is physically undefined in this case, and the velocities are meaningless.

---

## Interpolation-Based Initial Guess

The solver supports an **interpolation-based initial guess** using precomputed polynomial coefficients from the Fortran ivLam2 implementation. This is enabled by default and can be disabled with the `lightweight` Cargo feature.

The coefficient data covers the domain of typical zero-rev and multi-rev transfers using piecewise polynomial patches:
- **Zero-rev**: 347 patches in $(\tau, T/S)$ space, 8th-degree polynomials on a square grid
- **Multi-rev**: 917 patches in $(\tau, T/S, \ln N)$ space, 5th/7th-degree polynomials on a cube grid

The interpolation provides an initial $k$ guess that is typically within the convergence basin of Newton's method, reducing iteration counts by 1–3 compared to the analytical initial guess. When inputs fall outside the interpolation domain, the solver falls back to the analytical guess.

---

## Multi-Revolution Solutions

For $|N| > 0$ complete revolutions, the Lambert problem generally has **two solutions**:

- **Short-period** (positive $\tilde{N}$): Higher energy, shorter semi-major axis
- **Long-period** (negative $\tilde{N}$): Lower energy, larger semi-major axis

These correspond to two distinct values of $k$ in the elliptic domain $k \in (-\sqrt{2}, +\sqrt{2})$.

A **minimum time of flight** $T_{\min}(N)$ exists for each revolution count. If $T_* < T_{\min}(N)$, no solution exists and the solver returns `NoSolutionForRevolutions`.

The `solve_lambert_multi_rev` function returns both solutions in a `MultiRevSolution` struct containing `short_period` and `long_period` fields.

---

## Numerical Considerations

### Precision-Preserving Techniques

**Alternative $\tau$ formula near $\theta = \pi$:**

The standard formula for $\tau$ loses precision when $1 + \cos\theta \approx 0$. An alternative is used when $1 + \cos\theta < 10^{-8}$:

$$|\tau| = \frac{|\mathbf{r}_1 \times \mathbf{r}_2|}{(r_1 + r_2)\sqrt{(1 - \cos\theta) \, r_1 \, r_2}}$$

**Precomputed reciprocals:** $1/r_1$ and $1/r_2$ are stored in the `Geometry` struct to avoid repeated division.

**Log-form TOF for large flight times:** When $T_*/S > 10^4$, the solver uses:

$$\log\left[\sqrt{p}\left(\tau + p \, W\right)\right] = \log(T_*/S)$$

This prevents overflow and improves numerical conditioning.

### Error Conditions

| Error | Cause | Detection |
|-------|-------|-----------|
| `IdenticalPositions` | $\mathbf{r}_1 = \mathbf{r}_2$ | $\|\mathbf{r}_2 - \mathbf{r}_1\| < \varepsilon \|\mathbf{r}_1\|$ |
| `InvalidTimeOfFlight` | $T_* \leq 0$ | Direct check |
| `HalfRevolutionSingularity` | $\theta = \pi$ exactly | $|\tau| < 10^{-154}$ |
| `NoSolutionForRevolutions` | $T_* < T_{\min}(N)$ | Iteration exits bounds |
| `ConvergenceFailed` | Numerical issues | iterations $\geq 25$ and $|F| > 10^{-10}$ |

---

## References

1. Russell, R. P., **"On the Solution to Every Lambert Problem,"** Celestial Mechanics and Dynamical Astronomy, Vol. 131, Article 50, 2019, pp. 1–33. [doi:10.1007/s10569-019-9927-z](https://dx.doi.org/10.1007/s10569-019-9927-z)

2. Russell, R. P., **"Complete Lambert Solver Including Second-Order Sensitivities,"** Journal of Guidance, Control, and Dynamics, Vol. 45, No. 2, 2022, pp. 196–212. [doi:10.2514/1.G006089](https://doi.org/10.2514/1.G006089)

3. Arora, N., Russell, R. P., Strange, N., and Ottesen, D., **"Partial Derivatives of the Solution to the Lambert Boundary Value Problem,"** Journal of Guidance, Control, and Dynamics, Vol. 38, No. 9, 2015, pp. 1563–1572.
