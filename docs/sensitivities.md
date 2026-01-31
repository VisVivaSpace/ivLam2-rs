# Lambert Solver Sensitivities: First- and Second-Order Derivatives

This document describes the mathematical framework for computing analytical first-order (Jacobian) and second-order (Hessian) sensitivities of the Lambert solution. Based on Russell (2022) and Arora et al. (2015).

**Contents:**

1. [Overview](#overview)
2. [Generic Root-Solver Sensitivity Framework](#generic-root-solver-sensitivity-framework)
3. [Application to Lambert's Problem](#application-to-lamberts-problem)
4. [Jacobian: First-Order Sensitivities](#jacobian-first-order-sensitivities)
5. [Hessian: Second-Order Sensitivities](#hessian-second-order-sensitivities)
6. [Computational Cost](#computational-cost)
7. [Validation Strategy](#validation-strategy)
8. [References](#references)

---

## Overview

The Lambert solver outputs $\mathbf{z} = [\mathbf{v}_1, \mathbf{v}_2]$ (6 components) as a function of inputs $\mathbf{y} = [\mathbf{r}_1, \mathbf{r}_2, T_*]$ (7 components). The sensitivities are:

- **Jacobian** (first-order): $\frac{d\mathbf{z}}{d\mathbf{y}}$ — a $6 \times 7$ matrix
- **Hessian** (second-order): $\frac{d^2 z_i}{d\mathbf{y}^2}$ — six $7 \times 7$ symmetric matrices (one per output component)

These sensitivities are essential for:
- **Trajectory optimization** — gradient-based methods need $d\mathbf{v}/d\mathbf{r}$ and $d\mathbf{v}/dT$
- **Covariance propagation** — mapping position uncertainties to velocity uncertainties
- **Newton's method on outer problems** — Hessians enable second-order convergence in mission design

The key insight is that these derivatives can be computed **analytically** using the **implicit function theorem** (IFT), avoiding finite differences entirely. This gives exact gradients at a fraction of the computational cost.

---

## Generic Root-Solver Sensitivity Framework

The Lambert problem involves an implicit root solve: find $\hat{\alpha}$ such that $F(\hat{\alpha}, \mathbf{z}) = 0$, where $\mathbf{z}$ are problem parameters. The following framework applies to any such problem (Russell 2022, Eqs. 15–25).

### First-Order: Implicit Function Theorem

Given a root equation $F(\alpha, \mathbf{z}) = 0$ where $\hat{\alpha}(\mathbf{z})$ is the converged root, differentiate implicitly:

$$\frac{\partial F}{\partial \alpha}\frac{d\hat{\alpha}}{d\mathbf{z}} + \frac{\partial F}{\partial \mathbf{z}} = 0$$

Solving for the root sensitivity:

$$\frac{d\hat{\alpha}}{d\mathbf{z}} = -\left(\frac{\partial F}{\partial \alpha}\right)^{-1} \frac{\partial F}{\partial \mathbf{z}}$$

For an output quantity $\mathbf{p}(\alpha, \mathbf{z})$, the total derivative is:

$$\frac{d\mathbf{p}}{d\mathbf{z}} = \frac{\partial \mathbf{p}}{\partial \mathbf{z}} + \frac{\partial \mathbf{p}}{\partial \alpha}\frac{d\hat{\alpha}}{d\mathbf{z}}$$

This separates the **explicit** dependence (direct effect of $\mathbf{z}$ on $\mathbf{p}$) from the **implicit** dependence (effect through the root $\hat{\alpha}$).

### Second-Order: Hessian of the Root

Differentiating the IFT equation a second time gives the second derivative of the root (Russell 2022, Eq. 21):

$$\frac{d^2\hat{\alpha}}{dz_j \, dz_l} = -\left(\frac{\partial F}{\partial \alpha}\right)^{-1} \left[\frac{\partial^2 F}{\partial z_j \, \partial z_l} + \frac{\partial^2 F}{\partial \alpha \, \partial z_j}\frac{d\hat{\alpha}}{dz_l} + \frac{\partial^2 F}{\partial \alpha \, \partial z_l}\frac{d\hat{\alpha}}{dz_j} + \frac{\partial^2 F}{\partial \alpha^2}\frac{d\hat{\alpha}}{dz_j}\frac{d\hat{\alpha}}{dz_l}\right]$$

### Second-Order: Hessian of the Output

For each output component $p_i$, the full Hessian is (Russell 2022, Eq. 25):

$$\frac{d^2 p_i}{dz_j \, dz_l} = \underbrace{\frac{\partial^2 p_i}{\partial z_j \, \partial z_l}}_{\text{Term 1}} + \underbrace{\frac{\partial^2 p_i}{\partial \alpha \, \partial z_j}\frac{d\hat{\alpha}}{dz_l}}_{\text{Term 2}} + \underbrace{\frac{\partial^2 p_i}{\partial \alpha \, \partial z_l}\frac{d\hat{\alpha}}{dz_j}}_{\text{Term 3}} + \underbrace{\frac{\partial^2 p_i}{\partial \alpha^2}\frac{d\hat{\alpha}}{dz_j}\frac{d\hat{\alpha}}{dz_l}}_{\text{Term 4}} + \underbrace{\frac{\partial p_i}{\partial \alpha}\frac{d^2\hat{\alpha}}{dz_j \, dz_l}}_{\text{Term 5}}$$

This 5-term formula is the core of the Hessian computation. Each term requires specific partial derivatives that are computed analytically from the Lambert equation.

---

## Application to Lambert's Problem

### Variable Mapping

The generic framework maps to Lambert's problem as follows (Russell 2022, Eq. 27):

| Generic | Lambert | Description |
|---------|---------|-------------|
| $F$ | $T(k, \mathbf{r}_1, \mathbf{r}_2) - T_*$ | TOF residual |
| $\alpha$ | $k$ | Iteration variable |
| $\mathbf{z}$ | $[\mathbf{r}_1, \mathbf{r}_2, T_*]$ | 7 input parameters |
| $\mathbf{p}$ | $[\mathbf{v}_1, \mathbf{v}_2]$ | 6 output velocities |

### Required Partial Derivatives

The sensitivity computation requires (Russell 2022, Eq. 26):

**For the Jacobian:**
- $\partial F / \partial k$ — already computed during Newton iteration ($= F'$)
- $\partial F / \partial \mathbf{y}$ — partials of TOF residual w.r.t. inputs (through $\tau$, $S$, $T_*/S$)
- $\partial \mathbf{v} / \partial k$ — velocity partials w.r.t. $k$ (through Lagrange coefficients)
- $\partial \mathbf{v} / \partial \mathbf{y}|_{\text{explicit}}$ — direct dependence of velocities on inputs

**Additionally for the Hessian:**
- $\partial^2 F / \partial k^2$ — already available ($= F''$)
- $\partial^2 F / \partial k \, \partial \mathbf{y}$ — mixed second partials of TOF
- $\partial^2 F / \partial \mathbf{y}^2$ — second partials of TOF w.r.t. inputs
- $\partial^2 \mathbf{v} / \partial k^2$, $\partial^2 \mathbf{v} / \partial k \, \partial \mathbf{y}$, $\partial^2 \mathbf{v} / \partial \mathbf{y}^2$ — second partials of velocities

### Chain Rule Through Intermediate Quantities

The inputs $\mathbf{y}$ affect the Lambert equation through intermediate quantities:

$$\mathbf{y} = [\mathbf{r}_1, \mathbf{r}_2, T_*] \xrightarrow{\text{geometry}} (\tau, S, T_*/S) \xrightarrow{\text{solve}} k \xrightarrow{\text{velocity}} (\mathbf{v}_1, \mathbf{v}_2)$$

The geometry partials $\partial\tau/\partial\mathbf{y}$, $\partial S/\partial\mathbf{y}$ (and their second derivatives for Hessians) are computed analytically from:

$$\tau = d \cdot \frac{\sqrt{r_1 r_2 (1+\cos\theta)}}{r_1 + r_2}, \qquad S = \sqrt{\frac{(r_1+r_2)^3}{\mu}}$$

---

## Jacobian: First-Order Sensitivities

### Layout

The Jacobian is a $6 \times 7$ matrix:

```
              r1_x  r1_y  r1_z  r2_x  r2_y  r2_z  tof
           ┌                                           ┐
    v1_x   │  ·     ·     ·     ·     ·     ·     ·   │
    v1_y   │  ·     ·     ·     ·     ·     ·     ·   │
    v1_z   │  ·     ·     ·     ·     ·     ·     ·   │
    v2_x   │  ·     ·     ·     ·     ·     ·     ·   │
    v2_y   │  ·     ·     ·     ·     ·     ·     ·   │
    v2_z   │  ·     ·     ·     ·     ·     ·     ·   │
           └                                           ┘
```

Columns 0–2: $\partial / \partial \mathbf{r}_1$, Columns 3–5: $\partial / \partial \mathbf{r}_2$, Column 6: $\partial / \partial T_*$

### Accessor Methods

The `LambertSensitivities` struct provides accessors for common submatrices:

| Method | Returns | Shape |
|--------|---------|-------|
| `dv1_dr1()` | $\partial\mathbf{v}_1/\partial\mathbf{r}_1$ | $3 \times 3$ |
| `dv1_dr2()` | $\partial\mathbf{v}_1/\partial\mathbf{r}_2$ | $3 \times 3$ |
| `dv1_dtof()` | $\partial\mathbf{v}_1/\partial T_*$ | $3 \times 1$ |
| `dv2_dr1()` | $\partial\mathbf{v}_2/\partial\mathbf{r}_1$ | $3 \times 3$ |
| `dv2_dr2()` | $\partial\mathbf{v}_2/\partial\mathbf{r}_2$ | $3 \times 3$ |
| `dv2_dtof()` | $\partial\mathbf{v}_2/\partial T_*$ | $3 \times 1$ |
| `transpose()` | $J^T$ | $7 \times 6$ |

### Computational Cost

The Jacobian computation reuses all intermediate values from the core solve (stored in `SolverState` and `Geometry`). It adds roughly one Newton iteration's worth of computation — approximately $0.88\times$ the cost of a single Lambert solve (Russell 2022, Fig. 11), compared to $7$–$14\times$ for central finite differences.

---

## Hessian: Second-Order Sensitivities

### Layout

The Hessian consists of six $7 \times 7$ symmetric matrices, one for each output component $z_i$:

$$H_i[j][l] = \frac{d^2 z_i}{dy_j \, dy_l}, \qquad i = 0,\ldots,5; \quad j,l = 0,\ldots,6$$

where $z_0 = v_{1x}, \ldots, z_5 = v_{2z}$ and $y_0 = r_{1x}, \ldots, y_6 = T_*$.

Each matrix is symmetric: $H_i[j][l] = H_i[l][j]$, yielding 28 unique elements per output, or 168 total.

### Five-Term Assembly

Each element $H_i[j][l]$ is computed from the 5-term formula (Russell 2022, Eq. 25):

1. **Term 1**: $\partial^2 z_i / \partial y_j \partial y_l$ — explicit second derivatives of velocity w.r.t. inputs
2. **Term 2**: $(\partial^2 z_i / \partial k \partial y_j) \cdot (dk/dy_l)$ — cross terms
3. **Term 3**: $(\partial^2 z_i / \partial k \partial y_l) \cdot (dk/dy_j)$ — cross terms (symmetric with Term 2)
4. **Term 4**: $(\partial^2 z_i / \partial k^2) \cdot (dk/dy_j)(dk/dy_l)$ — second-order through $k$
5. **Term 5**: $(\partial z_i / \partial k) \cdot (d^2k/dy_j dy_l)$ — effect of Hessian of root

### Storage

The Hessian is stored as `Option<[[[f64; 7]; 7]; 6]>` in the `LambertSensitivities` struct. It is `None` when only the Jacobian was requested, and `Some(...)` when `solve_lambert_with_hessian` is called.

### Computational Cost

The Hessian computation adds significant work beyond the Jacobian — approximately $5.8\times$ the cost of a single Lambert solve (Russell 2022, Fig. 11). However, this is still much cheaper than finite differences, which would require $35$–$98\times$ (7 or 14 Lambert solves for central differences, each requiring its own Jacobian).

---

## Computational Cost

Summary of relative costs (Russell 2022, Fig. 11):

| Operation | Relative cost | Comparison |
|-----------|--------------|------------|
| Core Lambert solve | $1\times$ | Baseline |
| Jacobian (analytical) | $\sim 0.88\times$ | vs. $7$–$14\times$ for FD |
| Hessian (analytical) | $\sim 5.8\times$ | vs. $35$–$98\times$ for FD |
| **Solve + Jacobian** | **$\sim 1.88\times$** | |
| **Solve + Hessian** | **$\sim 6.8\times$** | |

---

## Validation Strategy

The sensitivity implementation is validated through three complementary approaches:

### 1. Finite-Difference Validation

Central finite differences provide a reference:

$$\frac{dz_i}{dy_j} \approx \frac{z_i(y_j + h) - z_i(y_j - h)}{2h}$$

$$\frac{d^2 z_i}{dy_j \, dy_l} \approx \frac{J_i(y_l + h) - J_i(y_l - h)}{2h}$$

where $J_i$ is the $i$-th row of the Jacobian. Step size $h \sim 10^{-7}$ balances truncation and roundoff error.

### 2. Symmetry Checks

The Hessian matrices must be symmetric: $H_i[j][l] = H_i[l][j]$. This is checked to machine precision and catches many implementation bugs.

### 3. Fortran Cross-Validation

When the `ivlam-ffi` feature is enabled, the Rust sensitivities are compared against the reference Fortran ivLam2 implementation, which computes the same derivatives using the same mathematical framework.

---

## References

1. Russell, R. P., **"Complete Lambert Solver Including Second-Order Sensitivities,"** Journal of Guidance, Control, and Dynamics, Vol. 45, No. 2, 2022, pp. 196–212. [doi:10.2514/1.G006089](https://doi.org/10.2514/1.G006089)

2. Arora, N., Russell, R. P., Strange, N., and Ottesen, D., **"Partial Derivatives of the Solution to the Lambert Boundary Value Problem,"** Journal of Guidance, Control, and Dynamics, Vol. 38, No. 9, 2015, pp. 1563–1572.

3. Russell, R. P., **"On the Solution to Every Lambert Problem,"** Celestial Mechanics and Dynamical Astronomy, Vol. 131, Article 50, 2019, pp. 1–33. [doi:10.1007/s10569-019-9927-z](https://dx.doi.org/10.1007/s10569-019-9927-z)
