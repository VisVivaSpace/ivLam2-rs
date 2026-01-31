# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**ivLam2-rs** is a Rust implementation of the vercosine Lambert solver based on Ryan P. Russell's ivLam2 algorithm. Given two position vectors, time of flight, and gravitational parameter (μ), it finds the connecting velocity vectors. The core Lambert equation is `T/S = √p · (τ + p·W(k))`, solved via Newton-Raphson iteration with higher-order corrections.

Single crate (`lambert_solver`), no workspace. Zero runtime dependencies; `approx` for tests only.

## Build Commands

```bash
cargo build                          # Build
cargo test                           # Run all tests
cargo test test_name                 # Run a single test
cargo test -- --nocapture            # Tests with stdout
cargo clippy                         # Lint
cargo test --features gooding-ffi     # Tests including C Gooding cross-validation
DYLD_LIBRARY_PATH=fortran cargo test --features gooding-ffi,ivlam-ffi  # All tests including Fortran ivLam
```

## Architecture

The solver pipeline flows through four modules in sequence:

1. **geometry.rs** — Precomputes all geometric quantities from r₁, r₂, and transfer direction. Produces `Geometry` struct with τ (normalized geometry), S (time scaling), cos(θ), and precomputed reciprocals. Handles near-π singularities with alternative formulas.

2. **stumpff.rs** — Computes W(k) and its derivatives. Region-based: Taylor series near k≈0 and k≈√2, arccos for elliptic (k²<2), acosh/log for hyperbolic (k²>2). The W function unifies all conic types into a single framework.

3. **solver.rs** — Main entry points `solve_lambert()`, `solve_lambert_multi_rev()`, and `solve_lambert_with_jacobian()`. Iterates on k using Newton-Raphson with up to 3rd-order corrections. Convergence: |F| < 1e-14 · max(T/S, 1), max 25 iterations. Multi-rev returns both short-period and long-period solutions.

4. **velocity.rs** — Computes v₁ and v₂ from Lagrange coefficients (f, g, ġ) using the solved k value.

5. **sensitivities.rs** — First-order Jacobian ∂[v₁,v₂]/∂[r₁,r₂,tof] via implicit function theorem: dk/dy = -(∂F/∂k)⁻¹·(∂F/∂y), total Dz/Dy = ∂z/∂y + (∂z/∂k)·(dk/dy). Validated against finite differences and Fortran ivLam reference.

**Key variable**: k is the iteration variable that determines conic type: k < √2 = ellipse, k = √2 = parabola, k > √2 = hyperbola.

## Skills

These skills are available:
- **rust-mastery** — idiomatic Rust, ownership, borrowing, lifetimes, traits, generics
- **space-mission-design** — reference frames, coordinate systems, time systems, orbital mechanics
- **aerospace-numerical-methods** — IEEE 754 pitfalls, tolerance tiers, stable formulas, testing strategies for scientific computing

## Notes Directory

- **ALGORITHM.md** — detailed algorithm description with mathematical notation
- **API.md** — API reference and design for derivative enhancement (`solve_lambert_with_jacobian`)
- **russell-2021-complete-lambert-solver-including-second-order-sensitivities.pdf** — Russell 2022 JGCD paper (PDF)
- **lamb.c & lamb.h** — C Gooding Lambert solver reference implementation, USE ONLY FOR TESTING. Entry point: `lambert(double gm, double r1[3], double r2[3], int nrev, double tdelt, double v1[3], double v2[3])`
- **ivLamRuntimeV2p50_739981p74401.f90** — Fortran implementation of the ivLam2 algorithm, USE ONLY FOR TESTING

## Dependencies

**Runtime:** None. Zero runtime dependencies.

**Dev/test only:**
- **approx** (dev-dependency) — floating-point comparison macros for tests
- **cc** (build-dependency, optional) — compiles `csrc/lamb.c` when the `gooding-ffi` feature is enabled for cross-validation tests

**Feature-gated FFI (testing only):**
- `gooding-ffi` — links C Gooding Lambert solver via `cc`
- `ivlam-ffi` — links pre-built Fortran `libivlam.dylib` (requires `DYLD_LIBRARY_PATH=fortran` at runtime)

## Key Patterns

1. **Error Handling**: Custom `LambertError` enum, `Result<T, LambertError>` returns.
2. **Tolerance rule:** Non-iterative, deterministic computations (closed-form formulas, algebraic identities, pure arithmetic) must assert within floating-point precision (1e-10 to 1e-14). Only use loose tolerances for inherently approximate operations (iterative solvers, interpolated ephemeris data, coordinate conversions through trig functions) — and document why.

## Workflow Instructions

1. Think through the problem, read the codebase for relevant files, and write a plan to `tasks/todo.md`.
2. The plan should have a list of todo items that you can check off as you complete them.
3. For each step, the plan should list parts of the code that SHOULD NOT be modified while working on that step.
4. Before you begin working, check in with me and I will verify the plan. Allow me to chat about the plan before we start.
5. Then, begin working on the todo items, marking them as complete as you go.
6. Please every step of the way just give me a high level explanation of what changes you made.
7. Make every task and code change you do as simple as possible. We want to avoid making any massive or complex changes. Every change should impact as little code as possible. Everything is about simplicity.
8. After you complete a major phase of the plan, add new tests as needed and commit the changes, then ask me to review your work. If user testing is needed, list what you want me to test.
9. Finally, add a review section to the `todo.md` file with a summary of the changes you made and any other relevant information.
10. DO NOT BE LAZY. NEVER BE LAZY. IF THERE IS A BUG FIND THE ROOT CAUSE AND FIX IT. NO TEMPORARY FIXES. YOU ARE A SENIOR DEVELOPER. NEVER BE LAZY
11. MAKE ALL FIXES AND CODE CHANGES AS SIMPLE AS HUMANLY POSSIBLE. THEY SHOULD ONLY IMPACT NECESSARY CODE RELEVANT TO THE TASK AND NOTHING ELSE. IT SHOULD IMPACT AS LITTLE CODE AS POSSIBLE. YOUR GOAL IS TO NOT INTRODUCE ANY BUGS. IT'S ALL ABOUT SIMPLICITY
