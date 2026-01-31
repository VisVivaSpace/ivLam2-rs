# ivLam2-rs: Expanded Testing, Benchmarks & Coefficient File Integration

## Phase 1: Expand Gooding Cross-Validation Tests

### 1.1 Add wide-range test cases to `tests/cross_validate_gooding.rs`
- [x] Angle sweep: 5°, 15°, 25°, 30°, 60°, 75°, 100°, 110°, 140°, 150°, 160°
- [x] Radius ratio sweep (90° transfer): r2/r1 = 0.3, 0.5, 0.7, 1.5, 2.0, 3.0, 5.0
- [x] TOF variation (fixed 90° geometry): TOF = 0.2, 0.5, 1.0, 2.0, 5.0, 10.0
- [x] 3D transfers: several out-of-plane cases (30°, 60° inclination, both-off-plane, different radii inclined)
- [x] Physical units: Earth-Mars, LEO-to-MEO, Sun-centered Venus
- [x] Large eccentricity: r2/r1 = 7, 10, and 6 at 60°

**DO NOT modify:** `src/` files, existing test cases, `CROSS_TOL` constant

### 1.2 Speed benchmark: Rust vercosine vs C Gooding
- [x] Create `benches/solver_benchmark.rs` using Criterion
- [x] Add criterion dependency + `[[bench]]` to `Cargo.toml`
- [x] Benchmark: single solve (90°), batch of 100 (angle sweep), LEO→GEO, hyperbolic

**DO NOT modify:** `src/` files

---

## Phase 2: Expand Fortran Cross-Validation & Fuzz Testing

### 2.1 Add wide-range test cases to Fortran cross-validation
- [x] Add cases to `tests/cross_validate_three_way.rs` (angle sweep, radius ratio sweep, TOF sweep, 3D, high eccentricity)
- [x] Add wider derivative cases to `tests/derivative_validation.rs` (10°, 160°, ratio 3, ratio 0.5, 3D inclined, long TOF, high eccentricity)
- [x] Add ivLam-only derivative tests (small angle, large angle, diff radii, hyperbolic, long TOF, 3D inclined)

### 2.2 Fuzz testing with randomized inputs
- [x] Create `tests/fuzz_validation.rs` (behind `ivlam-ffi` feature)
- [x] Simple inline seeded PRNG (xoshiro256**, no external dependency)
- [x] 1000 random velocity cases, Rust vs Fortran ivLam
- [x] 100 random Jacobian cases, Rust vs Fortran ivLam

---

## Phase 3: Coefficient File Integration

### 3.1 Binary parser tool
- [x] Create `tools/parse_coefficients.rs` — standalone binary that reads Fortran `.bin` and generates `src/generated_coefficients.rs`
- [x] Add `[[bin]]` entry to Cargo.toml
- [x] Run parser to generate `src/generated_coefficients.rs`

### 3.2 Interpolation module
- [x] Create `src/interpolation.rs` with zero-rev and multi-rev interpolation functions
- [x] Polynomial evaluators: `square_poly_eval_8`, `cube_poly_eval_5`, `square_poly_eval_7`
- [x] Transform functions: `get_x_from_tau`, `get_y_from_gamma_tp_nzero`, `get_x_bin_zero_rev`

### 3.3 Wire into solver
- [x] Modify `src/solver.rs`: add `#[cfg(not(feature = "lightweight"))]` branch in `compute_initial_guess`
- [x] Modify `src/lib.rs`: add conditional `mod interpolation` and `mod generated_coefficients`
- [x] Add `lightweight` feature to `Cargo.toml`

### 3.4 Validation tests
- [x] Create `tests/interpolation_validation.rs`
- [x] Verify solutions match existing tests
- [x] Verify iteration count reduction with interpolation

### 3.5 Final verification
- [x] `cargo test` passes (default, with interpolation) — 73 tests pass
- [x] `cargo test --features lightweight` passes (analytical fallback) — all tests pass

---

## Phase 4: Second-Order Sensitivities (Hessians)

### 4.1 Implement `compute_with_hessians()` in `sensitivities.rs`
- [x] Refactor: extract shared first-order intermediates into helper struct
- [x] Implement `compute_second_order_geometry` (d2tau, d2s)
- [x] Implement `compute_second_order_f_partials` (d2F/dk2, d2F/dkdy, d2F/dy2)
- [x] Implement `compute_second_order_ift` (d2k/dy2)
- [x] Implement `compute_second_order_velocity` (d2z/dk2, d2z/dkdy, d2z/dy2)
- [x] Assembly: fill in `compute_with_hessians()` using all helper functions
- [x] `cargo test` passes

**DO NOT modify:** `geometry.rs`, `solver.rs`, `velocity.rs`, `stumpff.rs`, `compute_first_order()`

### 4.2 Add `solve_lambert_with_hessian` API
- [x] Add `solve_lambert_with_hessian()` to `solver.rs`
- [x] Export from `lib.rs`
- [x] `cargo test` passes

**DO NOT modify:** `solve_lambert()`, `solve_lambert_with_jacobian()`

### 4.3 Extend Fortran FFI for second-order validation
- [x] Add `IvlamHessianResult` struct and `ivlam_with_hessians()` to `tests/common/ivlam_ffi.rs`

**DO NOT modify:** `ivlam_c_wrapper.f90`, Fortran build

### 4.4 Validation tests
- [x] Create `tests/hessian_validation.rs`
- [x] Finite-difference Hessians (central difference of Jacobians)
- [x] Symmetry check (H[j][l] == H[l][j])
- [x] Fortran FFI comparison (feature-gated `ivlam-ffi`)
- [x] Full test suite passes

### 4.5 Bug fixes in Hessian computation
- [x] Fix Bug #1: `d2_tof_by_s` parenthesization error (line 784)
- [x] Fix Bug #2: `d2sqrtp_jl` sign error (line 1122)
- [x] Remove debug instrumentation and test files
- [x] All 85 tests pass, clippy clean

---

## Phase 5: Create `docs/` Directory with Algorithm Documentation

### 5.1 Create `docs/algorithm.md`
- [x] Detailed algorithm walkthrough with LaTeX math blocks
- [x] Covers: problem statement, vercosine formulation, W function, root-finding, velocity recovery, interpolation, multi-rev

### 5.2 Create `docs/sensitivities.md`
- [x] First- and second-order sensitivity derivations
- [x] Generic IFT framework, application to Lambert's problem, Jacobian/Hessian layout

### 5.3 Create `docs/llm-context.md`
- [x] Compact technical reference for LLMs
- [x] Architecture, key variables, equation reference, code mapping, pitfalls

### 5.4 Update `README.md`
- [x] Add "Documentation" section with link to `docs/`
- [x] Remove outdated "Known Limitations" (Hessians and interpolation are now implemented)
- [x] Add `solve_lambert_with_hessian` to API section

### 5.5 Update `src/lib.rs`
- [x] Add doc comment pointing to `docs/` for algorithm details
- [x] Mention Hessian API

**DO NOT modify:** Any `.rs` files other than `src/lib.rs` doc comments, `notes/`, test files, `CLAUDE.md`

---

## Phase 6: Code Review Fixes

### 6.1 H1 + H2 + M1: Input validation and deduplication
- [x] Extract shared `validate_inputs()` helper that validates mu > 0, finite checks on r1/r2, TOF > 0, identical positions
- [x] Extract shared `solve_lambert_core()` helper that does validation + geometry + k-bounds + iteration
- [x] Extract shared `build_solution()` helper for building the LambertSolution
- [x] Refactor `solve_lambert()`, `solve_lambert_with_jacobian()`, `solve_lambert_with_hessian()` to use the helpers
- [x] Add tests for mu <= 0, mu = NaN, mu = Inf, r1 with NaN, r2 with Inf, zero position, NaN tof

**DO NOT modify:** `geometry.rs`, `stumpff.rs`, `velocity.rs`, `sensitivities.rs`, `interpolation.rs`

### 6.2 M3: Document convergence threshold relationship
- [x] Add comment explaining why the post-loop fallback threshold (1e-10) differs from the iteration convergence criterion (1e-14)

### 6.3 M4: Document geometry constants
- [x] Add inline comments to `NEAR_N_PI_REV_WARNING`, `ALTERNATE_TAU_THRESHOLD`, `SQRT_TINY`, `FOUR_RT_TINY`, `TOF_BY_S_HUGE_BOUNDARY`

### 6.4 M2: Document stumpff.rs coefficients
- [x] Add source comment to near-k=-sqrt(2) series coefficients in `compute_w_general()`
- [x] Add comment explaining `K_CLOSE_TO_M_SQRT2` narrower band width

### 6.5 M5: Document intentional clamping in interpolation
- [x] Add comment at y-clamp sites in `interpolate_initial_k_zero_rev()` and `interpolate_initial_k_multi_rev()`

---

## Phase 7: Documentation Fixes & Examples

### 7.1 Export `MultiRevSolution` from `lib.rs`
- [x] Add `MultiRevSolution` to the `pub use solver::{...}` line

### 7.2 Fix stale `docs/llm-context.md`
- [x] Code mapping: `iterate_to_convergence()` → `solve_lambert_core()` (internal)
- [x] Code mapping: add `validate_inputs()` row
- [x] Code mapping: `interpolate_zero_rev()` → `interpolate_initial_k_zero_rev()`, `interpolate_multi_rev()` → `interpolate_initial_k_multi_rev()`

### 7.3 Fix `docs/algorithm.md` error table
- [x] Add `InvalidInput` row to the Error Conditions table

### 7.4 Fix `README.md` staleness
- [x] Update `LambertSensitivities` struct display to include `hessians` field
- [x] Update `mu` parameter description to mention "must be finite and positive"

### 7.5 Add doc-test examples to `lib.rs`
- [x] Hessian (`solve_lambert_with_hessian`)
- [x] Multi-rev (`solve_lambert_multi_rev`)
- [x] Retrograde direction
- [x] Physical units (Earth orbit, km/s)
- [x] Error handling showing `InvalidInput`

### 7.6 Add examples to `README.md`
- [x] Hessian usage with accessing `sens.hessians`
- [x] Multi-rev showing both short/long period
- [x] Physical units example (LEO-to-GEO)

### 7.7 Verification
- [x] `cargo test` — 93 tests pass (41 unit + 45 integration + 7 doc-tests)
- [x] `cargo test --doc` — all 7 doc-tests pass
- [x] `cargo clippy` — no new warnings (pre-existing only)

---

## Review

### Phase 7 Summary — Documentation Fixes & Examples

**Bug fix:** Exported `MultiRevSolution` from `lib.rs` — it was defined in `solver.rs` and used in the public API (`solve_lambert_multi_rev` return type) but not re-exported, so downstream users couldn't name the type.

**Stale docs fixed:**
- `docs/llm-context.md`: Updated code mapping table — `iterate_to_convergence()` → `solve_lambert_core()`, added `validate_inputs()` row, fixed interpolation function names to `interpolate_initial_k_zero_rev()` / `interpolate_initial_k_multi_rev()`
- `docs/algorithm.md`: Added `InvalidInput` row to the Error Conditions table
- `README.md`: Updated `LambertSensitivities` struct to show `hessians` field, added "must be finite and positive" to `mu` parameter description

**New examples:**
- `src/lib.rs`: 5 new doc-test examples (Hessian, multi-rev, retrograde, physical units LEO-to-GEO, error handling with `InvalidInput`)
- `README.md`: 3 new code examples (Hessian with `sens.hessians`, multi-rev short/long period, physical units LEO-to-GEO)

**Files modified:** `src/lib.rs`, `README.md`, `docs/algorithm.md`, `docs/llm-context.md`, `tasks/todo.md`

**Test results:** 93 tests pass (41 unit + 45 integration + 7 doc-tests), clippy clean (pre-existing warnings only).

---

### Phase 6 Summary — Code Review Fixes

**Refactoring (H1+H2+M1):**
Extracted three shared helpers in `solver.rs`:
- `validate_inputs()` — validates mu > 0 and finite, position vectors finite and non-zero, TOF finite and positive, non-identical positions. This addresses H1 (missing mu validation) and H2 (no NaN/Inf checks).
- `solve_lambert_core()` — does validation + geometry + half-rev check + k-bounds + initial guess + Newton-Raphson iteration + convergence check. Returns `(SolverState, Geometry, residual)`.
- `build_solution()` — computes velocities, builds warning string, assembles `LambertSolution`.

The three public functions (`solve_lambert`, `solve_lambert_with_jacobian`, `solve_lambert_with_hessian`) are now 3-6 line thin wrappers that call `solve_lambert_core()`, optionally compute sensitivities, and call `build_solution()`. This eliminates the M1 duplication (~180 lines of repeated code removed).

**New tests:** 3 tests added — `test_error_invalid_mu` (mu=0, mu<0, mu=NaN, mu=Inf), `test_error_invalid_position_vectors` (r1 NaN, r2 Inf, r=0, tof NaN), `test_validation_applies_to_all_entry_points` (verifies all 3 public APIs reject bad mu).

**Documentation (M2-M5):**
- M2: Added source comment to stumpff.rs near-k=-sqrt(2) coefficients (Taylor expansion from Fortran ivLam reference)
- M3: Added comment to convergence fallback threshold explaining the 1e-10 vs 1e-14 relationship
- M4: Added inline comments to all 5 geometry threshold constants explaining their numerical motivation
- M5: Added comments at both y-clamp sites in interpolation.rs explaining why clamping is safe

**Files modified:** `src/solver.rs` (structural refactor + tests), `src/geometry.rs` (comments only), `src/stumpff.rs` (comments only), `src/interpolation.rs` (comments only), `tasks/todo.md`

**Test results:** 88 tests pass (41 unit + 47 integration), clippy clean (pre-existing warnings only).

---

### Phase 4.5 Summary — Hessian Bug Fixes

**Two bugs** in `compute_with_hessians()` in `src/sensitivities.rs` caused 24 of 294 Hessian elements to fail FD validation. Symmetry tests passed because both bugs produced symmetric errors.

**Bug #1 — `d2_tof_by_s` parenthesization (line 784):**
The formula `d²(T/S)/dy²` had incorrect operator precedence. The `*inv_s` was applied to the entire parenthesized expression instead of just the `d2s[j][l]` term:
- Wrong: `tof_by_s * (2.0 * ds[j] * ds[l] * inv_s_sq - d2s[j][l]) * inv_s`
- Correct: `tof_by_s * (2.0 * ds[j] * ds[l] * inv_s_sq - d2s[j][l] * inv_s)`
This caused ±0.695 errors in 16 elements where both j,l correspond to non-zero unit vector components.

**Bug #2 — `d2sqrtp_jl` sign error (line 1122):**
The second derivative `d²(√p)/dy_j·dy_l` had a sign error in the second term. The chain rule for `d/dy_l[−k·dτ_j/(2√p)]` produces a negative sign via `d(1/√p)/dy = −dp/(2p√p)`, but the code had a positive sign:
- Wrong: `−k·d²τ/(2√p) + k²·dτ_j·dτ_l/(4p√p)`
- Correct: `−k·d²τ/(2√p) − k²·dτ_j·dτ_l/(4p√p)`
This caused ±0.125 errors in the remaining 8 elements.

**Debugging approach:** Term-by-term decomposition of the 5-part Hessian assembly formula, with FD validation of each intermediate quantity (d2_tof_by_s, d2f_dydy, d2k_dy, d2g_lag, d2sqrtp).

**Cleanup:** Removed all temporary debug `eprintln!` statements from `compute_with_hessians()`, removed the internal `debug_fd_explicit_d2z` test, and deleted `tests/hessian_debug.rs` and `tests/hessian_term_debug.rs`.

**Validation:** All 12 hessian_validation tests pass (6 FD + 6 symmetry), 85 total tests pass, clippy clean.

**User testing needed:**
- `DYLD_LIBRARY_PATH=fortran cargo test --features ivlam-ffi` — Fortran ivLam Hessian cross-validation

---

### Phase 3 Summary

**New files:**
- `tools/parse_coefficients.rs` — Standalone binary that reads Fortran unformatted `.bin` coefficient file and generates Rust const arrays
- `src/generated_coefficients.rs` — ~3MB auto-generated file with all interpolation data as const arrays (zero-rev: 347 patches, multi-rev: 917 patches + kbot data)
- `src/interpolation.rs` — Zero-rev and multi-rev interpolation logic: coordinate transforms (tau→x, tof→y), custom bin lookup, polynomial evaluation
- `tests/interpolation_validation.rs` — 10 validation tests covering angle sweep, hyperbolic, retrograde, 3D, physical units, multi-rev, energy conservation

**Modified files:**
- `Cargo.toml` — Added `[[bin]]` for parser, `lightweight` feature flag
- `src/lib.rs` — Added conditional `mod generated_coefficients` and `mod interpolation` (excluded with `lightweight` feature)
- `src/solver.rs` — Added ~15-line `#[cfg(not(feature = "lightweight"))]` branch in `compute_initial_guess` that tries interpolation before falling back to analytical guess

**Key design decisions:**
- Interpolation is **on by default** — the `lightweight` feature excludes it
- Coefficients are compiled into the binary as const arrays — zero runtime file I/O
- Falls back gracefully to analytical initial guess when inputs are outside interpolation domain
- Multi-rev uses 3D interpolation (tau × tof × ln(N)) with kbot polynomial for minimum-TOF curve
- All polynomial evaluators are direct ports from Fortran with identical coefficient ordering

**User testing needed:**
- `DYLD_LIBRARY_PATH=fortran cargo test --features gooding-ffi` — Gooding cross-validation with interpolation
- `DYLD_LIBRARY_PATH=fortran cargo test --features gooding-ffi,ivlam-ffi` — Full three-way cross-validation
