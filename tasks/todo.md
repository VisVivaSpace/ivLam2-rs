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

## Review

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
