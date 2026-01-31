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

### 3.1 Analysis and decision
- [ ] Confirm approach: Option A+C (feature-gated const arrays)

### 3.2 Write Fortran .bin parser
- [ ] Create `tools/parse_coefficients.rs`

### 3.3 Implement interpolation-based initial guess
- [ ] Create `src/interpolation.rs`
- [ ] Create `src/coefficients.rs` (generated)
- [ ] Modify `src/solver.rs` with `#[cfg(feature = "interpolation")]` branch
- [ ] Update `src/lib.rs` and `Cargo.toml`

### 3.4 Validate interpolation accuracy
- [ ] Create `tests/interpolation_validation.rs`
