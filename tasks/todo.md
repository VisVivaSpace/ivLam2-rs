# ivLam2-rs: Testing & First Derivative Implementation

## Phase 0: Documentation & Attribution Updates
- [x] 0.1 Update CLAUDE.md — add PDF reference to Notes Directory
- [x] 0.2 Update README.md — add UT Austin link in Acknowledgments

## Phase 1: Fix Failing Tests & Improve Coverage
- [x] 1.1 Fix `test_w_ellipse_basic` — near-zero series coefficients A0–A8 were 2x too large; also fixed k^7 coefficient (-8/35 not -2/7)
- [x] 1.2 Fix `test_solve_lambert_hyperbolic` — clamp k_initial to [k_left, k_right] after bounds computed
- [x] 1.3 Tighten `test_solve_lambert_90_degree` tolerance from ±0.1 to 1e-10
- [x] 1.4 Add unit tests per module (stumpff: 8 new, geometry: 6 new, solver: 4 new, velocity: 2 new)
- [x] 1.5 Round-trip integration tests (tests/round_trip.rs with Kepler propagator — 8 tests)
- [ ] 1.6 Planetary data tests — deferred (requires SPICE kernels)

## Phase 2: C Gooding Solver Cross-Validation
- [ ] 2.1 Set up C FFI build (csrc/, build.rs, gooding-ffi feature)
- [ ] 2.2 Safe Rust wrapper (tests/common/gooding.rs)
- [ ] 2.3 Systematic cross-validation tests (tests/cross_validate_gooding.rs)

## Phase 3: Fortran ivLam Cross-Validation
- [ ] 3.1 Build Fortran shared library (fortran/, ivlam-ffi feature)
- [ ] 3.2 Safe Rust wrappers (tests/common/ivlam_ffi.rs)
- [ ] 3.3 Three-way cross-validation tests
- [ ] 3.4 Capture derivative reference values

## Phase 4: First Derivative Support
- [ ] 4.1 GeometryPartials struct in geometry.rs
- [ ] 4.2 dF_dtau helper in solver.rs
- [ ] 4.3 Implement compute_first_order in sensitivities.rs
- [ ] 4.4 Add solve_lambert_with_jacobian to solver.rs
- [ ] 4.5 Export new API from lib.rs
- [ ] 4.6 Validation tests (finite-diff, Fortran comparison, symmetry, energy)

## Phase 5: Documentation & README Update
- [ ] 5.1 Update README.md with new API docs
- [ ] 5.2 Update CLAUDE.md with new architecture
- [ ] 5.3 Update lib.rs doc comments
