//! Benchmarks comparing Rust vercosine Lambert solver vs C Gooding solver.
//!
//! Run with: cargo bench --features gooding-ffi

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lambert_solver::{solve_lambert, Direction};
use std::f64::consts::PI;

#[cfg(feature = "gooding-ffi")]
mod gooding_ffi {
    extern "C" {
        pub fn lambert(
            gm: f64,
            r1: *const f64,
            r2: *const f64,
            nrev: i32,
            dt: f64,
            v1: *mut f64,
            v2: *mut f64,
        ) -> i32;
    }

    pub fn gooding_solve(mu: f64, r1: &[f64; 3], r2: &[f64; 3], tof: f64) -> ([f64; 3], [f64; 3]) {
        let mut v1 = [0.0f64; 3];
        let mut v2 = [0.0f64; 3];
        unsafe {
            lambert(mu, r1.as_ptr(), r2.as_ptr(), 0, tof, v1.as_mut_ptr(), v2.as_mut_ptr());
        }
        (v1, v2)
    }
}

fn bench_single_solve_90_deg(c: &mut Criterion) {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    let tof = PI / 2.0;
    let mu = 1.0;

    let mut group = c.benchmark_group("single_90deg");

    group.bench_function("rust_vercosine", |b| {
        b.iter(|| {
            solve_lambert(
                black_box(&r1), black_box(&r2), black_box(tof), black_box(mu),
                Direction::Prograde, 0,
            ).unwrap()
        })
    });

    #[cfg(feature = "gooding-ffi")]
    group.bench_function("c_gooding", |b| {
        b.iter(|| {
            gooding_ffi::gooding_solve(black_box(mu), black_box(&r1), black_box(&r2), black_box(tof))
        })
    });

    group.finish();
}

fn bench_batch_angle_sweep(c: &mut Criterion) {
    let r1 = [1.0, 0.0, 0.0];
    let mu = 1.0;

    // Pre-compute 100 test cases: angle sweep from 5° to 170°
    let cases: Vec<([f64; 3], f64)> = (0..100)
        .map(|i| {
            let angle = (5.0 + 165.0 * i as f64 / 99.0).to_radians();
            let r2 = [angle.cos(), angle.sin(), 0.0];
            (r2, angle) // TOF = circular orbit time for the angle
        })
        .collect();

    let mut group = c.benchmark_group("batch_100_angles");

    group.bench_function("rust_vercosine", |b| {
        b.iter(|| {
            for (r2, tof) in &cases {
                let _ = solve_lambert(
                    black_box(&r1), black_box(r2), black_box(*tof), black_box(mu),
                    Direction::Prograde, 0,
                );
            }
        })
    });

    #[cfg(feature = "gooding-ffi")]
    group.bench_function("c_gooding", |b| {
        b.iter(|| {
            for (r2, tof) in &cases {
                gooding_ffi::gooding_solve(black_box(mu), black_box(&r1), black_box(r2), black_box(*tof));
            }
        })
    });

    group.finish();
}

fn bench_physical_units_leo_geo(c: &mut Criterion) {
    let mu = 398600.4418;
    let r1 = [6678.0, 0.0, 0.0];
    let r2 = [0.0, 42164.0, 0.0];
    let tof = 5.0 * 3600.0;

    let mut group = c.benchmark_group("leo_to_geo");

    group.bench_function("rust_vercosine", |b| {
        b.iter(|| {
            solve_lambert(
                black_box(&r1), black_box(&r2), black_box(tof), black_box(mu),
                Direction::Prograde, 0,
            ).unwrap()
        })
    });

    #[cfg(feature = "gooding-ffi")]
    group.bench_function("c_gooding", |b| {
        b.iter(|| {
            gooding_ffi::gooding_solve(black_box(mu), black_box(&r1), black_box(&r2), black_box(tof))
        })
    });

    group.finish();
}

fn bench_hyperbolic(c: &mut Criterion) {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    let tof = 0.1;
    let mu = 1.0;

    let mut group = c.benchmark_group("hyperbolic");

    group.bench_function("rust_vercosine", |b| {
        b.iter(|| {
            solve_lambert(
                black_box(&r1), black_box(&r2), black_box(tof), black_box(mu),
                Direction::Prograde, 0,
            ).unwrap()
        })
    });

    #[cfg(feature = "gooding-ffi")]
    group.bench_function("c_gooding", |b| {
        b.iter(|| {
            gooding_ffi::gooding_solve(black_box(mu), black_box(&r1), black_box(&r2), black_box(tof))
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_single_solve_90_deg,
    bench_batch_angle_sweep,
    bench_physical_units_leo_geo,
    bench_hyperbolic,
);
criterion_main!(benches);
