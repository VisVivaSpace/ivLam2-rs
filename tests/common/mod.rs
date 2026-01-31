//! Shared test utilities for Lambert solver integration tests.

pub mod kepler;

#[cfg(feature = "gooding-ffi")]
pub mod gooding;

#[cfg(feature = "ivlam-ffi")]
pub mod ivlam_ffi;
