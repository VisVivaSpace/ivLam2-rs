#!/bin/bash
# Build the ivLam2 Fortran shared library for FFI testing.
# Usage: cd fortran && ./build.sh

set -e

FFLAGS="-shared -fPIC -O2 -ffree-line-length-none"
SRC="ivLamRuntimeV2p50_739981p74401.f90"
WRAPPER="ivlam_c_wrapper.f90"
OUTPUT="libivlam.dylib"

echo "Building $OUTPUT..."
gfortran $FFLAGS -o $OUTPUT $SRC $WRAPPER

echo "Built $OUTPUT successfully."
ls -la $OUTPUT
