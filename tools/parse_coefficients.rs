//! Standalone tool to parse the Fortran unformatted binary coefficient file
//! and generate `src/generated_coefficients.rs` with all data as const arrays.
//!
//! Usage: cargo run --bin parse_coefficients
//!
//! Input:  fortran/ivLamTree_20210202_160219_i2d8.bin
//! Output: src/generated_coefficients.rs

use std::io::{self, Read};
use std::fs::File;
use std::path::Path;

/// Read a Fortran unformatted sequential record.
/// Each record is: [4-byte i32 len][data bytes][4-byte i32 len]
fn read_record(f: &mut File) -> io::Result<Vec<u8>> {
    let mut buf4 = [0u8; 4];

    f.read_exact(&mut buf4)?;
    let len_prefix = i32::from_le_bytes(buf4) as usize;

    let mut data = vec![0u8; len_prefix];
    f.read_exact(&mut data)?;

    f.read_exact(&mut buf4)?;
    let len_suffix = i32::from_le_bytes(buf4) as usize;

    assert_eq!(len_prefix, len_suffix, "Record framing mismatch");
    Ok(data)
}

fn read_i16_le(data: &[u8], offset: usize) -> i16 {
    i16::from_le_bytes([data[offset], data[offset + 1]])
}

fn read_i32_le(data: &[u8], offset: usize) -> i32 {
    i32::from_le_bytes([data[offset], data[offset + 1], data[offset + 2], data[offset + 3]])
}

fn read_f64_le(data: &[u8], offset: usize) -> f64 {
    f64::from_le_bytes([
        data[offset], data[offset + 1], data[offset + 2], data[offset + 3],
        data[offset + 4], data[offset + 5], data[offset + 6], data[offset + 7],
    ])
}

fn fmt_f64(v: f64) -> String {
    // Use enough precision to round-trip f64
    format!("{:.17e}", v)
}

fn main() {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")
        .unwrap_or_else(|_| ".".to_string());
    let bin_path = Path::new(&manifest_dir).join("fortran/ivLamTree_20210202_160219_i2d8.bin");
    let out_path = Path::new(&manifest_dir).join("src/generated_coefficients.rs");

    eprintln!("Reading: {}", bin_path.display());
    let mut f = File::open(&bin_path).expect("Cannot open .bin file");

    // ---- Record 1: precision markers ----
    let rec1 = read_record(&mut f).unwrap();
    assert_eq!(rec1.len(), 8, "Record 1 should be 2 x i32");
    let intknd = read_i32_le(&rec1, 0);
    let realknd = read_i32_le(&rec1, 4);
    eprintln!("Record 1: intknd={}, realknd={}", intknd, realknd);
    assert_eq!(intknd, 2, "Expected int kind = 2 (i16)");
    assert_eq!(realknd, 8, "Expected real kind = 8 (f64)");

    // ---- Record 2: zero-rev metadata ----
    let rec2 = read_record(&mut f).unwrap();
    assert_eq!(rec2.len(), 8, "Record 2 should be 4 x i16");
    let zrev_ncofs = read_i16_le(&rec2, 0) as usize;
    let zrev_ipatches = read_i16_le(&rec2, 2) as usize;
    let zrev_nx = read_i16_le(&rec2, 4) as usize;
    let zrev_ny = read_i16_le(&rec2, 6) as usize;
    eprintln!("Record 2: ncofs={}, ipatches={}, nx={}, ny={}", zrev_ncofs, zrev_ipatches, zrev_nx, zrev_ny);

    // ---- Record 3: zero-rev pointer array + patch data ----
    let rec3 = read_record(&mut f).unwrap();
    let expected_point_bytes = zrev_nx * zrev_ny * 2; // i16 each
    let expected_data_bytes = (4 + zrev_ncofs) * zrev_ipatches * 8; // f64 each
    assert_eq!(rec3.len(), expected_point_bytes + expected_data_bytes,
        "Record 3 size mismatch: got {}, expected {}", rec3.len(), expected_point_bytes + expected_data_bytes);

    // Parse pointer array (column-major: zRevPoint(0:nx-1, 0:ny-1))
    // Fortran column-major: fastest index is first (x varies fastest)
    let mut zrev_point = vec![0i16; zrev_nx * zrev_ny];
    for i in 0..(zrev_nx * zrev_ny) {
        zrev_point[i] = read_i16_le(&rec3, i * 2);
    }

    // Parse patch data: zRevData(4+ncofs, ipatches), column-major
    // Each patch q has: [xcorn_x, xcorn_y, inv_dx, inv_dy, coefs[ncofs]]
    let data_offset = expected_point_bytes;
    let patch_stride = 4 + zrev_ncofs; // number of f64 per patch
    let mut zrev_data: Vec<Vec<f64>> = Vec::with_capacity(zrev_ipatches);
    for q in 0..zrev_ipatches {
        let mut patch = Vec::with_capacity(patch_stride);
        for j in 0..patch_stride {
            // Column-major: element (j, q) is at offset (j + q * patch_stride) * 8
            let off = data_offset + (j + q * patch_stride) * 8;
            patch.push(read_f64_le(&rec3, off));
        }
        zrev_data.push(patch);
    }

    // ---- Record 4: zero-rev bounds ----
    let rec4 = read_record(&mut f).unwrap();
    let mut off = 0usize;

    // izoneBoundCustom(2,1) - 2 i16
    let _izbc_1 = read_i16_le(&rec4, off); off += 2;
    let _izbc_2 = read_i16_le(&rec4, off); off += 2;

    // dzoneBoundCustom(2,1) - 2 f64
    let dzbc_1 = read_f64_le(&rec4, off); off += 8;
    let dzbc_2 = read_f64_le(&rec4, off); off += 8;

    // smallDX, largeDX - 2 f64
    let small_dx = read_f64_le(&rec4, off); off += 8;
    let large_dx = read_f64_le(&rec4, off); off += 8;

    // numberCustomXbins, numberYbins - 2 i16
    let _nx_check = read_i16_le(&rec4, off); off += 2;
    let _ny_check = read_i16_le(&rec4, off); off += 2;

    // addChunk(2) - 2 i16
    let add_chunk_1 = read_i16_le(&rec4, off); off += 2;
    let add_chunk_2 = read_i16_le(&rec4, off); off += 2;

    // order, ncofs, ipatches - 3 i16
    let _order = read_i16_le(&rec4, off); off += 2;
    let _ncofs_check = read_i16_le(&rec4, off); off += 2;
    let _ipatches_check = read_i16_le(&rec4, off); off += 2;

    // Xlow(2) - 2 f64
    let zrev_xlow_1 = read_f64_le(&rec4, off); off += 8;
    let zrev_xlow_2 = read_f64_le(&rec4, off); off += 8;

    // oneByDeltX(2) - 2 f64
    let zrev_obdx_1 = read_f64_le(&rec4, off); off += 8;
    let zrev_obdx_2 = read_f64_le(&rec4, off); off += 8;

    // Xhi(2) - 2 f64
    let zrev_xhi_1 = read_f64_le(&rec4, off); off += 8;
    let zrev_xhi_2 = read_f64_le(&rec4, off); off += 8;

    // maxval - 1 f64
    let zrev_maxval = read_f64_le(&rec4, off);

    eprintln!("Zero-rev bounds: Xlow=[{}, {}], Xhi=[{}, {}], maxval={}",
        zrev_xlow_1, zrev_xlow_2, zrev_xhi_1, zrev_xhi_2, zrev_maxval);
    eprintln!("  smallDX={}, largeDX={}, addChunk=[{}, {}]", small_dx, large_dx, add_chunk_1, add_chunk_2);
    eprintln!("  dzoneBoundCustom=[{}, {}]", dzbc_1, dzbc_2);

    // ---- Record 5: multi-rev metadata ----
    let rec5 = read_record(&mut f).unwrap();
    // ncofs, ipatchesM, kbotNcofs, numEachDir(3) - total 6 i16
    assert_eq!(rec5.len(), 12, "Record 5 should be 6 x i16");
    let mrev_ncofs = read_i16_le(&rec5, 0) as usize;
    let mrev_ipatches = read_i16_le(&rec5, 2) as usize;
    let mrev_kbot_ncofs = read_i16_le(&rec5, 4) as usize;
    let mrev_n1 = read_i16_le(&rec5, 6) as usize;
    let mrev_n2 = read_i16_le(&rec5, 8) as usize;
    let mrev_n3 = read_i16_le(&rec5, 10) as usize;
    eprintln!("Record 5: mrev ncofs={}, ipatches={}, kbot_ncofs={}, numEachDir=[{},{},{}]",
        mrev_ncofs, mrev_ipatches, mrev_kbot_ncofs, mrev_n1, mrev_n2, mrev_n3);

    // ---- Record 6: multi-rev pointer + patch + kbot data ----
    let rec6 = read_record(&mut f).unwrap();
    let mrev_point_bytes = mrev_n1 * mrev_n2 * mrev_n3 * 2;
    let mrev_data_bytes = (6 + mrev_ncofs) * mrev_ipatches * 8;
    let mrev_kbot_bytes = mrev_kbot_ncofs * mrev_n1 * mrev_n3 * 8;
    assert_eq!(rec6.len(), mrev_point_bytes + mrev_data_bytes + mrev_kbot_bytes,
        "Record 6 size mismatch: got {}, expected {}",
        rec6.len(), mrev_point_bytes + mrev_data_bytes + mrev_kbot_bytes);

    // Parse mrev pointer array: mrevPoint(0:n1-1, 0:n2-1, 0:n3-1), column-major i16
    let mut mrev_point = vec![0i16; mrev_n1 * mrev_n2 * mrev_n3];
    for i in 0..(mrev_n1 * mrev_n2 * mrev_n3) {
        mrev_point[i] = read_i16_le(&rec6, i * 2);
    }

    // Parse mrev patch data: mrevData(6+ncofs, ipatchesM), column-major f64
    let mrev_patch_stride = 6 + mrev_ncofs;
    let mrev_data_start = mrev_point_bytes;
    let mut mrev_data: Vec<Vec<f64>> = Vec::with_capacity(mrev_ipatches);
    for q in 0..mrev_ipatches {
        let mut patch = Vec::with_capacity(mrev_patch_stride);
        for j in 0..mrev_patch_stride {
            let off = mrev_data_start + (j + q * mrev_patch_stride) * 8;
            patch.push(read_f64_le(&rec6, off));
        }
        mrev_data.push(patch);
    }

    // Parse mrev kbot data: mrevKbotData(kbotNcofs, 0:n1-1, 0:n3-1), column-major f64
    // Total elements: kbotNcofs * n1 * n3
    let kbot_start = mrev_data_start + mrev_data_bytes;
    let kbot_total = mrev_kbot_ncofs * mrev_n1 * mrev_n3;
    let mut mrev_kbot_flat: Vec<f64> = Vec::with_capacity(kbot_total);
    for i in 0..kbot_total {
        mrev_kbot_flat.push(read_f64_le(&rec6, kbot_start + i * 8));
    }

    // ---- Record 7: multi-rev bounds ----
    let rec7 = read_record(&mut f).unwrap();
    let mut off = 0usize;

    // numEachDir(3) - 3 i16
    let _n1_check = read_i16_le(&rec7, off); off += 2;
    let _n2_check = read_i16_le(&rec7, off); off += 2;
    let _n3_check = read_i16_le(&rec7, off); off += 2;

    // order, ncofs, ipatchesM - 3 i16
    let _mrev_order = read_i16_le(&rec7, off); off += 2;
    let _mrev_ncofs_check = read_i16_le(&rec7, off); off += 2;
    let _mrev_ipatches_check = read_i16_le(&rec7, off); off += 2;

    // xlow(3) - 3 f64
    let mrev_xlow = [
        read_f64_le(&rec7, off), { off += 8; read_f64_le(&rec7, off) }, { off += 8; read_f64_le(&rec7, off) },
    ];
    off += 8;

    // oneByDeltX(3) - 3 f64
    let mrev_obdx = [
        read_f64_le(&rec7, off), { off += 8; read_f64_le(&rec7, off) }, { off += 8; read_f64_le(&rec7, off) },
    ];
    off += 8;

    // lims(2,3) - 6 f64, column-major: lims(1,1), lims(2,1), lims(1,2), lims(2,2), lims(1,3), lims(2,3)
    let mut mrev_lims = [[0.0f64; 3]; 2]; // [lo/hi][dim]
    for dim in 0..3 {
        mrev_lims[0][dim] = read_f64_le(&rec7, off); off += 8;
        mrev_lims[1][dim] = read_f64_le(&rec7, off); off += 8;
    }

    // maxval - f64
    let _mrev_maxval = read_f64_le(&rec7, off); off += 8;

    // kbotOrd, kbotNcofs - 2 i16
    let _kbot_ord = read_i16_le(&rec7, off); off += 2;
    let _kbot_ncofs_check = read_i16_le(&rec7, off); off += 2;

    // kbotMaxval - f64
    let _kbot_maxval = read_f64_le(&rec7, off);

    eprintln!("Multi-rev bounds: xlow={:?}, lims_lo={:?}, lims_hi={:?}",
        mrev_xlow, mrev_lims[0], mrev_lims[1]);

    // ===== Generate Rust source =====
    let mut out = String::with_capacity(4 * 1024 * 1024);
    out.push_str("// Auto-generated from ivLamTree_20210202_160219_i2d8.bin\n");
    out.push_str("// DO NOT EDIT — regenerate with: cargo run --bin parse_coefficients\n\n");
    out.push_str("#![allow(clippy::excessive_precision)]\n");
    out.push_str("#![allow(clippy::unreadable_literal)]\n\n");

    // Zero-rev constants
    out.push_str(&format!("pub const ZREV_NCOFS: usize = {};\n", zrev_ncofs));
    out.push_str(&format!("pub const ZREV_IPATCHES: usize = {};\n", zrev_ipatches));
    out.push_str(&format!("pub const ZREV_NX_BINS: usize = {};\n", zrev_nx));
    out.push_str(&format!("pub const ZREV_NY_BINS: usize = {};\n", zrev_ny));
    out.push_str(&format!("pub const ZREV_XLOW: [f64; 2] = [{}, {}];\n", fmt_f64(zrev_xlow_1), fmt_f64(zrev_xlow_2)));
    out.push_str(&format!("pub const ZREV_ONE_BY_DELT_X: [f64; 2] = [{}, {}];\n", fmt_f64(zrev_obdx_1), fmt_f64(zrev_obdx_2)));
    out.push_str(&format!("pub const ZREV_XHI: [f64; 2] = [{}, {}];\n", fmt_f64(zrev_xhi_1), fmt_f64(zrev_xhi_2)));
    out.push_str(&format!("pub const ZREV_SMALL_DX: f64 = {};\n", fmt_f64(small_dx)));
    out.push_str(&format!("pub const ZREV_LARGE_DX: f64 = {};\n", fmt_f64(large_dx)));
    out.push_str(&format!("pub const ZREV_DZONE_BOUND_CUSTOM: [f64; 2] = [{}, {}];\n", fmt_f64(dzbc_1), fmt_f64(dzbc_2)));
    out.push_str(&format!("pub const ZREV_ADD_CHUNK: [i16; 2] = [{}, {}];\n", add_chunk_1, add_chunk_2));
    out.push_str("\n");

    // Zero-rev pointer array (flat, column-major)
    out.push_str(&format!("pub const ZREV_POINT: &[i16; {}] = &[\n", zrev_nx * zrev_ny));
    for (i, &v) in zrev_point.iter().enumerate() {
        if i % 20 == 0 { out.push_str("    "); }
        out.push_str(&format!("{},", v));
        if i % 20 == 19 || i == zrev_point.len() - 1 { out.push('\n'); }
    }
    out.push_str("];\n\n");

    // Zero-rev patch data: each patch has 4+ncofs f64
    out.push_str(&format!("pub const ZREV_DATA: &[[f64; {}]; {}] = &[\n", patch_stride, zrev_ipatches));
    for patch in &zrev_data {
        out.push_str("    [");
        for (j, &v) in patch.iter().enumerate() {
            out.push_str(&fmt_f64(v));
            if j < patch.len() - 1 { out.push_str(", "); }
        }
        out.push_str("],\n");
    }
    out.push_str("];\n\n");

    // Multi-rev constants
    out.push_str(&format!("pub const MREV_NCOFS: usize = {};\n", mrev_ncofs));
    out.push_str(&format!("pub const MREV_IPATCHES: usize = {};\n", mrev_ipatches));
    out.push_str(&format!("pub const MREV_KBOT_NCOFS: usize = {};\n", mrev_kbot_ncofs));
    out.push_str(&format!("pub const MREV_NUM_EACH_DIR: [usize; 3] = [{}, {}, {}];\n", mrev_n1, mrev_n2, mrev_n3));
    out.push_str(&format!("pub const MREV_XLOW: [f64; 3] = [{}, {}, {}];\n",
        fmt_f64(mrev_xlow[0]), fmt_f64(mrev_xlow[1]), fmt_f64(mrev_xlow[2])));
    out.push_str(&format!("pub const MREV_ONE_BY_DELT_X: [f64; 3] = [{}, {}, {}];\n",
        fmt_f64(mrev_obdx[0]), fmt_f64(mrev_obdx[1]), fmt_f64(mrev_obdx[2])));
    out.push_str(&format!("pub const MREV_LIMS: [[f64; 3]; 2] = [\n    [{}, {}, {}],\n    [{}, {}, {}],\n];\n",
        fmt_f64(mrev_lims[0][0]), fmt_f64(mrev_lims[0][1]), fmt_f64(mrev_lims[0][2]),
        fmt_f64(mrev_lims[1][0]), fmt_f64(mrev_lims[1][1]), fmt_f64(mrev_lims[1][2])));
    out.push_str("\n");

    // Multi-rev pointer array (flat, column-major)
    let mrev_point_len = mrev_n1 * mrev_n2 * mrev_n3;
    out.push_str(&format!("pub const MREV_POINT: &[i16; {}] = &[\n", mrev_point_len));
    for (i, &v) in mrev_point.iter().enumerate() {
        if i % 20 == 0 { out.push_str("    "); }
        out.push_str(&format!("{},", v));
        if i % 20 == 19 || i == mrev_point.len() - 1 { out.push('\n'); }
    }
    out.push_str("];\n\n");

    // Multi-rev patch data
    out.push_str(&format!("pub const MREV_DATA: &[[f64; {}]; {}] = &[\n", mrev_patch_stride, mrev_ipatches));
    for patch in &mrev_data {
        out.push_str("    [");
        for (j, &v) in patch.iter().enumerate() {
            out.push_str(&fmt_f64(v));
            if j < patch.len() - 1 { out.push_str(", "); }
        }
        out.push_str("],\n");
    }
    out.push_str("];\n\n");

    // Multi-rev kbot data: grouped as [kbot_ncofs] per (xi, zi)
    // Fortran layout: mrevKbotData(kbotNcofs, 0:n1-1, 0:n3-1), column-major
    // We store it flat as a 1D array, indexed by (xi * n3 + zi) * kbot_ncofs
    // But actually Fortran column-major means the layout is:
    //   for zi in 0..n3: for xi in 0..n1: for coef in 0..kbot_ncofs: data[coef + xi*kbot_ncofs + zi*n1*kbot_ncofs]
    // We'll store as &[[f64; kbot_ncofs]; n1 * n3] with row = xi + zi * n1 (Fortran column-major on outer two dims)
    let kbot_entries = mrev_n1 * mrev_n3;
    out.push_str(&format!("/// Kbot data indexed by (xi + zi * MREV_NUM_EACH_DIR[0]), Fortran column-major\n"));
    out.push_str(&format!("pub const MREV_KBOT_DATA: &[[f64; {}]; {}] = &[\n", mrev_kbot_ncofs, kbot_entries));
    for entry_idx in 0..kbot_entries {
        // In Fortran column-major with dims (kbotNcofs, n1, n3):
        // The linear index for element (coef, xi, zi) = coef + xi*kbotNcofs + zi*n1*kbotNcofs
        // entry_idx corresponds to (xi, zi) pair where entry_idx = xi + zi * n1
        let zi = entry_idx / mrev_n1;
        let xi = entry_idx % mrev_n1;
        out.push_str("    [");
        for coef in 0..mrev_kbot_ncofs {
            let flat_idx = coef + xi * mrev_kbot_ncofs + zi * mrev_n1 * mrev_kbot_ncofs;
            out.push_str(&fmt_f64(mrev_kbot_flat[flat_idx]));
            if coef < mrev_kbot_ncofs - 1 { out.push_str(", "); }
        }
        out.push_str("],\n");
    }
    out.push_str("];\n");

    // Write output
    std::fs::write(&out_path, &out).expect("Failed to write output file");
    eprintln!("Generated: {} ({} bytes)", out_path.display(), out.len());
}
