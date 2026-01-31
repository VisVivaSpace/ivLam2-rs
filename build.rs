fn main() {
    #[cfg(feature = "gooding-ffi")]
    {
        cc::Build::new()
            .file("csrc/lamb.c")
            .warnings(false)
            .compile("gooding_lambert");
    }

    if std::env::var("CARGO_FEATURE_IVLAM_FFI").is_ok() {
        // Link to the pre-built Fortran shared library
        let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
        let fortran_dir = format!("{}/fortran", manifest_dir);
        println!("cargo:rustc-link-search=native={}", fortran_dir);
        println!("cargo:rustc-link-lib=dylib=ivlam");

        // Find gfortran runtime library path
        let gfortran_path = std::process::Command::new("gfortran")
            .arg("--print-file-name=libgfortran.dylib")
            .output()
            .expect("gfortran not found — install via `brew install gcc`");
        let gfortran_file = String::from_utf8(gfortran_path.stdout).unwrap();
        let gfortran_dir = std::path::Path::new(gfortran_file.trim())
            .parent()
            .expect("cannot determine gfortran lib dir")
            .canonicalize()
            .expect("gfortran lib dir does not exist");
        println!("cargo:rustc-link-search=native={}", gfortran_dir.display());
        println!("cargo:rustc-link-lib=dylib=gfortran");
    }
}
