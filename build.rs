fn main() {
    #[cfg(feature = "gooding-ffi")]
    {
        cc::Build::new()
            .file("csrc/lamb.c")
            .warnings(false)
            .compile("gooding_lambert");
    }
}
