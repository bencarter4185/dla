use std::{env, fs::File, io::Write, path::Path};

fn write_str(var: Option<&str>, mut f: &File, name: &str) {
    let var: &str = var.unwrap_or("./config.ini");

    let out: &str = &format!("const {}: &str = \"{}\";", name, var)[..];
    let out_bytes: &[u8] = out.as_bytes();

    f.write_all(out_bytes).expect("Could not write file");
    println!("cargo:rerun-if-env-changed={}", name);
}

// Not used unless I need to define a usize at runtime.
fn _write_usize(var: Option<&str>, mut f: &File, name: &str, d_max: bool) {
    let msg: &str = &format!("Could not parse {}", name)[..]; 
    
    let mut var: usize = var
        .map_or(Ok(10_000), str::parse)
        .expect(msg);

    if d_max == true {
        var = 2*var + 1;
    }

    let out: &str = &format!("const {}: usize = {};", name, var)[..];
    let out_bytes: &[u8] = out.as_bytes();

    f.write_all(out_bytes).expect("Could not write file");
    println!("cargo:rerun-if-env-changed={}", name);
}

fn main() {
    let out_dir = env::var("OUT_DIR").expect("No out dir");
    let dest_path = Path::new(&out_dir).join("constants.rs");
    let f = File::create(&dest_path).expect("Could not create file");

    let cfg = option_env!("CFG");
    write_str(cfg, &f, "CFG");
}
