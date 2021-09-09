extern crate configparser;
use configparser::ini::Ini;

use std::fmt::Debug;

use std::error::Error;
use std::str::FromStr;

pub fn open_config(cfg_path: &str) -> Result<Ini, Box<dyn Error>> {
    // Open the `config.ini` file and extract the contents
    let mut config = Ini::new();
    config.load(cfg_path)?;

    Ok(config)
}

pub fn parse_config_array<T>(
    config: &Ini,
    section: &str,
    key: &str,
) -> Result<Vec<T>, Box<dyn Error>>
where
    T: std::str::FromStr,
    <T as FromStr>::Err: Debug,
{
    // Parse the config.ini file for its specified key value pair.
    // This code works for parsing a vector of values only.
    let config_entry = &config
        .get(section, key)
        .expect("Invalid section/key pair in `config.ini`.")[..]; // Convert to string literal

    let mut vals: Vec<T> = Vec::new();

    let values = config_entry.split(',');
    for v in values {
        let temp: T = v.trim().parse::<T>().unwrap();
        vals.push(temp);
    }

    Ok(vals)
}

pub fn parse_config_bool(config: &Ini, section: &str, key: &str) -> Result<bool, Box<dyn Error>> {
    // Parse the config.ini file for its specified key value pair.
    // This code works for parsing boolean values only.
    let config_entry = &config
        .get(section, key)
        .expect("Invalid section/key pair in config.ini.")[..]; // convert to string literal &str

    let aliases_true: Vec<&str> = vec!["true", "True", "t", "yes", "Yes", "y"];
    let aliases_false: Vec<&str> = vec!["false", "False", "f", "no", "No", "n"];

    let val: bool = match config_entry {
        config_entry if aliases_true.iter().any(|&i| i == config_entry) => true,
        config_entry if aliases_false.iter().any(|&i| i == config_entry) => false,
        _ => {
            println!(
                r"WARNING: The entry for `{}, {}` is malformed. Setting value to `false` by default.
            Please specify one of {:?} for `true`,
            or one of {:?} for `false`.",
                section, key, aliases_true, aliases_false
            );
            false
        }
    };

    Ok(val)
}

pub fn parse_config_int<T>(config: &Ini, section: &str, key: &str) -> Result<T, Box<dyn Error>>
where
    T: std::str::FromStr,
    <T as FromStr>::Err: Debug,
{
    // Parse the config.ini file for its specified key value pair.
    // This code works for parsing i32 values only.
    let config_entry = &config
        .get(section, key)
        .expect("Invalid section/key pair in config.ini.")[..]; // convert to string literal &str

    let val: T = config_entry.parse().unwrap();
    Ok(val)
}