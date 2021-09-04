/*
Library file for parsing config in `config.ini` and passing to main
*/

/*
Imports
*/

extern crate configparser;
use configparser::ini::Ini;

use std::error::Error;

/*
Structs
*/

pub struct SimulationParams {
    pub max_radius: i32,
    pub max_seed: i32,
    pub init_seed: i32,
    pub add_plot: bool,
    pub write_data: bool,
}

impl SimulationParams {
    pub fn new(
        total_params: InputParams,
        max_radius: i32,
        max_seed: i32,
    ) -> SimulationParams {
        let current_params: SimulationParams = SimulationParams {
            max_radius: max_radius,
            max_seed: max_seed,
            init_seed: total_params.init_seed,
            add_plot: total_params.add_plot,
            write_data: total_params.write_data,
        };
        current_params
    }
}

pub struct InputParams {
    pub substrate_radii: Vec<i32>,
    pub seeds: Vec<i32>,
    pub init_seed: i32,
    pub add_plot: bool,
    pub write_data: bool,
}

impl InputParams {
    pub fn new(displ_params: bool) -> InputParams {
        // Instantiate our parameters for the simulation
        let config = open_config().expect("Failed to open `config.ini`.");

        // Parse key/value pairs in our config
        let substrate_radii: Vec<i32> =
            parse_config_array(&config, "simulation_params", "substrate_radii")
                .expect("Failed to parse substrate radii.");

        let seeds: Vec<i32> = parse_config_array(&config, "simulation_params", "seeds")
            .expect("Failed to parse number of seeds.");

        let init_seed: i32 = parse_config_i32(&config, "options", "init_seed")
            .expect("Failed to parse initial random number seed.");

        let add_plot: bool = parse_config_bool(&config, "options", "add_plot")
            .expect("Failed to parse whether to add plot.");

        let write_data: bool = parse_config_bool(&config, "options", "write_data")
            .expect("Failed to parse whether to write data to file.");

        // Run-once to show the user the parameters about to be simulated
        if displ_params == true {
            println!(
                r"
        Parsed config from `./config/config_sim.ini`...
            Running simulations for parameters:
                Substrate Radii = {:?},
                Seeds = {:?},
                Initial Seed = {:?}
                Add plot: {:?},
                        ",
                substrate_radii, seeds, init_seed,  add_plot
            );
        };

        let params: InputParams = InputParams {
            substrate_radii: substrate_radii,
            seeds: seeds,
            init_seed: init_seed,
            add_plot: add_plot,
            write_data: write_data,
        };

        params
    }
}

pub fn open_config() -> Result<Ini, Box<dyn Error>> {
    // Open the `config.ini` file and extract the contents
    let mut config = Ini::new();
    config.load("./config.ini")?;

    Ok(config)
}

pub fn parse_config_array(
    config: &Ini,
    section: &str,
    key: &str,
) -> Result<Vec<i32>, Box<dyn Error>> {
    // Parse the config.ini file for its specified key value pair.
    // This code works for parsing a vector of values only.
    let config_entry = &config
        .get(section, key)
        .expect("Invalid section/key pair in `config.ini`.")[..]; // Convert to string literal

    let mut vals: Vec<i32> = Vec::new();

    let values = config_entry.split(',');
    for v in values {
        let temp: i32 = v.trim().parse().unwrap_or(0);
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

pub fn parse_config_i32(config: &Ini, section: &str, key: &str) -> Result<i32, Box<dyn Error>> {
    // Parse the config.ini file for its specified key value pair.
    // This code works for parsing i32 values only.
    let config_entry = &config
        .get(section, key)
        .expect("Invalid section/key pair in config.ini.")[..]; // convert to string literal &str

    let val: i32 = config_entry.parse().unwrap();
    Ok(val)
}
