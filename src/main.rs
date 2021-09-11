include!(concat!(env!("OUT_DIR"), "/constants.rs"));

/*
Imports
*/

mod sim;

mod ini_parser;
use std::time::Instant;

use ini_parser::*;

use itertools::iproduct;

/*
Structs
*/

pub struct InputParams {
    n_particles: Vec<usize>,
    array_sizes: Vec<usize>,
    d_max_vals: Vec<u8>,
    seeds: Vec<usize>,
    init_seed: usize,
    write_data: bool,
    plot_tree: bool,
}

impl InputParams {
    fn new() -> InputParams {
        let config = open_config(CFG).expect("Failed to open config file.");

        let params: InputParams = InputParams {
            n_particles: parse_config_array(&config, "simulation_params", "n_particles")
                .expect("Failed to parse number of particles."),
            array_sizes: parse_config_array(&config, "simulation_params", "array_sizes")
                .expect("Failed to parse array sizes."),
            d_max_vals: parse_config_array(&config, "simulation_params", "d_max_vals")
                .expect("Failed to parse d_max values."),
            seeds: parse_config_array(&config, "simulation_params", "seeds")
                .expect("Failed to parse seeds."),
            init_seed: parse_config_int(&config, "options", "init_seed")
                .expect("Failed to parse initial random number seed."),
            write_data: parse_config_bool(&config, "options", "write_data")
                .expect("Failed to parse whether to write data to disk."),
            plot_tree: parse_config_bool(&config, "options", "plot_tree")
                .expect("Failed to parse whether to plot tree."),
        };
        params
    }
}

fn main() {
    // Parse in the params from the `.ini` file.
    let params: InputParams = InputParams::new();

    // Iterate through the params provided and run simulations
    for (n, a, d_max, max_seed) in iproduct!(
        &params.n_particles,
        &params.array_sizes,
        &params.d_max_vals,
        &params.seeds
    ) {
        println!(
            "Running sim for n={} a={} d_max={} max_seed={}",
            n, a, d_max, max_seed
        );

        let now = Instant::now();
        sim::run(*n, *a, *d_max, *max_seed, &params);
        let new_now = Instant::now();
        println!(
            r"Done! Total time = {:?}
                ",
            new_now.duration_since(now)
        );
    }
}
