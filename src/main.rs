include!(concat!(env!("OUT_DIR"), "/constants.rs"));

/*
Imports
*/

mod sim;

mod ini_parser;
use std::time::Instant;

use ini_parser::*;

use itertools::iproduct;

use std::error::Error;

use crate::sim::writer::write_g_radii;

use rayon::prelude::*;

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
    write_tree: bool,
    write_g: bool,
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
            write_tree: parse_config_bool(&config, "options", "write_tree")
                .expect("Failed to parse whether to plot tree."),
            write_g: parse_config_bool(&config, "fractals", "write_g_radius")
                .expect("Failed to parse whether to write the gyration radius."),
        };
        params
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse in the params from the `.ini` file.
    let params: InputParams = InputParams::new();

    let mut sim_params: Vec<(usize, usize, u8, usize)> = Vec::new();

    for (n, a, d_max, max_seed) in iproduct!(
        &params.n_particles,
        &params.array_sizes,
        &params.d_max_vals,
        &params.seeds
    ) {
        sim_params.push((*n, *a, *d_max, *max_seed))
    }

    let results: Vec<(usize, f64, f64, usize, u8)> = sim_params
        .par_iter()
        .map(|(n, a, d_max, max_seed)| {
            let now: Instant = Instant::now();
            let (cpu_time, r_avg) = match sim::run(*n, *a, *d_max, *max_seed, &params) {
                Ok((cpu_time, r_avg)) => (cpu_time, r_avg),
                Err(_) => panic!("Performing the simulation has failed!"),
            };
            let new_now: Instant = Instant::now();
            println!(
                r"Sim for n={} a={} d_max={} max_seed={}, init_seed={},
Complete in {:?}
",
                n,
                a,
                d_max,
                max_seed,
                params.init_seed,
                new_now.duration_since(now)
            );

            (*n, cpu_time, r_avg, *max_seed, *d_max)
        })
        .collect();

    if params.write_g == true {
        write_g_radii(&results)?;
    }

    Ok(())
}
