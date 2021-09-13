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

    if params.write_g == false {
        setup_tree_sims(&params)?
    } else {
        setup_fractal_sims(&params)?
    }

    Ok(())
}

fn setup_fractal_sims(params: &InputParams) -> Result<(), Box<dyn Error>> {
    println!(r"Collecting gyration radii for particle numbers:
{:?}
", params.n_particles);
    
    // As we're collecting data across simulations, we have to collect the data here
    let n_particles: Vec<usize> = params.n_particles.clone();
    let mut g_radii: Vec<f64> = vec![0.0; n_particles.len()];
    
    // Hard-coded to run for d_max = 90, seeds = 1000, array_size = 40001.
    // Don't want to adjust this in config.
    let a: usize = 40001;
    let d_max: u8 = 90;
    let max_seed: usize = 1000;

    let mut g_radius: f64;
    let mut i: usize = 0;

    for n in n_particles.iter() {
        g_radius = sim::run_fractal(*n, a, d_max, max_seed, params)?;
        g_radii[i] = g_radius;
        i += 1;
    }

    // Now we've collected our data, we can save to disk
    write_g_radii(params, &n_particles, &g_radii, max_seed, d_max)?;

    Ok(())
}

fn setup_tree_sims(params: &InputParams) -> Result<(), Box<dyn Error>> {
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
        sim::run(*n, *a, *d_max, *max_seed, &params)?;
        let new_now = Instant::now();
        println!(
            r"Done! Total time = {:?}
                ",
            new_now.duration_since(now)
        );
    }

    Ok(())
}
