include!(concat!(env!("OUT_DIR"), "/constants.rs"));

/*
Imports
*/

mod sim;

mod ini_parser;
use std::time::Instant;

use ini_parser::*;

use itertools::iproduct;
use sim::writer::erase_g_radii;

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
    parallel: u8,
}

impl InputParams {
    fn new() -> InputParams {
        let config = open_config(CFG).expect("Failed to open config file.");

        let sim_parallel: u8 = parse_config_bool(&config, "parallelisation", "sim_parallel")
            .expect("Failed to parse whether to run sim in parallel.")
            as u8;
        let seed_parallel: u8 = parse_config_bool(&config, "parallelisation", "seed_parallel")
            .expect("Failed to parse whether to run seeds in parallel.")
            as u8
            * 2;

        let parallel = sim_parallel + seed_parallel;

        match parallel {
            0..=2 => (),
            3 => panic!("Sim parallel and seed parallel cannot both be true!"),
            _ => panic!("Sim parallel and seed parallel is malformed!"),
        };

        let params: InputParams = InputParams {
            n_particles: parse_config_array(&config, "simulation_params", "n_particles")
                .expect("Failed to parse number of particles."),
            array_sizes: parse_config_array(&config, "simulation_params", "array_sizes")
                .expect("Failed to parse array sizes."),
            d_max_vals: parse_config_array(&config, "simulation_params", "d_max_vals")
                .expect("Failed to parse d_max values."),
            seeds: parse_config_array(&config, "simulation_params", "seeds")
                .expect("Failed to parse seeds."),
            init_seed: parse_config_int(&config, "misc_options", "init_seed")
                .expect("Failed to parse initial random number seed."),
            write_data: parse_config_bool(&config, "misc_options", "write_data")
                .expect("Failed to parse whether to write data to disk."),
            write_tree: parse_config_bool(&config, "misc_options", "write_tree")
                .expect("Failed to parse whether to plot tree."),
            write_g: parse_config_bool(&config, "misc_options", "write_g_radius")
                .expect("Failed to parse whether to write the gyration radius."),
            parallel: parallel,
        };
        params
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse in the params from the `.ini` file.
    let params: InputParams = InputParams::new();

    // Erase g_radii if exists
    if params.write_g == true {
        erase_g_radii(&params)?
    }
    
    // Run parallel sims if asked to, else run single-threaded
    match params.parallel {
        1 => setup_parallel(&params)?,
        _ => setup(&params)?,
    }

    Ok(())
}

fn setup(params: &InputParams) -> Result<(), Box<dyn Error>> {
    let mut sim_params: Vec<(usize, usize, u8, usize)> = Vec::new();

    for (n, a, d_max, max_seed) in iproduct!(
        &params.n_particles,
        &params.array_sizes,
        &params.d_max_vals,
        &params.seeds
    ) {
        sim_params.push((*n, *a, *d_max, *max_seed))
    }

    if params.parallel == 2 {
        sim_params.iter().for_each(|(n, a, d_max, max_seed)| {
            let now: Instant = Instant::now();
            let (cpu_time, r_avg) = match sim::run_parallel(*n, *a, *d_max, *max_seed, &params) {
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

            let results = (*n, cpu_time, r_avg, *max_seed, *d_max);

            if params.write_g == true {
                match write_g_radii(&results, params) {
                    Ok(_) => (),
                    Err(_) => panic!("Writing results to csv has failed!"),
                };
            }
        });
    } else {
        sim_params.iter().for_each(|(n, a, d_max, max_seed)| {
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

            let results = (*n, cpu_time, r_avg, *max_seed, *d_max);

            if params.write_g == true {
                match write_g_radii(&results, params) {
                    Ok(_) => (),
                    Err(_) => panic!("Writing results to csv has failed!"),
                };
            }

            
        });
    }

    Ok(())
}

fn setup_parallel(params: &InputParams) -> Result<(), Box<dyn Error>> {
    let mut sim_params: Vec<(usize, usize, u8, usize)> = Vec::new();

    for (n, a, d_max, max_seed) in iproduct!(
        &params.n_particles,
        &params.array_sizes,
        &params.d_max_vals,
        &params.seeds
    ) {
        sim_params.push((*n, *a, *d_max, *max_seed))
    }

    // Doesn't check for the option to run the seeds in parallel,
    //  since we're already running the sims in parallel
    sim_params.par_iter().for_each(|(n, a, d_max, max_seed)| {
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

        let results = (*n, cpu_time, r_avg, *max_seed, *d_max);

        if params.write_g == true {
            match write_g_radii(&results, params) {
                Ok(_) => (),
                Err(_) => panic!("Writing results to csv has failed!"),
            };
        }
    });

    Ok(())
}
