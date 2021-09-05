/*
First attempt at writing some DLA code

Not going to be especially clever or optimised
*/

mod sim;
use sim::csv_writer;

mod ini_parser;
use ini_parser::{InputParams, SimulationParams};

use itertools::{Itertools, iproduct};

use std::{error::Error, time::Instant, vec};

use inline_python::python;

fn main() -> Result<(), Box<dyn Error>> {
    // Load the params in the `.ini` file
    let params = InputParams::new(true); 

    if params.benchmark == true { // Ignore input params and run a benchmark 
        let max_seed = 100; // Make use of multithreading but nothing crazy
        let max_radii = (16..=params.benchmark_max).collect_vec(); // Minimum radius = 16, get strange behaviour for small r

        let mut time: Vec<u128> = vec![0; max_radii.len()];
        let mut particles: Vec<f64> = vec![0.0; max_radii.len()];
        let mut i: usize = 0;

        for max_radius in max_radii.clone() {

            println!("Running for max radius: {}", max_radius);

            let params: InputParams = InputParams::new(false);
            
            let current_params: SimulationParams = SimulationParams::new(params, max_radius, max_seed);

            let now = Instant::now();
            particles[i] = sim::run(current_params)?;
            let new_now = Instant::now();

            time[i] = new_now.duration_since(now).as_millis();

            i += 1;
        }
        csv_writer::write_bench(&max_radii, &time, &particles, params.benchmark_max, max_seed)?;

        plot_bench(&max_radii, &time, &particles);
    } else {
        for (max_radius, max_seed) in iproduct!(params.substrate_radii, params.seeds) {
            // Re-parse the params so we're not trying to use a moved value every loop
            //  Wasteful? Probably. But insignificant compared to the runtime of the simulations
            let params: InputParams = InputParams::new(false);
    
            let current_params = SimulationParams::new(params, max_radius, max_seed);
            
            let now = Instant::now();
            sim::run(current_params)?;
            let new_now = Instant::now();
            println!(
                r"Done! Took {:?}
                    ",
                new_now.duration_since(now)
            );
        }
    }
    
    Ok(())
}

fn plot_bench(max_radii: &Vec<i32>, time: &Vec<u128>, particles: &Vec<f64>) {
    python! {
        import matplotlib.pyplot as plt
        import numpy as np

        // Use seaborn for pretty graphs
        plt.style.use("seaborn")

        // Create figure
        fig1 = plt.figure(1)
        fig1.set_size_inches(12.5, 7.5)
        ax1 = plt.subplot(111)

        fig2 = plt.figure(2)
        fig2.set_size_inches(12.5, 7.5)
        ax2 = plt.subplot(111)

        ax1.set_title("Benchmark of total runtime for Max Seed = 100, as function of Max Radius")
        ax1.set_xlabel("".format(["ln(", r"$r_{max}$", ")"]))
        ax1.set_xlabel("".format(["ln(", r"$t$", ")"]))

        ax2.set_title("Benchmark of total runtime for Max Seed = 100, as function of Particles Added")
        ax2.set_xlabel("".format(["ln(", r"$N_p$", ")"]))
        ax2.set_xlabel("".format(["ln(", r"$t$", ")"]))

        ax1.plot(np.log('max_radii), np.log('time))
        ax2.scatter(np.log('particles), np.log('time))

        plt.show()
    }
}