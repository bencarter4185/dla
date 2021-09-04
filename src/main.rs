/*
First attempt at writing some DLA code

Not going to be especially clever or optimised
*/

mod sim;

mod ini_parser;
use ini_parser::{InputParams, SimulationParams};

use itertools::iproduct;

use std::{error::Error, time::Instant};

fn main() -> Result<(), Box<dyn Error>> {
    // Load the params in the `.ini` file
    let params = InputParams::new(true); 

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
    
    Ok(())
}
