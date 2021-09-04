/*
Library file used for writing csv files to disk.
*/

extern crate csv;
use std::error::Error;
use std::{env, fs};

use super::{SimulationParams, SimulationResults};

pub fn check_folder_exists(folder: &String) -> Result<bool, Box<dyn Error>> {
    // Generate path to folder
    let mut path = env::current_dir()?;
    path.push(folder);

    // Check if metadata for the folder exists. If it fails, return false
    let metadata = fs::metadata(path);
    let folder_exists = match metadata {
        Ok(_) => true,
        Err(_) => false,
    };

    Ok(folder_exists)
}

pub fn write_csv(params: &SimulationParams, results: &SimulationResults) -> Result<(), Box<dyn Error>> {
    // Unpack struct of params
    let r_max: f64 = params.max_radius as f64;
    let max_seed: i32 = params.max_seed;
    let init_seed: i32 = params.init_seed;
    let _add_plot: bool = params.add_plot;
    
    // Check if the `data` folder exists. If not, create it
    let folder: String = String::from("data");
    if !(check_folder_exists(&folder)?) {
        fs::create_dir(&folder)?
    }

    // Create a filename based on input parameters
    let filepath = format!("./{}/r{}_seeds{}_iseed{}.csv", folder, r_max, max_seed, init_seed);

    let mut wtr = csv::Writer::from_path(filepath)?;

    for i in 0..results.radii.len() {
        wtr.write_record(&[&results.radii[i].to_string(), &results.n_tree[i].to_string()])?;
    }

    Ok(())
}
