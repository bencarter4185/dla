extern crate csv;
use std::error::Error;
use std::fs::OpenOptions;
use std::{env, fs};

use crate::InputParams;

use super::Data;

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

pub fn erase_g_radii(params: &InputParams) -> Result<(), Box<dyn Error>> {
    // Check if the `data` folder exists. If not, create it
    let folder: String = String::from("data");
    if !(check_folder_exists(&folder)?) {
        fs::create_dir(&folder)?
    }

    // Get max and min n_particles
    let mut n_min: usize = usize::MAX;
    let mut n_max: usize = 0;

    for radius in &params.n_particles {
        n_min = if *radius < n_min { *radius } else { n_min };
        n_max = if *radius > n_max { *radius } else { n_max };
    }

    // Create filename
    let filepath = format!("./{}/g_radii_nmin{}_nmax{}.csv", folder, n_min, n_max);

    match fs::remove_file(&filepath) {
        Ok(_) => (),
        Err(_) => (),
    };

    Ok(())
}

pub fn write_g_radii(
    results: &(usize, f64, f64, usize, u8),
    params: &InputParams,
) -> Result<(), Box<dyn Error>> {
    // Check if the `data` folder exists. If not, create it
    let folder: String = String::from("data");
    if !(check_folder_exists(&folder)?) {
        fs::create_dir(&folder)?
    }

    // Get max and min n_particles
    let mut n_min: usize = usize::MAX;
    let mut n_max: usize = 0;

    for radius in &params.n_particles {
        n_min = if *radius < n_min { *radius } else { n_min };
        n_max = if *radius > n_max { *radius } else { n_max };
    }

    // Create filename
    let filepath = format!("./{}/g_radii_nmin{}_nmax{}.csv", folder, n_min, n_max);

    let mut wtr = match OpenOptions::new().write(true).append(true).open(&filepath) {
        // If ok, file exists, and we can append to file
        //  else, create a new file with contents
        Ok(file) => csv::Writer::from_writer(file),
        Err(_) => csv::Writer::from_path(&filepath)?,
    };

    // Iterate through results and add to file
    let n: String = results.0.to_string();
    let cpu_time = results.1.to_string();
    let g_radius: String = results.2.to_string();
    let max_seed: String = results.3.to_string();
    let d_max: String = results.4.to_string();

    wtr.write_record([n, cpu_time, g_radius, max_seed, d_max])?;

    Ok(())
}

pub fn write_data(
    n: usize,
    params: &InputParams,
    radii: Vec<f64>,
    n_tree: Vec<f64>,
    max_seed: usize,
    d_max: u8,
) -> Result<(), Box<dyn Error>> {
    // Check if the `data` folder exists. If not, create it
    let folder: String = String::from("data");
    if !(check_folder_exists(&folder)?) {
        fs::create_dir(&folder)?
    }

    // Create a filename based on input parameters
    let filepath = format!(
        "./{}/data_n{}_dmax{}_seeds{}_iseed{}.csv",
        folder, n, d_max, max_seed, params.init_seed
    );

    let mut wtr = csv::Writer::from_path(filepath)?;

    // Iterate through data and write to file
    for i in 0..radii.len() {
        let radius: String = radii[i].to_string();
        let n: String = n_tree[i].to_string();

        wtr.write_record([radius, n])?;
    }

    Ok(())
}

pub fn write_tree(data: &Data, params: &InputParams, seed: usize) -> Result<(), Box<dyn Error>> {
    // Check if the `data` folder exists. If not, create it
    let folder: String = String::from("data");
    if !(check_folder_exists(&folder)?) {
        fs::create_dir(&folder)?
    }

    // Create a filename based on input parameters
    let filepath = format!(
        "./{}/tree_n{}_dmax{}_seed{}_iseed{}.csv",
        folder, data.n, data.d_max, seed, params.init_seed
    );

    let mut wtr = csv::Writer::from_path(filepath)?;

    // Iterate through Omega
    for item in &data.omega {
        let x: String = item.1 .0.to_string();
        let y: String = item.1 .1.to_string();
        wtr.write_record([x, y])?;
    }

    Ok(())
}
