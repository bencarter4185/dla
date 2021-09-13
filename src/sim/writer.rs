extern crate csv;
use std::error::Error;
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

pub fn write_g_radii(
    params: &InputParams,
    n_particles: &Vec<usize>,
    g_radii: &Vec<f64>,
    max_seed: usize,
    d_max: u8
) -> Result<(), Box<dyn Error>> {
    // Check if the `data` folder exists. If not, create it
    let folder: String = String::from("data");
    if !(check_folder_exists(&folder)?) {
        fs::create_dir(&folder)?
    }

    let n_min: usize = *n_particles.iter().min().unwrap();
    let n_max: usize = *n_particles.iter().max().unwrap();

    // Create filename
    let filepath = format!(
        "./{}/g_radii_nmin{}_nmax{}_dmax{}_max_seed{}_iseed{}.csv", folder, n_min, n_max, d_max, max_seed, params.init_seed,
    );

    let mut wtr = csv::Writer::from_path(filepath)?;

    // Iterate through results and add to file
    for i in 0..n_particles.len() {
        let n: String = n_particles[i].to_string();
        let g_radius: String = g_radii[i].to_string();

        wtr.write_record([n, g_radius])?;
    }

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
        let x: String = item.1.0.to_string();
        let y: String = item.1.1.to_string();
        wtr.write_record([x, y])?;
    }

    Ok(())
}

// pub fn write_csv(
//     params: &SimulationParams,
//     results: &SimulationResults,
// ) -> Result<(), Box<dyn Error>> {
//     // Unpack struct of params
//     let r_max: f64 = params.max_radius as f64;
//     let max_seed: i32 = params.max_seed;
//     let init_seed: i32 = params.init_seed;

//     // Check if the `data` folder exists. If not, create it
//     let folder: String = String::from("data");
//     if !(check_folder_exists(&folder)?) {
//         fs::create_dir(&folder)?
//     }

//     // Create a filename based on input parameters
//     let filepath = format!(
//         "./{}/r{}_seeds{}_iseed{}.csv",
//         folder, r_max, max_seed, init_seed
//     );

//     let mut wtr = csv::Writer::from_path(filepath)?;

//     for i in 0..results.radii.len() {
//         wtr.write_record(&[
//             &results.radii[i].to_string(),
//             &results.n_tree[i].to_string(),
//         ])?;
//     }

//     Ok(())
// }

// pub fn write_bench(max_radii: &Vec<i32>, time: &Vec<u128>, particles: &Vec<f64>, r_max: i32, max_seed: i32) -> Result<(), Box<dyn Error>> {
//     // Check if the `data` folder exists. If not, create it
//     let folder: String = String::from("data");
//     if !(check_folder_exists(&folder)?) {
//         fs::create_dir(&folder)?
//     }

//     let filepath = format!("./{}/bench_r{}_seeds{}.csv", folder, r_max, max_seed);

//     let mut wtr = csv::Writer::from_path(filepath)?;

//     for i in 0..max_radii.len() {
//         wtr.write_record(&[
//             &max_radii[i].to_string(),
//             &time[i].to_string(),
//             &particles[i].to_string(),
//         ])?;
//     }

//     Ok(())
// }
