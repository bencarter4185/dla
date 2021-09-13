/*
Library file used to run DLA simulations.
*/

/*
Imports
*/

use itertools::iproduct;

use rayon::prelude::*;

use std::error::Error;

pub(crate) mod writer;

use std::time::Instant;

use itertools::Itertools;
use ndarray::{Array, ArrayBase, Dim, OwnedRepr};

use math::round::floor;

use kurbo::common::solve_quadratic as quad;

use crate::{
    sim::{
        random::Ran2Generator,
        writer::{write_data, write_tree},
    },
    InputParams,
};

mod random;

use hashbrown::HashMap;

/*
Constants
*/

const PI_2: f32 = std::f32::consts::PI * 2.0;

/*
Structs
*/

pub struct Data {
    // Input Parameters
    n: usize, // number of particles
    a: usize, // array size
    d_max: u8,
    // Simulation data points
    rp: f32,      // particle radii
    l_min: f32,   // minimum step size before collision
    r_max: f32,   // Current maximum radius
    r_g: f32,     // Generation radius
    r_kill: f32,  // Kill radius
    x0: f32,      // Coordinates of origin
    y0: f32,      // Coordinates of origin
    ix0_a: usize, // Index of origin for square of side a
    iy0_a: usize, // Index of origin for square of side a
    ix0_s: usize, // Index of origin for square of side s
    iy0_s: usize, // Index of origin for square of side s
    c: u8,        // Collision condition
    seed: usize,
    // Arrays
    /* Omega: Off-lattice array of size n (number of particles) x 2 */
    omega: HashMap<usize, (f32, f32)>,
    /* Theta: Vicinity grid */
    theta: ArrayBase<OwnedRepr<u8>, Dim<[usize; 2]>>,
    /* Psi: On-lattice distance grid */
    psi: ArrayBase<OwnedRepr<u8>, Dim<[usize; 2]>>,
    // Random number generator
    rng: Ran2Generator,
}

impl Data {
    fn new(n: usize, a: usize, d_max: u8, seed: usize, params: &InputParams) -> Data {
        let idum: i32 = -1 * (params.init_seed + seed) as i32;

        let mut rng: Ran2Generator = Ran2Generator::new(idum);

        // Warm up random number generator
        for _ in 0..100 {
            rng.next();
        }

        // Define (x,y) coordinates of origin
        let x0: f32 = 0.0;
        let y0: f32 = 0.0;

        // Array indices of origin
        let (ix0_a, iy0_a) = get_array_index(x0, y0, a);

        let s: usize = 2 * d_max as usize + 1;

        let (ix0_s, iy0_s) = get_array_index(x0, y0, s);

        // Hard-code some variables
        let rp: f32 = 1.0;
        let l_min: f32 = 1.0;

        // Generate collision condition
        let c: u8 = (2.0 * rp + l_min + 1.0) as u8;

        // Kill radius
        let r_max: f32 = 1.0; // Initial value
        let r_kill: f32 = 5.0 * r_max;

        let data: Data = Data {
            n: n,
            a: a,
            d_max: d_max,
            rp: 1.0,    // hard-coded for now
            l_min: 1.0, // hard-coded for now
            r_max: r_max,
            r_g: 3.0, // Initial value
            r_kill: r_kill,
            x0: x0,
            y0: y0,
            ix0_a: ix0_a,
            iy0_a: iy0_a,
            ix0_s: ix0_s,
            iy0_s: iy0_s,
            c: c,
            seed: seed,
            omega: HashMap::with_capacity(n),
            theta: generate_theta(d_max),
            psi: Array::from_elem((a, a), d_max),
            rng: rng,
        };

        data
    }
}

pub fn run_fractal(
    n: usize,
    a: usize,
    d_max: u8,
    max_seed: usize,
    params: &InputParams,
) -> Result<f64, Box<dyn Error>> {
    let seeds: Vec<usize> = (0..max_seed).collect_vec();
    let mut g_radius: f64 = 0.0;

    if params.run_parallel == true {
        let results: Vec<f64> = seeds
            .par_iter()
            .map(|seed| {
                let g_radius = match do_sim(n, a, d_max, *seed, params) {
                    Ok((_, _, _, g_radius)) => g_radius,
                    Err(_) => panic!("Performing the simulation has failed!"),
                };
                g_radius
            })
            .collect();

        // Average over the ensemble
        for i in 0..max_seed {
            g_radius += results[i] / max_seed as f64;
        }
    } else {
        for seed in seeds {
            let g_r = match do_sim(n, a, d_max, seed, params) {
                Ok((_, _, _, g_radius)) => g_radius,
                Err(_) => panic!("Performing the simulation has failed!"),
            };

            g_radius += g_r / max_seed as f64;
        }
    }

    Ok(g_radius)
}

pub fn run(
    n: usize,
    a: usize,
    d_max: u8,
    max_seed: usize,
    params: &InputParams,
) -> Result<(), Box<dyn Error>> {
    let seeds = (0..max_seed).collect_vec();

    let results_len = generate_timestep(a);

    let results: Vec<(Vec<f64>, Vec<f64>)> = seeds
        .par_iter()
        .map(|seed| {
            let (radii, n_tree) = match do_sim(n, a, d_max, *seed, params) {
                Ok((radii, n_tree, _, _)) => (radii, n_tree),
                Err(_) => panic!("Performing the simulation has failed!"),
            };
            (radii, n_tree)
        })
        .collect();

    // Need to work out ensemble averages
    // Generate vectors to store ensemble averages
    let mut radii: Vec<f64> = vec![0.0; results_len];
    let mut n_tree: Vec<f64> = vec![0.0; results_len];

    // Don't need to average radius, so can do it in its own loop
    for j in 0..results_len {
        radii[j] += results[0].0[j];
    }

    // Average across the ensemble
    for (i, j) in iproduct!(0..max_seed as usize, 0..results_len) {
        n_tree[j] += results[i].1[j] / max_seed as f64;
    }

    if params.write_data == true {
        write_data(n, params, radii, n_tree, max_seed, d_max)?;
    }

    Ok(())
}

fn do_sim(
    n: usize,
    a: usize,
    d_max: u8,
    seed: usize,
    params: &InputParams,
) -> Result<(Vec<f64>, Vec<f64>, u128, f64), Box<dyn Error>> {
    // Create new struct of the data in this simulation
    let mut data: Data = Data::new(n, a, d_max, seed, params);

    // Temporarily store the origin as local variables
    let x0: f32 = data.x0;
    let y0: f32 = data.y0;

    // Add seed particle at origin
    data.omega.insert(cantor(data.ix0_a, data.iy0_a), (x0, y0));
    overlap_psi_theta(&mut data, x0, y0);

    // Calculate `cpu time` (kind of)
    let now = Instant::now();

    /* We have generated the arrays, so can now perform a simulation */
    launch_particles(&mut data);

    let cpu_time = Instant::now().duration_since(now).as_millis();

    // If asked, write this tree to disk
    if params.write_tree == true {
        write_tree(&data, params, seed)?;
    }

    let (radii, n_tree, r_avg) = count_tree(&data);

    Ok((radii, n_tree, cpu_time, r_avg))
}

fn generate_timestep(a: usize) -> usize {
    // Generate empty vectors to store results
    let max_radius: f64 = (a as f64 - 1.0) / 2.0;

    let mut r: f64 = 8.0;
    let mut r_points: usize = 0;
    while r < max_radius {
        // Perform a logarithmic timestep
        let n: f64 = r / 100.0 + 0.1;
        r += n as f64;
        r_points += 1;
    }

    r_points
}

fn count_tree(data: &Data) -> (Vec<f64>, Vec<f64>, f64) {
    // Generate empty vectors to store results
    let max_radius: f64 = (data.a as f64 - 1.0) / 2.0;

    let mut radii: Vec<f64> = Vec::new();

    let mut r: f64 = 8.0;
    let mut r_points: usize = 0;
    while r < max_radius {
        // Perform a logarithmic timestep
        let n: f64 = r / 100.0 + 0.1;
        r += n as f64;
        r_points += 1;
        radii.push(r)
    }

    let mut n_tree: Vec<f64> = vec![0.0; r_points];

    let mut x: f64;
    let mut y: f64;
    let mut r: f64;

    let mut r_tot: f64 = 0.0;

    for particle in &data.omega {
        x = particle.1 .0 as f64;
        y = particle.1 .1 as f64;
        r = (x.powi(2) + y.powi(2)).sqrt();

        r_tot += r;

        let mut i: usize = r_points - 1;

        loop {
            if r < radii[i] {
                n_tree[i] += 1.0;
            } else {
                break;
            }

            if i == 0 {
                break;
            }

            i -= 1;
        }
    }

    let r_avg: f64 = r_tot / data.n as f64;

    (radii, n_tree, r_avg)
}

fn calc_rg(data: &mut Data) -> f32 {
    // Calculate the new generation radius
    data.r_max * 2.0 as f32 + data.rp * 2.0
}

fn launch_particles(data: &mut Data) {
    /*
    Seed particle is already placed at (x0, y0). Psi has been recalculated for the seed particle.
    Launch a walker from from the entrance radius with a small step, L_min
    */

    let mut i: usize = 2; // First particle to stick will be i = 2
    while i <= data.n {
        // Generate new particle
        let mut alpha: f32 = next_angle(&mut data.rng);
        let mut x: f32 = data.r_g * alpha.cos(); // x coordinate of particle
        let mut y: f32 = data.r_g * alpha.sin(); // y coordinate of particle

        let mut d_wc: u8;

        // Size of next (non-colliding) step
        let mut step: f32;

        // Radius of particle
        let mut r: f32;

        // Perform a random walk
        loop {
            // Get current particle radius
            r = (x.powi(2) + y.powi(2)).sqrt();

            // Break if particle is larger than r_kill
            if r > data.r_kill {
                break;
            }

            // Generate a random angle
            alpha = next_angle(&mut data.rng);

            // Get the current distance d_wc
            let (xi, yi) = get_array_index(x, y, data.a);
            // Catch index out of range exceptions by defaulting to d_max
            d_wc = match data.psi.get([xi, yi]) {
                Some(val) => *val,
                None => data.d_max,
            };

            // Decide which situation the particle is currently in
            match d_wc {
                d_wc if d_wc <= data.c => {
                    // Possibility of collision on next step
                    // Get all the particles in upsilon in a square of size (2*c + 1) around the current x and y
                    let particles: Vec<(f32, f32)> = find_particles(data, x, y);
                    let lh = check_collisions(data, particles, &mut x, &mut y, alpha, i);

                    if lh != f32::MAX {
                        // Do the collision
                        do_step(&mut x, &mut y, lh, alpha);

                        // If the radius of this new particle is greater than the current highest radius:
                        //  Generate new r_gen
                        r = (x.powi(2) + y.powi(2)).sqrt();
                        if r > data.r_max {
                            data.r_max = r;
                            data.r_kill = 5.0 * r; // temp, set back to 5.0 later
                            data.r_g = calc_rg(data);
                        }

                        // Update arrays
                        let (xi, yi) = get_array_index(x, y, data.a);
                        data.omega.insert(cantor(xi, yi), (x, y));
                        overlap_psi_theta(data, x, y);

                        // Increment i
                        i += 1;

                        // Particle deposited. Break from the infinite loop
                        break;
                    }

                    // If no collision, do a step of size l_min
                    do_step(&mut x, &mut y, data.l_min, alpha);
                }

                d_wc if d_wc == data.d_max => {
                    // Calculate distance of particle from origin
                    let d_w0: f32 = (x.powi(2) + y.powi(2)).sqrt();

                    // Calculate new step size
                    step = max_f32(
                        d_w0 - data.r_max - data.c as f32,
                        data.d_max as f32 - data.c as f32,
                    );

                    // Do a big step without the need to check for collision
                    do_step(&mut x, &mut y, step, alpha);
                }

                _ => {
                    // Can do a step of size d_wc - c without needing to check for collision
                    do_step(&mut x, &mut y, d_wc as f32 - data.c as f32, alpha);
                }
            }
        }
    }
}

fn find_particles(data: &mut Data, x: f32, y: f32) -> Vec<(f32, f32)> {
    // Assume size means radius, i.e. the square sizer is actually of side (2*data.c + 1)
    let side: usize = (2 * data.c + 1) as usize;

    // Get the indices in psi in which the particle currently resides
    let (ixp, iyp) = get_array_index(x, y, data.a);

    // Generate an empty square of side c
    let square: ArrayBase<OwnedRepr<i32>, Dim<[usize; 2]>> = Array::from_elem((side, side), 0);

    // Get the indices of the centre of square of size c
    let (ixt0, iyt0) = get_array_index(0.0, 0.0, side);

    // Calculate the offset in indices between psi and square
    let x_offset = ixp - ixt0;
    let y_offset = iyp - iyt0;

    let mut ixp: usize; // x index of cell in Psi
    let mut iyp: usize; // y index of cell in Psi
    let mut val: u8;

    let mut particles: Vec<(f32, f32)> = Vec::new();

    let mut x: f32;
    let mut y: f32;

    for ((ixt, iyt), _) in square.indexed_iter() {
        ixp = ixt + x_offset;
        iyp = iyt + y_offset;

        // Get the value in the cell
        val = match data.psi.get([ixp, iyp]) {
            Some(val) => *val,
            // Handle out of range exceptions
            None => panic!("Trying to add a particle out of bounds! Set the array size to be larger and try again."),
        };

        // If there is a particle in the cell, add it to particles
        if val == 0 {
            // Get the `label` of this grid position
            let can: usize = cantor(ixp, iyp);
            x = data.omega[&can].0;
            y = data.omega[&can].1;

            particles.push((x, y))
        }
    }

    particles
}

fn overlap_psi_theta(data: &mut Data, x: f32, y: f32) {
    /*
    Overlap Psi and Theta. Replace Psi with the minimum of the two matrices.
    */

    // Get the indices of the pivot point of psi -- the point at which we've added a particle
    let (ixp, iyp) = get_array_index(x, y, data.a);

    // Calculate the offset in indices between psi and theta
    let ixo = ixp - data.ix0_s;
    let iyo = iyp - data.iy0_s;

    let mut ixp: usize; // x index of cell in Psi
    let mut iyp: usize; // y index of cell in Psi

    // Overlap cells in theta and psi
    for ((ixt, iyt), val_t) in data.theta.indexed_iter() {
        ixp = ixt + ixo;
        iyp = iyt + iyo;

        let val_psi = match data.psi.get([ixp, iyp]) {
            Some(val) => *val,
            None => panic!("Calculating new distances has caused an out of bounds error. Try again with a larger array or less particles."),
        };

        // Replace the cell with the smallest of the two values
        data.psi[[ixp, iyp]] = min(val_psi, *val_t);
    }
}

fn next_angle(rng: &mut Ran2Generator) -> f32 {
    rng.next() as f32 * PI_2
}

fn check_collisions(
    data: &mut Data,
    particles: Vec<(f32, f32)>,
    x: &mut f32,
    y: &mut f32,
    alpha: f32,
    i: usize,
) -> f32 {
    // Type convert all our variables to f64 to increase precision
    let x: f64 = *x as f64;
    let y: f64 = *y as f64;
    let alpha: f64 = alpha as f64;

    // Define variables for the quadratic solver
    let mut a: f64;
    let mut b: f64;
    let mut c: f64;

    // Calculate the diameter of the particle
    let dp: f64 = 2.0 * data.rp as f64;

    // Keep a note of the smallest jump distance
    let mut lh: f32 = f32::MAX;

    for particle in particles {
        let xp = particle.0 as f64;
        let yp = particle.1 as f64;

        // Define our variables for the quadratic solver
        a = 1.0;
        b = 2.0 * (alpha.cos() * (x - xp) + alpha.sin() * (y - yp));
        c = (xp - x).powi(2) + (yp - y).powi(2) - dp.powi(2);

        // Solve the quadratic equation and get the roots
        let roots = quad(c, b, a);

        /*
        5 situations:
            1 - No roots. Walker does not collide with particle and step of size L can be taken.
            2 - Lh = 0 (ans contains 0). Particle has already stuck and we error
            3 - Lh < 0. Walker does not collide with the particle because it would need to travel in the
                opposite direction. Take a step of size L.
            4 - Lh > L. Walker does not collide with particle because the normal step size is too short.
                Take a step of size L.
            5 - 0 < Lh < L. Walker collides with particle. Take a step of size Lh.
        */

        // Extract the smaller root
        // Can now safely convert back to f32
        let root: f32;
        match roots.get(0) {
            Some(val) => {
                root = *val as f32;
            }
            // Situation 1: Particle doesn't touch. Take a step of size L.
            None => {
                continue;
            }
        };

        // Situation 2: If one of the roots is 0.0, panic
        if root == 0.0 {
            panic!(
                r"Error! This particle has already stuck!
Error occured during simulation for n = {}, seed = {}.
Error occured when adding particle i = {}",
                data.n, data.seed, i
            )
        }
        // Situation 3 and 4, wrong direction or too far:
        else if root < 0.0 || root > data.l_min {
            continue;
        }
        // Situation 5, walker collides with particle
        else {
            if root < lh {
                lh = root
            }
        }
    }

    // All potential collisions checked. If any, choose the collision of lowest lh
    if lh != f32::MAX {}

    lh
}

fn do_step(x: &mut f32, y: &mut f32, step: f32, alpha: f32) {
    // Get new values of x and y
    *x += step * alpha.cos();
    *y += step * alpha.sin();
}

fn get_array_index(x: f32, y: f32, side_length: usize) -> (usize, usize) {
    /*
    Get array index of cell given (x, y) coordinates, and side length of array
    */
    let xi = ((side_length as f32 + 1.0) / 2.0 - x.round()) as usize;
    let yi = ((side_length as f32 + 1.0) / 2.0 - y.round()) as usize;
    (xi, yi)
}

fn generate_theta(d_max: u8) -> ArrayBase<OwnedRepr<u8>, Dim<[usize; 2]>> {
    // Define side length and of 2d array and generate
    let side: usize = 2 * d_max as usize + 1;
    let mut theta: ArrayBase<OwnedRepr<u8>, Dim<[usize; 2]>> = Array::zeros((side, side));

    // Define array indices of the origin (x0, y0)
    let (ix0, iy0) = get_array_index(0.0, 0.0, side);

    // Convert to i64 so we don't have to worry about the order of subtraction
    let ix0: i64 = ix0 as i64;
    let iy0: i64 = iy0 as i64;

    for ((ix, iy), value) in theta.indexed_iter_mut() {
        // Type convert out indices
        let ix: i64 = ix as i64;
        let iy: i64 = iy as i64;

        // Calculate the distance between the indices as f32, then round down
        let d: usize = floor(
            (((ix - ix0).abs().pow(2) + (iy - iy0).abs().pow(2)) as f64).sqrt(),
            0,
        ) as usize;

        // Round down to the nearest integer and store value in the array if less than d_max
        *value = if d <= d_max as usize {
            d as u8
        } else {
            d_max as u8
        };
    }

    theta
}

// Return the smallest of two values
fn min<T: Ord>(a: T, b: T) -> T {
    if a <= b {
        a
    } else {
        b
    }
}

// Return the largest of two f32 values
fn max_f32(a: f32, b: f32) -> f32 {
    if a >= b {
        a
    } else {
        b
    }
}

fn cantor(k1: usize, k2: usize) -> usize {
    // Return a unique integer based on the two numbers
    ((k1 + k2) * (k1 + k2 + 1)) / 2 + k2
}
