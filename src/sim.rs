/*
Library file for running DLA simulations
*/

/*
Imports
*/
mod random;
use random::Ran2Generator;

use inline_python::python;

use math::round::ceil;

use std::error::Error;

use super::ini_parser::SimulationParams;

mod csv_writer;
use csv_writer::write_csv;

/*
Structs
*/

pub struct SimulationResults {
    radii: Vec<i64>,
    n_tree: Vec<i64>,
}

impl SimulationResults {
    pub fn new(radii: Vec<i64>, n_tree: Vec<i64>) -> SimulationResults {
        let results: SimulationResults = SimulationResults {
            radii: radii,
            n_tree: n_tree,
        };
        results
    }
}

/*
Functions
*/

// Conversions between array indices and coordinates
fn to_array_index(val: i64, r_max: f64) -> usize {
    (val + r_max as i64) as usize
}

fn pbc(index: usize, s_len: usize) -> (usize, usize) {
    // Deal with periodic boundary conditions
    let index_bl: usize; // The index of the site below or left
    let index_ar: usize; // The index of the site above or right

    if index as i64 - 1 < 0 {
        index_bl = s_len - 1
    } else {
        index_bl = index - 1
    }

    if index + 1 > s_len - 1 {
        index_ar = 0
    } else {
        index_ar = index + 1
    }

    (index_bl, index_ar)
}

/*
Constants
*/

// Pi
const PI: f64 = std::f64::consts::PI;

pub fn run(params: SimulationParams) -> Result<(), Box<dyn Error>> {
    // Unpack input parameters
    let r_max: f64 = params.max_radius as f64;
    let max_seed: i32 = params.max_seed;
    let init_seed: i32 = params.init_seed;
    let add_plot: bool = params.add_plot;
    let write_data: bool = params.write_data;

    println!(
        r"Running simulation for
        Substrate Radius = {}, Max Seed = {},
        Initial Seed = {}, Add Plot = {}...",
        r_max, max_seed, init_seed, add_plot
    );

    // Create a new random number generator based on the ran2() function
    let mut rng: Ran2Generator = Ran2Generator::new(0); // Initial seed = 0

    // Warm up the random number generator (not necessary, I just want to)
    for _ in 0..100 {
        rng.next();
    }

    // Generate empty 2D array of side 2*r_max (for now)
    let s_len = 2 * r_max as usize;
    let mut s: Vec<Vec<bool>> = vec![vec![false; s_len]; s_len];

    // Instantiate our initial particle at (0,0)
    let xi: usize = 0 + r_max as usize;
    let yi: usize = 0 + r_max as usize;
    s[xi][yi] = true;

    // Generate some angles
    let mut theta: f64;
    let mut x: i64;
    let mut y: i64;
    let mut xi: usize;
    let mut yi: usize;

    // Start at a small radius r = 2 (ish)
    let mut r: f64 = 2.0;
    let r_kill_ratio: f64 = 1.5;
    let mut r_kill: f64 = r * r_kill_ratio;

    let mut _n_particles: usize = 0; // Might use it later

    loop {
        // Generate random angle
        theta = rng.next() * 2.0 * PI;

        // Calculate associated x and y indices
        // Generates points with a radius just slightly less than r_max
        x = (r * theta.cos()) as i64;
        y = (r * theta.sin()) as i64;

        // Get associated array indices
        xi = to_array_index(x, r_max);
        yi = to_array_index(y, r_max);

        // Saturation/exit condition: If point is already occupied: break
        if ((x * x + y * y) as f64).sqrt() > r_max {
            break;
        }

        // // Print value
        // println!("Angle generated: {}", theta);
        // println!("Associated (x,y) coordinates: ({},{})", x, y);
        // println!("Array (xi,yi) indices: ({},{})", xi, yi);
        // println!(
        //     "Radius: {}",
        //     ((x as f64).powi(2) + (y as f64).powi(2)).sqrt()
        // );
        // println!("");

        // Temp values for generated numbers
        let mut x_temp: f64;
        let mut y_temp: f64;

        // Do a random walk
        loop {
            // Perform the random walk
            x_temp = rng.next();
            y_temp = rng.next();

            // Equal chance of moving up/moving down/not moving
            if x_temp <= 0.5 {
                if xi > 0 {
                    xi -= 1;
                    x -= 1;
                } else {
                    xi = s_len - 1;
                    x = s_len as i64 / 2;
                }
            } else {
                if xi < s_len - 1 {
                    xi += 1;
                    x += 1;
                } else {
                    xi = 0;
                    x = -(s_len as i64) / 2;
                }
            }

            if y_temp <= 0.5 {
                if yi > 0 {
                    yi -= 1;
                    y -= 1;
                } else {
                    yi = s_len - 1;
                    y = s_len as i64 / 2;
                }
            } else {
                if yi < s_len - 1 {
                    yi += 1;
                    y += 1;
                } else {
                    yi = 0;
                    y = -(s_len as i64) / 2;
                }
            }

            // Generate indices for the sites, below, above, right, and left
            let (x_l, x_r) = pbc(xi, s_len);
            let (y_b, y_a) = pbc(yi, s_len);

            // Exit condition: is a nearest neighbour occupied?
            if [s[x_l][yi], s[x_r][yi], s[xi][y_b], s[xi][y_a]].contains(&true) {
                // Deposit at the current site
                s[xi][yi] = true;

                /*
                Need to generate new radius for next particle.
                Get radius of current particle: perform a ceiling function, add 2 (for safety) and store for next
                */

                // Need to generate new radius for next particle.
                // Get radius of deposited particle
                let r_current = ceil(((x * x + y * y) as f64).sqrt(), 0) + 2.0;
                if r_current > r {
                    r = r_current;
                    r_kill = r * r_kill_ratio;
                }

                //             println!(
                //                 r"It stuck!
                // (x,y) = ({},{})
                // (xi,yi) = ({},{})
                // i = {}",
                //                 x, y, xi, yi, i
                //             );
                break;
            }

            // Safety valve
            if ((x * x + y * y) as f64).sqrt() > r_kill {
                break;
            }
        }

        _n_particles += 1;
    }

    if add_plot == true {
        // Create a clone of s that we don't mind losing
        let s_py = s.clone();

        python! {
            import matplotlib.pyplot as plt
            from matplotlib import colors
            import numpy as np

            // Use seaborn for pretty graphs
            plt.style.use("default")

            theta = np.linspace(0, 2*np.pi, 100)

            r = 's_len / 2

            x1 = r*np.cos(theta) + r
            x2 = r*np.sin(theta) + r


            s = np.array('s_py)

            fig = plt.figure()
            fig.set_size_inches(12.5, 12.5)

            // Create discrete colormap
            cmap = colors.ListedColormap(["white", "black"])
            bounds = [0,0.5,1] # white if 0, black if 1
            norm = colors.BoundaryNorm(bounds, cmap.N)

            // Create a subplot
            ax = plt.subplot(111)

            // Plot the empty cells and the blocks
            ax.imshow(s, cmap=cmap, norm=norm)

            ax.plot(x1, x2)

            // Set the tick parameters so we don't see any ticks or labels
            plt.tick_params(
                axis="both",          # changes apply to the x-axis
                which="both",      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                left=False,
                right=False,
                labelleft=False,
                labelbottom=False) # labels along the bottom edge are off

            // // Create figure
            // fig = plt.figure()
            // fig.set_size_inches(12.5, 7.5)
            // ax = plt.subplot(111)

            // // plt.scatter('s_py)
        }
    }

    // Vectors of results
    let radii: Vec<i64> = (1..r_max as i64).collect(); // Start from a radius of 1

    let mut n_tree_vec: Vec<i64> = vec![0; radii.len()];
    let mut i: usize = 0;

    // We've now added blocks 'to saturation'
    for r in radii.iter() {
        let mut x: i64 = 0 - r;
        let mut xi: usize = to_array_index(x, r_max);
        let mut y: i64;
        let mut yi: usize;

        let mut n_tree: i64 = 0; // Number of points in the circle part of the tree

        while x < *r {
            // Reset y
            y = 0 - r;
            yi = to_array_index(y, r_max);

            while y < *r {
                if x * x + y * y < r * r {
                    // Now we're in the circle, check if there is a tree or empty space
                    if s[xi][yi] == true {
                        n_tree += 1;
                    }
                }

                y += 1;
                yi += 1;
            }
            x += 1;
            xi += 1;
        }

        //         println!(r"Checked radius {}
        // Number of points in circle: {}
        // Number of points in tree: {}", r, n_circle, n_tree);

        n_tree_vec[i] = n_tree;

        i += 1;
    }

    if add_plot == true {
        // Create clones of the data that we don't mind losing
        let radii_py: Vec<i64> = radii.clone();
        let n_tree_vec_py: Vec<i64> = n_tree_vec.clone();

        python! {
            import matplotlib.pyplot as plt
            import numpy as np

            // Use seaborn for pretty graphs
            plt.style.use("seaborn")

            // Create figure
            fig = plt.figure()
            fig.set_size_inches(12.5, 7.5)
            ax = plt.subplot(111)

            plt.title("Number of Points in Brownian Tree as function of Radius", fontsize = 20)
            plt.xlabel(r"$ln(r)$")
            plt.ylabel(r"$ln(N)$")

            plt.plot(np.log('radii_py), np.log('n_tree_vec_py))

            plt.show()
        }
    }

    if write_data == true {
        // We want to write the results of this simulation to disk
        let results: SimulationResults = SimulationResults::new(radii, n_tree_vec);

        write_csv(&params, &results)?;
    }

    Ok(())
}
