use std::vec;

use random::Ran2Generator;

use inline_python::python;

/*
First attempt at writing some DLA code

Not going to be especially clever or optimised
*/

/*
Imports
*/
mod random;

/*
Functions
*/

// Conversions between array indices and coordinates
fn to_array_index(val: i64) -> usize {
    (val + R_GEN as i64) as usize
}
fn _to_coordinate(val: usize) -> i64 {
    val as i64 - R_GEN as i64
}

/*
Constants
*/

// Radius of 'generation circle'
const R_GEN: f64 = 128.0;
const PI: f64 = std::f64::consts::PI;

// x = x_ind - R_GEN
// x_ind = x + R_GEN

fn main() {
    // Create a new random number generator based on the ran2() function
    let mut rng: Ran2Generator = Ran2Generator::new(0); // Initial seed = 0

    // Warm up the random number generator (not necessary, I just want to)
    for _ in 0..100 {
        rng.next();
    }

    // Generate empty 2D array of side 2*R_GEN (for now)
    let s_len = 2 * R_GEN as usize;
    let mut s: Vec<Vec<bool>> = vec![vec![false; s_len]; s_len];
    // let mut s: ArrayBase<OwnedRepr<bool>, Dim<[usize; 2]>> =
    //     Array::from_shape_fn((s_len, s_len), |(_, _)| false);

    // Instantiate our initial particle at (0,0)
    let xi: usize = 0 + R_GEN as usize;
    let yi: usize = 0 + R_GEN as usize;
    s[xi][yi] = true;  //[[xi, yi]] = true;

    // Generate some angles
    let mut theta: f64;
    let mut x: i64;
    let mut y: i64;
    let mut xi: usize;
    let mut yi: usize;

    loop {
        // Generate random angle
        theta = rng.next() * 2.0 * PI;

        // Calculate associated x and y indices
        // Generates points with a radius just slightly less than R_GEN
        x = (R_GEN * theta.cos()) as i64;
        y = (R_GEN * theta.sin()) as i64;

        // Get associated array indices
        xi = to_array_index(x);
        yi = to_array_index(y);

        // If point is already occupied: break
        if s[xi][yi] == true { break }

        // // Print value
        // println!("Angle generated: {}", theta);
        // println!("Associated (x,y) coordinates: ({},{})", x, y);
        // println!("Array (xi,yi) indices: ({},{})", xi, yi);
        // println!(
        //     "Radius: {}",
        //     ((x as f64).powi(2) + (y as f64).powi(2)).sqrt()
        // );
        // println!("");

        // Iteration counter
        let mut i: usize = 0;

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
            if i >= 10000000 {
                // println!("Didn't stick");
                break;
            }

            // Increment iterations
            i += 1;
        }
    }

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

    // Vectors of results
    let radii: Vec<i64> = (1..R_GEN as i64).collect(); // Start from a radius of 1

    let mut results: Vec<i64> = vec![0; radii.len()];
    let mut i: usize = 0;

    // We've now added blocks 'to saturation'
    for r in radii.iter() {



        let mut x: i64 = 0 - r;
        let mut xi: usize = to_array_index(x);
        let mut y: i64;
        let mut yi: usize;

        let mut n_circle: i64 = 0; // Number of points in the circle
        let mut n_tree: i64 = 0; // Number of points in the circle part of the tree

        while x < *r {
            // Reset y
            y = 0 - r;
            yi = to_array_index(y);

            while y < *r {
                if x * x + y * y < r * r {
                    // We're in the circle we want to check.
                    if s[xi][yi] == true {
                        n_tree += 1;
                    }

                    // println!("In circle: ({},{})", x, y);
                    n_circle += 1;
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

    results[i] = n_tree;

    i+= 1;
    }

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

        plt.plot(np.log('radii), np.log('results))

        plt.show()
    }


}

fn pbc(index: usize, s_len: usize) -> (usize, usize) {
    // Deal with periodic boundary conditions in x
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
