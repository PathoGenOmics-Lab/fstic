pub mod common;
pub mod fst;
#[cfg(test)]
mod tests;
pub mod gst;
pub mod nei;
pub mod chord;
pub mod bray_curtis;
pub mod jost_d;
pub mod reynolds;
pub mod rogers;

pub fn generate_sample_pairs(samples: &[String]) -> Vec<(usize, usize)> {
    let n = samples.len();
    if n < 2 {
        return Vec::new();
    }
    let mut pairs = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            pairs.push((i, j));
        }
    }
    pairs
}

pub fn create_distance_matrix(results: &[((usize, usize), f64)], num_samples: usize) -> Vec<Vec<f64>> {
    let mut matrix = vec![vec![0.0; num_samples]; num_samples];
    for &((i, j), dist) in results {
        matrix[i][j] = dist;
        matrix[j][i] = dist;
    }
    matrix
}
