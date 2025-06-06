pub mod fst;
pub mod gst;
pub mod nei;
pub mod chord;
pub mod bray_curtis;
pub mod jost_d;
pub mod reynolds;
pub mod rogers;

pub fn generate_sample_pairs(samples: &[String]) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for i in 0..samples.len() {
        for j in (i + 1)..samples.len() {
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
