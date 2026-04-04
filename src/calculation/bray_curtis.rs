use super::common::get_all_freqs_at_pos;
use crate::types::PositionalData;
use rayon::prelude::*;
use std::collections::HashSet;

/// Calculates Bray-Curtis dissimilarity: BC = 0.5 Σ |p_i − q_i|.
///
/// Equivalent to the Manhattan distance on allele frequency profiles.
/// With `--normalize` the sum is divided by the number of loci.
pub fn calculate_bray_curtis_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
    normalize: bool,
    num_loci: usize,
) -> f64 {
    let sum_dist: f64 = all_positions
        .par_iter()
        .map(|&pos| {
            let psf = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));
            let abs_diff: f64 = psf
                .all_alleles
                .iter()
                .map(|a| {
                    let p = psf.freqs1.get(a).copied().unwrap_or(0.0);
                    let q = psf.freqs2.get(a).copied().unwrap_or(0.0);
                    (p - q).abs()
                })
                .sum();
            0.5 * abs_diff
        })
        .sum();

    if normalize && num_loci > 0 {
        sum_dist / num_loci as f64
    } else {
        sum_dist
    }
}
