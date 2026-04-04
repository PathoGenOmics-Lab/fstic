use crate::types::PositionalData;
use super::common::get_all_freqs_at_pos;
use rayon::prelude::*;
use std::collections::HashSet;

/// Calculates the Bray-Curtis dissimilarity, which for frequency data is
/// equivalent to the Absolute or Manhattan distance: D = 0.5 * sum(|pi - qi|).
/// This implementation sums the per-locus distances.
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
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            let sum_abs_diff: f64 = all_alleles
                .iter()
                .map(|allele| {
                    let p_i = *freqs1.get(allele).unwrap_or(&0.0);
                    let q_i = *freqs2.get(allele).unwrap_or(&0.0);
                    (p_i - q_i).abs()
                })
                .sum();

            0.5 * sum_abs_diff
        })
        .sum();

    if normalize && num_loci > 0 {
        sum_dist / num_loci as f64
    } else {
        sum_dist
    }
}
