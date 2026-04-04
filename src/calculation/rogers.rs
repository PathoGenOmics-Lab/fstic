use crate::types::PositionalData;
use super::common::get_all_freqs_at_pos;
use rayon::prelude::*;
use std::collections::HashSet;

/// Calculates Rogers' distance (1972), a geometric distance based on Euclidean distance.
/// D = sqrt( (1/2L) * sum_loci(sum_alleles((pi - qi)^2)) )
/// This is a global metric, already normalized by L, so --normalize does not apply.
pub fn calculate_rogers_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
    num_loci: usize,
) -> f64 {
    if num_loci == 0 {
        return 0.0;
    }

    let total_sum_sq_diff: f64 = all_positions
        .par_iter()
        .map(|&pos| {
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            all_alleles
                .iter()
                .map(|allele| {
                    let p_i = *freqs1.get(allele).unwrap_or(&0.0);
                    let q_i = *freqs2.get(allele).unwrap_or(&0.0);
                    (p_i - q_i).powi(2)
                })
                .sum::<f64>()
        })
        .sum();

    (total_sum_sq_diff / (2.0 * num_loci as f64)).sqrt()
}
