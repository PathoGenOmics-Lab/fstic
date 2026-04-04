use super::common::get_all_freqs_at_pos;
use crate::types::PositionalData;
use rayon::prelude::*;
use std::collections::HashSet;

/// Calculates Rogers' distance (1972): D_R = √(Σ(p_i − q_i)² / 2L).
///
/// Already normalized by L, so `--normalize` does not apply.
pub fn calculate_rogers_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
    num_loci: usize,
) -> f64 {
    if num_loci == 0 {
        return 0.0;
    }

    let total_sum_sq: f64 = all_positions
        .par_iter()
        .map(|&pos| {
            let psf = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));
            psf.all_alleles
                .iter()
                .map(|a| {
                    let p = psf.freqs1.get(a).copied().unwrap_or(0.0);
                    let q = psf.freqs2.get(a).copied().unwrap_or(0.0);
                    (p - q).powi(2)
                })
                .sum::<f64>()
        })
        .sum();

    (total_sum_sq / (2.0 * num_loci as f64)).sqrt()
}
