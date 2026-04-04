use super::common::get_all_freqs_at_pos;
use crate::types::{GenomicPos, PositionalData};
use rayon::prelude::*;

/// Calculates Cavalli-Sforza & Edwards' chord distance.
///
/// Per-locus: D_ch = √(2 (1 − Σ √(p_i q_i))).
/// With `--normalize` the sum is divided by the number of loci.
pub fn calculate_chord_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &[GenomicPos],
    normalize: bool,
    num_loci: usize,
) -> f64 {
    let sum_dist: f64 = all_positions
        .par_iter()
        .map(|pos| {
            let psf = get_all_freqs_at_pos(data1.get(pos), data2.get(pos));
            let sum_sqrt: f64 = psf.all_alleles.iter().map(|a| {
                let p = psf.freqs1.get(a).copied().unwrap_or(0.0);
                let q = psf.freqs2.get(a).copied().unwrap_or(0.0);
                (p * q).sqrt()
            }).sum();
            let term = 2.0 * (1.0 - sum_sqrt);
            if term > 0.0 { term.sqrt() } else { 0.0 }
        })
        .sum();

    if normalize && num_loci > 0 { sum_dist / num_loci as f64 } else { sum_dist }
}
