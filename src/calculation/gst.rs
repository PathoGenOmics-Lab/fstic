use crate::types::PositionalData;
use super::common::get_all_freqs_at_pos;
use rayon::prelude::*;
use std::collections::HashSet;

/// Calculates Nei's global GST, as the ratio of sums over loci.
/// GST = sum(Ht - Hs) / sum(Ht). This is not affected by the --normalize flag.
pub fn calculate_gst_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
) -> f64 {
    let (total_numerator, total_denominator) = all_positions
        .par_iter()
        .map(|&pos| {
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            let h1: f64 = 1.0 - freqs1.values().map(|p| p.powi(2)).sum::<f64>();
            let h2: f64 = 1.0 - freqs2.values().map(|p| p.powi(2)).sum::<f64>();
            let h_s = (h1 + h2) / 2.0;

            let h_t: f64 = 1.0 - all_alleles
                .iter()
                .map(|allele| {
                    let p_i = freqs1.get(allele).unwrap_or(&0.0);
                    let q_i = freqs2.get(allele).unwrap_or(&0.0);
                    let mean_freq = (p_i + q_i) / 2.0;
                    mean_freq.powi(2)
                })
                .sum::<f64>();

            (h_t - h_s, h_t)
        })
        .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));

    if total_denominator > 0.0 {
        total_numerator / total_denominator
    } else {
        0.0
    }
}
