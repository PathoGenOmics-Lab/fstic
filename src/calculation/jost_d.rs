use crate::types::PositionalData;
use super::common::get_all_freqs_at_pos;
use rayon::prelude::*;
use std::collections::HashSet;

/// Calculates Jost's D, a measure of true population differentiation.
/// For two populations: D = (n/(n-1)) * ((Ht - Hs) / (1 - Hs))
/// where n=2, Ht = total heterozygosity, Hs = mean within-population heterozygosity.
pub fn calculate_jost_d_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
    normalize: bool,
    num_loci: usize,
) -> f64 {
    let n = 2.0_f64;

    let sum_dist: f64 = all_positions
        .par_iter()
        .map(|&pos| {
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            let h1: f64 = 1.0 - freqs1.values().map(|p| p.powi(2)).sum::<f64>();
            let h2: f64 = 1.0 - freqs2.values().map(|p| p.powi(2)).sum::<f64>();
            let h_s = (h1 + h2) / n;

            let h_t: f64 = 1.0 - all_alleles
                .iter()
                .map(|allele| {
                    let p_i = *freqs1.get(allele).unwrap_or(&0.0);
                    let q_i = *freqs2.get(allele).unwrap_or(&0.0);
                    let mean_freq = (p_i + q_i) / n;
                    mean_freq.powi(2)
                })
                .sum::<f64>();

            if (1.0 - h_s) > 0.0 {
                (n / (n - 1.0)) * (h_t - h_s) / (1.0 - h_s)
            } else {
                0.0
            }
        })
        .sum();

    if normalize && num_loci > 0 {
        sum_dist / num_loci as f64
    } else {
        sum_dist
    }
}
