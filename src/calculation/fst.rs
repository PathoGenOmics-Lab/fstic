use super::common::{get_all_freqs_at_pos, heterozygosity, pooled_heterozygosity};
use crate::types::PositionalData;
use rayon::prelude::*;
use std::collections::HashSet;

/// Calculates FST as sum-of-per-site Nei's GST: Σ (Ht - Hs) / Ht.
///
/// Note: this is Nei's (1973) approach, *not* the Weir & Cockerham (1984) θ estimator.
/// With `--normalize` the result is divided by the number of loci.
pub fn calculate_fst_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
    normalize: bool,
    num_loci: usize,
) -> f64 {
    let sum_fsts: f64 = all_positions
        .par_iter()
        .map(|&pos| {
            let psf = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));
            let h_s = (heterozygosity(&psf.freqs1) + heterozygosity(&psf.freqs2)) / 2.0;
            let h_t = pooled_heterozygosity(&psf);
            if h_t > 0.0 { (h_t - h_s) / h_t } else { 0.0 }
        })
        .sum();

    if normalize && num_loci > 0 {
        sum_fsts / num_loci as f64
    } else {
        sum_fsts
    }
}
