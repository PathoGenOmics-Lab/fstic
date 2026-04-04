use super::common::{get_all_freqs_at_pos, heterozygosity, pooled_heterozygosity};
use crate::types::{GenomicPos, PositionalData};
use rayon::prelude::*;

const N_POPS: f64 = 2.0;

/// Calculates Jost's D (2008): D = (n/(n-1)) × (Ht - Hs) / (1 - Hs), for n=2 populations.
///
/// With `--normalize` the sum is divided by the number of loci.
pub fn calculate_jost_d_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &[GenomicPos],
    normalize: bool,
    num_loci: usize,
) -> f64 {
    let correction = N_POPS / (N_POPS - 1.0);

    let sum_d: f64 = all_positions
        .par_iter()
        .map(|pos| {
            let psf = get_all_freqs_at_pos(data1.get(pos), data2.get(pos));
            let h_s = (heterozygosity(&psf.freqs1) + heterozygosity(&psf.freqs2)) / 2.0;
            let h_t = pooled_heterozygosity(&psf);
            let denom = 1.0 - h_s;
            if denom > 0.0 { correction * (h_t - h_s) / denom } else { 0.0 }
        })
        .sum();

    if normalize && num_loci > 0 { sum_d / num_loci as f64 } else { sum_d }
}
