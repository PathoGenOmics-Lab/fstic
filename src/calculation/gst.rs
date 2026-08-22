use super::common::{get_all_freqs_at_pos, heterozygosity, pooled_heterozygosity, sum_per_locus_n};
use crate::types::{GenomicPos, PositionalData};

/// Calculates Nei's global GST as ratio-of-sums: Σ(Ht - Hs) / Σ(Ht).
///
/// Not affected by `--normalize`.
pub fn calculate_gst_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &[GenomicPos],
) -> f64 {
    let [total_num, total_den] = sum_per_locus_n(all_positions, |pos| {
        let psf = get_all_freqs_at_pos(data1.get(pos), data2.get(pos));
        let h_s = (heterozygosity(&psf.freqs1) + heterozygosity(&psf.freqs2)) / 2.0;
        let h_t = pooled_heterozygosity(&psf);
        [h_t - h_s, h_t]
    });

    if total_den > 0.0 {
        total_num / total_den
    } else {
        0.0
    }
}
