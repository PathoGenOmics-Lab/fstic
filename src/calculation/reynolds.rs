use super::common::get_all_freqs_at_pos;
use crate::types::{GenomicPos, PositionalData};
use rayon::prelude::*;

/// Calculates Reynolds, Weir & Cockerham's (1983) distance: `D_R = -ln(1 - theta)`.
///
/// The coancestry coefficient is the ratio of sums
/// `theta = sum_l [ sum_i (p_i - q_i)^2 / 2 ] / sum_l [ 1 - sum_i p_i q_i ]`,
/// which is not Nei's GST: the two disagree at any locus that is not completely
/// fixed. Substituting GST gave numbers that could not be compared with adegenet,
/// Arlequin or GenAlEx output carrying the same name.
///
/// Not affected by `--normalize`.
pub fn calculate_reynolds_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &[GenomicPos],
) -> f64 {
    let (total_num, total_den) = all_positions
        .par_iter()
        .map(|pos| {
            let psf = get_all_freqs_at_pos(data1.get(pos), data2.get(pos));
            let mut sum_sq = 0.0_f64;
            let mut sum_prod = 0.0_f64;
            for allele in &psf.all_alleles {
                let p = psf.freqs1.get(allele).copied().unwrap_or(0.0);
                let q = psf.freqs2.get(allele).copied().unwrap_or(0.0);
                sum_sq += (p - q).powi(2);
                sum_prod += p * q;
            }
            (sum_sq / 2.0, 1.0 - sum_prod)
        })
        .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));

    // A zero denominator means every locus is monomorphic and shared, i.e. the two
    // samples are identical.
    if total_den <= 0.0 {
        return 0.0;
    }

    let theta = total_num / total_den;
    if theta < 1.0 {
        -(1.0 - theta).ln()
    } else {
        f64::INFINITY
    }
}
