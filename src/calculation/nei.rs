use super::common::get_all_freqs_at_pos;
use crate::types::PositionalData;
use rayon::prelude::*;
use std::collections::HashSet;

/// Calculates Nei's standard genetic distance: D = −ln(I).
///
/// Genetic identity I = Jxy / √(Jx · Jy), where J = Σ p_i · q_i.
/// Not affected by `--normalize`.
pub fn calculate_nei_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
) -> f64 {
    let (j_xy, j_x, j_y) = all_positions
        .par_iter()
        .map(|&pos| {
            let psf = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));
            let mut sxy = 0.0_f64;
            let mut sx = 0.0_f64;
            let mut sy = 0.0_f64;
            for allele in &psf.all_alleles {
                let p = psf.freqs1.get(allele).copied().unwrap_or(0.0);
                let q = psf.freqs2.get(allele).copied().unwrap_or(0.0);
                sxy += p * q;
                sx += p * p;
                sy += q * q;
            }
            (sxy, sx, sy)
        })
        .reduce(|| (0.0, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));

    let denom = (j_x * j_y).sqrt();
    if denom > 0.0 {
        let identity = j_xy / denom;
        if identity > 0.0 {
            -identity.ln()
        } else {
            f64::INFINITY
        }
    } else {
        f64::INFINITY
    }
}
