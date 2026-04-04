use crate::types::PositionalData;
use super::common::get_all_freqs_at_pos;
use rayon::prelude::*;
use std::collections::HashSet;

/// Calculates Nei's standard genetic distance D = -ln(I).
/// This is a global metric and is not affected by the --normalize flag.
pub fn calculate_nei_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
) -> f64 {
    let (j_xy, j_x, j_y) = all_positions
        .par_iter()
        .map(|&pos| {
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            let mut site_j_xy = 0.0;
            let mut site_j_x = 0.0;
            let mut site_j_y = 0.0;

            for allele in all_alleles {
                let p_i = *freqs1.get(&allele).unwrap_or(&0.0);
                let q_i = *freqs2.get(&allele).unwrap_or(&0.0);
                site_j_xy += p_i * q_i;
                site_j_x += p_i.powi(2);
                site_j_y += q_i.powi(2);
            }
            (site_j_xy, site_j_x, site_j_y)
        })
        .reduce(|| (0.0, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));

    let denominator = (j_x * j_y).sqrt();
    if denominator > 0.0 {
        let identity = j_xy / denominator;
        if identity > 0.0 {
            -identity.ln()
        } else {
            f64::INFINITY
        }
    } else {
        f64::INFINITY
    }
}
