use crate::types::{AlleleFrequencies, PositionalData, SiteData};
use rayon::prelude::*;
use std::collections::{HashSet};

/// Calculates Rogers' distance (1972), a geometric distance based on Euclidean distance.
/// D = sqrt( (1/2L) * sum_loci(sum_alleles((pi - qi)^2)) )
/// This is a global metric, already normalized by L, so --normalize does not apply.
pub fn calculate_rogers_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
    num_loci: usize,
) -> f64 {
    if num_loci == 0 {
        return 0.0;
    }

    let total_sum_sq_diff: f64 = all_positions
        .par_iter()
        .map(|&pos| {
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            let site_sum_sq_diff: f64 = all_alleles
                .iter()
                .map(|allele| {
                    let p_i = *freqs1.get(allele).unwrap_or(&0.0);
                    let q_i = *freqs2.get(allele).unwrap_or(&0.0);
                    (p_i - q_i).powi(2)
                })
                .sum();
            
            site_sum_sq_diff
        })
        .sum();

    // Final calculation for Rogers' distance
    (total_sum_sq_diff / (2.0 * num_loci as f64)).sqrt()
}

/// Helper to get complete frequency maps for both samples at a site.
fn get_all_freqs_at_pos<'a>(
    site1_data: Option<&'a SiteData>,
    site2_data: Option<&'a SiteData>,
) -> (HashSet<String>, AlleleFrequencies, AlleleFrequencies) {
    let mut all_alleles: HashSet<String> = HashSet::new();
    let mut freqs1 = AlleleFrequencies::new();
    let mut freqs2 = AlleleFrequencies::new();
    
    let mut ref_allele = String::new();
    if let Some(data) = site1_data {
        ref_allele = data.reference_allele.clone();
    } else if let Some(data) = site2_data {
        ref_allele = data.reference_allele.clone();
    }

    if let Some(data) = site1_data {
        freqs1 = data.freqs.clone();
        all_alleles.extend(data.freqs.keys().cloned());
    }

    if let Some(data) = site2_data {
        freqs2 = data.freqs.clone();
        all_alleles.extend(data.freqs.keys().cloned());
    }
    
    if !ref_allele.is_empty() {
        all_alleles.insert(ref_allele.clone());
    }

    let sum1: f64 = freqs1.values().sum();
    if !ref_allele.is_empty() {
        freqs1.insert(ref_allele.clone(), 1.0 - sum1);
    }

    let sum2: f64 = freqs2.values().sum();
    if !ref_allele.is_empty() {
        freqs2.insert(ref_allele, 1.0 - sum2);
    }

    (all_alleles, freqs1, freqs2)
}
