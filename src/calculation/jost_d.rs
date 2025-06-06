use crate::types::{AlleleFrequencies, PositionalData, SiteData};
use rayon::prelude::*;
use std::collections::{HashSet};

/// Calculates Jost's D, a measure of population differentiation.
/// This implementation sums the per-locus D values.
pub fn calculate_jost_d_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
    normalize: bool,
    num_loci: usize,
) -> f64 {
    let sum_dist: f64 = all_positions
        .par_iter()
        .map(|&pos| {
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            // Homozygosity = sum of squared allele frequencies
            let j1: f64 = freqs1.values().map(|p| p.powi(2)).sum();
            let j2: f64 = freqs2.values().map(|p| p.powi(2)).sum();

            // Mean homozygosity within populations
            let j_s = (j1 + j2) / 2.0;

            // Homozygosity of the mean allele frequencies
            let j_t: f64 = all_alleles
                .iter()
                .map(|allele| {
                    let p_i = *freqs1.get(allele).unwrap_or(&0.0);
                    let q_i = *freqs2.get(allele).unwrap_or(&0.0);
                    let mean_freq = (p_i + q_i) / 2.0;
                    mean_freq.powi(2)
                })
                .sum();
            
            // Per-locus Jost's D
            // Formula is D = (Jt - Js) / (1 - Js)
            if (1.0 - j_s) > 0.0 {
                (j_t - j_s) / (1.0 - j_s)
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
