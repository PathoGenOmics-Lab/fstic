use crate::types::{AlleleFrequencies, PositionalData, SiteData};
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
    let n = 2.0_f64; // number of populations

    let sum_dist: f64 = all_positions
        .par_iter()
        .map(|&pos| {
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            // Heterozygosity within each population: H = 1 - sum(p_i^2)
            let h1: f64 = 1.0 - freqs1.values().map(|p| p.powi(2)).sum::<f64>();
            let h2: f64 = 1.0 - freqs2.values().map(|p| p.powi(2)).sum::<f64>();

            // Mean within-population heterozygosity
            let h_s = (h1 + h2) / n;

            // Total heterozygosity from pooled allele frequencies
            let h_t: f64 = 1.0 - all_alleles
                .iter()
                .map(|allele| {
                    let p_i = *freqs1.get(allele).unwrap_or(&0.0);
                    let q_i = *freqs2.get(allele).unwrap_or(&0.0);
                    let mean_freq = (p_i + q_i) / n;
                    mean_freq.powi(2)
                })
                .sum::<f64>();

            // Per-locus Jost's D = (n/(n-1)) * ((Ht - Hs) / (1 - Hs))
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
