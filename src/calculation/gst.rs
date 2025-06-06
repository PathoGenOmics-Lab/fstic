use crate::types::{AlleleFrequencies, PositionalData, SiteData};
use rayon::prelude::*;
use std::collections::{HashSet};

/// Calculates Nei's global GST, as the ratio of sums over loci.
/// GST = sum(Ht - Hs) / sum(Ht). This is not affected by the --normalize flag.
pub fn calculate_gst_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
) -> f64 {
    // Calculate the numerator (Ht - Hs) and denominator (Ht) for each site,
    // then sum them up in parallel.
    let (total_numerator, total_denominator) = all_positions
        .par_iter()
        .map(|&pos| {
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            // Heterozygosity within each population
            let h1: f64 = 1.0 - freqs1.values().map(|p| p.powi(2)).sum::<f64>();
            let h2: f64 = 1.0 - freqs2.values().map(|p| p.powi(2)).sum::<f64>();

            // Average heterozygosity within populations
            let h_s = (h1 + h2) / 2.0;

            // Heterozygosity in the total "pooled" population
            let h_t: f64 = 1.0 - all_alleles
                .iter()
                .map(|allele| {
                    let p_i = freqs1.get(allele).unwrap_or(&0.0);
                    let q_i = freqs2.get(allele).unwrap_or(&0.0);
                    let mean_freq = (p_i + q_i) / 2.0;
                    mean_freq.powi(2)
                })
                .sum::<f64>();
            
            // Return a tuple of (numerator, denominator) for this site
            (h_t - h_s, h_t)
        })
        // Reduce the tuples by summing them element-wise
        .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
    
    // Perform the final division only once
    if total_denominator > 0.0 {
        total_numerator / total_denominator
    } else {
        0.0
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
