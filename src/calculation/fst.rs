use crate::types::{AlleleFrequencies, PositionalData, SiteData};
use rayon::prelude::*;
use std::collections::{HashSet};

/// Calculates a cumulative FST-like metric by summing per-site Nei's GST values.
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
            calculate_nei_gst_at_site(data1.get(&pos), data2.get(&pos))
        })
        .sum();
    
    if normalize && num_loci > 0 {
        sum_fsts / num_loci as f64
    } else {
        sum_fsts
    }
}

/// Calculates Nei's GST for a single site: (Ht - Hs) / Ht.
/// This now correctly handles all alleles, including indels.
fn calculate_nei_gst_at_site(
    site1_data: Option<&SiteData>,
    site2_data: Option<&SiteData>,
) -> f64 {
    let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(site1_data, site2_data);

    // Heterozygosity (gene diversity) within each population
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
    
    if h_t > 0.0 {
        (h_t - h_s) / h_t
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
